* One source per fresh Stata process: no baseline/candidate switching in Mata.
args source fixture output mode repetitions precision processors control_count
version 13
clear all
set more off
if "`processors'"=="" local processors 1
capture set processors `processors'
if _rc | c(processors)!=`processors' {
    display as error "Requested processors `processors' unavailable; actual " c(processors)
    exit 198
}
display "REDBENCH_START"
adopath ++ "`source'/src"
findfile ppmlhdfe.ado
display "SOURCE_ADO `r(fn)'"
findfile ppmlhdfe.mata
display "SOURCE_MATA `r(fn)'"
capture noisily which ftools
capture noisily which reghdfe
capture findfile reghdfe5.mata
if !_rc display "DEPENDENCY_REGHDFE5 `r(fn)'"
capture findfile ftools.mata
if !_rc display "DEPENDENCY_FTOOLS `r(fn)'"
display "STATA_VERSION " c(stata_version) " FLAVOR " c(flavor) " PROCESSORS " c(processors)
use "`fixture'", clear
local design : char _dta[design]
local absorb "worker firm"
if inlist("`design'","communities_3","akm_low","akm_medium","akm_high") local absorb "worker firm time"
if "`design'"=="imbalanced_4" local absorb "worker firm time product"
if "`design'"=="balanced_4" local absorb "balanced1 balanced2 balanced3 balanced4"
local controls "x1 x2 x3 x4"
if "`design'"=="balanced_4" local controls ""
if "`control_count'"!="" {
    local controls ""
    if `control_count'>0 {
        forvalues column=1/`control_count' {
            confirm variable x`column'
            local controls "`controls' x`column'"
        }
    }
}
local actual_controls : word count `controls'
display "BENCHMARK_CONTROLS `actual_controls'"
local weights ""
if "`design'"=="communities_3" local weights "[pw=pweight]"
if "`design'"=="imbalanced_4" local weights "[fw=fweight]"
local extra ""
if "`design'"=="imbalanced_4" local extra "offset(offset)"
if "`precision'"=="" local precision "tol(1e-10)"
local option ""
if "`mode'"=="candidate" local option "use_schur(1)"
if "`mode'"=="default" local option "use_schur(0)"
local first = 0
if "`mode'"=="default" {
    local first = 1
    local repetitions = 1
}
if "`mode'"=="tight" {
    local first = 1
    local repetitions = 1
    local precision "tol(1e-12) use_exact_partial(1) start_inner_tol(1e-12) itol(1e-12)"
}
display "BENCHMARK_PRECISION `precision'"

tempname results
local diagnostics "N deviance converged ic ic2 schur_calls schur_iterations schur_fallbacks schur_col_calls schur_col_iters schur_col_fails num_separated num_singletons"
postfile `results' str24 design str10 mode str10 phase int rep int rc double seconds processors control_count `diagnostics' using "`output'/timings.dta", replace
tempname csv
file open `csv' using "`output'/timings.csv", write replace
local diagnostics_csv : subinstr local diagnostics " " ",", all
file write `csv' "design,mode,phase,rep,rc,seconds,processors,control_count,`diagnostics_csv'" _n
file close `csv'
forvalues rep = `first'/`repetitions' {
    local phase = cond(`rep'==0,"warmup","timed")
    if inlist("`mode'","tight","default") local phase "validation"
    tempname stage
    file open `stage' using "`output'/stage.txt", write replace
    file write `stage' "BEGIN `phase' `rep'" _n
    file close `stage'
    capture drop __d
    timer clear 1
    if "`phase'"!="validation" timer on 1
    capture quietly ppmlhdfe y `controls' `weights', absorb(`absorb') `extra' `precision' `option' d(__d)
    local rc = _rc
    scalar bench_seconds = .
    if "`phase'"!="validation" {
        timer off 1
        quietly timer list 1
        scalar bench_seconds = r(t1)
    }
    local posted ""
    foreach field of local diagnostics {
        scalar bench_`field' = .
        if !`rc' capture scalar bench_`field' = e(`field')
        local posted "`posted' (bench_`field')"
    }
    post `results' ("`design'") ("`mode'") ("`phase'") (`rep') (`rc') (bench_seconds) (c(processors)) (`actual_controls') `posted'
    * Close after each row so completed fits survive a later fit timeout.
    file open `csv' using "`output'/timings.csv", write append
    file write `csv' "`design',`mode',`phase',`rep',`rc'," (bench_seconds) "," (c(processors)) ",`actual_controls'"
    foreach field of local diagnostics {
        file write `csv' "," (bench_`field')
    }
    file write `csv' _n
    file close `csv'
    display "FIT_RESULT `design' `mode' `phase' `rep' rc=`rc' seconds=" bench_seconds
    if `rc' continue, break
    if `rep'>0 {
        estimates save "`output'/fit_`rep'.ster", replace
        preserve
        gen byte sample = e(sample)
        predict double mu if sample, mu
        gen double comparison_weight = 1
        if "`design'"=="communities_3" replace comparison_weight=pweight
        if "`design'"=="imbalanced_4" replace comparison_weight=fweight
        keep rowid sample mu comparison_weight
        sort rowid
        save "`output'/fit_`rep'.dta", replace
        restore
    }
}
postclose `results'
use "`output'/timings.dta", clear
export delimited using "`output'/timings.csv", replace
tempname stage
file open `stage' using "`output'/stage.txt", write replace
file write `stage' "DONE" _n
file close `stage'
display "REDBENCH_FITS_OK"
