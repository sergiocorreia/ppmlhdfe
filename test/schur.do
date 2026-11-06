* End-to-end compatibility of the optional categorical Schur projector.
* From the checkout root: do test/schur.do .
args root
if `"`root'"' == "" local root ..
clear all
set more off
set linesize 180
capture set processors 1
adopath ++ "`root'/src"
set seed 380623

program define schur_compare
    syntax [varlist(default=none)], Absorb(string) [Options(string) Weight(string) Eligible(integer 1) Mustseparate(integer 0) Mustaccelerate(integer 0) Comparetolerance(real 1e-11)]
    tempname b V dev ll n sep dofa rank
    tempvar samp mu0 mu1 d0 d1
    quietly ppmlhdfe y `varlist' `weight', a(`absorb') `options' d(`d0') tol(`comparetolerance') use_schur(0)
    matrix `b'=e(b)
    matrix `V'=e(V)
    local cn0 : colfullnames `b'
    scalar `dev'=e(deviance)
    scalar `ll'=e(ll)
    scalar `n'=e(N)
    scalar `sep'=e(num_separated)
    scalar `dofa'=e(df_a)
    scalar `rank'=e(rank)
    local singletons=e(num_singletons)
    _ms_omit_info e(b)
    tempname omitted
    matrix `omitted'=r(omit)
    gen byte `samp'=e(sample)
    predict double `mu0' if e(sample), mu
    quietly ppmlhdfe y `varlist' `weight', a(`absorb') `options' d(`d1') tol(`comparetolerance') use_schur(1)
    local cn1 : colfullnames e(b)
    assert "`cn0'"=="`cn1'"
    assert e(N)==`n'
    assert e(num_separated)==`sep'
    assert e(df_a)==`dofa'
    assert e(rank)==`rank'
    assert e(num_singletons)==`singletons'
    _ms_omit_info e(b)
    assert mreldif(`omitted',r(omit))==0
    assert `samp'==e(sample)
    assert mreldif(`b',e(b))<1e-7
    assert mreldif(`V',e(V))<1e-7
    assert reldif(`dev',e(deviance))<1e-9
    assert reldif(`ll',e(ll))<1e-9
    predict double `mu1' if e(sample), mu
    assert abs(`mu1'-`mu0')/(1+abs(`mu0'))<1e-5 if `samp'
    assert (e(schur_calls)>0)==`eligible'
    assert e(schur_fallbacks)<=e(schur_calls)
    assert e(schur_col_fails)<=e(schur_col_calls)
    if `mustaccelerate' {
        assert e(schur_calls)>e(schur_fallbacks)
        assert e(schur_col_calls)>e(schur_col_fails)
    }
    if `mustseparate' assert e(num_separated)>0
    di as result "SCHUR_PASS absorb=`absorb' vars=`varlist' options=`options' weight=`weight' calls=" e(schur_calls) " sep=" e(num_separated)
end

set obs 2400
gen long id=ceil(_n/6)
gen int firm=1+mod(floor(17*runiform())+_n,31)
gen byte year=1+mod(_n,6)
gen byte region=mod(id,4)
gen double x=rnormal()
gen double z=rnormal()
gen double off=.4*sin(_n)
gen double exposure=exp(off)
gen double w=exp(3*runiform()-1.5)
gen int fw=1+mod(_n,5)
gen double sumfe=sin(id)+cos(firm)+.3*year
gen double y=exp(.2*x-.1*z+.15*sumfe+.2*rnormal())
gen double duplicate=2*x
gen double zero=0
gen byte constant=1
gen double sortkey=runiform()
sort sortkey

schur_compare x z, absorb(id firm) mustaccelerate(1)
schur_compare x z, absorb(id firm year) mustaccelerate(1)
schur_compare x z, absorb(year firm id region) options("offset(off) vce(cluster id year)") weight("[pw=w]")
schur_compare x z, absorb(firm id year) options("exposure(exposure) guess(ols)") weight("[fw=fw]")
schur_compare x z, absorb(id firm year) options("use_exact_solver(1) use_exact_partial(1) use_heuristic_tol(0)")
schur_compare x z, absorb(id firm year) options("accel(lsmr)")
schur_compare x z, absorb(id firm year) options("accel(none) transform(sym)")
schur_compare x z, absorb(id firm year) options("standardize_data(0)")
schur_compare x z duplicate zero constant sumfe, absorb(id firm year)
schur_compare sumfe zero constant, absorb(id firm year)
schur_compare, absorb(id firm year)
schur_compare x z, absorb(id firm region)
schur_compare x z, absorb(id id firm year)
schur_compare x z, absorb(region#year firm id)
schur_compare x z, absorb(id firm year##c.off) eligible(0)
schur_compare x z, absorb(id firm year#c.off) eligible(0)
schur_compare x z, absorb(id) eligible(0)
schur_compare x z, absorb(id firm, prune) eligible(0)
schur_compare x z, absorb(id firm year) options("remove_collinear_variables(0)")

* First IRLS pass is the important path when initial check is disabled.
gen double nearly=sumfe+1e-7*rnormal()
schur_compare x z nearly, absorb(id firm year) options("remove_collinear_variables(0)")
replace nearly=sumfe+1e-4*rnormal()
schur_compare x z nearly, absorb(id firm year) options("remove_collinear_variables(0)")

* This checks off/on compatibility at globally tiny weights, not full
* estimator scale invariance. Upstream's deviance floor (.1) and inner
* stopping rules are scale sensitive; even two agreeing paths can differ
* from a normal-weight fit. Use tighter convergence for BOTH paths and keep
* the same coefficient/covariance/prediction gates. The independent dummy
* oracle in schur_projection.do separately tests kernel weight invariance.
* With exact partialling, upstream may tighten HDFE tolerance below its usual
* 1e-12 target floor. This is existing behavior, not a Schur change.
replace w=w*1e-12
schur_compare x z, absorb(id firm year) weight("[pw=w]") comparetolerance(1e-23) options("use_exact_partial(1) start_inner_tol(1e-12) itol(1e-12)")
replace w=w*1e24
schur_compare x z, absorb(id firm year) weight("[pw=w]")
replace w=w*1e-12

* FE storage is identified only up to normalization; check total saved FE.
quietly ppmlhdfe y x z, absorb(feid=id fefirm=firm feyear=year) use_schur(1) d(fetotal) tol(1e-11)
gen double recovered=feid+fefirm+feyear
assert abs(recovered-fetotal)<1e-6 if e(sample)
drop feid fefirm feyear fetotal recovered

* Singleton removal changes the factor/sample alignment before initialization.
preserve
replace id=10000+_n if _n>2395
schur_compare x z, absorb(id firm year)
restore

* Explicit mu-only separation forces post-IRLS sample rebuild.
replace y=0 if id<=3
schur_compare x z, absorb(id firm year) options("separation(mu)") mustseparate(1)

* ReLU and simplex trimming followed by Schur use on rebuilt factors.
gen byte sep=(id<=15 & year==1)
replace y=0 if sep
schur_compare x z sep, absorb(id firm year) options("separation(relu)") mustseparate(1)
schur_compare x z sep, absorb(id firm year) options("separation(simplex mu)") mustseparate(1)

* Invalid option values must fail and must not corrupt the next fit.
capture quietly ppmlhdfe y x, absorb(id firm year) use_schur(2)
assert _rc!=0
schur_compare x z sep, absorb(id firm year)

di as result "SCHUR_REGRESSION_OK"

program drop schur_compare
