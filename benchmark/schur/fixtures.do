* Deterministic fixtures; run once, outside every estimation timer.
* args output_directory requested_N max_controls selected_cases
* N is a multiple of 5000, at least 10000; AKM cases require exactly 200000.
* An omitted case list generates the four retained 100,000-row design families.
args output nobs max_controls selected_cases
version 13
clear all
set more off
if "`max_controls'"=="" local max_controls 4
if "`selected_cases'"=="" local selected_cases "chain_2 communities_3 imbalanced_4 balanced_4"
display "REDBENCH_START"
local workers = `nobs'/10
local firms = `workers'/20
foreach name in chain_2 communities_3 imbalanced_4 balanced_4 {
    local design = cond("`name'"=="chain_2",2,cond("`name'"=="communities_3",5,cond("`name'"=="imbalanced_4",7,8)))
    if "`selected_cases'"!="" & !strpos(" `selected_cases' "," `name' ") continue
    clear
    set seed 98241
    set obs `workers'
    gen long worker = _n
    gen long mover_index = ceil(worker/5)
    gen byte mover = mod(worker,5)==0
    gen int home = 1+mod(mover_index-1,`firms')
    gen int away = 1+mod(home,`firms')
    replace away = 1+floor(runiform()*`firms') if mover_index>`firms'
    replace away = 1+mod(home,`firms') if away==home

    * Chain: all mobility stays between adjacent firms.
    if `design'==2 replace away = cond(home==`firms',home-1,home+1)

    * Communities have 25 firms, connected by a sparse ring of bridges.
    if `design'==5 {
        replace away = 25*floor((home-1)/25)+1+mod(home,25)
        replace away = 25*floor((home-1)/25)+1+floor(runiform()*25) if mover_index>`firms'
        replace away = 25*floor((home-1)/25)+1+mod(home,25) if away==home
        if `design'==5 replace away = 1+mod(home,`firms') if mover_index<=`firms' & mod(home,25)==0
    }
    expand 10
    bysort worker: gen byte time = _n
    gen int firm = cond(mover & time>5,away,home)
    gen int region = ceil(firm/25)
    gen int product = 1+floor(runiform()*40)
    gen long rowid = _n
    gen double offset = .1*rnormal()
    gen double pweight = exp(rnormal())
    gen int fweight = 1+mod(rowid,4)
    gen double x1 = rnormal()
    gen double x2 = rnormal()
    gen double x3 = rnormal()
    gen double x4 = rnormal()
    gen double eta = .2*sin(worker)+.3*cos(firm)+.1*(x1+x2+x3+x4)
    if inlist(`design',5,7) replace eta = eta+.1*cos(time)
    if `design'==7 replace eta = eta+.1*sin(product)+offset
    gen double y = exp(eta+.3*rnormal())
    if `design'==7 replace y = rpoisson(exp(eta))

    * Balanced four-way outcome; explicit scaling controls are zero-effect noise.
    if `design'==8 {
        forvalues block=1/4 {
            gen int balanced`block' = 1+floor(runiform()*`firms')
        }
        replace eta = .1*sin(balanced1)+.1*cos(balanced2)+.1*sin(balanced3)+.1*cos(balanced4)
        replace y = rpoisson(exp(eta))
    }
    * Add zero-effect noise controls after ALL original outcome/factor draws.
    * This independent stream preserves x1-x4, y and topology exactly.
    if `max_controls'>4 {
        set seed 183271
        forvalues column=5/`max_controls' {
            gen double x`column' = rnormal()
        }
    }
    char _dta[max_controls] "`max_controls'"
    char _dta[extra_controls_seed] "183271"
    char _dta[design] "`name'"
    char _dta[seed] "98241"
    char _dta[requested_N] "`nobs'"
    compress
    save "`output'/`name'.dta", replace
    display "FIXTURE_OK `name' N=" _N " workers=`workers' firms=`firms'"
}


mata:
real scalar ppmlbench_akm_components()
{
    real matrix links
    real colvector parent
    real scalar i, a, b
    links = st_data(.,("previous_firm","firm"),"moved")
    parent = (1::1000)
    for (i=1; i<=rows(links); i++) {
        a = links[i,1]
        b = links[i,2]
        while (parent[a]!=a) {
            parent[a] = parent[parent[a]]
            a = parent[a]
        }
        while (parent[b]!=b) {
            parent[b] = parent[parent[b]]
            b = parent[b]
        }
        parent[a] = b
    }
    return(sum(parent:==(1::1000)))
}
end

* Matched employer/employee panels. All three mobility rates share the same
* initial employers, move/destination draws, covariates and outcome shocks.
* Mobility is a job change to a DIFFERENT firm in each of nine transitions.
foreach name in akm_low akm_medium akm_high {
    if "`selected_cases'"!="" & !strpos(" `selected_cases' "," `name' ") continue
    assert `nobs'==200000
    local mobility = cond("`name'"=="akm_low",.02,cond("`name'"=="akm_medium",.10,.30))
    clear
    set seed 73192
    set obs 20000
    gen long worker = _n
    gen int home = 1+mod(worker-1,1000)
    expand 10
    bysort worker: gen byte time = _n
    sort worker time
    gen long rowid = _n
    gen double move_draw = runiform()
    gen double destination_draw = runiform()
    gen byte moved = time>1 & move_draw<`mobility'
    gen int firm = home
    forvalues year = 2/10 {
        by worker (time): replace firm = cond(moved, ///
            1+mod(firm[_n-1]+floor(999*destination_draw),1000),firm[_n-1]) if time==`year'
    }
    by worker (time): gen int previous_firm = firm[_n-1] if time>1
    by worker: egen byte ever_moved = max(moved)
    forvalues column = 1/4 {
        gen double x`column' = rnormal()
    }
    gen double outcome_shock = rnormal()
    gen double eta = .2*sin(worker)+.3*cos(firm)+.1*cos(time)+.1*(x1+x2+x3+x4)
    * Positive wage-like outcome with E(y | worker, firm, time, X) = exp(eta).
    gen double y = exp(eta+.3*outcome_shock-.045)
    if `max_controls'>4 {
        set seed 183271
        forvalues column = 5/`max_controls' {
            gen double x`column' = rnormal()
        }
    }
    quietly count if moved
    local moves = r(N)
    quietly count if time==1 & ever_moved
    local movers = r(N)
    * Connected components of the firm mobility graph, including isolated firms.
    * This is a fixture diagnostic, not part of the estimator or its timer.
    mata: st_numscalar("akm_components",ppmlbench_akm_components())
    char _dta[design] "`name'"
    char _dta[requested_N] "200000"
    char _dta[workers] "20000"
    char _dta[firms] "1000"
    char _dta[years] "10"
    char _dta[annual_mobility] "`mobility'"
    char _dta[moves] "`moves'"
    char _dta[movers] "`movers'"
    char _dta[components] "`=scalar(akm_components)'"
    char _dta[max_controls] "`max_controls'"
    char _dta[seed] "73192"
    char _dta[extra_controls_seed] "183271"
    compress
    save "`output'/`name'.dta", replace
    display "FIXTURE_OK `name' N=" _N " workers=20000 firms=1000 years=10 annual_mobility=`mobility' moves=`moves' movers=`movers' components=" scalar(akm_components)
}
display "REDBENCH_FIXTURES_OK"
