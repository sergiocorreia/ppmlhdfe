* Deterministic fixtures; run once, outside every estimation timer.
* args output_directory requested_N (a multiple of 5000, at least 10000)
args output nobs
version 13
clear all
set more off
display "REDBENCH_START"
local workers = `nobs'/10
local firms = `workers'/20
forvalues design = 1/8 {
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

    * Communities have 25 firms. Case 3 never crosses a community boundary;
    * case 5 connects consecutive communities with a sparse ring of bridges.
    if inlist(`design',3,5) {
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
    if inlist(`design',4,5,7) replace eta = eta+.1*cos(time)
    if `design'==7 replace eta = eta+.1*sin(product)+offset
    gen double y = exp(eta+.3*rnormal())
    if inlist(`design',4,7) replace y = rpoisson(exp(eta))
    if `design'==4 replace y = 0 if worker<=2

    * Balanced, highly connected four-way case; no covariates are estimated.
    if `design'==8 {
        forvalues block=1/4 {
            gen int balanced`block' = 1+floor(runiform()*`firms')
        }
        replace eta = .1*sin(balanced1)+.1*cos(balanced2)+.1*sin(balanced3)+.1*cos(balanced4)
        replace y = rpoisson(exp(eta))
    }
    local names "random_2 chain_2 disconnected_2 random_3 communities_3 redundant_3 imbalanced_4 balanced_4"
    local name : word `design' of `names'
    char _dta[design] "`name'"
    char _dta[seed] "98241"
    char _dta[requested_N] "`nobs'"
    compress
    save "`output'/`name'.dta", replace
    display "FIXTURE_OK `name' N=" _N " workers=`workers' firms=`firms'"
}
display "REDBENCH_FIXTURES_OK"
