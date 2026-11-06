* args benchmark_directory new_external_output_directory
args benchmark output
version 13
set more off
capture mkdir "`output'"
do "`benchmark'/fixtures.do" "`output'" 200000 4 "akm_low akm_medium akm_high"
foreach name in akm_low akm_medium akm_high {
    use "`output'/`name'.dta", clear
    assert _N==200000
    isid worker time
    assert inrange(worker,1,20000) & inrange(firm,1,1000) & inrange(time,1,10)
    bysort worker (time): assert _N==10 & time==_n
    assert home==1+mod(worker-1,1000)
    assert firm==home & moved==0 if time==1
    by worker (time): assert moved==(firm!=firm[_n-1]) if time>1
    by worker (time): assert previous_firm==firm[_n-1] if time>1
    local rate = cond("`name'"=="akm_low",.02,cond("`name'"=="akm_medium",.10,.30))
    assert moved==(time>1 & move_draw<`rate')
    quietly count if moved
    local observed_moves = r(N)
    assert abs(`observed_moves'/180000-`rate')<.002
    local stored_moves : char _dta[moves]
    assert `observed_moves'==`stored_moves'
    local components : char _dta[components]
    assert inrange(`components',1,1000)
    assert y>0 & y<. & !missing(x1,x2,x3,x4)
    assert reldif(y,exp(eta+.3*outcome_shock-.045))<1e-14
    bysort firm: gen byte first_firm = _n==1
    quietly count if first_firm
    assert r(N)==1000
    sort worker time
    cf worker home time rowid move_draw destination_draw x1 x2 x3 x4 outcome_shock ///
        using "`output'/akm_low.dta"
    display "AKM_TEST_OK `name' moves=`observed_moves' components=`components'"
}
display "AKM_FIXTURE_TESTS_OK"
