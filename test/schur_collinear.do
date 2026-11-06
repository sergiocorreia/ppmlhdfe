* Omission decisions under exact and near FE/regressor collinearity.
* From the checkout root: do test/schur_collinear.do .
args root
if `"`root'"' == "" local root ..
clear all
set more off
capture set processors 1
adopath ++ "`root'/src"
set seed 843911
set obs 3000
gen long worker=ceil(_n/5)
gen int firm=1+floor(40*runiform())
gen byte year=mod(_n,7)
gen double x=rnormal()
gen double z=rnormal()
gen double duplicate=2*x
gen double mixed=3*x-2*z
gen double absorbed=sin(worker)+cos(firm)+.2*year
gen double noise=rnormal()
gen double w=1+mod(_n,4)
gen double off=.1*rnormal()
gen double y=exp(.2*x-.1*z+.1*absorbed+.1*rnormal())

program define schur_collinear_compare
    syntax varlist, [Extra(string) Weight(string) Fallback(integer 0)]
    tempname b V omit dev
    tempvar d0 d1 sample mu0 mu1
    quietly ppmlhdfe y `varlist' `weight', a(worker firm year) `extra' d(`d0') tol(1e-10)
    matrix `b'=e(b)
    matrix `V'=e(V)
    local names0 : colfullnames e(b)
    _ms_omit_info e(b)
    matrix `omit'=r(omit)
    scalar `dev'=e(deviance)
    gen byte `sample'=e(sample)
    predict double `mu0' if e(sample), mu
    quietly ppmlhdfe y `varlist' `weight', a(worker firm year) `extra' d(`d1') tol(1e-10) use_schur(1)
    local names1 : colfullnames e(b)
    _ms_omit_info e(b)
    assert "`names0'"=="`names1'"
    assert mreldif(`omit',r(omit))==0
    assert mreldif(`b',e(b))<1e-7
    assert mreldif(`V',e(V))<1e-7
    assert reldif(`dev',e(deviance))<1e-9
    assert `sample'==e(sample)
    predict double `mu1' if e(sample), mu
    assert abs(`mu1'-`mu0')/(1+abs(`mu0'))<1e-5 if `sample'
    assert e(schur_col_calls)==1
    if `fallback' assert e(schur_col_fails)==1
    di "COLLINEAR_PASS `varlist'; `weight'; `extra'"
end

schur_collinear_compare x z duplicate mixed absorbed, fallback(1)
schur_collinear_compare absorbed mixed duplicate z x, weight("[pw=w]") extra("offset(off) vce(cluster worker)") fallback(1)
schur_collinear_compare x z duplicate mixed absorbed, weight("[fw=w]") fallback(1)
schur_collinear_compare absorbed, fallback(1)
gen double near=absorbed+1e-5*noise
schur_collinear_compare x z near, fallback(1)
replace near=absorbed+1e-8*noise
schur_collinear_compare x z near, fallback(1)
replace near=absorbed+.01*noise
schur_collinear_compare x z near
replace near=x+1e-7*noise
schur_collinear_compare x z near, fallback(1)
* A correlation guard must also protect the separate small invsym rank decision.
replace near=x+1e-5*noise
schur_collinear_compare x z near, fallback(1)
gen double shuffle=runiform()
sort shuffle
schur_collinear_compare x z duplicate mixed absorbed, fallback(1)
replace firm=firm+100*(worker>300)
schur_collinear_compare x z duplicate mixed absorbed, fallback(1)
program drop schur_collinear_compare
display "SCHUR_COLLINEAR_OK"
