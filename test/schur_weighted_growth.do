* From the checkout root: do test/schur_weighted_growth.do .
* Weighted projection can increase the unweighted norm legitimately.
args root
if `"`root'"' == "" local root ..
clear all
set more off
capture set processors 1
adopath ++ "`root'/src"
set obs 800
gen int worker=ceil(_n/10)
gen byte region=mod(worker,4)
gen byte x=mod(_n-1,10)==0
gen double w=cond(x,100,1)
gen double y=exp(.3*x+.1*sin(worker)+.02*cos(_n))

* This nested-effect model has an exact within-worker projection. Giving the
* head observation more weight spreads its value across the other nine rows,
* increasing the regressor's unweighted squared norm by more than sevenfold.
local accuracy "tol(1e-11) itol(1e-12) start_inner_tol(1e-12) use_exact_partial(1) use_heuristic_tol(0)"
quietly ppmlhdfe y x [pw=w], absorb(worker region) d(d0) `accuracy'
matrix b0=e(b)
matrix V0=e(V)
scalar dev0=e(deviance)
gen byte sample0=e(sample)
predict double mu0 if sample0, mu
quietly ppmlhdfe y x [pw=w], absorb(worker region) d(d1) `accuracy' use_schur(1)
assert e(schur_col_calls)==1
assert e(schur_col_fails)==0
assert e(schur_calls)>e(schur_fallbacks)
assert sample0==e(sample)
assert mreldif(b0,e(b))<1e-7
assert mreldif(V0,e(V))<1e-7
assert reldif(dev0,e(deviance))<1e-9
predict double mu1 if sample0, mu
assert abs(mu1-mu0)/(1+abs(mu0))<1e-6 if sample0
display "WEIGHTED_GROWTH_ESTIMATION_PASS"

* A deliberately misleading intermediate norm: four identical worker blocks
* with one high-weight firm-1 observation and nine low-weight firm-2 rows.
* A loose moment tolerance permits a within-worker starting transformation.
* The final residual must survive the original norm criterion, although using
* the larger intermediate norm would incorrectly zero it.
clear
set obs 40
gen byte worker=ceil(_n/10)
gen byte pos=mod(_n-1,10)
gen byte firm=cond(pos==0,1,2)
gen double w=cond(pos==0,1,1e-12)
gen byte noise=cond(pos==0,0,cond(pos<=4,-1,cond(pos==5,0,1)))
gen double input=(pos==0)+sqrt(.5e-6)*noise
clonevar worker_copy=worker
mata: mata clear
findfile ppmlhdfe.mata
include "`r(fn)'"

mata:
mata set matastrict on
void schur_growth_elimination(string scalar effects, real scalar expected)
{
    class FixedEffects scalar H
    class PpmlSchur scalar s
    real scalar g, selected_count
    H=fixed_effects(effects,"","","",0,0)
    s.init(H)
    assert(s.enabled)
    assert(s.eliminated==expected)
    selected_count=H.factors[s.eliminated].num_levels
    for (g=1;g<=H.G;g++) {
        assert(selected_count>=H.factors[g].num_levels)
        if (g<s.eliminated) assert(selected_count>H.factors[g].num_levels)
    }
}

void schur_growth_false_zero()
{
    class FixedEffects scalar H, R
    class PpmlSchur scalar s
    real colvector w, original, candidate, expected, internal, actual
    real scalar norm0, norm1, cutoff
    H=fixed_effects("worker firm","","","",0,0)
    R=fixed_effects("worker firm","","","",0,0)
    w=st_data(.,"w")
    original=st_data(.,"input")
    H.load_weights("aweight","w",w,0)
    R.load_weights("aweight","w",w,0)
    H.tolerance=R.tolerance=1e-4
    H.maxiter=R.maxiter=1000
    H.verbose=R.verbose=-1
    cutoff=min((1e-6,H.tolerance/10))
    norm0=quadcross(original,original)

    expected=original
    R._partial_out(expected,0,0,0,1)
    assert(R.kept==1)
    assert(quadcross(expected,expected)/norm0>cutoff)

    s.init(H)
    candidate=original
    assert(s.partial_out(candidate,w,H.tolerance,H.maxiter))
    assert(s.iterations==0)
    norm1=quadcross(candidate,candidate)
    assert(norm1>8*norm0)
    assert(quadcross(expected,expected)/norm1<cutoff)

    internal=candidate
    H._partial_out(internal,0,0,0,1)
    assert(H.kept==0)
    assert(all(internal:==0))

    s.init(H)
    actual=original
    s.project(H,actual,w)
    assert(s.calls==1)
    assert(s.fallbacks==1)
    assert(H.kept==R.kept)
    assert(mreldif(H.means,R.means)<1e-14)
    assert(mreldif(actual,expected)<1e-14)
    printf("WEIGHTED_GROWTH_FALSE_ZERO_PASS original_ratio=%g candidate_ratio=%g growth=%g\n",quadcross(expected,expected)/norm0,quadcross(expected,expected)/norm1,norm1/norm0)
}
schur_growth_elimination("worker firm",1)
schur_growth_elimination("firm worker",2)
schur_growth_elimination("worker_copy worker firm",1)
schur_growth_false_zero()
end
display "SCHUR_WEIGHTED_GROWTH_OK"
