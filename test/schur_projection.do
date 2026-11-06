* From the checkout root: do test/schur_projection.do .
* Independent weighted-dummy oracle and numerical failure-safety tests.
args root
if `"`root'"' == "" local root ..
clear all
set more off
set linesize 200
capture set processors 1
adopath ++ "`root'/src"
set seed 710536
set obs 240
gen int f1=ceil(_n/6)
gen int f2=1+mod(_n+f1,13)
gen int f3=1+mod(floor(_n/3)+2*f1,11)
gen int f4=1+mod(_n+floor(_n/7),9)
gen int f5=1+mod(_n+floor(_n/11),7)
gen double x=rnormal()
gen double y=exp(.2*x+.1*rnormal())
quietly ppmlhdfe y x, absorb(f1 f2)
mata: mata clear
findfile ppmlhdfe.mata
include "`r(fn)'"

mata:
mata set matastrict on
real matrix schur_adversarial_design(class FixedEffects scalar H)
{
    real matrix D
    real scalar width, g, offset, i
    width=0
    for (g=1;g<=H.G;g++) width=width+H.factors[g].num_levels
    D=J(rows(H.sample),width,0)
    offset=0
    for (g=1;g<=H.G;g++) {
        for (i=1;i<=rows(H.sample);i++) D[i,offset+H.factors[g].levels[i]]=1
        offset=offset+H.factors[g].num_levels
    }
    return(D)
}

void schur_adversarial_case(string scalar label, string scalar effects, real scalar allow_fallback)
{
    class FixedEffects scalar H, R
    class PpmlSchur scalar s
    real matrix D, Y, oracle, candidate, reference, tighter, qroracle
    real matrix expected, polished, before
    real colvector seq, basew, w, wn
    real scalar N, scenario, ok, err, moment, tighter_ok, denom
    real scalar upstream_error, polished_error
    H=fixed_effects(effects,"","","",0,0)
    s.init(H)
    assert(s.enabled)
    D=schur_adversarial_design(H)
    N=rows(D)
    seq=(1::N)
    Y=(sin(seq),cos(seq),sin(seq:/7),J(N,1,0),D*(sin((1::cols(D)))))
    // One zero RHS, one fully absorbed RHS, and three independent RHS.
    basew=exp(3:*sin(seq:/11))
    reference=.
    for (scenario=1;scenario<=5;scenario++) {
        if (scenario==1) w=J(N,1,1)
        if (scenario==2) w=basew
        if (scenario==3) w=basew:*1e-150
        if (scenario==4) w=basew:*1e150
        if (scenario==5) w=basew:*1e200
        wn=w:/max(w)
        oracle=Y-D*pinv(sqrt(wn):*D)*(sqrt(wn):*Y)
        candidate=Y
        ok=s.partial_out(candidate,w,1e-12,5000)
        if (!ok) {
            assert(candidate==Y)
            assert(allow_fallback)
            printf("KERNEL_FALLBACK %s weights=%g\n",label,scenario)
            continue
        }
        err=max(abs(candidate-oracle))
        moment=max(abs(quadcross(D,wn,candidate)))
        assert(err<2e-8)
        assert(moment<2e-9)
        if (scenario==2) reference=candidate
        if (scenario>2) assert(max(abs(candidate-reference))<2e-8)
        printf("ORACLE_PASS %s weights=%g iter=%g err=%g moment=%g\n",label,scenario,s.iterations,err,moment)
    }

    // A scalar weight has the same projection as unit weights.
    oracle=Y-D*pinv(D)*Y
    candidate=Y
    assert(s.partial_out(candidate,1e200,1e-12,5000))
    assert(max(abs(candidate-oracle))<2e-8)

    // Strongly unequal weights exercise a different preconditioner scale.
    w=exp(18:*sin(seq:/11))
    wn=w:/max(w)
    oracle=Y-D*pinv(sqrt(wn):*D)*(sqrt(wn):*Y)
    candidate=Y
    ok=s.partial_out(candidate,w,1e-12,5000)
    if (ok) {
        err=sqrt(sum(wn:*(candidate-oracle):^2))/sqrt(sum(wn:*Y:^2))
        printf("WEIGHT_RANGE_RESULT %s err=%g iter=%g\n",label,err,s.iterations)
        if (err>=2e-8) {
            tighter=Y
            tighter_ok=s.partial_out(tighter,w,1e-14,5000)
            assert(tighter_ok)
            assert(sqrt(sum(wn:*(tighter-oracle):^2))/sqrt(sum(wn:*Y:^2))<2e-8)
            qroracle=Y-D*qrsolve(sqrt(wn):*D,sqrt(wn):*Y)
            printf("WEIGHT_RANGE_DIAG %s tightok=%g tighterr=%g qrerr=%g oracleqr=%g\n",label,tighter_ok,sqrt(sum(wn:*(tighter-oracle):^2))/sqrt(sum(wn:*Y:^2)),sqrt(sum(wn:*(candidate-qroracle):^2))/sqrt(sum(wn:*Y:^2)),sqrt(sum(wn:*(oracle-qroracle):^2))/sqrt(sum(wn:*Y:^2)))
        }
        // Moment tolerances alone do not bound solution error under poor
        // conditioning. Assess the production wrapper against upstream too.
        R=fixed_effects(effects,"","","",0,0)
        H.load_weights("aweight","redteam_weight",w,0)
        R.load_weights("aweight","redteam_weight",w,0)
        H.tolerance=R.tolerance=1e-12
        H.maxiter=R.maxiter=20000
        H.verbose=R.verbose=-1
        expected=Y[.,1..3]
        R._partial_out(expected,0,0,0,1)
        polished=Y[.,1..3]
        s.project(H,polished,w)
        denom=sqrt(sum(wn:*Y[.,1..3]:^2))
        upstream_error=sqrt(sum(wn:*(expected-oracle[.,1..3]):^2))/denom
        polished_error=sqrt(sum(wn:*(polished-oracle[.,1..3]):^2))/denom
        printf("WEIGHT_RANGE_POLISH %s upstream=%g polished=%g fallbacks=%g\n",label,upstream_error,polished_error,s.fallbacks)
        assert(H.kept==R.kept)
        assert(mreldif(H.means,R.means)<1e-12)
        // At this extreme weight ratio the original absorption can be
        // inaccurate too. Keep its error visible and require no regression
        // beyond the same 2e-8 absolute weighted-oracle allowance.
        assert(polished_error<=upstream_error+2e-8)
        if (s.fallbacks) assert(mreldif(polished,expected)<2e-8)

        // Moment tolerance does not bound solution error under poor
        // conditioning. Tightening both solvers must achieve the unchanged
        // 2e-8 independent-oracle gate on an accepted accelerated projection.
        H.tolerance=R.tolerance=1e-14
        expected=Y[.,1..3]
        R._partial_out(expected,0,0,0,1)
        polished=Y[.,1..3]
        s.init(H)
        s.project(H,polished,w)
        upstream_error=sqrt(sum(wn:*(expected-oracle[.,1..3]):^2))/denom
        polished_error=sqrt(sum(wn:*(polished-oracle[.,1..3]):^2))/denom
        assert(s.fallbacks==0)
        assert(H.kept==R.kept)
        assert(mreldif(H.means,R.means)<1e-12)
        assert(polished_error<2e-8)
        printf("WEIGHT_RANGE_TIGHT %s upstream=%g polished=%g\n",label,upstream_error,polished_error)
        printf("WEIGHT_RANGE_PASS %s kernel_err=%g\n",label,err)
    }
    else {
        assert(candidate==Y)
        printf("WEIGHT_RANGE_SAFE_FALLBACK %s\n",label)
    }

    // All rejected inputs and exhausted solves must be transactional.
    w=basew
    candidate=Y
    ok=s.partial_out(candidate,w,1e-15,0)
    if (label=="regular_G2") assert(!ok)
    if (!ok) assert(candidate==Y)
    candidate=Y
    w[1]=0
    assert(!s.partial_out(candidate,w,1e-12,1000))
    assert(candidate==Y)
    w[1]=-1
    assert(!s.partial_out(candidate,w,1e-12,1000))
    assert(candidate==Y)
    w[1]=.
    assert(!s.partial_out(candidate,w,1e-12,1000))
    assert(candidate==Y)
    w=basew
    assert(!s.partial_out(candidate,w[1::(N-1)],1e-12,1000))
    assert(candidate==Y)
    assert(!s.partial_out(candidate,w,0,1000))
    assert(candidate==Y)
    assert(!s.partial_out(candidate,w,.,1000))
    assert(candidate==Y)
    assert(!s.partial_out(candidate,w,1e-12,-1))
    assert(candidate==Y)
    assert(!s.partial_out(candidate,w,1e-12,.))
    assert(candidate==Y)
    candidate=Y
    candidate[1,1]=.
    before=candidate
    assert(!s.partial_out(candidate,w,1e-12,1000))
    assert(candidate==before)
    printf("IMMUTABILITY_PASS %s\n",label)
}

schur_adversarial_case("regular_G2","f1 f2",0)
schur_adversarial_case("regular_G3","f1 f2 f3",0)
schur_adversarial_case("regular_G4","f1 f2 f3 f4",0)
schur_adversarial_case("regular_G5","f1 f2 f3 f4 f5",0)
schur_adversarial_case("reordered_G5","f3 f5 f1 f4 f2",0)
end

preserve
replace f4=f2
replace f5=f3
mata: schur_adversarial_case("duplicate_retained_G5","f1 f2 f3 f4 f5",0)
restore

preserve
replace f3=1+mod(f1,3)
replace f4=ceil(f1/5)
replace f5=1+mod(f1,7)
mata: schur_adversarial_case("nested_retained_G5","f1 f2 f3 f4 f5",0)
restore

preserve
replace f2=f1
replace f3=f1
replace f4=f1
replace f5=f1
mata: schur_adversarial_case("all_duplicate_G2","f1 f2",0)
mata: schur_adversarial_case("all_duplicate_G3","f1 f2 f3",0)
mata: schur_adversarial_case("all_duplicate_G4","f1 f2 f3 f4",0)
mata: schur_adversarial_case("all_duplicate_G5","f1 f2 f3 f4 f5",0)
restore

preserve
forvalues g=1/5 {
    replace f`g'=f`g'+100*(_n>120)
}
mata: schur_adversarial_case("disconnected_G2","f1 f2",0)
mata: schur_adversarial_case("disconnected_G3","f1 f2 f3",0)
mata: schur_adversarial_case("disconnected_G4","f1 f2 f3 f4",0)
mata: schur_adversarial_case("disconnected_G5","f1 f2 f3 f4 f5",0)
gen double shuffle=runiform()
sort shuffle
mata: schur_adversarial_case("disconnected_permuted_G5","f5 f2 f1 f4 f3",0)
restore



mata:
void schur_adversarial_polish(string scalar effects)
{
    class FixedEffects scalar H, R
    class PpmlSchur scalar s
    real matrix D, Y, expected, candidate
    real colvector seq, w, f, r
    real scalar N, cutoff, C, F, Q, epsilon, scenario
    H=fixed_effects(effects,"","","",0,0)
    R=fixed_effects(effects,"","","",0,0)
    D=schur_adversarial_design(H)
    N=rows(D)
    seq=(1::N)
    w=exp(3:*sin(seq:/11))
    H.load_weights("aweight","redteam_weight",w,0)
    R.load_weights("aweight","redteam_weight",w,0)
    H.tolerance=R.tolerance=1e-12
    H.maxiter=R.maxiter=5000
    H.verbose=R.verbose=-1
    f=D*sin((1::cols(D)))
    r=sin(seq)-D*pinv(sqrt(w):*D)*(sqrt(w):*sin(seq))
    cutoff=min((1e-6,H.tolerance/10))
    C=quadcross(f,r)
    F=quadcross(f,f)
    Q=quadcross(r,r)
    epsilon=(cutoff*C+sqrt(cutoff^2*C^2+(1-cutoff)*Q*cutoff*F))/((1-cutoff)*Q)
    for (scenario=1;scenario<=6;scenario++) {
        if (scenario==1) Y=(sin(seq),cos(seq),f,J(N,1,0))
        if (scenario==2) Y=(sin(seq),f+epsilon*.999*r,f+epsilon*1.001*r)
        if (scenario==3) Y=(f+epsilon*.01*r,f+epsilon*100*r,sin(seq))
        if (scenario==4) Y=(sin(seq)*1e-8,cos(seq)*1e8,f+1e-8*r)
        if (scenario>=5) {
            Y=(sin(seq),cos(seq),sin(seq:/7))
            w=scenario==5 ? J(N,1,1) : exp(.2:*sin(seq:/11))
            H.load_weights("aweight","redteam_weight",w,0)
            R.load_weights("aweight","redteam_weight",w,0)
        }
        expected=Y
        R._partial_out(expected,0,0,0,1)
        candidate=Y
        s.init(H)
        s.project(H,candidate,w)
        assert(H.kept==R.kept)
        assert(mreldif(H.means,R.means)<1e-12)
        assert(mreldif(candidate,expected)<2e-8)
        if (scenario>=5) assert(s.fallbacks==0)
        else assert(s.fallbacks==1)
        printf("POLISH_PASS effects=%s scenario=%g calls=%g fallback=%g\n",effects,scenario,s.calls,s.fallbacks)
    }
}
schur_adversarial_polish("f1 f2")
schur_adversarial_polish("f1 f2 f3")
schur_adversarial_polish("f1 f2 f3 f4 f5")
end
display "SCHUR_PROJECTION_OK"
