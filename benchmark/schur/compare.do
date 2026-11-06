* Untimed numerical comparison in a third fresh process.
args reference target output require_exact
version 13
if "`require_exact'"=="" local require_exact 0
clear all
set more off
display "REDBENCH_START"
estimates use "`reference'.ster"
matrix b0=e(b)
matrix V0=e(V)
local names0 : colfullnames e(b)
_ms_omit_info e(b)
matrix omitted0=r(omit)
scalar dev0=e(deviance)
scalar sep0=e(num_separated)
scalar singleton0=e(num_singletons)
scalar N0=e(N)
estimates use "`target'.ster"
matrix b1=e(b)
matrix V1=e(V)
local names1 : colfullnames e(b)
_ms_omit_info e(b)
matrix omitted1=r(omit)
scalar stripe_equal = (`"`names0'"'==`"`names1'"')
scalar omitted_equal=0
scalar b_scaled=.
scalar b_se=.
scalar V_scaled=.
scalar V_cor_scaled=.
scalar V_diag_relative=.
scalar se_relative=.
if stripe_equal {
    scalar omitted_equal=mreldif(omitted0,omitted1)==0
    scalar b_scaled=mreldif(b0,b1)
    scalar V_scaled=mreldif(V0,V1)
    mata:
        v0=st_matrix("V0"); v1=st_matrix("V1")
        diag0=diagonal(v0); diag1=diagonal(v1)
        keep=selectindex(diag0:>0)
        if (length(keep)) {
            sd=sqrt(diag0[keep])
            bdiff=st_matrix("b1")'-st_matrix("b0")'
            st_numscalar("b_se",max(abs(bdiff[keep]):/sd))
            st_numscalar("V_cor_scaled",max(abs(v1[keep,keep]-v0[keep,keep]):/(sd*sd')))
            st_numscalar("V_diag_relative",max(abs(diag1[keep]-diag0[keep]):/diag0[keep]))
            st_numscalar("se_relative",max(abs(sqrt(diag1[keep])-sd):/sd))
        }
        else {
            st_numscalar("b_se",max(abs(st_matrix("b1")-st_matrix("b0"))))
            st_numscalar("V_cor_scaled",max(abs(v1-v0)))
            st_numscalar("V_diag_relative",max(abs(v1-v0)))
            st_numscalar("se_relative",max(abs(v1-v0)))
        }
    end
}
scalar dev_scaled=reldif(dev0,e(deviance))
scalar counts_equal=(sep0==e(num_separated) & singleton0==e(num_singletons) & N0==e(N))
use "`reference'.dta", clear
rename (sample mu comparison_weight) (sample0 mu0 weight0)
merge 1:1 rowid using "`target'.dta"
gen byte sample_mismatch=_merge!=3 | sample0!=sample
quietly count if sample_mismatch
scalar sample_equal = (r(N)==0)
quietly count if sample & (missing(mu,mu0,weight0,comparison_weight) | mu<=0 | mu0<=0 | weight0!=comparison_weight)
scalar predictions_valid = (r(N)==0)
gen double mu_error=abs(mu-mu0)/(1+abs(mu0)) if sample & sample0
quietly summarize mu_error, meanonly
scalar mu_scaled=r(max)
gen double eta_error=abs(ln(mu)-ln(mu0)) if sample & sample0
quietly summarize eta_error, meanonly
scalar eta_max=r(max)
mata:
    use=selectindex(st_data(.,"sample"):==1 :& st_data(.,"sample0"):==1)
    m0=st_data(use,"mu0"); m1=st_data(use,"mu"); w=st_data(use,"weight0")
    st_numscalar("mu_weighted_relative",sqrt(quadcross(m1-m0,w,m1-m0)/quadcross(m0,w,m0)))
end
scalar exact_results=stripe_equal & omitted_equal & counts_equal & sample_equal & predictions_valid & b_scaled==0 & V_scaled==0 & dev_scaled==0 & mu_scaled==0
scalar passed=stripe_equal & omitted_equal & counts_equal & sample_equal & predictions_valid & b_scaled<1e-7 & b_se<1e-4 & V_scaled<1e-7 & V_cor_scaled<1e-4 & V_diag_relative<1e-4 & se_relative<1e-4 & dev_scaled<1e-9 & mu_scaled<1e-4
if `require_exact' scalar passed=exact_results
clear
set obs 1
foreach item in stripe_equal omitted_equal counts_equal sample_equal predictions_valid b_scaled b_se V_scaled V_cor_scaled V_diag_relative se_relative dev_scaled mu_scaled eta_max mu_weighted_relative exact_results passed {
    gen double `item'=scalar(`item')
}
export delimited using "`output'", replace
display "REDBENCH_COMPARE_OK passed=" passed[1]
if !passed[1] exit 459
