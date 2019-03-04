noi cscript "ppmlhdfe with offset()" adofile ppmlhdfe
set linesize 150

* Dataset
	sysuse auto
	gen int idx = _n
	gen byte c = 1
	bys turn: gen t = _n
	xtset turn t
	sort turn t idx
	*clonevar y = price


* Comparison list
	loc list_poisson ///
		scalar: N converged ll /// chi2 p ???
		matrix: trim_b trim_V ///
		macros: wexp wtype

	loc list_xtpoisson ///
		scalar: N converged df_m rank /// ll
		matrix: trim_b trim_V ///
		macros: wexp wtype


* [TEST] No absorb
	loc lhs price
	loc rhs weight gear length

	* 1. Run benchmark
	poisson `lhs' `rhs', vce(robust) offset(mpg) tol(1e-12) nrtol(1e-12)
	notrim
 	storedresults save benchmark e()

	* 2. Run ppmlhdfe
	ppmlhdfe `lhs' `rhs', offset(mpg) tol(1e-12) v(-1) noa maxiter(100)
	notrim
	storedresults compare benchmark e(), tol(1e-9) include(`list_poisson')
	assert inrange(e(ic), 1, 1000)
	storedresults drop benchmark


* [TEST] No absorb and frequency weights
	loc lhs price
	loc rhs weight gear length

	gen double W = foreign + 1
	loc fw "[fw=W]"

	* 1. Run benchmark
	poisson `lhs' `rhs' `fw', vce(robust) offset(mpg) tol(1e-14) nrtol(1e-12)
	notrim
 	storedresults save benchmark e()

	* 2. Run ppmlhdfe
	ppmlhdfe `lhs' `rhs' `fw', offset(mpg) tol(1e-12)  maxiter(100) // v(-1)
	notrim
	storedresults compare benchmark e(), tol(1e-8) include(`list_poisson') // TODO: change tol to 1e-9
	assert inrange(e(ic), 1, 1000)
	storedresults drop benchmark


* [TEST] Absorb turn
	loc lhs price
	loc rhs weight gear length
	local K : list sizeof rhs
	loc absvars turn

	* 1. Run benchmark
	poisson `lhs' `rhs' i.`absvars', offset(mpg) vce(robust) tol(1e-10)
	trim_cons `K'
 	storedresults save benchmark e()

	* 2. Run ppmlhdfe
	ppmlhdfe `lhs' `rhs', absorb(`absvars') offset(mpg) tol(1e-10) keepsingletons v(1) // v(-1)
	trim_cons
	storedresults compare benchmark e(), tol(1e-8) include(`list_poisson')
	assert inrange(e(ic), 1, 1000)

exit
