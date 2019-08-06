set more off
/*
adopath ++ "..\github\ppml_hdfe_demo-master\src"
adopath ++ "..\github\ftools-master\src"
adopath ++ "..\github\reghdfe-master\src"
adopath ++ "..\datasets\"
*/


clear

which reghdfe
which ftools
which ppmlhdfe


** Example 1
sjlog using example1_a, replace
input y x1 x2 x3
0	1	2	1
0	0	0	2
0	2	3	3
1	1	2	4
2	2	4	5
3	1	2	6
end
sjlog close, replace

sjlog using example1_b, replace
poisson y x1 x2 x3, nolog
sjlog close, replace

sjlog using example1_c, replace
ppml y x1 x2 x3
sjlog close, replace

sjlog using example1_d, replace
ppmlhdfe y x1 x2 x3, nolog
sjlog close, replace

** Example 1
*webuse ships, clear
*use ships, clear
*sjlog using input1a, replace
webuse ships, clear
*xtpoisson acc op_75_79 co_65_69 co_70_74 co_75_79, exp(service) irr fe nolog
*sjlog close, replace

** Example 2
sjlog using input1b, replace
ppmlhdfe acc op_75_79 co_65_69 co_70_74 co_75_79, a(ship) exp(service) irr nolog
sjlog close, replace

** Example 3
sjlog using input1c, replace
ppmlhdfe acc op_75_79 co_65_69, a(ship co_70_74 co_75_79) exp(service) irr nolog
sjlog close, replace

* Example 3

* same as
* ppml_panel_sg trade fta if category == "TOTAL", ex(isoexp) im(isoimp) y(year)

*use example_trade_fta_data if category == "TOTAL", clear
sjlog using input1d, replace
use http://fmwww.bc.edu/RePEc/bocode/e/EXAMPLE_TRADE_FTA_DATA if category=="TOTAL", clear
cap egen imp=group(isoimp)
cap egen exp=group(isoexp)
ppmlhdfe trade fta, a(imp#year exp#year imp#exp) cluster(imp#exp) nolog
sjlog close, replace

* Example 4

/*
sjlog using input1, replace
ppmlhdfe trade fta if category == "TOTAL", absorb(i.expid#i.year i.impid#i.year i.expid#i.impid i.year) vce(cluster i.impid#i.expid) 
sjlog close, replace
*/

use citations_example, clear
estimates clear
ppmlhdfe cit nbaut, a(issn type jel2 pubyear)
eststo
ppmlhdfe cit nbaut, a(issn#c.year type jel2 pubyear)
eststo
ppmlhdfe cit nbaut, a(issn#year type jel2 pubyear)
eststo
estfe *, labels(type "Article type FE" jel2 "JEL code FE" pubyear "Publication year FE" issn "ISSN FE" issn#c.year "Year trend by ISSN" issn#year "ISSN-Year FE")
sjlog using input1e, replace
esttab, indicate(`r(indicate_fe)', labels("Yes" "")) b(3) se(3) varwidth(25) label ///
	stat(N ll, fmt(%12.0fc %13.1fc)) se starlevels(* 0.1 ** 0.05 *** 0.01) compress
sjlog close, replace
exit
