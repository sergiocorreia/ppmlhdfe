* Focused regression check: original x1-x4/y/topology unchanged by extra draws.
args original_fixture extended_fixture
version 13
clear all
set more off
use "`original_fixture'", clear
ds
local original_variables "`r(varlist)'"
use "`extended_fixture'", clear
cf `original_variables' using "`original_fixture'", all
confirm variable x64
assert !missing(x5,x64)
display "REDBENCH_FIXTURE_PREFIX_OK"
