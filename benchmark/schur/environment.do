* Untimed environment/processor preflight, without license or user information.
args output processors source
version 13
clear all
set more off
adopath ++ "`source'/src"
capture set processors `processors'
if _rc | c(processors)!=`processors' {
    display as error "Requested processors `processors' unavailable; actual " c(processors)
    exit 198
}
tempname metadata dependencies
file open `metadata' using "`output'/environment.tsv", write replace
file write `metadata' "key" _tab "value" _n
file write `metadata' "stata_version" _tab (c(stata_version)) _n
file write `metadata' "processors" _tab (c(processors)) _n
file write `metadata' "os" _tab "`c(os)'" _n
file write `metadata' "machine_type" _tab "`c(machine_type)'" _n
file close `metadata'
file open `dependencies' using "`output'/dependencies.tsv", write replace
file write `dependencies' "name" _tab "path" _n
local required "ftools.ado ftools.mata ftools_type_aliases.mata ftools_common.mata ftools_main.mata ftools_hash1.mata ftools_plugin.mata fcollapse_functions.mata reghdfe.ado reghdfe5.mata reghdfe5_parse.ado"
foreach name of local required {
    findfile `name'
    file write `dependencies' "`name'" _tab "`r(fn)'" _n
}
local optional "lftools.mlib lmoremata.mlib lreghdfe.mlib ms_get_version.ado ms_compile_mata.ado ms_parse_absvars.ado ms_fvstrip.ado ms_expand_varlist.ado"
foreach name of local optional {
    capture findfile `name'
    if !_rc file write `dependencies' "`name'" _tab "`r(fn)'" _n
}
file close `dependencies'
display "REDBENCH_ENVIRONMENT_OK"
