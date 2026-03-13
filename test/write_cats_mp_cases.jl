# Read CATS test case and write mutiple single-scenario matpower cases.

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Test
using JuMP
using Gurobi, Ipopt
using PowerModels
using Markdown, DataFrames

# TODO: Sort the data properly
# TODO: Tests
filename = "CaliforniaTestSystem"
num_scenarios = 10
output_dir = "input/"

if isdir(output_dir)
    if !isempty(output_dir)
        rm(output_dir, recursive = true)
    end
end
mkdir(output_dir)

params = TEP.Parameters()

# Read scenarios and solve coresponding DCP PMs
TEP.log(params, "Build cases", true)
mp_cases = TEP.build_cats_cases(params, num_scenarios)

for (i, mpc) in enumerate(mp_cases)
    fname = output_dir * filename * "Scen#$i.m"
    # @warn mpc["shunt"]
    # ["bus", "source_type", "name", "dcline", "source_version", "gen", 
    #  "branch", "storage", "switch", "baseMVA", "per_unit", "shunt", "load"]
    PowerModels.export_matpower(fname, mpc)
end