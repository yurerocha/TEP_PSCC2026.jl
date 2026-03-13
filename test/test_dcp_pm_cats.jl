module TestDCPPMCATS

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Test
using JuMP
using Gurobi, Ipopt
using PowerModels
using Markdown, DataFrames

file = "../external/CATS-CaliforniaTestSystem/MATPOWER/CaliforniaTestSystem.m"
num_scenarios = 10

params = TEP.Parameters()
TEP.config_dcp_pm_tests!(params)

eps = 0.5

# Some objective function values (OVF) differ, with the TEP model OFV larger, 
# but if the solution of the DCP PM, in terms of generation, is enforced as 
# constraints within the TEP model, the OVFs equal. This is probably due to the 
# quadractic terms in the objective function

@testset "[DCP PM] California Test System" begin
    # Read scenarios and solve coresponding DCP PMs
    TEP.log(params, "Build instances and solve corresponding DCP PMs", true)
    inst, sols = TEP.build_cats_instance(params, num_scenarios)

    for scen in 1:num_scenarios
        TEP.log(params, "Test $scen $file")
        params.log_file = "log/" * TEP.get_inst_name(file) * "_$scen.txt"
        
        # Run TEP
        mip = TEP.build_mip(inst, params, scen)
        # TEP.enforce_sol(inst, mip, sols[scen]["solution"])
        results = TEP.solve!(params, mip)
        
        @info results[7], sols[scen]["objective"]
        @test abs(results[7] - sols[scen]["objective"]) < eps
    end
end

end # module