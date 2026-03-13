module TestDCPPMPGLibOPF

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Test
using JuMP
using Gurobi, Ipopt
using PowerModels
using Markdown, DataFrames

path = "../external/pglib-opf/"
num_tests = 40
files = TEP.select_files(path, num_tests)

params = TEP.Parameters()
TEP.config_dcp_pm_tests!(params)

eps = 1e-3

# BASELINE.md solution costs do not have precision

@testset "[DCP PM] PG Lib OPF" begin
    for (i, file) in enumerate(files)
        TEP.log(params, "Test $i $file")
        
        # params.log_file = "log/" * TEP.get_inst_name(file) * ".txt"

        filepath = path * file

        mp_data = PowerModels.parse_file(filepath)
        # pm = instantiate_model(mp_data, DCPPowerModel, PowerModels.build_opf)
        # TEP.print_constrs(pm.model, "TEP.jl/model2.lp")

        # Run DC-OPF
        dc_opf = PowerModels.solve_opf(mp_data, 
                                       DCPPowerModel, 
                                       # DCMPPowerModel, 
                                       params.model.optimizer)
        
        # Run TEP
        inst = TEP.build_instance(params, filepath)
        mip = TEP.build_mip(inst, params)
        # force_solution(inst, mip, dc_opf["solution"], mp_data)
        # TEP.print_constrs(mip.jump_model, "TEP.jl/model1.lp")
        tep = TEP.solve!(inst, params, mip)

        TEP.log(params, "Test $i $file")
        TEP.log(params, "$(tep["objective"]), $(dc_opf["objective"])", true)
        @test abs(tep["objective"] - dc_opf["objective"]) < eps
        # readline()
    end
end

end # module
