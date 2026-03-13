module TestPTDFModel

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Test
using JuMP
using Gurobi
using PowerModels

path = "../external/pglib-opf/"
num_tests = 26
files = TEP.select_files(path, num_tests)

params = TEP.Parameters()
params.log_level = 0

eps = 1e-1

@testset "PTDF Model" begin
    for (i, file) in enumerate(files)
        TEP.log(params, "Test $i $file")

        # params.log_file = "log/" * TEP.get_inst_name(file) * ".txt"

        inst = TEP.build_instance(params, path * file)
        
        ptdf = TEP.build_ptdf(inst, params, T = Float64)
        JuMP.optimize!(ptdf.jump_model)
        
        lp = TEP.build_lp(inst, params)
        JuMP.optimize!(lp.jump_model)

        ptdf_obj_val = objective_value(ptdf.jump_model)
        lp_obj_val = objective_value(lp.jump_model)

        TEP.log(params, "Test $i $file")
        TEP.log(params, "$ptdf_obj_val, $lp_obj_val", true)

        @test abs(ptdf_obj_val - lp_obj_val) < eps
    end
end

end # module