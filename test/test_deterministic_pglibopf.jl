module TestDeterministicPGLibOPF

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using JuMP
using PowerModels
using Random
using Serialization
using Test

params = TEP.Parameters()

start_file = 2
end_file = 2
dir = "input/"
# num_tests = 10
# files = TEP.select_files(path, num_tests)
log_dir = "test/det"
log_file = "$log_dir/log.md"

try
    TEP.rm_dir(log_dir)
catch e
    @warn e
end

# rng = Random.MersenneTwister(123)

TEP.log_header(log_file)

files = TEP.select_files(dir, end_file)
# Sort files so that the smallest instances are solved first
sort!(files, by=x->parse(Int, match(r"\d+", x).match))
skip = []

@testset "[Deterministic] PG Lib OPF" begin
    # TEP.log(params, "Test $file")
    for (i, file) in enumerate(files[start_file:end_file])
        if file in skip
            TEP.log(params, "Skipping instance $file")
            continue
        end
        TEP.log(params, "Test $file num $(start_file + i - 1)", true)

        filepath = "$dir/$file"
    
        params.log_dir = log_dir
        params.log_file = "$log_dir/$file"

        inst = open(filepath, "r") do io
            deserialize(io)
        end

        try
            results = TEP.init_results()

            mip, _ = TEP.build_deterministic(inst, params)
            JuMP.set_attribute(mip.jump_model, 
                               MOI.RawOptimizerAttribute("TimeLimit"), 
                               params.solver.time_limit)
            start = time()
            JuMP.optimize!(mip.jump_model)
            results["ph_time"] = time() - start
            results["ub"] = JuMP.result_count(mip.jump_model) > 0 ?
                        JuMP.objective_value(mip.jump_model) : const_infinite

            TEP.log_instance(log_file, file, inst, results)
        catch e
            @warn e
            TEP.log_instance(log_file, "<s>" * file * "</s>", inst, Dict())
        end
    end
end

end # module