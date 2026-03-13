module TestBeamSearchCATS

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using JuMP
using PowerModels
using Random
using Test

params = TEP.Parameters()
params.instance.cost_mult = 1e3

log_file = ARGS[1]
log_dir = ARGS[2]
dir = ARGS[3]
file = ARGS[4]

filepath = "$dir/$file"
params.log_file = "$log_dir/$file"

inst = TEP.build_instance(params, filepath)

results = TEP.init_results()

try
    results = TEP.run_parallel_bs!(inst, params)
    TEP.log_instance(log_file, file, inst, results)
catch e
    @warn e
    TEP.log_instance(log_file, "<s>" * file * "</s>", inst, Dict())
end

end # module