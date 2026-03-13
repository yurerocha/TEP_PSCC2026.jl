module TestBeamSearchPGLibOPF

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Test
using JuMP
using PowerModels
using Random

params = TEP.Parameters()

log_file = ARGS[1]
start_file = parse(Int64, ARGS[2])
end_file = parse(Int64, ARGS[3])
log_dir = ARGS[4]
dir = ARGS[5]
file = ARGS[6]
scen = 1

filepath = "$dir/$file"

params.log_dir = log_dir
params.log_file = "$log_dir/$file"

try
    results = TEP.init_results()

    inst = TEP.build_instance(params, filepath)

    lp = TEP.build_lp(inst, params, scen)
    c = TEP.Cache(0, 0)

    results = TEP.run_parallel_bs!(inst, params, scen, lp, false, c)

    TEP.log_instance(log_file, file, inst, results)
catch e
    @warn e
    TEP.log_instance(log_file, "<s>" * file * "</s>", inst, Dict())
end

end # module