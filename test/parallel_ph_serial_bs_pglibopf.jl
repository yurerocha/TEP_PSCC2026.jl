module TestProgressiveHedgingBeamSearchPGLibOPF

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Serialization

params = TEP.Parameters()

log_file = ARGS[1]
log_dir = ARGS[2]
dir = ARGS[3]
file = ARGS[4]
index = parse(Int64, ARGS[5])

filepath = "$dir/$file"

params.log_dir = log_dir
params.log_file = "$log_dir/$file"

# if index <= 4
#     @warn "update time limits for binary and beam search"
#     params.binary_search.time_limit = 300.0
#     params.beam_search.time_limit = 300.0
# end

inst = open(filepath, "r") do io
    deserialize(io)
end

try
    results = TEP.init_results()

    el, cost, is_feas, lb_cost, ub_cost, inserted = 
                                    TEP.run_parallel_ph_serial_bs!(inst, params)
    results["add_rat"] = length(inserted) / inst.num_K
    results["time"] = el
    results["best"] = cost
    results["lb"] = lb_cost
    results["ub"] = ub_cost
    results["is_feas"] = is_feas
    
    TEP.log_instance(log_file, file, inst, results)
catch e
    @warn e
    TEP.log_instance(log_file, "<s>" * file * "</s>", inst, Dict())
end

end # module