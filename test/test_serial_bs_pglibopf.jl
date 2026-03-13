module TestBeamSearchPGLibOPF

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using MPI
using Random
using PowerModels
# using Profile, PProf
using Serialization

# start_file = 38
start_file = 38
end_file = 42 # 61
dir = "external/pglib-opf"
scen = 1

# try
#     TEP.rm_dir(log_dir)
# catch e
#     @warn e
# end

# rng = Random.MersenneTwister(123)

files = TEP.select_files(dir, end_file)
# Sort files so that the smallest instances are solved first
sort!(files, by=x->parse(Int, match(r"\d+", x).match))
# Run solver with binary decision variables
skip = ["pglib_opf_case8387_pegase.m"]
# "pglib_opf_case30_ieee.m" Infeasible

project_dir = dirname(Base.active_project())
parallel_script = joinpath(@__DIR__, "parallel_bs_pglibopf.jl")

# test_params = [[it, it_wo_impr] for it in [5, 10, 15], it_wo_impr in [1]]
# test_params = [[cands_per_batch, rcl_rat] for cands_per_batch in [5e-4, 1e-3, 5e-3], rcl_rat in [1.0]]
# test_params = [[cands_per_batch, rcl_rat] for cands_per_batch in [0.25e-2], rcl_rat in [1.0]]

for tp in [20, 5]
    log_dir = "test/tuning/bs/test$tp"
    log_file = "$log_dir/log$tp.md"
    try
        TEP.rm_dir(log_dir)
    catch e
        @warn e
    end
    TEP.log_header(log_file)

    # params.binary_search.max_it = tp[1]
    # params.binary_search.num_max_it_wo_impr = tp[2]

    # params.beam_search.candidates_per_batch_mult = tp[1]
    # params.beam_search.restricted_list_ratio = tp[2]

    for (i, file) in enumerate(files[start_file:end_file])
        params = TEP.Parameters()
        params.beam_search.num_max_it_wo_impr = tp
        if file in skip
            TEP.log(params, "Skipping instance $file")
            continue
        end
        TEP.log(params, "Test $file num $(start_file + i - 1)", true)

        filepath = "$dir/$file"
        
        params.log_dir = log_dir
        params.log_file = "$log_dir/$file"
        
        inst = TEP.build_instance(params, filepath)
        # try
            lp = TEP.build_lp(inst, params, scen)
            cache = TEP.WorkerCache(TEP.Cache(inst, params))
            inserted = Set{TEP.CandType}(keys(inst.K))
            removed = Set{TEP.CandType}()
            start_time = time()
            TEP.update_lp!(inst, params, lp, inserted)
            init_cost, _ = TEP.comp_penalized_cost(inst, params, scen, 
                                                    lp, cache, inserted)
        
            # Profile.clear()
            # @profile inserted = TEP.run_serial_bs!(inst, params, scen, lp, cache, 
            #                                   inserted, removed, cost, start_time)
            # pprof()
            el, inserted = TEP.run_serial_bs!(inst, params, scen, lp, cache, 
                                    inserted, removed, init_cost, start_time)
            cost, _ = TEP.comp_penalized_cost(inst, params, scen, 
                                                lp, cache, inserted)

            results = TEP.init_results()
            
            st = TEP.Status("rr", inst.num_K - length(inserted), inst.num_K, 
                            cost, init_cost, start_time)

            results["rnr_time"] = el
            results["rnr_rm_rat"] = st.rm_ratio
            results["rnr_impr_rat"] = st.impr_ratio
            results["ub"] = cost


            # data = Profile.fetch()
            # open("profile_data.jls", "w") do io
            #     serialize(io, data)
            # end

            TEP.log_instance(log_file, file, inst, results)
        # catch e
        #     @warn e
        #     TEP.log_instance(log_file, "<s>" * file * "</s>", inst, Dict())
        # end
    end
end

end # module