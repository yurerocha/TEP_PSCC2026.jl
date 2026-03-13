module TestBeamSearchPGLibOPF

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using JuMP
using MPI
# using Random
using PowerModels
# using Profile, PProf
using Serialization

start_file = 1
end_file = 6 # 61
scen = 1
log_dir = "test/det"
is_heur_en = true

log_dir *= is_heur_en ? "_ia" : "_mip"
log_file = "$log_dir/log.md"

params = TEP.Parameters()
params.beam_search.candidates_per_batch_mult = 0.1e-3

try
    TEP.rm_dir(log_dir)
catch e
    @warn e
end

files = [
         "pglib_opf_case3012wp_k.m",
         "pglib_opf_case6495_rte.m",
         "pglib_opf_case7336_epigrids.m",
         "CaliforniaTestSystem.m",
         "pglib_opf_case9591_goc.m",
         "pglib_opf_case10000_goc.m"
        ]

# files = TEP.select_files(dir, end_file)
# Sort files so that the smallest instances are solved first
# sort!(files, by=x->parse(Int, match(r"\d+", x).match))
# Run solver with binary decision variables
skip = ["pglib_opf_case8387_pegase.m"]

project_dir = dirname(Base.active_project())
parallel_script = joinpath(@__DIR__, "parallel_bs_pglibopf.jl")

TEP.log_header(log_file)

for (i, file) in enumerate(files[start_file:end_file])
    if file in skip
        TEP.log(params, "Skipping instance $file")
        continue
    end
    TEP.log(params, "Test $file num $(start_file + i - 1)", true)

    dir = "external/pglib-opf"
    if i == 4
        dir = "external/CATS-CaliforniaTestSystem/MATPOWER/"
    end

    filepath = "$dir/$file"
    
    params.log_dir = log_dir
    params.log_file = "$log_dir/$file"
    
    inst = TEP.build_instance(params, filepath)
    # try
        TEP.log(params, "Test $file num $(start_file + i - 1)", true)
        lp = TEP.build_lp(inst, params, scen)
        mip = TEP.build_mip(inst, params, scen)
        TEP.set_state!(inst, mip)
        cache = TEP.WorkerCache(TEP.Cache(inst, params), TEP.run_method)
        inserted = Set{TEP.CandType}(keys(inst.K))
        removed = Set{TEP.CandType}()

        # Compute initial cost for bin status report
        TEP.fix_s_vars!(lp)
        TEP.update_lp!(inst, params, lp, inserted)
        start_cost, _ = TEP.comp_penalized_cost(inst, params, scen, 
                                                lp, cache, inserted)
        TEP.unfix_s_vars!(lp)

        # Run serial D&R-BS approach
        start_time = time()
        heur_runtime = 0.0
        if is_heur_en
            inserted, _, _ = 
                TEP.run_serial_bs!(inst, params, scen, lp, cache, 
                                   inserted, removed, start_cost, start_time)
            heur_runtime = time() - start_time
        end

        TEP.fix_s_vars!(lp)
        TEP.update_lp!(inst, params, lp, inserted)
        start_cost, _ = TEP.comp_penalized_cost(inst, params, scen, 
                                                lp, cache, inserted)

        @info start_cost
        solver_runtime = 0.0
        opt_lb = 0.0
        opt_ub = 0.0
        opt_gap = 0.0
        cost = 0.0
        if TEP.fix_start!(inst, params, scen, mip, inserted)
            el = time() - start_time
            tl = max(params.beam_search.time_limit - el, 0.0)
            JuMP.set_attribute(mip.jump_model, "TimeLimit", tl)
            solver_runtime = @elapsed JuMP.optimize!(mip.jump_model)
            inserted = TEP.get_inserted_candidates(inst, mip)
            opt_lb = JuMP.objective_bound(mip.jump_model)
            opt_ub = JuMP.objective_value(mip.jump_model)
            opt_gap = JuMP.relative_gap(mip.jump_model)
        end

        TEP.update_lp!(inst, params, lp, inserted)
        cost, _ = TEP.comp_penalized_cost(inst, params, scen, 
                                          lp, cache, inserted)

        @warn opt_ub, cost
        @assert TEP.iseq(opt_ub, cost, 1e-1) "diff obj values $opt_ub $cost"

        results = TEP.init_results()

        results["start"] = start_cost
        results["best"] = cost
        results["lb"] = opt_lb
        results["ub"] = cost
        results["gap"] = 100.0 * opt_gap
        results["add_rat"] = length(inserted) / inst.num_K
        results["heur_time"] = heur_runtime
        results["solver_time"] = solver_runtime
        results["time"] = heur_runtime + solver_runtime
        results["is_feas"] = true

        TEP.log_instance(log_file, file, inst, results)
    # catch e
    #     @warn e
    #     TEP.log_instance(log_file, "<s>" * file * "</s>", inst, Dict())
    # end
end

end # module