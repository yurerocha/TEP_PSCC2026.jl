"""
    run(logname::String, 
        gl_strategy::Int64, 
        gl_ins::Float64, 
        rnr_percent::Float64, 
        rnr_delta::Float64, 
        rnr_time_limit::Float64)

Solve all instances.
"""
# function run(logname::String = "log.md", 
#              gl_strategy::Int64 = 2, 
#              gl_ins::Float64 = 0.1, 
#              rnr_percent::Float64 = 0.8, 
#              rnr_delta::Float64 = 0.1, 
#              rnr_time_limit::Float64 = 5.0,
#              is_symmetry_en::Bool = false)
function run(logname::String = "log.md")
    params = Parameters()

    # Alterar logfile, start e end files
    # tun1: 0.4, 0.1
    # tun2: 0.5, 0.1
    logfile = "tep_wo_bs.md"
    start_file = 1 # 40
    end_file = 1 # 62
    is_heur_en = true
    log_dir = "log_wo_bs/"
    
    # Number of seconds since the Unix epoch
    # seed = Int(floor(datetime2unix(now())))
    # Random.seed!(seed)
    # rng = Random.MersenneTwister(123)

    log_header(logfile)

    dir = "external/pglib-opf"
    files = select_files(dir, end_file)
    # Sort files so that the smallest instances are solved first
    sort!(files, by=x->parse(Int, match(r"\d+", x).match))
    # Run solver with binary decision variables
    skip = []
    # skip = [
    #         "pglib_opf_case2853_sdet.m",    # Numerical trouble encountered
    #         "pglib_opf_case9241_pegase.m",  # Numerical trouble encountered
    #         "pglib_opf_case13659_pegase.m"  # Numerical trouble encountered
    # ]
    # files = ["pglib_opf_case4020_goc.m", "pglib_opf_case4661_sdet.m", 
    #          "pglib_opf_case6470_rte.m", "pglib_opf_case6515_rte.m", 
    #          "pglib_opf_case10000_goc.m"]
    counter = start_file
    files = ["pglib_opf_case8387_pegase.m"]
    for file in files[start_file:end_file]
        if file in skip
            log(params, "Skipping instance $file num $counter")
            continue
        end
        log(params, "Processing $file num $counter", true)

        filepath = "$(dir)/$file"
        logsolver = "$(dir)/log/$file"

        params.log_file = log_dir * file

        inst = build_instance(params, filepath)

        mip = nothing
        build_time = 0.0
        start_time = 0.0
        results = init_results()

        try
            log(params, "Build full model", true)
            build_time = @elapsed (mip = build_mip(inst, params))

            # TODO: Add parameter to indicate if an initial solution will be 
            # used
            log(params, "Build heuristic solution", true)
            (start, status) = binary_search(inst, params, is_heur_en)

            log(params, "Fix the start of the model", true)
            start_time = 
                    @elapsed (obj = fix_start!(inst, params, mip, start))
            results["objective"] = obj

            params.solver.time_limit -= (start_time + status.time)

            # TODO: Add tuning flag
            log(params, "Solve the model", true)
            results = solve!(inst, params, mip)

            results["build_time"] = build_time
            if isa(results["objective"], Number)
                results["build_obj_rat"] = comp_build_obj_rat(inst, 
                                                        results["objective"], 
                                                        start.inserted)
            end

            if is_heur_en
                results["rnr_impr_rat"] = status.impr_ratio
                results["rnr_time"] = status.time
                results["fix_start_time"] = start_time
            end
            
            log_instance(logfile, file, inst, results)
        catch e
            # @warn e
            log_instance(logfile, "<s>" * file * "</s>", inst, Dict())
        end

        # readline()
        counter += 1
    end
end

function tuning()
    # gl_strategy = [1, 2, 3, 4]
    # gl_ins = [0.1, 0.2, 0.4]
    # count_exp = 1
    # for s in gl_strategy, i in gl_ins
    #     @warn "Tuning experiment gl_strategy:$s gl_ins:$i"
    #     run("exp$(count_exp).md", s, i, 0.8, 0.1, 5.0)
    #     count_exp += 1
    # end
    # Best: 2, 0.1
    # New best: 4, 0.2
    # rnr_percent_delta = [[0.8, 0.2], [0.9, 0.1]]
    # count_exp = 13
    # for v in rnr_percent_delta
    #     @warn "Tuning experiment rnr_percent:$(v[1]) rnr_delta:$(v[2])"
    #     run("exp$(count_exp).md", 4, 0.2, v[1], v[2], 5.0)
    #     count_exp += 1
    # end
    # Best: 0.8, 0.1
    # New best: 0.8, 0.2
    count_exp = 15
    time_limit = [7.5, 10, 12.5, 15.0]
    for tl in time_limit
        @warn "Tuning experiment time_limit:$tl"
        run("exp$(count_exp).md", 4, 0.2, 0.9, 0.1, tl)
        count_exp += 1
    end
    # Best: 15.0
    # count_exp = 22
    # is_symmetry_en = [false, true]
    # for is_en in is_symmetry_en
    #     run("exp$(count_exp).md", 2, 0.1, 0.8, 0.1, 15.0, is_en)
    #     count_exp += 1
    # end
end

function solve(filepath::String, num_scenarios::Int64 = 1)
    params = Parameters()
    params.log_file *= "/" * get_inst_name(filepath) * ".txt"

    params.solution_strategy = Deterministic()

    inst = nothing
    if occursin("CATS-CaliforniaTestSystem", filepath)
        # Read the CATS instance, with multiple scenarios
        inst, _ = build_cats_instance(params, num_scenarios)
    else
        # Read the pglib-opf instances with single scenarios
        inst = build_instance(params, filepath)
        build_scenarios!(inst, num_scenarios - 1, 0.25)
    end

    if params.solution_strategy isa Deterministic
        solve_deterministic!(inst, params)
    elseif params.solution_strategy isa Serial
        run_serial_ph!(inst, params)
    elseif params.solution_strategy isa Parallel
        run_parallel_ph!(inst, params)
    end

    return nothing
end

function run_ptdf(filepath::String)
    params = Parameters()
    params.log_file *= "/" * get_inst_name(filepath) * ".txt"

    inst = build_instance(params, filepath)

    # build_ptdf(inst, params, 1, T = Float64)
    build_ptdf_system(inst, params, 1, T = Float64)

    return nothing
end