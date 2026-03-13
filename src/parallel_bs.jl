function run_parallel_bs!(inst::Instance, 
                          params::Parameters, 
                          scen::Int64, 
                          lp::LPModel, 
                          is_ph::Bool = false, 
                          cache::Cache = Cache())
    @info "Parallel beam search heuristic"

    # Config JobQueueMPI
    JQM.mpi_init()

    start_time = time()
    status = nothing
    inserted = nothing
    candidates = nothing
    if JQM.is_controller_process()
        # Build initial solution
        _, status, inserted, candidates = 
                        build_solution!(inst, params, scen, lp, is_ph, cache)
        fix_s_vars!(inst, lp)
    end

    JQM.mpi_barrier()

    @warn JQM.is_worker_process(), JQM.num_workers()
    if JQM.is_worker_process()
        params.solver.num_threads = 1
        bs_workers_loop(inst, params)
        return nothing
    end

    # Initialization
    controller = JQM.Controller(JQM.num_workers())

    params.solver.num_threads = 8

    remaining_time = params.beam_search.time_limit - status.time

    K = collect(candidates)

    num_ins_start = length(inserted)
    num_ins = 0
    obj = 0.0
    best_obj = 0.0
    obj_start = 0.0

    # num_candidates_per_batch = params.beam_search.num_candidates_per_batch

    it = 1
    scen = 1
    for bs_it in 1:params.beam_search.num_max_it
        update_lp!(inst, params, lp, inserted)
        obj = comp_penalized_cost(inst, params, scen, lp, inserted)
        # params.beam_search.num_candidates_per_batch = 
        #                                             num_candidates_per_batch

        num_it_wo_impr = 0
        if bs_it == 1
            best_obj = obj
            obj_start = obj
        end

        root = Node(obj, get_values(lp.f), 0.0, collect(inserted), K, [])
        Q = [root]
        while true
            if isg(time() - start_time, remaining_time)
                @info "BS Time limit reached"
                break
            end
            LB = []

            # Send jobs
            for (i, node) in enumerate(Q)
                lines = select_lines(inst, params, lp, 
                                     node, node.inserted)
                for k in lines
                    # Estimates the number of nodes required per level so that 
                    # the number of threads per node is set accordingly
                    num_threads = div(JQM.num_workers(), 
                                (params.beam_search.num_children_per_parent * 
                                 length(Q)))
                    msg = BSControllerMessage(it, i, node, 
                                              best_obj, k, num_threads)
                    JQM.add_job_to_queue!(controller, msg)
                end
            end

            has_impr = false
            # Recover results
            while !has_finished_all_jobs(controller)
                if !JQM.is_job_queue_empty(controller)
                    JQM.send_jobs_to_any_available_workers(controller)
                end
                if JQM.any_pending_jobs(controller)
                    job_answer = JQM.check_for_job_answers(controller)
                    if !isnothing(job_answer)
                        msg = JQM.get_message(job_answer)
                        node = Q[msg.node_idx]

                        add_node!(params, LB, best_obj, msg, node)

                        if msg.is_feas && isl(msg.obj, best_obj)
                            has_impr = true
                            best_obj = msg.obj
                            inserted = node.inserted
                            num_ins = length(LB[end].inserted)

                            # Log info
                            st = Status("bs", num_ins_start - num_ins, 
                                        inst.num_K, obj, obj_start, start_time)
                            log(params, st)
                        end
                    end
                end
            end
            params.log_level == 4 && flush(io)

            if has_impr
                num_it_wo_impr = 0
            else
                num_it_wo_impr += 1

                if num_it_wo_impr >= params.beam_search.num_max_it_wo_impr
                    @info "Max it wo impr reached"
                    break
                end
            end

            if length(LB) == 0
                @info "No new nodes to investigate"
                break
            end
            Q = select!(params, LB, K)
            
            it += 1
        end

        if bs_it >= params.beam_search.num_max_it || 
           isg(time() - start_time, remaining_time)
            break
        else
            @info "Change shuffle strategy"
            params.beam_search.is_shuffle_en = ~params.beam_search.is_shuffle_en
        end
    end
    elapsed_time = time() - start_time

    JQM.send_termination_message()
    JQM.mpi_finalize()

    @info it, 
          num_ins_start, num_ins, 
          num_ins / num_ins_start, 
          elapsed_time

    @info "Build full model"
    build_time = @elapsed (mip = build_mip(inst, params))

    update_lp!(inst, params, lp, inserted)
    # g_bus, bus_inj, viol_update = get_data(inst, lp, ptdf)

    f = get_values(lp.f)
    g = get_values(lp.g)
    start = Start(Set(inserted), f, g)

    results = init_results()

    @info "Fix the start of the model"
    fix_start_time = @elapsed (start_ub = fix_start!(inst, params, mip, start))
    
    params.beam_search.time_limit = max(remaining_time - elapsed_time, 0.0)
    
    # @info "Solve the model"
    # results = solve!(inst, params, mip)
    
    results["rnr_time"] = status.time
    results["rnr_rm_rat"] = status.rm_ratio
    results["rnr_impr_rat"] = status.impr_ratio
    results["fix_start_time"] = fix_start_time
    results["bs_time"] = elapsed_time
    results["start_ub"] = start_ub

    results["build_time"] = build_time
    # results["build_obj_rat"] = comp_build_obj_rat(inst, 
    #                                               results["objective"], 
    #                                               start.inserted)

    @info "Obj $(JuMP.objective_value(mip.jump_model))"

    params.log_level == 4 && flush(io)
    params.log_level == 4 && close(io)
    # params.log_level == 4 && redirect_stdout(stdout)

    return results
end
    
function bs_workers_loop(inst::Instance, params::Parameters)
    if JQM.is_worker_process()
        worker = JQM.Worker()
        # Build the model for the first scenario
        # Solve subscenarios with a single thread
        num_threads = params.solver.num_threads
        params.solver.num_threads = 1
        lp = build_lp(inst, params, scen)
        params.solver.num_threads = num_threads

        fix_s_vars!(inst, lp)
        scen = 1
        while true
            job = JQM.receive_job(worker)
            msg = JQM.get_message(job)
            if msg == JQM.TerminationMessage()
                break
            end

            set_attribute(lp.jump_model, 
                          MOI.RawOptimizerAttribute("Threads"), 
                          msg.num_threads)

            # TODO: Reduzir os dados em msg para o mínimo necessário
            ins_candidates = setdiff(msg.node.inserted, msg.k)
            update_lp!(inst, params, lp, ins_candidates)

            obj = const_infinite
            is_feas = false
            viol = 0.0
            # The flow values for the parent node
            f = msg.node.f
            # if JuMP.has_values(lp.jump_model)
            if JuMP.termination_status(lp.jump_model) == MOI.OPTIMAL
                obj = comp_penalized_cost(inst, params, scen, lp, ins_candidates)
                f = get_values(lp.f)
                is_feas = true
                # viol = comp_viol(lp)
            end
            
            ret_msg = BSWorkerMessage(msg.node_idx, msg.k, is_feas, 
                                      obj, f, viol)

            JQM.send_job_answer_to_controller(worker, ret_msg)
        end
        exit(0)
    end

    return nothing
end