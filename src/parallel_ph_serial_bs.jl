"""
    run_parallel_ph_serial_bs!(inst::Instance, params::Parameters)

Implementation of the sequential progressive hedging algorithm. 

Associated paper: https://link.springer.com/article/10.1007/s10107-016-1000-z
"""
# TODO: Consider non-binary first stage decision variables
function run_parallel_ph_serial_bs!(inst::Instance, params::Parameters)
    # @info "parallel solution strategy"

    # Config JobQueueMPI
    JQM.mpi_init()
    JQM.mpi_barrier()

    if JQM.is_worker_process()
        ph_serial_bs_workers_loop(inst, params)
        JQM.mpi_barrier()
        return nothing
    end

    # Initialization
    start_time = time()
    cache = Cache(inst, params)

    controller = JQM.Controller(JQM.num_workers())

    # Config logging
    io = open(get_log_filename(inst, params, 0), "w+")
    Logging.global_logger(ConsoleLogger(io))

    ph_cost = const_infinite
    lb_best_cost = const_infinite
    ub_best_cost = const_infinite
    is_global_feas = false
    elapsed_time = time() - start_time
    it = 1

    @info "-------------------- it 0 --------------------"
    jqm_comp_costs!(inst, params, cache, it, controller, start_time)
    ph_cost, is_global_feas, lb_best_cost, ub_best_cost = 
        update_cache_start_and_best_sols!(inst, params, cache, const_infinite, 
                                          false, const_infinite, const_infinite)
    @info "time to compute initial costs(s):$(time() - start_time)"
    flush(io)

    # return elapsed_time, ph_cost, is_global_feas, 
    #         lb_best_cost, ub_best_cost, cache.best_sol 
    
    while true
        for scen in 1:inst.num_scenarios
            tl = comp_bs_time_limit(inst, params, time() - start_time)
            
            wcache = WorkerCache(cache, run_method)
            msg = ControllerMessage(wcache, it, scen, tl)
            JQM.add_job_to_queue!(controller, msg)
        end
        while !has_finished_all_jobs(controller)
            if !JQM.is_job_queue_empty(controller)
                JQM.send_jobs_to_any_available_workers(controller)
            end
            if JQM.any_pending_jobs(controller)
                job_answer = JQM.check_for_job_answers(controller)
                if !isnothing(job_answer)
                    msg = JQM.get_message(job_answer)
                    update_cache_incumbent!(cache, msg)
                    # update_cache_sol_costs_and_viols!(cache, msg)
                    update_cache_status!(cache, msg)
                end
            end
        end

        @info "-------------------- it $it --------------------"

        # Update SEP-rho parameters
        update_cache_sep_rho_x_min_max!(inst, cache)
        if it == 1
            if params.progressive_hedging.is_sep_rho_en
                update_cache_sep_rho!(inst, cache)
            else
                update_cache_cost_proportional_rho!(inst, cache)
            end
        end
        
        # Aggregation
        update_cache_x_hat!(inst, cache)
        
        # Price update
        update_cache_omega!(inst, params, cache)
        
        # Termination criterion
        update_cache_x_average!(inst, params, cache)

        update_cache_sols!(inst, params, cache)

        update_cache_best_convergence_delta!(inst, params, cache, it)

        # TODO Implementar método que, ao final da iteração, já calcula os 
        # custos e as violações das soluções lb e ub. Caso sejam inviáveis, 
        # também já repara as soluções. Para isso, utilizar os workers.
        # Implementar tudo como um único método. Dessa forma, não será mais 
        # preciso duas iterações para calcular os custos. Na próxima it, já 
        # passar os dados da melhor solução para utilizar como mipstart.

        jqm_repair_sols!(inst, params, cache, it, controller, start_time)

        jqm_comp_costs!(inst, params, cache, it, controller, start_time)

        log_status(inst, cache)
        @info "ph time(s):$(time() - start_time)"

        ph_cost, is_global_feas, lb_best_cost, ub_best_cost = 
            update_cache_start_and_best_sols!(inst, params, cache, ph_cost, 
                                    is_global_feas, lb_best_cost, ub_best_cost)

        elapsed_time = time() - start_time
        if isl(cache.best_convergence_delta, 
               params.progressive_hedging.convergence_eps)
            @info "convergence reached: $(cache.best_convergence_delta)"
            break
            # break
        elseif isg(elapsed_time, params.progressive_hedging.time_limit)
            @info "ph time limit reached $(round(elapsed_time, digits = 2))"
            break
        elseif it == params.progressive_hedging.max_it
            @info "max num it reached"
            break
        end

        # update_cache_detect_cycles!(inst, cache)

        it += 1
        flush(io)
    end

    # for scen in 1:inst.num_scenarios
    #     println("Scen#$(scen): $(cache.scenarios[scen].state)")
    # end

    # JQM.mpi_barrier()
    JQM.send_termination_message()
    JQM.mpi_finalize()

    # ph_cost = comp_ph_cost(inst, params, cache)
    @info "obj:$ph_cost elapsed_time:$elapsed_time"
    if params.debugging_level == 1
        mip, subproblems = build_deterministic(inst, params)
        @info "fix the start of the model"
        det_cost = fix_start!(inst, mip, subproblems, cache.best_sol)
        @info "$ph_cost $det_cost"
        @assert iseq(ph_cost, det_cost, 1e-3) "ph diff costs $ph_cost $det_cost"
    end
    close(io)

    return elapsed_time, ph_cost, is_global_feas, 
            lb_best_cost, ub_best_cost, cache.best_sol 
end

function ph_serial_bs_workers_loop(inst::Instance, params::Parameters)
    if JQM.is_worker_process()
        worker = JQM.Worker()
        # Build models for the first scenario
        current_model_scen = 1

        # Set the number of threads for the LP models
        num_threads = params.solver.num_threads
        params.solver.num_threads = 1

        mip = build_mip(inst, params, current_model_scen)
        set_state!(inst, mip)

        lp = build_lp(inst, params, current_model_scen)

        # Reset the number of threads to the default value
        params.solver.num_threads = num_threads

        while true
            job = JQM.receive_job(worker)
            msg = JQM.get_message(job)
            if msg == JQM.TerminationMessage()
                break
            end

            state_values = Float64[]
            sol_info_lb = SolutionInfo(0.0, 0.0, 0.0, Set{CandType}())
            sol_info_ub = SolutionInfo(0.0, 0.0, 0.0, Set{CandType}())
            start_time = time()
            bs_runtime = 0.0
            solver_runtime = 0.0
            bin_rm_rat = 0.0
            bs_rm_rat = 0.0
            repair_st = (false, false)
            reinsert_st = (false, false, 0)

            io = open(get_log_filename(inst, params, msg.scen), "a")
            # Logging.with_logger(ConsoleLogger(io[msg.scen])) do
            Logging.with_logger(ConsoleLogger(io)) do
                # Update the model according to the current scenario
                if msg.scen != current_model_scen
                    # update_model!(inst, params, msg.scen, mip)
                    update_model!(inst, params, msg.cache, msg.scen, lp)
                    current_model_scen = msg.scen
                end

                if msg.cache.option == repair_sols
                    # Repair both lower bound and upper bound solutions, if 
                    # needed
                    sol_info_lb.reinsert, lb_rep_st, lb_rein_st = 
                                        repair!(inst, params, msg.cache, 
                                                msg.scen, lp, msg.cache.sol_lb)

                    sol_info_ub.reinsert, ub_rep_st, ub_rein_st = 
                                        repair!(inst, params, msg.cache, 
                                                msg.scen, lp, msg.cache.sol_ub)

                    # Update tuples with results from both lb and ub repairs
                    repair_st = map((a, b) -> a + b, lb_rep_st, ub_rep_st)
                    reinsert_st = map((a, b) -> a + b, lb_rein_st, ub_rein_st)


                elseif msg.cache.option == comp_g_costs
                    # Compute generation costs 
                    sol_info_lb, sol_info_ub = 
                                    comp_sol_info_lb_ub!(inst, params, msg.scen, 
                                                            lp, msg.cache)


                elseif msg.cache.option == run_method
                    @info "-------------------- it $(msg.it) " * 
                            "--------------------"
                    # Used in utils:bs.jl:comp_penalized_cost
                    params.progressive_hedging.is_en = msg.it > 1
                    
                    params.beam_search.time_limit = msg.time_limit

                    inserted = msg.cache.inserted
                    # TODO Compute once in controller
                    removed = setdiff(Set{CandType}(keys(inst.K)), inserted)
                    cost = msg.cache.costs[msg.scen]

                    bs_start_time = time()
                    inserted, bin_rm_rat, bs_rm_rat = 
                        run_serial_bs!(inst, params, msg.scen, lp, msg.cache, 
                                        inserted, removed, cost, start_time)

                    bs_runtime = time() - bs_start_time
                    solver_runtime = 0.0
                    
                    tl = max(msg.time_limit - (time() - start_time), 0.0)

                    has_mip_state_vals = false
                    if isg(tl, 0.0)
                        update_model!(inst, params, msg.cache, msg.scen, mip)

                        # If feasible, continue the execution
                        if fix_start!(inst, params, msg.scen, mip, inserted)
                            el = time() - start_time
                            tl = max(msg.time_limit - el, 0.0)

                            set_attribute(mip.jump_model, "TimeLimit", tl)

                            solver_runtime = 
                                        @elapsed JuMP.optimize!(mip.jump_model)

                            has_mip_state_vals = 
                                        JuMP.result_count(mip.jump_model) > 0
                        end
                    end

                    if has_mip_state_vals
                        state_values = get_state_values(mip)
                    else
                        if params.debugging_level == 1
                            is_opt = JuMP.termination_status(lp.jump_model) == 
                                        MOI.OPTIMAL
                            @assert is_opt "Not opt scen#$(msg.scen)"
                        end
                        state_values = get_state_values(inst, inserted)
                    end
                end
                # flush(io[msg.scen])
                flush(io)
            end
            close(io)

            status = ScenarioStatus(bs_runtime, solver_runtime, bin_rm_rat, 
                                    bs_rm_rat, repair_st, reinsert_st)
            ret_msg = WorkerMessage(state_values, msg.it, msg.scen, sol_info_lb, 
                                    sol_info_ub, status)

            JQM.send_job_answer_to_controller(worker, ret_msg)
        end
        # Close log files
        # close.(io)
        exit(0)
    end

    return nothing
end