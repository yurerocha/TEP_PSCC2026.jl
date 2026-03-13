"""
    run_parallel_ph!(inst::Instance, params::Parameters)

Implementation of the sequential progressive hedging algorithm. 

Associated paper: https://link.springer.com/article/10.1007/s10107-016-1000-z
"""
# TODO: Consider non-binary first stage decision variables
function run_parallel_ph!(inst::Instance, params::Parameters)
    log(params, "Parallel solution strategy", true)

    # Config JobQueueMPI
    JQM.mpi_init()
    JQM.mpi_barrier()

    if JQM.is_worker_process()
        ph_workers_loop(inst, params)
        JQM.mpi_barrier()
        return nothing
    end

    # Initialization
    cache = Cache(inst.num_scenarios, inst.num_K)
    models = Vector{MIPModel}(undef, inst.num_scenarios)

    controller = JQM.Controller(JQM.num_workers())

    for it in 1:params.progressive_hedging.max_it
        # TODO: Guardar os modelos por cenário
        # TODO: Cenários como um vector, pois pode ser mais difícil causar erros
        # Imagine se o usuário dá um id que não corresponde a um índice válido
        for scen in 1:inst.num_scenarios
            # MPI does no copy state variables
            msg = ControllerMessage(cache, it, scen)
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
                end
            end
        end
        # Aggregation
        update_cache_x_hat!(inst, cache)

        # Price update
        update_cache_omega!(inst, params, cache)

        # Termination criterion
        update_cache_x_average!(inst, cache)

        update_cache_best_convergence_delta!(inst, cache, it)

        if isl(cache.best_convergence_delta, 
               params.progressive_hedging.convergence_eps)
            log(params, "Convergence reached: $(cache.best_convergence_delta)")
            break
        end
        # readline()
    end

    for scen in 1:inst.num_scenarios
        println("Scen#$(scen): $(cache.scenarios[scen].state)")
    end

    # JQM.mpi_barrier()
    JQM.send_termination_message()
    JQM.mpi_finalize()

    return cache
end

function ph_workers_loop(inst::Instance, params::Parameters)
    if JQM.is_worker_process()
        worker = JQM.Worker()
        # Build model for the first scenario
        current_model_scen = 1
        mip = build_mip(inst, params, current_model_scen)
        set_state!(inst, mip)
        while true
            job = JQM.receive_job(worker)
            msg = JQM.get_message(job)
            if msg == JQM.TerminationMessage()
                break
            end

            # Update the model according to the current scenario
            if msg.scen != current_model_scen
                update_model!(inst, params, msg.scen, mip)
                current_model_scen = msg.scen
            end

            if msg.it > 1
                update_model_obj!(params, msg.cache, msg.scen, mip)
            end

            solve!(inst, params, mip)
            
            ret_msg = WorkerMessage(get_state_values(mip), msg.it, msg.scen)

            JQM.send_job_answer_to_controller(worker, ret_msg)
        end
        exit(0)
    end

    return nothing
end