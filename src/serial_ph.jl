"""
    run_serial_ph!(inst::Instance, params::Parameters)

Implementation of the sequential progressive hedging algorithm. 

Associated paper: https://link.springer.com/article/10.1007/s10107-016-1000-z
"""
# TODO: Consider non-binary first stage decision variables
function run_serial_ph!(inst::Instance, params::Parameters)
    log(params, "Serial solution strategy", true)

    # Initialization
    cache = Cache(inst.num_scenarios, inst.num_K)
    models = Vector{MIPModel}(undef, inst.num_scenarios)

    for it in 1:params.progressive_hedging.max_it
        # TODO: Guardar os modelos por cenário
        # TODO: Cenários como um vector, pois pode ser mais difícil causar erros
        # Imagine se o usuário dá um id que não corresponde a um índice válido
        for scen in 1:inst.num_scenarios
            mip = nothing
            if it == 1
                # TODO: Change LP objective as well
                # TODO: Run heuristic in every it
                mip = build_mip(inst, params, scen)
                set_state!(inst, mip)
                
                # (start, _) = build_solution(inst, params, scen)
                # fix_start!(inst, params, scen, mip, start)

                # Store model at first it
                models[scen] = mip
            else
                mip = models[scen]

                update_model_obj!(params, cache, scen, mip)
            end

            solve!(params, mip)
            
            update_cache_incumbent!(cache, scen, mip)
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
    end

    for scen in 1:inst.num_scenarios
        println("Scen#$(scen): $(cache.scenarios[scen].state)")
    end

    return cache
end