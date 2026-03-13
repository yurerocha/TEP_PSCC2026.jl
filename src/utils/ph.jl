# --------------------------- Parallel and Serial PH ---------------------------
# function update_cache_src_obj!(cache::Cache, scen::Int64, mip::MIPModel)
#     cache.scenarios[scen].src_obj = JuMP.objective_function(mip.jump_model)
#     return nothing
# end

function update_cache_incumbent!(cache::Cache, scen::Int64, mip::MIPModel)
    cache.scenarios[scen].state = get_state_values(mip)
    return nothing
end

function update_cache_x_hat!(inst::Instance, cache::Cache)
    cache.x_hat = sum(inst.scenarios[scen].p * cache.scenarios[scen].state 
                      for scen in 1:inst.num_scenarios)
    return nothing
end

function update_cache_omega!(inst::Instance, params::Parameters, cache::Cache)
    for scen in 1:inst.num_scenarios
        delta = cache.scenarios[scen].state - cache.x_hat
        cache.scenarios[scen].omega += params.progressive_hedging.rho * delta
    end
    return nothing
end

function update_cache_x_average!(inst::Instance, 
                                 params::Parameters, 
                                 cache::Cache)
    cache.x_average = sum(inst.scenarios[scen].p * cache.scenarios[scen].state 
                          for scen in 1:inst.num_scenarios) / 
                    sum(inst.scenarios[scen].p for scen in 1:inst.num_scenarios)

    return nothing
end

function update_cache_sols!(inst::Instance, 
                            params::Parameters, 
                            cache::Cache)
    th = params.progressive_hedging.lb_threshold
    cache.sol_lb.insert = Set{CandType}([k for (i, k) in enumerate(keys(inst.K)) 
                                            if isg(cache.x_average[i], th)])
    cache.sol_ub.insert = Set{CandType}([k for (i, k) in enumerate(keys(inst.K)) 
                                            if !iseq(cache.x_average[i], 0.0)])

    LoggingExtras.withlevel(Info; verbosity = params.log_level) do
        n = inst.num_K - length(cache.sol_ub.insert)
        t = cache.sol_lb.count_use + cache.sol_ub.count_use
        # Avoid NaN in init_sol(%) inter and union
        t = max(1, t)
        @infov 1 "candidates(%) discard:$(roundp(n, inst.num_K)) " * 
                 "lb:$(roundp(cache.sol_lb.insert, inst.num_K)) " * 
                 "ub:$(roundp(cache.sol_ub.insert, inst.num_K))"
        @infov 2 "init_sol(%) " * 
                 "lb:$(roundp(cache.sol_lb.count_use, t)) " * 
                 "ub:$(roundp(cache.sol_ub.count_use, t))"
    end
end

function update_cache_sol_costs_and_viols!(cache::Cache, msg::WorkerMessage)
    cache.sol_lb.costs[msg.scen] = msg.sol_info_lb.cost
    cache.sol_ub.costs[msg.scen] = msg.sol_info_ub.cost

    cache.sol_lb.g_costs[msg.scen] = msg.sol_info_lb.g_cost
    cache.sol_ub.g_costs[msg.scen] = msg.sol_info_ub.g_cost

    cache.sol_lb.viols[msg.scen] = msg.sol_info_lb.viol
    cache.sol_ub.viols[msg.scen] = msg.sol_info_ub.viol

    return nothing
end

function update_cache_status!(cache::Cache, msg::WorkerMessage)
    cache.status[msg.scen] = msg.status
    return nothing
end

function log_status(inst::Instance, cache::Cache)
    bs_avg_rt = 0.0
    solver_avg_rt = 0.0
    bin_avg_rm_rat = 0.0
    bin_min_rm_rat = 1.0
    bin_min_rm_scen = 0
    bs_avg_rm_rat = 0.0
    bs_min_rm_rat = 1.0
    bs_min_rm_scen = 0
    repair_count_use = 0
    repair_count_success = 0
    reinsert_count_use = 0
    reinsert_count_success = 0
    reinsert_avg_it = 0
    reinsert_max_it = 0
    reinsert_max_it_scen = 0

    for scen in eachindex(inst.scenarios)
        bs_avg_rt += cache.status[scen].beam_search_runtime
        solver_avg_rt += cache.status[scen].solver_runtime

        bin_avg_rm_rat += cache.status[scen].bin_search_rm_ratio
        if isl(cache.status[scen].bin_search_rm_ratio, bin_min_rm_rat)
            bin_min_rm_rat = cache.status[scen].bin_search_rm_ratio
            bin_min_rm_scen = scen
        end

        bs_avg_rm_rat += cache.status[scen].beam_search_rm_ratio
        if isl(cache.status[scen].beam_search_rm_ratio, bs_min_rm_rat)
            bs_min_rm_rat = cache.status[scen].beam_search_rm_ratio
            bs_min_rm_scen = scen
        end
        
        repair_count_use += cache.status[scen].repair[1]
        repair_count_success += cache.status[scen].repair[2]

        reinsert_count_use += cache.status[scen].reinsert[1]
        reinsert_count_success +=  cache.status[scen].reinsert[2]
        reinsert_avg_it += cache.status[scen].reinsert[3]
        if cache.status[scen].reinsert[3] > reinsert_max_it
            reinsert_max_it = cache.status[scen].reinsert[3]
            reinsert_max_it_scen = scen
        end
    end

    bs_avg_rt /= inst.num_scenarios
    solver_avg_rt /= inst.num_scenarios
    total_avg_rt = bs_avg_rt + solver_avg_rt

    @info "avg runtime(%) bs:$(roundp(bs_avg_rt, total_avg_rt))" *
            " solver:$(roundp(solver_avg_rt, total_avg_rt))"

    bin_avg_rm_rat = roundp(bin_avg_rm_rat / inst.num_scenarios, 1.0)
    bin_min_rm_rat = roundp(bin_min_rm_rat, 1.0)
    bs_avg_rm_rat = roundp(bs_avg_rm_rat / inst.num_scenarios, 1.0)
    bs_min_rm_rat = roundp(bs_min_rm_rat, 1.0)

    @info "avg rm rat(%) bin:$bin_avg_rm_rat bs:$bs_avg_rm_rat"
    @info "min rm rat(%) bin:$bin_min_rm_rat(scen#$bin_min_rm_scen) " * 
            "bs:$bs_min_rm_rat(scen#$bs_min_rm_scen)"

    if repair_count_use > 0
        @info "avg repair rat(%):" * 
                            "$(roundp(repair_count_success, repair_count_use))"
    end

    if reinsert_count_use > 0
        reinsert_avg_it = 
                        round(reinsert_avg_it / reinsert_count_use, digits = 2)
        @info "avg reinsert rat(%):" * 
                "$(roundp(reinsert_count_success, reinsert_count_use)) " * 
                "it:$reinsert_avg_it max_it:$reinsert_max_it" * 
                "(scen#$reinsert_max_it_scen)"
    end


    return nothing
end

"""
    update_cache_best_convergence_delta!(inst::Instance, 
                                         params::Parameters, 
                                         cache::Cache, 
                                         it::Int64)

Compute and update the best convergence delta as in 
https://doi.org/10.1007/s10287-010-0125-4
"""
function update_cache_best_convergence_delta!(inst::Instance, 
                                              params::Parameters, 
                                              cache::Cache, 
                                              it::Int64)
    conv_delta = 0.0
    x_avg = cache.x_average
    for scen in eachindex(inst.scenarios)
        x =  cache.scenarios[scen].state
        # cache.deltas[scen] = maximum(abs, x - cache.x_average)
        cache.deltas[scen] = 
            sum(abs(x[i] - x_avg[i]) / x_avg[i] for i in 1:inst.num_K 
                                            if isg(x_avg[i], 0.0); init = 0.0) 
    end
    # cache.t_deviation = maximum(cache.deltas)
    # Compute the average per-scenario deviation from the average solution
    cache.t_deviation = sum(cache.deltas) / inst.num_scenarios
    cache.quality_deviation = comp_quality_deviation(inst, cache)

    LoggingExtras.withlevel(Info; verbosity = params.log_level) do
        @infov 1 "ph td:$(round(cache.t_deviation, digits = 2)) " * 
                 "best_td:$(round(cache.best_convergence_delta, digits = 2)) " * 
                 "qd:$(round(cache.quality_deviation, digits = 2))"
        # @infov 2 join(round.(cache.deltas, digits = 2), " ")
    end

    if isl(cache.t_deviation, cache.best_convergence_delta)
        cache.best_convergence_delta = cache.t_deviation
        cache.best_it = it
        # cache.best_sol = cache.sol_ub.insert

        return true
    end

    return false
end

"""
    update_model_obj!(params::Parameters, 
                      cache::T, 
                      scen::Int64, 
                      tep::TEPModel) where T <: Union{Cache, WorkerCache}

Update the model objective function according to the progressive hedging 
algorithm.
"""
function update_model_obj!(params::Parameters, 
                           cache::T, 
                           scen::Int64, 
                           tep::TEPModel) where T <: Union{Cache, WorkerCache}
    delta_obj = comp_delta_obj(params, cache, scen, tep)
    @objective(tep.jump_model, Min, tep.obj + delta_obj)
    return nothing
end

"""
    comp_delta_obj(params::Parameters, 
                   cache::T, 
                   scen::Int64, 
                   tep::TEPModel) where T <: Union{Cache, WorkerCache}

Compute the delta objective associated with the last progressive hedging 
iteration to incorporate in the new objective function.
"""
function comp_delta_obj(params::Parameters, 
                        cache::T, 
                        scen::Int64, 
                        tep::TEPModel) where T <: Union{Cache, WorkerCache}
    x = tep.jump_model.ext[:state]
    x_hat = cache.x_hat 
    # Piece-wise linear function for the squared two-norm:
    #     (a - b)² = a² - 2ab + b² (binary first stage variables)
    # TODO: decide which linearization to use based on the type of the variable
    squared_twonorm = x - 2.0 * x .* x_hat + x_hat .^ 2

    penalty::Float64 = params.progressive_hedging.rho / 2.0

    return cache.scenarios[scen].omega' * x + penalty * sum(squared_twonorm)
end

function get_state(params::Parameters, mip::MIPModel)
    if :state in keys(mip.jump_model.ext)
        return mip.jump_model.ext[:state]
    else
        log(params, "Error: symbol :state not found in model", true)
    end

    return nothing
end

"""
    update_cache_sep_rho_x_min_max!(inst::Instance, cache::Cache)

Update SEP-rho x_min and x_max.
"""
function update_cache_sep_rho_x_min_max!(inst::Instance, cache::Cache)
    cache.sep_rho_x_min = cache.scenarios[1].state
    cache.sep_rho_x_max = cache.scenarios[1].state

    for scen in 2:inst.num_scenarios
        cache.sep_rho_x_min = 
                    min.(cache.sep_rho_x_min, cache.scenarios[scen].state)
        cache.sep_rho_x_max = 
                    max.(cache.sep_rho_x_max, cache.scenarios[scen].state)
    end

    return nothing
end

"""
    update_cache_sep_rho!(inst::Instance, cache::Cache)

Update rho according to the SEP-rho heuristic.
"""
function update_cache_sep_rho!(inst::Instance, cache::Cache)
    for k in keys(inst.K)
        i = inst.key_to_index[k]
        cache.rho[i] = inst.K[k].cost / 
                        (cache.sep_rho_x_max[i] - cache.sep_rho_x_min[i] + 1)
    end

    return nothing
end

"""
    update_cache_cost_proportional_rho!(inst::Instance, 
                                        cache::Cache, 
                                        cost_multiplier::Float64 = 1.0)

Set rho(k) equal to a multiplier of the candidate unit cost C^{inv}_k.
"""
function update_cache_cost_proportional_rho!(inst::Instance, 
                                             cache::Cache, 
                                             cost_multiplier::Float64 = 1.0)
    for k in keys(inst.K)
        i = inst.key_to_index[k]
        cache.rho[i] = cost_multiplier * inst.K[k].cost
    end

    return nothing
end

"""
    comp_ph_cost(inst::Instance, params::Parameters, cache::Cache)

Compute multi-scenario problem objective value.
"""
function comp_ph_cost(inst::Instance, params::Parameters, cache::Cache)
    build = cache.sol_ub.insert
    cost = comp_build_cost(inst, build)

    for scen in 1:inst.num_scenarios
        lp = build_lp(inst, params, scen)
        update_lp!(inst, params, lp, build)
        cost += inst.scenarios[scen].p * comp_g_cost(inst, params, scen, lp)
    end

    return cost
end

function comp_ph_penalized_cost(inst::Instance, 
                                params::Parameters, 
                                cache::Cache)
    lb_build = comp_build_cost(inst, cache.sol_lb.insert)
    ub_build = comp_build_cost(inst, cache.sol_ub.insert)

    lb_cost = lb_build
    ub_cost = ub_build
    lb_viol = 0.0
    ub_viol = 0.0

    g_lb = 0.0
    g_ub = 0.0
    for scen in eachindex(inst.scenarios)
        lb_cost += inst.scenarios[scen].p * cache.sol_lb.g_costs[scen]
        ub_cost += inst.scenarios[scen].p * cache.sol_ub.g_costs[scen]
        g_lb += inst.scenarios[scen].p * cache.sol_lb.g_costs[scen]
        g_ub += inst.scenarios[scen].p * cache.sol_ub.g_costs[scen]

        lb_viol += cache.sol_lb.viols[scen]
        ub_viol += cache.sol_ub.viols[scen]
    end
    @info "costs lb:$(round(lb_cost, digits = 2)) " * 
          "ub:$(round(ub_cost, digits = 2))"

    @info "build lb:$(round(lb_build, digits = 2)) " * 
          "ub:$(round(ub_build, digits = 2))"

    g_lb = round(g_lb, digits = 2)
    g_ub = round(g_ub, digits = 2)
    @info "g lb:$g_lb ub:$g_ub"

    @info "build(%) lb:$(roundp(lb_build, lb_cost)) " * 
          "ub:$(roundp(ub_build, ub_cost))"

    pen = params.progressive_hedging.penalty_mult
    lb_cost += pen * lb_viol
    ub_cost += pen * ub_viol

    @info "viols lb:$(round(lb_viol, digits = 2)) " * 
          "ub:$(round(ub_viol, digits = 2))"

    @info "pcosts lb:$(round(lb_cost, digits = 2)) " * 
          "ub:$(round(ub_cost, digits = 2))"

    lb_is_feas = iseq(lb_viol, 0.0)
    ub_is_feas = iseq(ub_viol, 0.0)

    best_cost = lb_cost
    is_feas = lb_is_feas || ub_is_feas
    best_sol = cache.sol_lb.insert
    best_costs = cache.sol_lb.costs
    cache.sol_lb.count_use += 1

    if (lb_is_feas == ub_is_feas && isl(ub_cost, lb_cost)) || 
            (!lb_is_feas && ub_is_feas)
        best_cost = ub_cost
        best_sol = cache.sol_ub.insert
        best_costs = cache.sol_ub.costs
        cache.sol_lb.count_use -= 1
        cache.sol_ub.count_use += 1
    end

    return best_cost, is_feas, lb_cost, ub_cost, copy(best_sol), best_costs
end

"""
    update_cache_start_and_best_sols!(inst::Instance, params::Parameters, 
                                cache::Cache, best_cost::Float64, 
                                is_global_feas::Bool, 
                                lb_best_cost::Float64, ub_best_cost::Float64)

Update the cache best solution and return the new costs.
"""
function update_cache_start_and_best_sols!(inst::Instance, params::Parameters, 
                                cache::Cache, best_cost::Float64, 
                                is_global_feas::Bool, 
                                lb_best_cost::Float64, ub_best_cost::Float64)
    cost, is_feas, lb_cost, ub_cost, sol, costs = 
                                    comp_ph_penalized_cost(inst, params, cache)

    # 0 0 && isl
    # 0 1 true
    # 1 0 false
    # 1 1 && isl
    @info "status g_feas:$is_global_feas feas:$is_feas"
    if !is_global_feas && is_feas || 
        (is_global_feas == is_feas && isl(cost, best_cost))
        best_cost = cost
        cache.best_sol = copy(sol)
        if is_feas
            is_global_feas = true
        end
    end
    if isl(lb_cost, lb_best_cost)
        lb_best_cost = lb_cost
    end
    if isl(ub_cost, ub_best_cost)
        ub_best_cost = ub_cost
    end

    if is_feas
        cache.start_sol = sol
        cache.start_costs = costs
    end

    @info "best cost:$(round(best_cost, digits = 2))"

    return best_cost, is_global_feas, lb_best_cost, ub_best_cost
end

function update_cache_sol_ub_feas!(inst::Instance, cache::Cache)
    cache.sol_ub.feas_insert = cache.sol_ub.insert
    return nothing
end

"""
    comp_bs_time_limit(inst::Instance, 
                       params::Parameters, 
                       elapsed_time::Float64)

Compute beam search' time limit based on the elapsed time and progressive 
hedging's time limit.
"""
function comp_bs_time_limit(inst::Instance, 
                            params::Parameters, 
                            elapsed_time::Float64)
    den = 1.0
    workers_num_th = params.progressive_hedging.num_threads - 1
    # When the number of threads is less thant the number of scenarios, the time
    # limit has to be divided so that all scenarios can be solved
    if workers_num_th < inst.num_scenarios
        den = ceil(Int64, inst.num_scenarios / workers_num_th)
    end

    tl = max(params.progressive_hedging.time_limit - elapsed_time, 0.0)
    return min(tl / den, params.beam_search.time_limit)
end

function comp_quality_deviation(inst::Instance, cache::Cache)
    vmin = inst.costs' *  cache.sep_rho_x_min
    if iseq(vmin, 0.0)
        return 100.0
    end

    return 100.0 * (inst.costs' *  cache.sep_rho_x_max) / vmin
end

function comp_hash_values(inst::Instance, cache::Cache)
    hash = Vector{Float64}(undef, inst.num_K)
    for k in cache.sol_ub.insert
        i = inst.key_to_index[k]
        hash[i] = sum(cache.hash_weights[scen] * cache.scenarios[scen].omega[i] 
                        for scen in eachindex(inst.scenarios))
    end

    return hash
end

function update_cache_detect_cycles!(inst::Instance, cache::Cache)
    hash_values = comp_hash_values(inst, cache)
    for k in cache.sol_ub.insert
        i = inst.key_to_index[k]
        # TODO: Modificar tipo de cache.fixed_x_variables para Set
        if iseq(cache.hash_values[i], hash_values[i]) 
            if iseq(cache.sep_rho_x_max[i], 1.0) && 
                !(k in cache.fixed_x_variables)
                if cache.count_cycle_it[i] == 5
                    push!(cache.fixed_x_variables, k)
                    cache.count_cycle_it[i] = 0
                else
                    cache.count_cycle_it[i] += 1
                end
            else
                cache.count_cycle_it[i] = 0
            end
        end
    end
    cache.hash_values = hash_values

    return nothing
end

"""
select_best_warm_start!(inst::Instance, params::Parameters, 
                        cache::WorkerCache, lp::LPModel, 
                        lb_cost::Float64, lb_viol::Float64, 
                        ub_cost::Float64, ub_viol::Float64)

Select the best between the lower bound and upper bound solutions based on their
associated costs and violations.
"""
function select_best_warm_start!(inst::Instance, params::Parameters, 
                                 cache::WorkerCache, lp::LPModel, 
                                 lb_cost::Float64, lb_viol::Float64, 
                                 ub_cost::Float64, ub_viol::Float64)
    cost = const_infinite
    inserted =  Set{CandType}()
    viol = 0.0
    lb_count_use = 0
    ub_count_use = 0

    # 0 0 smallest cost sol
    # 0 1 ub sol
    # 1 0 lb sol
    # 1 1 smallest cost sol
    lb_is_feas = iseq(lb_viol, 0.0)
    ub_is_feas = iseq(ub_viol, 0.0)
    if (lb_is_feas && !ub_is_feas) || 
            (lb_is_feas == ub_is_feas && isl(lb_cost, ub_cost))
        cost = lb_cost
        inserted = copy(cache.sol_lb)
        viol = lb_viol
        lb_count_use = 1
        @info "selected:lb"
    else
        cost = ub_cost
        inserted = copy(cache.sol_ub)
        viol = ub_viol
        ub_count_use = 1
        @info "selected:ub"
        if iseq(ub_viol, 0.0)
            update_lp!(inst, params, lp, cache.sol_ub)
        end
    end

    return cost, inserted, viol, lb_count_use, ub_count_use
end

"""
    repair!(inst::Instance, params::Parameters, cache::WorkerCache, 
            scen::Int64, lp::LPModel, inserted::Set{CandType})

Repair a solution with the binary search repair operator, if possible, and 
return the new solution; otherwise, get the last feasible upper bound solution.
"""
function repair!(inst::Instance, params::Parameters, cache::WorkerCache, 
                 scen::Int64, lp::LPModel, inserted::Set{CandType})
    unfix_s_vars!(lp)

    update_lp!(inst, params, lp, inserted)
    viol = comp_viol(lp)

    # TODO COnverter todos para inteiro como tava antes
    repair_count = 0
    repair_success = 0
    reinsert_st = (0, 0, 0)

    reinsert = Set{CandType}()
    if isg(viol, 0.0)
        repair_count = 1
        v = viol

        reinsert = setdiff(Set{CandType}(keys(inst.K)), inserted)
        removed = copy(reinsert)
        viol, _ = repair!(inst, params, scen, lp, removed, viol)

        if iseq(viol, 0.0)
            repair_success = 1
            @info "repaired: $v -> 0.0"
            setdiff!(reinsert, removed)
            # @info "reinserted: $v -> $viol"
        else
            # v = viol
            # viol, reinsert_st = 
            #         reinsert!(inst, params, scen, lp, inserted, removed, viol)
            # setdiff!(reinsert, removed)
            @info "viol:$v -> $viol"
            # If still infeasible, get the last feasible ub solution
            # reinsert = setdiff(cache.inserted, inserted)
            setdiff!(reinsert, removed)
        end
    end
    repair_st = (repair_count, repair_success)
    
    return reinsert, repair_st, reinsert_st
end

function reinsert!(inst::Instance, 
                   params::Parameters, 
                   scen::Int64, 
                   lp::LPModel, 
                   inserted::Set{CandType}, 
                   removed::Set{CandType}, 
                   best_viol::Float64)
    reinsert = Set{CandType}()
    b = comp_candidates_per_batch_mult(inst, params, removed) * length(removed)
    b = max(floor(Int64, b), 1)
    @info "b:$b"
    rm_vec = collect(removed)
    samples = disjoint_samples(rm_vec, b, true)

    # Sort samples in non-descending order of sum of costs
    sort!(samples, by = x -> sum([inst.K[k].cost for k in x]), rev = false)

    it = 0
    for s in samples
        if iseq(best_viol, 0.0)
            break
        end

        new_inserted = union(inserted, s)
        update_lp!(inst, params, lp, new_inserted)
        viol = comp_viol(lp)
        @info it viol, best_viol

        if isl(viol, best_viol)
            best_viol = viol
            union!(reinsert, s)
            inserted = new_inserted
        end
        it += 1
    end

    l = length(rm_vec)
    setdiff!(removed, reinsert)
    # best_viol = repair!(inst, params, scen, lp, removed, best_viol)
    @info "reinforce it:$it removed:$l -> $(length(removed))"

    num_applied = 1
    num_successes = iseq(best_viol, 0.0) ? 1 : 0

    return best_viol, (num_applied, num_successes, it)
end

function jqm_repair_sols!(inst::Instance, 
                          params::Parameters, 
                          cache::Cache, it::Int64, 
                          controller, 
                          start_time::Float64)
    num_lb_cands = length(cache.sol_lb.insert)
    num_ub_cands = length(cache.sol_ub.insert)
    for scen in 1:inst.num_scenarios
        tl = comp_bs_time_limit(inst, params, time() - start_time)
        
        wcache = WorkerCache(cache, repair_sols)
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
                update_cache_repair!(cache, msg)
            end
        end
    end
    @info "size total:$(inst.num_K) " * 
            "lb:$(num_lb_cands) -> $(length(cache.sol_lb.insert)) " * 
            "ub:$(num_ub_cands) -> $(length(cache.sol_ub.insert))"

    return nothing
end

function update_cache_repair!(cache::Cache, msg::WorkerMessage)
    cache.status[msg.scen].repair = msg.status.repair
    cache.status[msg.scen].reinsert = msg.status.reinsert

    union!(cache.sol_lb.insert, msg.sol_info_lb.reinsert)
    union!(cache.sol_ub.insert, msg.sol_info_ub.reinsert)

    return nothing
end

function jqm_comp_costs!(inst::Instance, 
                         params::Parameters, 
                         cache::Cache, it::Int64, 
                         controller, 
                         start_time::Float64)
    for scen in 1:inst.num_scenarios
        tl = comp_bs_time_limit(inst, params, time() - start_time)
        
        wcache = WorkerCache(cache, comp_g_costs)
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
                update_cache_sol_costs_and_viols!(cache, msg)
            end
        end
    end

    return nothing
end

function comp_costs_and_viol!(inst::Instance, 
                              params::Parameters,
                              scen::Int64, 
                              lp::LPModel, 
                              cache::WorkerCache, 
                              inserted::Set{CandType})
    update_lp!(inst, params, lp, inserted)
    cost, g_cost = comp_penalized_cost(inst, params, scen, lp, cache, inserted)

    return cost, g_cost, comp_viol(lp)
end

"""
    comp_sol_info_lb_ub!(inst::Instance, 
                         params::Parameters,
                         scen::Int64, 
                         lp::LPModel, 
                         cache::WorkerCache)

Compute generation costs of both lower and upper bound solutions.
"""
function comp_sol_info_lb_ub!(inst::Instance, 
                              params::Parameters,
                              scen::Int64, 
                              lp::LPModel, 
                              cache::WorkerCache)
    unfix_s_vars!(lp)

    cost, g_cost, viol = comp_costs_and_viol!(inst, params, scen, lp, 
                                                cache, cache.sol_lb)
    sol_info_lb = SolutionInfo(cost, g_cost, viol, Set{CandType}())
    isg(viol, 0.0) && @info "lb viol:$viol"
    
    cost, g_cost, viol = comp_costs_and_viol!(inst, params, scen, lp, 
                                                cache, cache.sol_ub)
    sol_info_ub = SolutionInfo(cost, g_cost, viol, Set{CandType}())
    isg(viol, 0.0) && @info "ub viol:$viol"

    return sol_info_lb, sol_info_ub
end

# --------------------------------- Parallel PH --------------------------------
function has_finished_all_jobs(controller::JobQueueMPI.Controller)
    return JobQueueMPI.is_job_queue_empty(controller) && 
           !JobQueueMPI.any_pending_jobs(controller)
end

function update_cache_incumbent!(cache::Cache, msg::WorkerMessage)
    cache.scenarios[msg.scen].state = msg.state_values
    # cache.sol_lb.count_use += msg.sol_info_lb.count_use
    # cache.sol_ub.count_use += msg.sol_info_ub.count_use

    return nothing
end

function update_cache_sol_average_perturb!(cache::Cache)

end