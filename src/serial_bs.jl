function run_serial_bs!(inst::Instance, 
                        params::Parameters, 
                        scen::Int64, 
                        lp::LPModel, 
                        cache::WorkerCache, 
                        inserted::Set{CandType}, 
                        removed::Set{CandType}, 
                        cost::Float64, 
                        start_time::Float64)
    params.binary_search.time_limit = 
                        params.beam_search.time_limit - (time() - start_time)
    num_ins_start = length(inserted)

    # Compute BS improvement wrt the initial cost
    start_cost = cost

    cost = binary_search!(inst, params, scen, lp, cache, 
                            inserted, removed, cost, start_time)
    bin_rm_rat = (num_ins_start - length(inserted)) / inst.num_K

    # params.solver.log_level = 2
    # lp = build_lp(inst, params, scen)
    # @warn "aqui"

    # readline()
    
    # Update time limit
    time_limit = params.beam_search.time_limit - (time() - start_time)

    fix_s_vars!(lp)

    num_ins_start = length(inserted)
    num_ins = length(inserted)
    best_cost = cost

    cache_in, cache_rm = init_cache_in_rm(inst)
    # Set the inserted lines according to cache_in and cache_rm
    update_lp!(inst, params, lp, cache_in, false)

    cands_per_batch_m = comp_candidates_per_batch_mult(inst, params, inserted)

    has_reached_tl = false
    for bs_it in 1:params.beam_search.num_max_it

        it = 1
        num_it_wo_impr = 0
        root = Node(cost, 0.0, inserted, removed, Set{CandType}())
        Q = [root]
        while true # Evaluate levels in the tree
            # println("bs level $it") 

            UB = Vector{Node}()
            has_impr = false
            for (i, node) in enumerate(Q) # Evaluate nodes in the same level
                batches = 
                    select_batches!(inst, params, lp, node, cands_per_batch_m)
                for lines in batches
                    el = time() - start_time
                    if isg(el, time_limit)
                        has_reached_tl = true
                        break
                    end
                    in_cands = setdiff(node.inserted, lines)
                    # Set max time limit according to elapsed time
                    # set_attribute(lp.jump_model, "TimeLimit", time_limit - el)
                    set_time_limit!(params, lp, start_time, 
                                    params.beam_search.time_limit)
                    update_lp!(inst, params, lp, cache_in, cache_rm, in_cands)

                    cost = const_infinite
                    is_feas = false
                    viol = 0.0
                    if JuMP.termination_status(lp.jump_model) == MOI.OPTIMAL
                        cost, _ = comp_penalized_cost(inst, params, scen, 
                                                        lp, cache, in_cands)
                        is_feas = true
                        # The model is either feasible or infeasible as s
                        # variables are fixed at zero
                        # viol = comp_viol(lp)
                    end

                    msg = BSWorkerMessage(i, lines, is_feas, cost, viol)

                    add_node!(params, UB, msg, node)

                    if is_feas && isl(cost, best_cost)
                        # log(params, "Inserted update", true)
                        has_impr = true
                        best_cost = cost
                        inserted = UB[end].inserted
                        num_ins = length(UB[end].inserted)

                        # Log info
                        LoggingExtras.withlevel(Info; 
                                                verbosity = params.log_level) do
                            st = Status("bs level:$it", num_ins_start - num_ins, 
                                        num_ins_start, 
                                        cost, start_cost, start_time)
                            @infov 2 log(st)
                        end
                    end
                end
                if has_reached_tl
                    break
                end
            end
            if has_impr
                num_it_wo_impr = 0
            else
                num_it_wo_impr += 1
            end

            if has_reached_tl
                break
            elseif num_it_wo_impr >= params.beam_search.num_max_it_wo_impr
                st = Status("bs", num_ins_start - num_ins, 
                            num_ins_start, best_cost, start_cost, start_time)
                @info log(st)
                @info "max it wo impr reached"
                break
            elseif length(UB) == 0
                @info "no new nodes to investigate"
                break
            end
            
            Q = select!(params, UB)
            it += 1
        end
        if has_reached_tl
            @info "bs tl reached $(round(time() - start_time, digits = 2))"
            break
        else
            @info "enable shuffle strategy"
            params.beam_search.is_shuffle_en = ~params.beam_search.is_shuffle_en
        end
    end
    # union!(inserted, Set(cache.fixed_x_variables))
    JuMP.set_attribute(lp.jump_model, "TimeLimit", GRB_INFINITY)
    update_lp!(inst, params, lp, cache_in, cache_rm, inserted)

    @info "bs best cost:$best_cost"
    
    bs_rm_rat = (num_ins_start - length(inserted)) / inst.num_K

    return inserted, bin_rm_rat, bs_rm_rat 
end