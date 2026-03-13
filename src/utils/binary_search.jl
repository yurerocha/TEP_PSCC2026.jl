function has_reached_stop(params::Parameters, it::Int64, it_wo_impr::Int64, 
                         num_cands_prev_it::Int64, rm_cands::Set{CandType}, 
                         start_time::Float64, is_log_en::Bool = false)
    c1 = it > params.binary_search.max_it
    c2 = it_wo_impr >= params.binary_search.num_max_it_wo_impr
    c3 = length(rm_cands) == num_cands_prev_it
    c4 = isempty(rm_cands)
    c5 = isg(time() - start_time, params.binary_search.time_limit)

    if is_log_en
        @info "bin $c1 $c2 $c3 $c4 $c5"
    end

    return c1 || c2 || c3 || c4 || c5
end

"""
    sort_by_residual_flows(inst::Instance, 
                           lp::LPModel, 
                           f::Dict{Any, Float64}, 
                           inserted::Set{CandType}, 
                           rev = true)

Sort the inserted lines in non-ascending order of residual flows.
"""
function sort_by_residual_flows(inst::Instance, 
                                lp::LPModel, 
                                f::Dict{Any, Float64}, 
                                inserted::Set{CandType}, 
                                rev = true)
    res_flows = Vector{Tuple{CandType, Float64}}()

    for k in inserted
        delta = inst.K[k].f_bar - abs(f[k])
        r = delta / inst.K[k].f_bar
        push!(res_flows, (k, r))
    end

    # Sort lines in non-ascending order of residual flows
    sort!(res_flows, by = x->x[2], rev = rev)

    return Vector{CandType}([res_flows[i][1] for i in eachindex(res_flows)])
end

"""
    divide_into_rm_in(lines::Vector{CandType}, rm_ratio::Float64)

Divide lines into removal and insertion sets based on the remove ratio.
"""
function divide_into_rm_in(lines::Vector{CandType}, rm_ratio::Float64)
    n = floor(Int64, rm_ratio * length(lines))
    return Set{CandType}(lines[1:n]), Set{CandType}(lines[n + 1:end])
end

"""
    select_with_viol(inst::Instance, 
                     params::Parameters, 
                     s::Dict{Any, Float64}, 
                     lines::Set{CandType})

Select candidate lines whose associated existing line has capacity violation.
"""
function select_with_viol(inst::Instance, 
                          params::Parameters, 
                          s::Dict{Any, Float64}, 
                          lines::Set{CandType})
    viols = Vector{Tuple{Any, Float64}}()
    for k in lines
        j = get_existing_line(k)
        if isg(s[k], 0.0)
            push!(viols, (k, s[k]))
        elseif isg(s[j], 0.0)
            push!(viols, (k, s[j]))
        end
    end
    # Sort lines in non-ascending order of residuals
    sort!(viols, by = x->x[2], rev = true)

    return Set{CandType}([v[1] for v in viols])
end

# ------------------------- Rm and add lines functions -------------------------
"""
    rm_lines!(inst::Instance, 
              params::Parameters, 
              lp::LPModel,  
              lines::Set{CandType}, 
              is_opt_req::Bool = false)

Remove lines from the model by setting the susceptances to a small value.
"""
function rm_lines!(inst::Instance, 
                   params::Parameters, 
                   lp::LPModel,  
                   lines::Set{CandType}, 
                   is_opt_req::Bool = false)
    for k in lines
        if !lp.has_fixed_s_vars && !is_fixed(lp.s[k])
            fix(lp.s[k], 0.0; force = true)
        end

        # set_normalized_coefficient([lp.f_cons[k], lp.f_cons[k]], 
        #                            [lp.theta[k[1][2]], lp.theta[k[1][3]]], 
        #                            [0, 0])
        # if !iseq(normalized_coefficient(lp.f_cons[k], lp.Dtheta[k[1][2:3]]), 0)
        set_normalized_coefficient(lp.f_cons[k], lp.Dtheta[k[1][2:3]], 0)
        # end
        # fix(lp.f[k], 0.0; force = true)
    end

    if is_opt_req
        optimize!(lp.jump_model)
    end

    return nothing
end

"""
    add_lines!(inst::Instance, 
               lp::LPModel, 
               lines::lines::Set{CandType}, 
               is_opt_req::Bool = true)

Insert lines in the model by setting the diagonal terms of the susceptance.
"""
function add_lines!(inst::Instance, 
                    params::Parameters, 
                    lp::LPModel,
                    lines::Set{CandType}, 
                    is_opt_req::Bool = true)
    for k in lines
        # if params.model.is_lp_model_s_var_set_req && is_fixed(lp.s[k])
        # Attention! Beam search requires the s variables to be always fixed at 
        # zero
        if !lp.has_fixed_s_vars && is_fixed(lp.s[k])
            unfix(lp.s[k])
            set_lower_bound(lp.s[k], 0.0)
        end

        # if is_fixed(lp.f[k])
        #     unfix(lp.f[k])
        # end

        # set_normalized_coefficient([lp.f_cons[k], lp.f_cons[k]], 
        #                            [lp.theta[k[1][2]], lp.theta[k[1][3]]], 
        #                            [-inst.K[k].gamma, inst.K[k].gamma])
        # if iseq(normalized_coefficient(lp.f_cons[k], lp.Dtheta[k[1][2:3]]), 
        #         -inst.K[k].gamma)
        set_normalized_coefficient(lp.f_cons[k], 
                                   lp.Dtheta[k[1][2:3]], 
                                   -inst.K[k].gamma)
        # end
        # unfix(lp.f[k])
    end

    if is_opt_req
        optimize!(lp.jump_model)
    end

    return nothing
end

"""
    update_lp!(inst::Instance, 
               params::Parameters, 
               lp::LPModel, 
               inserted::Set{CandType}, 
               is_opt_req::Bool = true)
    
Remove all candidate lines, next insert candidates from a set and reoptimize.
"""
function update_lp!(inst::Instance, 
                    params::Parameters, 
                    lp::LPModel, 
                    inserted::Set{CandType}, 
                    is_opt_req::Bool = true)
    # log(params, "It update $it", true)
    rm_lines!(inst, params, lp, Set{CandType}(keys(inst.K)), false)
    add_lines!(inst, params, lp, inserted, is_opt_req)

    # return termination_status(lp.jump_model) == MOI.OPTIMAL
    return nothing
end

"""
    update_lp!(inst::Instance, 
               params::Parameters, 
               lp::LPModel, 
               cache_rm::Set{CandType}, 
               insert::Set{CandType})
    
Remove all candidate lines, next insert candidates from a set and reoptimize.

For multiple iterations, it can be more efficient by storing previously removed 
candidates lines in an initially empty set.
"""
function update_lp!(inst::Instance, 
                    params::Parameters, 
                    lp::LPModel, 
                    cache_in::Set{CandType}, 
                    cache_rm::Set{CandType}, 
                    insert::Set{CandType})
    # Remove inserted lines out of the new insert set
    rm = setdiff(cache_in, insert)
    rm_lines!(inst, params, lp, rm, false)

    # Insert lines not previously inserted
    ins = setdiff(insert, cache_in)
    add_lines!(inst, params, lp, ins, true)

    # Update the set of removed lines
    empty!(cache_in)
    union!(cache_in, insert)
    union!(cache_rm, rm)
    setdiff!(cache_rm, ins)

    if params.debugging_level == 1
        @assert length(cache_in) + length(cache_rm) == inst.num_K
    end

    return nothing
end