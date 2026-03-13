
"""
    config!(inst::Instance, params::Parameters, scen::Int64, tep::TEPModel)

Configure the solver parameters.
"""
function config!(inst::Instance, params::Parameters, scen::Int64, tep::TEPModel)
    set_attribute(tep.jump_model, "Threads", params.solver.num_threads)
    # set_attribute(tep.jump_model, "DualReductions", 0)
    # TODO: Investigate reason for numerical trouble with barrier algorithm for 
    # deterministic systems
    # if tep isa LPModel
    #     set_attribute(tep.jump_model, "Method", 2)
    #     set_attribute(tep.jump_model, "Crossover", 0)
    # end
    # set_attribute(tep.jump_model, "Method", 6)
    # set_attribute(tep.jump_model, "GURO_PAR_PDHGGPU", 1)
    config_log!(inst, params, scen, tep)

    return nothing
end

function config_log!(inst::Instance, 
                     params::Parameters, 
                     scen::Int64, 
                     tep::TEPModel)
    if params.model.optimizer == Gurobi.Optimizer
        # if tep isa LPModel
        if tep isa LPModel || params.solver.log_level == 0
            JuMP.set_silent(tep.jump_model)
            # set_attribute(tep.jump_model, "OutputFlag", 0)
        elseif params.solver.log_level == 1 || params.solver.log_level == 3
            set_attribute(tep.jump_model, "LogFile", 
                            get_log_filename(inst, params, scen))
            set_attribute(tep.jump_model, "LogToConsole", 0)
        elseif params.solver.log_level == 2 || params.solver.log_level == 3
            set_attribute(tep.jump_model, "LogToConsole", 1)
        end
    end

    return nothing
end

function add_vars!(inst::Instance, 
                   params::Parameters, 
                   scen::Int64, 
                   mip::MIPModel)
    for j in keys(inst.J)
        mip.f[j] = @variable(mip.jump_model, base_name = "f[$j]")
    end
    for k in keys(inst.K)
        mip.f[k] = @variable(mip.jump_model, base_name = "f[$k]")
        mip.x[k] = @variable(mip.jump_model, binary = true, base_name = "x[$k]")
    end

    for i in inst.I
        mip.theta[i] = @variable(mip.jump_model, base_name = "theta[$i]")
    end

    return nothing
end

function add_vars!(inst::Instance, params::Parameters, scen::Int64, lp::LPModel)
    for j in keys(inst.J)
        lp.f[j] = @variable(lp.jump_model, base_name = "f[$j]")
        # if params.model.is_lp_model_s_var_set_req
            lp.s[j] = @variable(lp.jump_model, 
                                lower_bound = 0.0, 
                                base_name = "s[$j]")
        # end
    end
    for k in keys(inst.K)
        lp.f[k] = @variable(lp.jump_model, base_name = "f[$k]")
        # if params.model.is_lp_model_s_var_set_req
            lp.s[k] = @variable(lp.jump_model, 
                                lower_bound = 0.0, 
                                base_name = "s[$k]")
        # end
    end

    for i in inst.I
        lp.theta[i] = @variable(lp.jump_model, base_name = "theta[$i]")
    end

    return nothing
end

function add_Dtheta_vars_cons!(inst::Instance, tep::TEPModel)
    for j in keys(inst.J)
        br = j[2:3]
        if !haskey(tep.Dtheta, br)
            tep.Dtheta[br] = @variable(tep.jump_model, 
                                       base_name = "Dtheta[$br]")
            @constraint(tep.jump_model, 
                        tep.Dtheta[br] == tep.theta[j[2]] - tep.theta[j[3]])
        end
    end
end

"""
    add_g_vars!(inst::Instance, 
                params::Parameters, 
                scen::Int64, 
                tep::TEPModel)

Add generation variables.
"""
function add_g_vars!(inst::Instance, 
                     params::Parameters, 
                     scen::Int64, 
                     tep::TEPModel)
    # Some buses may not have load or generation
    for k in keys(inst.scenarios[scen].G)
        bus = inst.scenarios[scen].G[k].bus
        lb = inst.scenarios[scen].G[k].lower_bound
        ub = inst.scenarios[scen].G[k].upper_bound

        if isl(ub, 0.0)
            log(params, "Negative g bounds: $lb, $ub, $k", true)
        end

        tep.g[k] = @variable(tep.jump_model, 
                             lower_bound = lb, 
                             upper_bound = ub, 
                             base_name = "g[$k]")
        if bus in keys(tep.g_bus)
            tep.g_bus[bus] += tep.g[k]
        else
            tep.g_bus[bus] = tep.g[k]
        end
    end

    return nothing
end

"""
    update_g_vars!(inst::Instance, scen::Int64, tep::TEPModel)

Update the model g variables according to the scenario.
"""
function update_g_vars!(inst::Instance, scen::Int64, tep::TEPModel)
    for k in keys(inst.scenarios[scen].G)
        lb = inst.scenarios[scen].G[k].lower_bound
        ub = inst.scenarios[scen].G[k].upper_bound
        JuMP.set_lower_bound(tep.g[k], lb)
        JuMP.set_upper_bound(tep.g[k], ub)
    end
end

"""
    add_symmetry_cons!(inst::Instance, params::Parameters, mip::MIPModel)

Add constraints to help breaking symmetry.

For candidates k, k + 1, ..., k + n with the same gamma, x_k >= x_{k + 1} >= ...
>= x_{k + n}.
"""
function add_symmetry_cons!(inst::Instance, params::Parameters, mip::MIPModel)
    for j in keys(inst.J)
        for l in 1:params.instance.num_candidates
            k = (j, l)
            k_plus = (j, l + 1)
            if k in keys(inst.K) && 
               iseq(inst.K[k].gamma, inst.K[k_plus].gamma)
                @constraint(mip.jump_model, 
                            mip.x[k] >= mip.x[k_plus], 
                            base_name = "sym[$k,$kplus]")
            end
        end
    end

    return nothing
end

function add_symmetry_cons!(inst::Instance, params::Parameters, lp::LPModel)
    return nothing
end

"""
    add_thermal_limits_cons!(inst::Instance, tep::TEPModel)

Add flow variables to the model and thermal limits constraints.
"""
function add_thermal_limits_cons!(inst::Instance, 
                                  params::Parameters, 
                                  tep::TEPModel)
    # Existing lines
    for j in keys(inst.J)
        JuMP.set_lower_bound(tep.f[j], -inst.J[j].f_bar)
        JuMP.set_upper_bound(tep.f[j],  inst.J[j].f_bar)
    end
    # Candidates
    for k in keys(inst.K)
        @constraint(tep.jump_model, -tep.f[k] <= inst.K[k].f_bar * tep.x[k])
        @constraint(tep.jump_model,  tep.f[k] <= inst.K[k].f_bar * tep.x[k])
    end

    return nothing
end

"""
    add_thermal_limits_cons!(inst::Instance, lp::LPModel)

Add flow variables to the model and thermal limits constraints considering the
slack variables.
"""
function add_thermal_limits_cons!(inst::Instance, 
                                  params::Parameters, 
                                  lp::LPModel)
    # Existing lines
    for j in keys(inst.J)
        # if params.model.is_lp_model_s_var_set_req
            @constraint(lp.jump_model, -lp.f[j] <= lp.s[j] + inst.J[j].f_bar)
            @constraint(lp.jump_model,  lp.f[j] <= lp.s[j] + inst.J[j].f_bar)
        # else
        #     JuMP.set_lower_bound(lp.f[j], -inst.J[j].f_bar)
        #     JuMP.set_upper_bound(lp.f[j],  inst.J[j].f_bar)
        # end
    end
    # Candidates
    for k in keys(inst.K)
        # if params.model.is_lp_model_s_var_set_req
            @constraint(lp.jump_model, -lp.f[k] <= lp.s[k] + inst.K[k].f_bar)
            @constraint(lp.jump_model,  lp.f[k] <= lp.s[k] + inst.K[k].f_bar)
        # else
        #     JuMP.set_lower_bound(lp.f[k], -inst.K[k].f_bar)
        #     JuMP.set_upper_bound(lp.f[k],  inst.K[k].f_bar)
        # end
    end

    return nothing
end

function add_ohms_law_cons!(inst::Instance, tep::TEPModel)
    # Ohm's law for existing circuits
    for j in keys(inst.J)
        # Dtheta = tep.theta[j[2]] - tep.theta[j[3]]
        @constraint(tep.jump_model, 
                    tep.f[j] == inst.J[j].gamma * tep.Dtheta[j[2:3]],
                    # tep.f[j] == inst.J[j].gamma * Dtheta,
                    base_name = "ol[$j]")
    end
    # Ohm's law for candidate circuits
    for k in keys(inst.K)
        # Dtheta = tep.theta[k[1][2]] - tep.theta[k[1][3]]
        # y = tep.f[k] - inst.K[k].gamma * Dtheta
        y = tep.f[k] - inst.K[k].gamma * tep.Dtheta[k[1][2:3]]

        bigM = comp_bigM(inst, k)

        @constraint(tep.jump_model, 
                    -y <= bigM * (1.0 - tep.x[k]), 
                    base_name = "ol.lb[$k]")
        @constraint(tep.jump_model, 
                    y <= bigM * (1.0 - tep.x[k]), 
                    base_name = "ol.ub[$k]")
    end

    return nothing
end

function add_ohms_law_cons!(inst::Instance, lp::LPModel)
    # Ohm's law for existing circuits
    for j in keys(inst.J)
        # Dtheta = lp.theta[j[2]] - lp.theta[j[3]]
        lp.f_cons[j] = @constraint(lp.jump_model, 
                                lp.f[j] == inst.J[j].gamma * lp.Dtheta[j[2:3]],
                                # lp.f[j] == inst.J[j].gamma * Dtheta,
                                base_name = "ol.j[$j]")
    end
    # Ohm's law for candidate circuits
    for k in keys(inst.K)
        # Dtheta = lp.theta[k[1][2]] - lp.theta[k[1][3]]
        lp.f_cons[k] = @constraint(lp.jump_model, 
                        lp.f[k] == inst.K[k].gamma * lp.Dtheta[k[1][2:3]],
                        # lp.f[k] == inst.K[k].gamma * Dtheta,
                        base_name = "ol.k[$k]")
    end

    return nothing
end

function add_Dtheta_bounds_cons!(inst::Instance, mip::MIPModel)
    bounds = Dict{Tuple{Int64, Int64}, Tuple{Float64, Float64}}()
    for j in keys(inst.J)
        if iseq(inst.J[j].Dtheta_bounds[1], 0.0) &&
           iseq(inst.J[j].Dtheta_bounds[2], 0.0)
            continue
        end
        k = (j[2], j[3])
        b = inst.J[j].Dtheta_bounds
        if k in keys(bounds)
            lb = max(bounds[k][1], b[1])
            ub = min(bounds[k][2], b[2])
            bounds[k] = (lb, ub)
        else
            bounds[k] = b
        end
    end
    for (k, b) in bounds
        Dtheta = mip.theta[k[1]] - mip.theta[k[2]]
        @constraint(mip.jump_model, Dtheta >= b[1])
        @constraint(mip.jump_model, Dtheta <= b[2])
    end

    return nothing
end

function add_Dtheta_bounds_cons!(inst::Instance, lp::LPModel)
    return nothing
end

"""
    add_fkl_cons!(inst::Instance, 
                  scen::Int64, 
                  tep::TEPModel)

Add first Kirchhoff law constraints.
"""
function add_fkl_cons!(inst::Instance, 
                       scen::Int64, 
                       tep::TEPModel)
    for i in inst.I
        e = comp_candidate_incident_flows(inst, tep, i)
        # e += comp_existing_incident_flows(inst, f, i)
        add_to_expression!(e, comp_existing_incident_flows(inst, tep, i))
        
        # Some buses may not have load or generation
        d = i in keys(inst.scenarios[scen].D) ? inst.scenarios[scen].D[i] : 0.0
        g = i in keys(tep.g_bus) ? tep.g_bus[i] : AffExpr(0.0)

        cons = @constraint(tep.jump_model, e + g == d, base_name = "fkl[$i]")

        # if tep isa MIPModel
            tep.fkl_cons[i] = cons
        # end
    end

    return nothing
end

"""
    update_fkl_cons_rhs!(inst::Instance, 
                         scen::Int64, 
                         tep::TEPModel)

Update the right-hand side coefficients of the first Kirchhoff law constraints.
"""
function update_fkl_cons_rhs!(inst::Instance, 
                              scen::Int64, 
                              tep::TEPModel)
    for i in inst.I
        d = i in keys(inst.scenarios[scen].D) ? inst.scenarios[scen].D[i] : 0.0
        JuMP.set_normalized_rhs(tep.fkl_cons[i], d)
    end
end

"""
    add_ref_bus_cons!(inst::Instance, tep::TEPModel)

Add constraint on the reference bus theta value.
"""
function add_ref_bus_cons!(inst::Instance, tep::TEPModel)
    @constraint(tep.jump_model, tep.theta[inst.ref_bus] == 0.0)
    
    return nothing
end

"""
    set_obj!(inst::Instance, 
             params::Parameters, 
             scen::Int64, 
             is_build_obj_req::Bool, 
             tep::TEPModel)

Set the objective to minimize the costs of expanding the network.
"""
function set_obj!(inst::Instance, 
                  params::Parameters, 
                  scen::Int64, 
                  is_build_obj_req::Bool, 
                  tep::TEPModel)
    # Cost of building new candidate lines
    if (is_build_obj_req)
        for k in keys(inst.K)
            add_to_expression!(tep.obj, inst.K[k].cost, tep.x[k])
        end
    end
    # Generation costs
    for k in keys(tep.g)
        add_to_expression!(tep.obj, 
                comp_g_obj(params, tep.g[k], inst.scenarios[scen].G[k].costs))
    end

    @objective(tep.jump_model, Min, tep.obj)
    
    return nothing
end

"""
    set_obj!(inst::Instance, 
             params::Parameters, 
             scen::Int64, 
             lp::LPModel)

Set the objective to minimize the slacks.
"""
function set_obj!(inst::Instance, 
                  params::Parameters, 
                  scen::Int64, 
                  is_build_obj_req::Bool, 
                  lp::LPModel)
    # Generation costs
    pen = 0.0
    for k in keys(lp.g)
        g_info = inst.scenarios[scen].G[k]
        add_to_expression!(lp.obj, comp_g_obj(params, lp.g[k], g_info.costs))
        # pen = max(pen, comp_largest_g(params, c, 
        #                               inst.scenarios[scen].G[k].upper_bound))
        pen = max(pen, comp_max_g_obj_value(params, g_info))
    end
    pen *= params.model.penalty
    # if params.model.is_lp_model_s_var_set_req
        add_to_expression!(lp.obj, sum([pen * s for s in values(lp.s)]))
    # end

    @objective(lp.jump_model, Min, lp.obj)
    
    return nothing
end

"""
   comp_g_obj(params::Parameters, 
                   g::JuMP.VariableRef, 
                   costs::Vector{Float64})

Compute g in the objective function.
"""
function comp_g_obj(params::Parameters, 
                         g::T, 
                         costs::Vector{Float64}) where 
                                        T <: Union{Float64, JuMP.VariableRef}
    l = length(costs)
    if l == 0
        return 0.0
    elseif l == 1
        return costs[1]
    elseif l == 2 || !params.model.is_dcp_power_model_en
        # When is_dcp_power_model_en == false, the objective function is, at 
        # most, linear in terms of g.
        return costs[1] + costs[2]*g
    else
        if l > 3
            log(params, 
                "Length of costs > 3, but quadratic polynomial considered", 
                true)
        end
        return costs[1] + costs[2]*g + costs[3]*g^2
    end
end

function comp_max_g_obj_value(params::Parameters, g_info::GeneratorInfo)
    l = length(g_info.costs)
    if l == 0
        return 0.0
    elseif l == 1
        return g_info.costs[1]
    elseif l == 2 || !params.model.is_dcp_power_model_en
        # When is_dcp_power_model_en == false, the objective function is, at 
        # most, linear in terms of g.
        return max(g_info.costs[1], g_info.costs[2]*g_info.upper_bound)
    else
        if l > 3
            log(params, 
                "Length of costs > 3, but quadratic polynomial considered", 
                true)
        end
        return maximum([g_info.costs[1], 
                        g_info.costs[2]*g_info.upper_bound, 
                        g_info.costs[3]*g_info.upper_bound^2])
    end
end

# """
#     comp_largest_g(params::Parameters, 
#                    costs::Vector{Float64}, 
#                    upper_bound::Float64)
# Compute the largest value between C^{mai} and C^{ope} * G^{max}_g
# """
# function comp_largest_g(params::Parameters, 
#                         costs::Vector{Float64}, 
#                         upper_bound::Float64)
    
#     if length(costs) == 1
#         return costs[1]
#     elseif length(costs) == 2
#         @info costs[1], costs[2], upper_bound, costs[2]*upper_bound
#         return max(costs[1], costs[2]*upper_bound)
#     else
#         return 0.0
#     end
# end

"""
    update_model!(inst::Instance, 
                  params::Parameters, 
                  cache::WorkerCache, 
                  scen::Int64, 
                  tep::TEPModel)

Update the model g variables and the right-hand side coefficients of the first 
Kirchhoff law constraints, according to scenario.
"""
function update_model!(inst::Instance, 
                       params::Parameters, 
                       cache::WorkerCache, 
                       scen::Int64, 
                       tep::TEPModel)
    update_g_vars!(inst, scen, tep)
    update_fkl_cons_rhs!(inst, scen, tep)

    config_log!(inst, params, scen, tep)

    if params.progressive_hedging.is_en && tep isa MIPModel
        update_model_obj!(params, cache, scen, tep)
    end

    return nothing
end

"""
    print_cons(model::JuMP.Model, filename::String)

Print model constraints according to their types.
"""
function print_cons(model::JuMP.Model, filename::String)
    open(filename, "w") do f
        for (F, S) in JuMP.list_of_constraint_types(model)
            for cref in JuMP.all_constraints(model, F, S)
                println(f, cref)
            end
        end
    end
end

"""
    enforce_sol(inst::Instance, 
                tep::TEPModel, 
                sol::Dict{String, Any})

Enforce solution as constraints within the model.
"""
function enforce_sol(inst::Instance, 
                     tep::TEPModel, 
                     sol::Dict{String, Any})
    # for b in sol["bus"]
    #     i = parse(Int64, b[1])
    #     @constraint(tep.jump_model, tep.theta[i] == -b[2]["va"])
    # end
    for g in sol["gen"]
        i = parse(Int64, g[1])
        @constraint(tep.jump_model, tep.g[i] == g[2]["pg"])
    end

    return nothing
end

"""
    comp_viol(lp::LPModel)

Compute violation as the sum of values of slack variables.
"""
function comp_viol(lp::LPModel)
    if JuMP.termination_status(lp.jump_model) == MOI.OPTIMAL
        sum(JuMP.value(s) for s in values(lp.s))
    else
        return const_infinite
    end
end

"""
    comp_viol_and_max(lp::LPModel)

Compute violation as the sum of the slack variables, and the maximum violation.
"""
function comp_viol_and_max(lp::LPModel)
    if JuMP.termination_status(lp.jump_model) == MOI.OPTIMAL
        vals = [JuMP.value(s) for s in values(lp.s)]
        return sum(vals), maximum(vals)
    else
        return const_infinite, const_infinite
    end
end

function check_sol(inst::Instance, tep::TEPModel, md)
    is_feas = true
    filename = "check_sol.txt"

    x = get_values(tep.x)
    f = get_values(tep.f)
    th = get_values(tep.theta)

    open(filename, "w") do file
        for k in keys(inst.J)
            br_info = inst.J[k]
            if isg(abs(f[k]), br_info.f_bar)
                println(file, "F viol k:$k f:$(f[k]) f_bar:$(br_info.f_bar)")
                is_feas = false
            end
        end

        for k in keys(inst.K)
            br_info = inst.K[k]
            if iseq(x[k], 1.0)
                if isg(abs(f[k]), br_info.f_bar)
                    println(file, 
                            "F viol k:$k f:$(f[k]) f_bar:$(br_info.f_bar)")
                    is_feas = false
                end

                fr = k[1][2]
                to = k[1][3]
                if !iseq(f[k], br_info.gamma * (th[fr] - th[to]))
                    println(file, 
                            "T error k:$k f:$(f[k]) Gamma: $(br_info.gamma) " * 
                            " t_fr:$(th[fr]) t_to:$(th[to])")
                    is_feas = false
                end
            elseif !iseq(f[k], 0.0)
                println(file, "F viol k:$k f:$(f[k]) f_bar:0")
                is_feas = false
            end
        end

        # Check if a constraint is infeasible based on its attribute slack
        for i in inst.I
            a = MOI.get(md, 
                        Gurobi.ConstraintAttribute("Slack"), 
                        index(tep.fkl_cons[i]))
            s = JuMP.value(a)
            if isg(s, 0.0)
                println(file, "C viol $i")
                is_feas = false
            end
        end

        g = get_values(tep.g)
        for k in keys(g)
            g_info = inst.scenarios[1].G[k]
            lb = g_info.lower_bound
            ub = g_info.upper_bound
            if isl(g[k], lb) || isg(g[k], ub)
                println(file, "G viol k:$k g:$(g[k]) lb:$(lb) ub:$(ub)")
                is_feas = false
            end
        end
    end

    return is_feas
end

"""
    fix_s_vars!(lp::LPModel)

Fix the s variables.
"""
function fix_s_vars!(lp::LPModel)
    lp.has_fixed_s_vars = true
    for k in keys(lp.s)
        if !JuMP.is_fixed(lp.s[k])
            JuMP.fix(lp.s[k], 0.0; force = true)
        end
    end
    return nothing
end

"""
    unfix_s_vars!(lp::LPModel)

Unfix and set the lower bounds of the s variables.
"""
function unfix_s_vars!(lp::LPModel)
    lp.has_fixed_s_vars = false
    for k in keys(lp.s)
        if JuMP.is_fixed(lp.s[k])
            JuMP.unfix(lp.s[k])
            JuMP.set_lower_bound(lp.s[k], 0.0)
        end
    end
    return nothing
end

"""
    fix_start!(inst::Instance, 
               params::Parameters, 
               mip::MIPModel, 
               scen::Int64, 
               start::Start)

Fix start the model with the lines, generation, and flows of the start struct.
"""
function fix_start!(inst::Instance, 
                    params::Parameters, 
                    scen::Int64, 
                    mip::MIPModel, 
                    inserted::Set{CandType})
    JuMP.set_attribute(mip.jump_model, "SolutionLimit", 1)
    # JuMP.set_attribute(mip.jump_model, 
    #                    MOI.RawOptimizerAttribute("FeasibilityTol"), 
    #                    1e-3)

    if params.model.is_symmetry_en
        fix_for_symmetry_contrs(inst, params, mip, start)
    end

    for k in keys(inst.K)
        JuMP.fix(mip.x[k], 0.0; force = true)
    end
    for k in inserted
        JuMP.fix(mip.x[k], 1.0; force = true)
    end
    # all_keys = vcat(collect(keys(inst.J)), collect(keys(inst.K)))
    # for l in all_keys
    #     JuMP.fix(mip.f[l], start.f[l], force = true)
    # end
    # for k in keys(inst.scenarios[scen].G)
    #     JuMP.fix(mip.g[k], start.g[k]; force = true)
    #     # fix(md.theta[i], start.theta[i])
    # end

    optimize!(mip.jump_model)
    # To compare the objective value of this solution with the objective value
    # of the mip start solution in the next solver call, the fix constraints on
    # the flows and generation have to be removed

    model = mip.jump_model
    status = JuMP.termination_status(model)
    cost = const_infinite
    is_feas = true
    if status == MOI.OPTIMAL
        cost = JuMP.objective_value(model)
    elseif status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
        @info "infeasible model"
        if params.solver.log_level >= 1
            JuMP.compute_conflict!(model)
            if JuMP.get_attribute(model, MOI.ConflictStatus()) == 
               MOI.CONFLICT_FOUND
                iis_model, _ = copy_conflict(model)
                print(iis_model)
                print_cons(iis_model, "model.iis")
            end
        end
        is_feas = false
    end
    if params.debugging_level == 1
        @assert is_feas "scen#$scen infeasible sol"
    end

    for k in keys(inst.K)
        JuMP.unfix(mip.x[k])
        # set_lower_bound(mip.x[k], 0.0)
        # set_upper_bound(mip.x[k], 1.0)
    end
    # for l in all_keys
    #     JuMP.unfix(mip.f[l])
    # end
    # # Some buses may not have generation
    # for k in keys(inst.scenarios[scen].G)
    #     JuMP.unfix(mip.g[k])
    #     lb = inst.scenarios[scen].G[k].lower_bound
    #     ub = inst.scenarios[scen].G[k].upper_bound
    #     JuMP.set_lower_bound(mip.g[k], lb)
    #     JuMP.set_upper_bound(mip.g[k], ub)
    #     # unfix(md.theta[i])
    # end

    set_attribute(mip.jump_model, "SolutionLimit", MAXINT)

    return is_feas
end

function solve!(inst::Instance, params::Parameters, mip::MIPModel)
    model = mip.jump_model

    if params.model.optimizer == Gurobi.Optimizer
        set_attribute(model, MOI.RawOptimizerAttribute("SolutionLimit"), MAXINT)
        set_attribute(model, 
                      MOI.RawOptimizerAttribute("TimeLimit"), 
                      params.solver.time_limit)
        # terminates the optimization process as soon as pre-processing and the 
        # computation of the root relaxation is finished
        # set_attribute(model, MOI.RawOptimizerAttribute("NodeLimit"), 0)
        # set_attribute(model, MOI.RawOptimizerAttribute("IntFeasTol"), 1e-9)
        # set_attribute(model, MOI.RawOptimizerAttribute("FeasibilityTol"), 1e-9)
        # set_attribute(model, MOI.RawOptimizerAttribute("OptimalityTol"), 1e-9)
    end

    rt_runtime = 0.0
    incumbent_time = params.solver.time_limit
    rt_best_bound = 0.0
    mip_start_gap = Inf64
    has_solved_rt_relaxation = false
    function root_node_callback(cb_data, cb_where::Cint)
        # if cb_where == GRB_CB_MESSAGE
        #     buff = Vector{Cchar}(undef, 100)
        #     GRBcbget(cb_data, cb_where, GRB_CB_MSG_STRING, buff)
        #     buff = [Char(abs(ch)) for ch in buff]
        #     @warn String(buff)
        #     readline()
        # if has_solved_rt_relaxation
        #     runtime = Ref{Cdouble}()
        #     GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, runtime)
        #     # Terminate as soon as pre-processing and the computation of the 
        #     # root relaxation is finished and the beam search time limit is 
        #     # reached
        #     if isg(runtime[], params.beam_search.time_limit)
        #         GRBterminate(backend(model).optimizer.model.inner)
        #     end
        # end

        if cb_where == GRB_CB_MIPNODE
            node_count = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, node_count)
            if iseq(node_count[], 0.0)
                runtime = Ref{Cdouble}()
                root_bound = Ref{Cdouble}()
                GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, runtime)
                GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBND, root_bound)

                rt_runtime = runtime[]
                rt_best_bound = root_bound[]
                has_solved_rt_relaxation = true
                # GRBterminate(backend(model).optimizer.model.inner)
            end
        elseif cb_where == GRB_CB_MIPSOL && 
               iseq(incumbent_time, params.solver.time_limit)
            # Prevents incumbent_time from being updated after the first 
            # incumbent solution is found
            runtime = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, runtime)
            incumbent_time = runtime[]
            # if isl(runtime[], 0.5)
            #     obj_bst = Ref{Cdouble}()
            #     obj_bnd = Ref{Cdouble}()
            #     GRBcbget(model, cb_where, GRB_CB_MIPSOL_OBJBST, obj_bst)
            #     GRBcbget(model, cb_where, GRB_CB_MIPSOL_OBJBND, obj_bnd)
            #     # https://support.gurobi.com/hc/en-us/articles/360025036151-Why-are-there-large-or-increasing-MIP-gap-values
            #     mip_start_gap = (obj_bnd[] - obj_bst[]) / obj_bst[]
            #     @warn obj_bnd[], obj_bst[], mip_start_gap
            #     readline()
            # end
        end
    end
    if params.model.optimizer == Gurobi.Optimizer
        set_attribute(model, Gurobi.CallbackFunction(), root_node_callback)
    end
    
    JuMP.optimize!(model)

    status = JuMP.termination_status(model)
    
    upper_bound = const_infinite
    obj = const_infinite
    gap = const_infinite
    build_obj_rat = const_infinite

    # If the solver found a solution
    if JuMP.has_values(model)
        upper_bound = JuMP.objective_value(model)
        obj = JuMP.objective_value(model)
        if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
            gap = 0.0
        else # if status == MOI.TIME_LIMIT
            gap = 100.0 * JuMP.relative_gap(model)
        end

        # grb_model = JuMP.backend(model)
        # is_feas = check_sol(inst, mip, grb_model)
        # if !is_feas 
        #     status = "MODEL_ERROR"
        # end
    # elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
    #     grb_model = backend(model)
    #     # relaxobjtype = 2: the objective of the feasibility relaxation is to 
    #     # minimize the total number of bound and constraint violations
    #     relaxobjtype = 1
    #     lbpen = Cdouble[10000.0]
    #     ubpen = Cdouble[10000.0]
    #     rhspen = Cdouble[1.0]
    #     feasobjP = Ref{Cdouble}()
    #     GRBfeasrelax(grb_model, relaxobjtype, 1, 
    #                  lbpen, ubpen, rhspen, feasobjP)
    #     # set_attribute(model, 
    #     #               MOI.RawOptimizerAttribute("DualReductions"), 
    #     #               0)
    #     # GRBreset(model, 0)
    #     # optimize!(model)
    #     GRBoptimize(grb_model)
    #     check_sol(inst, mip, grb_model)
    elseif status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
        # grb_model = backend(model)
        # GRBreset(grb_model, 1)
        # GRPoptimize(grb_model)
        # set_attribute(grb_model, 
        #               MOI.RawOptimizerAttribute("DualReductions"), 
        #               0)
        # TODO: Add param to compute conflict when infeasible
        # https://jump.dev/JuMP.jl/stable/manual/solutions/#Conflicts
        compute_conflict!(model)
        if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
            iis_model, _ = copy_conflict(model)
            print_cons(iis_model, "model.iis")
        end
    end

    results = Dict(
        # "incumbent_time" => incumbent_time, 
        # "solve_time" => solve_time(model), 
        # "status" => status, 
        # "root_best_bound" => rt_best_bound, 
        "lb" => upper_bound, 
        "ub" => upper_bound, 
        "best" => obj, 
        "gap" => gap
    )
end

"""
    get_g(inst, model_dt)

Get g values from the model.
"""
function get_g(inst::Instance, model_dt::MIPModel)
    g = Dict{Tuple{Int64, Int64}, Float64}
    for k in keys(inst.G)
        g[i] = value(model_dt.g[i])
    end
    return g
end