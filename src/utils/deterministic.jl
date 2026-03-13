"""
    solve_deterministic(inst::Instance, params::Parameters)

Implementation of the deterministic equivalent stochastic model.
"""
function solve_deterministic!(inst::Instance, params::Parameters)
    log(params, "Deterministic solution strategy", true)

    # Initialization
    cache = Cache(inst.num_scenarios, inst.num_K)

    mip, subproblems = build_deterministic(inst, params)
    
    results = solve!(inst, params, mip)

    if has_values(mip.jump_model)
        for scen in 1:inst.num_scenarios
            update_cache_incumbent!(cache, scen, subproblems[scen])
            x = cache.scenarios[scen].state
            println("Scen#$(scen): $(round(Int64, sum(x)))/$(length(x))")
        end
    end

    return cache, results
end

"""
    add_subproblem!(inst::Instance, 
                    model::JuMP.Model, 
                    scen::Int64, 
                    subproblem::JuMP.Model)

Add a subproblem as a block to the deterministic problem.
"""
function add_subproblem!(inst::Instance, 
                         model::JuMP.Model, 
                         scen::Int64, 
                         subproblem::JuMP.Model)
    # Associate every variable in the subproblem with a variable in the model
    sp_vars = vcat(JuMP.all_variables(subproblem))
    md_vars = @variable(model, [v = 1:length(sp_vars)])
    var_in_md = Dict{JuMP.VariableRef, JuMP.VariableRef}()
    # Naming the new model variables according to their subproblem names
    for (src, dest) in zip(sp_vars, md_vars)
        var_in_md[src] = dest
        name = JuMP.name(src)
        JuMP.set_name(dest, string(name, "_scen#", scen))
    end
    # Build constraints in the model according to the constraints in the 
    # subproblem
    for (F, S) in JuMP.list_of_constraint_types(subproblem)
        for cref in JuMP.all_constraints(subproblem, F, S)
            cobj = JuMP.constraint_object(cref)
            c = build_expr(cobj.func, var_in_md)
            @constraint(model, c in cobj.set)
        end
    end
    # Update the objective function to consider the new variables
    incumbent = JuMP.objective_function(model)
    delta_obj =  build_expr(JuMP.objective_function(subproblem), var_in_md)

    @objective(model, Min, incumbent + inst.scenarios[scen].p * delta_obj)

    # Se ligar nisso na hora de recuperar a função objetivo do det. equivalente
    # Atualizar as variáveis de estado dos subproblemas para que sejam as novas
    # variáveis x associadas ao problema principal
    update_state_vars!(model, subproblem, var_in_md)

    subproblem.ext[:var_in_md] = var_in_md

    return nothing
end

"""
    build_expr(cons::AffQuadExpr, 
               var_in_model::Dict{JuMP.VariableRef, JuMP.VariableRef})

Build a new constraint equal to the one in the subproblem, but involving the 
variables in the model.
"""
function build_expr(cons::T, 
                var_in_model::Dict{JuMP.VariableRef, JuMP.VariableRef}) where 
                                    T <: Union{AffQuadExpr, JuMP.VariableRef}
    if cons isa JuMP.AffExpr
        return JuMP.AffExpr(cons.constant, 
                            OrderedDict(var_in_model[var] => coef 
                                        for (var, coef) in cons.terms))
    elseif cons isa JuMP.QuadExpr
        # Build the aff expr with the vars in the model
        e = AffExpr()
        for (var, coef) in cons.aff.terms
            add_to_expression!(e, coef, var_in_model[var])
        end
        terms = OrderedDict{UnorderedPair{JuMP.VariableRef}, Float64}()
        for (var, coef) in cons.terms
            p = UnorderedPair{JuMP.VariableRef}(var_in_model[var.a], 
                                                var_in_model[var.b])
            terms[p] = coef
        end
        return JuMP.QuadExpr(e, terms)
    else
        return var_in_model[var]
    end
end

"""
    build_expr(var::JuMP.VariableRef, 
               var_in_model::Dict{JuMP.VariableRef, JuMP.VariableRef})

Build a new constraint equal to the one in the subproblem, but involving the 
variables in the model.
"""
function build_expr(var::JuMP.VariableRef, 
                    var_in_model::Dict{JuMP.VariableRef, JuMP.VariableRef})
    return var_in_model[var]
end

"""
    update_state_vars!(model::JuMP.Model, 
                            subproblem::JuMP.Model, 
                            var_in_model::Dict{JuMP.VariableRef, 
                                               JuMP.VariableRef})

Update the state variables of the subproblem with the state variables of the
model.
"""
function update_state_vars!(model::JuMP.Model, 
                            subproblem::JuMP.Model, 
                            var_in_model::Dict{JuMP.VariableRef, 
                                               JuMP.VariableRef})
    x = [build_expr(subproblem.ext[:state][k], 
                    var_in_model) for k in keys(inst.K)]
    subproblem.ext[:state] = x
    return nothing
end

"""
    add_non_anticipativity_cons!(inst::Instance, 
                                 model::JuMP.Model, 
                                 subproblems::Vector{MIPModel})
                                    
Add constraints to make the first-stage decisions for all subproblems to be the
same.
"""
function add_non_anticipativity_cons!(inst::Instance, 
                                      model::JuMP.Model, 
                                      subproblems::Vector{MIPModel})
    for scen in 2:inst.num_scenarios
        @constraint(model, subproblems[1].jump_model.ext[:state] .== 
                           subproblems[scen].jump_model.ext[:state])
    end
    return nothing
end

"""
    add_obj_build_costs!(inst::Instance, 
                         model::JuMP.Model, 
                         subproblems::Vector{MIPModel})

Update the objective function with the costs ob building candidates.
"""
function add_obj_build_costs!(inst::Instance, 
                              model::JuMP.Model, 
                              subproblems::Vector{MIPModel})
    incumbent = JuMP.objective_function(model)

    x = subproblems[1].jump_model.ext[:state]
    build = sum(inst.K[k].cost * x[i] for (i, k) in enumerate(keys(inst.K)))

    @objective(model, Min, build + incumbent)

    return nothing
end

"""
    fix_start!(inst::Instance, 
               mip::MIPModel, 
               subproblems::Vector{MIPModel}, 
               inserted::Set{CandType})

Fix the x variables according to a list and return the corresponding objective 
value.
"""
function fix_start!(inst::Instance, 
                    mip::MIPModel, 
                    subproblems::Vector{MIPModel}, 
                    inserted::Set{CandType})
    if !in(:var_in_md, keys(subproblems[1].jump_model.ext))
        return nothing
    end

    var_in_md = subproblems[1].jump_model.ext[:var_in_md]
    for k in keys(inst.K)
        JuMP.fix(var_in_md[subproblems[1].x[k]], 0.0; force = true)
    end
    for k in inserted
        JuMP.fix(var_in_md[subproblems[1].x[k]], 1.0; force = true)
    end

    JuMP.set_attribute(mip.jump_model, 
                       MOI.RawOptimizerAttribute("SolutionLimit"), 
                       1)

    optimize!(mip.jump_model)
    cost = JuMP.result_count(mip.jump_model) > 0 ? 
                                JuMP.objective_value(mip.jump_model) : nothing

    JuMP.set_attribute(mip.jump_model, 
                       MOI.RawOptimizerAttribute("SolutionLimit"), 
                       MAXINT)

    for k in keys(inst.K)
        JuMP.unfix(var_in_md[subproblems[1].x[k]])
    end

    return cost
end