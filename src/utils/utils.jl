"""
    isl(a::Float64, b::Float64)

Compute if a is less than b.
"""
function isl(a::Float64, b::Float64, eps::Float64 = const_eps)
    return a < b - eps
end

"""
    isg(a::Float64, b::Float64)

Compute if a is greater than b.
"""
function isg(a::Float64, b::Float64, eps::Float64 = const_eps)
    return a > b + eps
end

"""
    iseq(a::Float64, b::Float64)

Compute if a is equal to b.
"""
function iseq(a::Float64, b::Float64, eps::Float64 = const_eps)
    return abs(a - b) < eps
end

"""
    iseq(A::Matrix{T}, B::Matrix{T}) where T <: AbstractFloat

Return true if matrices A and B are equal.
"""
function iseq(A::Matrix{T}, B::Matrix{S}, eps::Float64 = const_eps) where 
                                        {T <: AbstractFloat, S <: AbstractFloat}
    return norm(A - B) < eps
end

function get_num(s::Vector{String}, i::Int64)
    return parse(Int, split(s[i], ":")[2])
end

function get_cand(inst::Instance, l::Int64)
    for i in 1:inst.num_K
        if inst.K[i].fr == l.fr && inst.K[i].to == l.to
            return i
        end
    end
end

"""
Compute incidence matrix with existing and candidate circuits.
"""
function comp_incidence_matrix(inst::Instance, 
                               f::Vector{JuMP.VariableRef}, 
                               i::Int64)
    e = 0
    for j in 1:inst.num_J
        circ = inst.J[j]
        if circ.to == i
            e += f[j]
        elseif circ.fr == i
            e -= f[j]
        end
    end
    for k in 1:inst.num_K
        circ = inst.K[k]
        if circ.to == i
            e += f[inst.num_J + k]
        elseif circ.fr == i
            e -= f[inst.num_J + k]
        end
    end
    return e
end

"""
    comp_incidence_matrix(inst, f, t, i)

Compute incidence matrix with existing and candidate circuits.
"""
function comp_incidence_matrix(inst::Instance,
                               md::JuMP.Model, 
                               f::Vector{JuMP.VariableRef}, 
                               i::Int64, 
                               is_add_cons::Bool,
                               is_cand_en::Bool,
                               lines::Vector{Int64})
    e = 0
    # Candidate lines
    is_cons_update_req = false
    if is_cand_en
        for l in lines
            circ = inst.K[l - inst.num_J]

            # l = inst.num_J + k

            if (circ.to == i || circ.fr == i) && !isdefined(f, l)
                is_cons_update_req = true
                f[l] = @variable(md, base_name = "f")
            end

            if circ.to == i
                e += f[l]
            elseif circ.fr == i
                e -= f[l]
            end
        end
    end

    # Existing lines
    if is_add_cons || is_cons_update_req
        for j in 1:inst.num_J
            circ = inst.J[j]
            if circ.to == i
                e += f[j]
            elseif circ.fr == i
                e -= f[j]
            end
        end
    end

    return e, is_cons_update_req
end

"""
    comp_existing_incident_flows(inst::Instance,
                                 tep::TEPModel, 
                                 i::Int64)

Compute the incident flow at node i for the existing lines.
"""
function comp_existing_incident_flows(inst::Instance, tep::TEPModel, i::Int64)
    e = AffExpr(0.0)
    for j in keys(inst.J)
        if j[2] == i
            add_to_expression!(e, -1.0, tep.f[j])
        elseif j[3] == i
            add_to_expression!(e, 1.0, tep.f[j])
        end
    end
    return e
end

"""
    comp_candidate_incident_flows(inst::Instance, 
                                  tep::TEPModel, 
                                  i::Int64)

Compute the incident flow at node i for the candidate lines.
"""
function comp_candidate_incident_flows(inst::Instance, tep::TEPModel, i::Int64)
    e = AffExpr(0.0)
    for k in keys(inst.K)
        if k[1][2] == i
            add_to_expression!(e, -1.0, tep.f[k])
        elseif k[1][3] == i
            add_to_expression!(e, 1.0, tep.f[k])
        end
    end
    return e
end

"""
    log(file::String, msg::String)

Log message to file.
"""
function log(file::String, msg::String)
    open(file, "a") do f
        write(f, msg)
    end
    return nothing
end

function log_header(outputfile::String)
    # outstr = "| Instance | L | N | L/N | Build/Obj (%) | Build (s) " *
    #          "| Incumbent (s) | Solve (s) | Status | Rt best bound " *
    #          "| Rt solve (s) | LB | UB | Gap (%) | Start (s) | RNRBS (s) " * 
    #          "| Rm | RNRBS impr | BS (s) | Start UB | PH (s) | \n"
    # outstr *= "|:---"^21 * "| \n"
    s = "| Instance | L | N | L/N | Start | Best | LB | UB | Gap | Add  " * 
        "| Heur | Solver | Time | Feas | \n"
    s *= "|:---"^14 * "| \n"
    log(outputfile, s)

    return nothing
end

function get_keys_results()
    # return ["build_obj_rat", "build_time", "incumbent_time", "solve_time", 
    #         "status", "root_best_bound", "root_time", "lb", "ub", "gap", 
    #         "fix_start_time", "rnr_time", "rm_rat", "rnr_impr_rat", 
    #         "bs_time", "start_ub", "ph_time"]
    return ["start", "best", "lb", "ub", "gap", "add_rat", 
            "heur_time", "solver_time", "time", "is_feas"]
end

"""
    log_instance(outputfile::String, 
                 inst_name::String, 
                 inst::Instance, 
                 results::Tuple)
"""
function log_instance(outputfile::String, 
                      inst_name::String, 
                      inst::Instance, 
                      results::Dict)
    N = inst.num_I
    L = (inst.num_K + inst.num_J)
    s = "| $inst_name | $L | $N | $(round(L / N, digits=2)) |"

    for k in get_keys_results()
        r = "-"
        if k in keys(results)
            r = results[k]
            if typeof(r) == Float64
                r = round(r, digits=2)
            end
        end
        s *= " $r |"
    end
    s *= "\n"

    log(outputfile, s)

    return nothing
end

function init_results()
    d = Dict{String, Any}()
    for k in get_keys_results()
        d[k] = const_infinite
    end
    return d
end

"""
    comp_bigM(inst, k::CandType)

Compute the big-M value for the model. 

The big-M is computed as discussed in Ghita's thesis, as there is at least one
existing line conecting every two buses.
"""
function comp_bigM(inst::Instance, k::CandType)
    bigM = const_infinite
    for j in keys(inst.J)
        if j[2] == k[1][2] && j[3] == k[1][3]
            tmp = inst.K[k].gamma * (inst.J[j].f_bar / inst.J[j].gamma)
            bigM = min(bigM, tmp)
        end
    end
    return abs(bigM)
end

"""
    is_one(I::Matrix{T}) where T <: AbstractFloat

Check if the matrix is the identity matrix.
"""
function is_one(I::Matrix{T}) where T <: AbstractFloat
    _, n = size(I)
    for i in 1:n, j in 1:n
        if i == j
            if !iseq(I[i, j], 1.0)
                return false
            end
        elseif !iseq(I[i, j], 0.0)
            return false
        end
    end
    return true
end

"""
    log(params::Parameters, msg::String, is_info::Bool = false)

Log message to console.
"""
function log(params::Parameters, msg::String, is_info::Bool = false)
    if params.log_level >= 1
        if is_info
            @info msg
        else
            println(msg)
        end
        # @info msg
    end
    return nothing
end

function log(st::Status)
    d = (st.rm_ratio, st.impr_ratio, st.time)
    d = round.(d, digits = 5)
    return "$(st.name) rm:$(d[1]) impr:$(d[2]) t:$(d[3])"
end

"""
    shift_to_existing_line(inst::Instance, params::Parameters, k::Int64)

Map candidate line k to corresponding existing circuit.
"""
function map_to_existing_line(inst::Instance, params::Parameters, k::Int64)
    return div(k - inst.num_J + params.instance.num_candidates - 1, 
               params.instance.num_candidates)
end

"""
    log_neigh_run(inst::Instance, 
                  params::Parameters, 
                  best_val::Float64, 
                  new_val::Float64, 
                  rm_ins::Set{Any}, 
                  runtime::Float64,
                  msg::String = "cost")

Log neighborhood run.
"""
function log_neigh_run(inst::Instance, 
                       params::Parameters, 
                       best_val::Float64, 
                       new_val::Float64, 
                       rm_ins::Set{Any}, 
                       runtime::Float64,
                       msg::String = "cost")
    log(params, "best_$msg:" * string(round(best_val, digits = 2)) *
                " new_$msg:" * string(round(new_val, digits = 2)) *
                " delta_perc:" * string(round(length(rm_ins) / 
                                            inst.num_K, digits = 2)) * 
                " runtime:" * string(round(runtime, digits = 2)))
    return nothing
end

"""
    comp_new_cost(inst::Instance, 
                  old_cost::Float64, 
                  removed::Vector{CandType})

Compute the new cost after removing some candidate circuits.
"""
function comp_new_cost(inst::Instance, 
                       old_cost::Float64, 
                       removed::Vector{CandType})
    new_cost = old_cost
    for k in removed
        new_cost -= inst.K[k].cost
    end
    return new_cost
end

"""
    get_values(vars::Dict{T, JuMP.VariableRef}) where T <: Union{Int64, Any}

Get values of vars variables.
"""
function get_values(vars::Dict{T, JuMP.VariableRef}) where 
                                                        T <: Union{Int64, Any}
    vals = Dict{T, Float64}()

    for (k, v) in vars
        vals[k] = value(v)
    end

    return vals
end

"""
    comp_gd_ratio(inst::Instance)

Compute ratio of generation over demand.
"""
function comp_gd_ratio(inst::Instance)
    sum_g = 0.0
    sum_d = 0.0
    for i in inst.I
        # Some buses may not have load or generation
        g = i in keys(inst.G) ? inst.G[i] : 0.0
        d = inst.D[i]

        sum_g += g
        sum_d += d
    end
    return sum_g / sum_d
end

"""
    fix_for_symmetry_contrs!(inst::Instance, 
                             params::Parameters, 
                             mip::MIPModel, 
                             start::Start)

Fix the start solution considering the symmetry constraints.
"""
function fix_for_symmetry_contrs!(inst::Instance, 
                                  params::Parameters, 
                                  mip::MIPModel, 
                                  start::Start)
    if params.model.is_symmetry_en
        for j in keys(inst.J)
            for l in params.instance.num_candidates:-1:1
                k = (j, l)
                k_plus = (j, l + 1)
                if iseq(inst.K[key].gamma, inst.K[k_plus].gamma) &&
                   !in(k, start.inserted) && in(kplus, start.inserted)
                    start.f[k], start.f[k_plus] = start.f[k_plus], start.f[k]
                    delete!(start.inserted, k_plus)
                    push!(start.inserted, k)
                end
            end
        end
    end
    return nothing
end

function set_state!(inst::Instance, 
                    mip::MIPModel)
    # In case the x is a single variable instead of a vector
    # if x isa JuMP.VariableRef
    #     x = [x]
    # end
    # if y isa JuMP.VariableRef
    #     y = [y]
    # end
    # if !(:state in keys(md.ext))
    #     md.ext[:state] = []
    # end
    mip.jump_model.ext[:state] = [mip.x[k] for k in keys(inst.K)]

    return nothing
end

function set_state!(inst::Instance, lp::LPModel)
    return nothing
end

function get_state_values(mip::MIPModel)
    return JuMP.value.(mip.jump_model.ext[:state])
end

function get_state_values(inst::Instance, inserted::Set{CandType})
    x = zeros(Float64, inst.num_K)

    for k in inserted
        i = inst.key_to_index[k]
        x[i] = 1.0
    end

    return x
end

function get_inserted_candidates(inst::Instance, mip::MIPModel)
    x = JuMP.value.(mip.jump_model.ext[:state])

    inserted = Set{CandType}()
    for (i, k) in enumerate(keys(inst.K))
        if iseq(x[i], 1.0)
            push!(inserted, k)
        end
    end

    return inserted
end

"""
    config_dcp_pm_tests!(params::Parameters)

Configure the parameters for tests against the DCP Power Models.
"""
function config_dcp_pm_tests!(params::Parameters)
    params.instance.num_candidates = 0
    params.instance.load_gen_mult = 1.0
    params.model.is_dcp_power_model_en = true
    params.log_level = 0
    return nothing
end    

"""
    comp_sparsity(A)

Compute the sparsity of a matrix.
"""
function comp_sparsity(A)
    return 1.0 - count(!iszero, A) / length(A)
end

"""
    comp_build_obj_rat(inst::Instance, 
                       obj::Float64, 
                       inserted::Vector{Any})

Compute the percentage ratio of the build cost over the total objective value.
"""
function comp_build_obj_rat(inst::Instance, 
                            obj::Float64, 
                            inserted)
    return 100.0 * comp_build_cost(inst, inserted) / obj
end


"""
    get_log_filename(inst::Instance, params::Parameters, scen::Int64)

Get the name of the log file.
"""
function get_log_filename(inst::Instance, params::Parameters, scen::Int64)
    return "$(params.log_dir)/$(inst.name)/#$scen.log"
end

function init_cache_in_rm(inst::Instance)
    return Set{CandType}(keys(inst.K)), Set{CandType}()
end

function roundp(sol::Set{CandType}, den::Int64, digits::Int64 = 2)
    return round(100.0 * length(sol) / den, digits = digits)
end

function roundp(num, den, digits::Int64 = 2)
    return round(100.0 * num / den, digits = digits)
end

function set_time_limit!(params::Parameters, 
                         lp::LPModel, 
                         start_time::Float64, 
                         time_limit::Float64)
    el = time() - start_time
    tl = max(time_limit - el, 0.0)
    # tl = min(tl, params.solver.lp_time_limit)
    JuMP.set_attribute(lp.jump_model, "TimeLimit", tl)
end