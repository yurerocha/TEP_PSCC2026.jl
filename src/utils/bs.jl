"""
    comp_f_residuals()

Compute the residuals of the line flows.
"""
function comp_f_residuals(inst::Instance, 
                          ptdf::PTDFSystem, 
                          f::Dict{Any, Float64}, 
                          removed)
    f_residuals = Vector{Tuple{Any, Float64}}()

    for k in removed
        # Shift to the existing lines
        # j = map_to_existing_line(inst, k)
        diff = inst.K[k].f_bar - abs(f[k[1]])
        # if !isl(diff, 0.0) # diff >= 0.0
        r = diff / inst.K[k].f_bar
        push!(f_residuals, (k, r))
        # end
    end

    # Sort lines in non-descending order of residuals
    sort!(f_residuals, by = x->x[2], rev = false)

    return [f_residuals[i][1] for i in 1:length(f_residuals)]
end

function get_g_bus_values(inst::Instance, 
                          ptdf::PTDFSystem, 
                          tep::TEPModel, 
                          scen::Int64)
    g_bus = zeros(length(inst.I))

    for k in keys(inst.scenarios[scen].G)
        bus = ptdf.bus_to_idx[inst.scenarios[scen].G[k].bus]
        g_bus[bus] += JuMP.value(tep.g[k])
    end

    return g_bus
end

"""
    select_batches!(inst::Instance, 
                    params::Parameters, 
                    lp::LPModel, 
                    node::Node, 
                    candidates_per_batch_mult::Float64)

Compute the residuals of the line flows.
"""
function select_batches!(inst::Instance, 
                         params::Parameters, 
                         lp::LPModel, 
                         node::Node, 
                         candidates_per_batch_mult::Float64)
    # Disregard nodes that lead to an infeasible solution
    setdiff!(node.inserted, node.ignore)
    K = collect(node.inserted)

    # Disregard lines with the largest costs
    l = round(Int64, params.beam_search.restricted_list_ratio * length(K))
    partialsort!(K, 1:l, by = k -> inst.K[k].cost, rev = true)
    K = K[1:l]

    w = params.beam_search.num_children_per_parent
    b = candidates_per_batch_mult * length(K)
    b = max(floor(Int64, b), 1)

    samples = disjoint_samples(K, b, params.beam_search.is_shuffle_en)

    return select_best_child_nodes!(inst, samples, w)
end


"""
    disjoint_samples(lines::Vector{CandType}, b::Int, is_shuffle_en::Bool)

Divide the lines into disjoint samples of length b.
"""
function disjoint_samples(lines::Vector{CandType}, b::Int, is_shuffle_en::Bool)
    # Shuffle the indices, if required
    indices = is_shuffle_en ? randperm(length(lines)) : 1:length(lines)
    selected = lines[indices]

    # Divide the elements into batches
    samples = Vector{Set{CandType}}()
    # n = trunc(Int64, length(selected) / b)
    for i in 1:b:length(selected)
        j = i + b - 1
        push!(samples, Set{CandType}(selected[i:min(j, end)]))
    end

    return samples
end

"""
    select_best_child_nodes!(inst::Instance, 
                             samples::Vector{Vector{Any}}, 
                             w::Int64)

Select the w samples whose elements sum the largest costs.
"""
function select_best_child_nodes!(inst::Instance, 
                                  samples::Vector{Set{CandType}}, 
                                  w::Int64)
    sort!(samples, by = x -> sum([inst.K[k].cost for k in x]), rev = true)

    return samples[1:min(w, end)]
end

function comp_viol(inst::Instance, 
                   ptdf::T, 
                   beta::Matrix{S}, 
                   bus_inj::Vector{S},) where 
                                            {T <: Union{PTDFModel, PTDFSystem}, 
                                             S <: AbstractFloat}
    f = beta * bus_inj
    viol = 0.0
    i = 1

    for (_, line) in inst.J
        viol += max(abs(f[i]) - line.f_bar, 0.0)
        i += 1
    end
    for (_, line) in inst.K
        viol += max(abs(f[i]) - line.f_bar, 0.0)
        i += 1
    end

    return viol, f
end


function select!(params::Parameters, UB)
    sort!(UB, by = x -> x.cost)
    # Select N the l best upper bound solutions 
    N = params.beam_search.num_children_per_level
    l = floor(Int64, (1.0 + params.beam_search.num_children_per_level_mult) * N)
    l = min(l, length(UB))
    N = min(N, length(UB))
    if l > N
        UB = Random.shuffle(UB[1:l])
    end
    # Return, randomly, N samples among the best l UBs
    return UB[1:N]
end

function delete_value!(vec::Vector, val)
    idx = findfirst(==(val), vec)
    if idx !== nothing
        deleteat!(vec, idx)
        return true
    end

    return false
end

function add_node!(params::Parameters, 
                   UB::Vector{Node}, 
                   msg::BSWorkerMessage, 
                   node::Node)
    inserted = Set{CandType}()
    removed = Set{CandType}()
    ignore = copy(node.ignore)
    # TODO: Se o custo for menor do que o custo anterior para o mesmo nó, 
    # remover. Só remover quando melhorar o best obj, pode fazer com que um 
    # mesmo batch seja reavaliado múltiplas vezes
    # TODO: Se não melhorar e !params.beam_search.is_shuffle_en, colocar na
    # lista de ignore. Fazer testes apenas com essa alteração. Em seguida, fazer
    # testes trocando os ifs a seguir.
    if isl(msg.cost, node.cost)
        inserted = setdiff(node.inserted, msg.lines)
        removed = union(node.removed, msg.lines)
    else
        inserted = copy(node.inserted)
        removed = copy(node.removed)
        # msg.obj = node.obj # Parent's obj
        if !params.beam_search.is_shuffle_en
            # If infeasible, place the lines in the ignore list
            ignore = union(node.ignore, msg.lines)
        end
    end
    push!(UB, Node(msg.cost, msg.viol, inserted, removed, ignore))

    return nothing
end

function get_data(inst::Instance, 
                  lp::LPModel, 
                  ptdf::PTDFSystem)
    g_bus = get_g_bus_values(inst, ptdf, lp)
    bus_inj = comp_bus_inj(ptdf.d, g_bus)
    viol = comp_viol(lp)
    
    return g_bus, bus_inj, viol
end

"""
    comp_penalized_cost(inst::Instance, 
                        params::Parameters, 
                        scen::Int64, 
                        lp::LPModel, 
                        cache::WorkerCache, 
                        inserted::Set{CandType})

Compute the penalized cost considering inserted candidate lines and the cache 
memory of the progressive hedging algorithm, if in the progressive hedging.
"""
function comp_penalized_cost(inst::Instance, 
                             params::Parameters, 
                             scen::Int64, 
                             lp::LPModel, 
                             cache::WorkerCache, 
                             inserted::Set{CandType})
    cost = params.progressive_hedging.penalty_mult * comp_viol(lp)
    g_cost = comp_g_cost(inst, params, scen, lp)

    if params.progressive_hedging.is_en
        sq2norm = 0.0
        acc_omega = 0.0
        for k in inserted
            cost += inst.K[k].cost
            # Progressive hedging data
            i = inst.key_to_index[k]
            # sq2norm += 1.0 - 2.0 * cache.x_hat[i] + cache.x_hat[i] ^ 2
            sq2norm += 
                cache.rho[i] * (1.0 - 2.0 * cache.x_hat[i] + cache.x_hat[i] ^ 2)
            acc_omega += cache.scenarios[scen].omega[i]
        end
        # cost += acc_omega + 
        #        (params.progressive_hedging.rho / 2.0) * sq2norm
        cost += acc_omega + sq2norm
    else
        cost += comp_build_cost(inst, inserted)
    end

    return cost + g_cost, g_cost
end

"""
    comp_build_cost(inst::Instance, inserted::Set{CandType})

Compute the cost of the solution.
"""
function comp_build_cost(inst::Instance, inserted::Set{CandType})
    return sum(inst.K[k].cost for k in inserted; init = 0.0)
end

function comp_g_cost(inst::Instance, 
                     params::Parameters, 
                     scen::Int64, 
                     tep::TEPModel)
    if JuMP.result_count(tep.jump_model) == 0
        return const_infinite
    end

    g = JuMP.value.(values(tep.g))
    cost = 0.0
    i = 1
    for k in keys(tep.g)
        cost += comp_g_obj(params, g[i], inst.scenarios[scen].G[k].costs)
        i += 1
    end

    return cost
    # return sum(comp_g_obj(params, g[i], 
    #     inst.scenarios[scen].G[k].costs) for (i, k) in enumerate(keys(tep.g)))
    # return sum(comp_g_obj(params, JuMP.value(tep.g[k]), 
    #                     inst.scenarios[scen].G[k].costs) for k in keys(tep.g))
end

function comp_candidates_per_batch_mult(inst::Instance, 
                                        params::Parameters, 
                                        inserted::Set{CandType})
    n1 = length(inserted) / inst.num_K
    n2 = inst.num_I / 1000.0
    return params.beam_search.candidates_per_batch_mult * n2 * n1
end