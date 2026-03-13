function build_clusters(inst::Instance)
    Gamma = fill(const_infinite, inst.num_I, inst.num_I)
    for k in keys(inst.J)
        b1 = inst.bus_to_index[k[2]]
        b2 = inst.bus_to_index[k[3]]
        Gamma[b1, b2] = inst.J[k].gamma
    end

    # Project into 2 dimensions
    start = time()
    coords = MultivariateStats.classical_mds(Gamma, 2)

    R = Clustering.kmeans(coords, 100)

    println(Clustering.assignments(R))
    println("Elapsed time: $(time() - start)")

    return nothing
end

function clusters_reinsert!(inst::Instance, 
                            params::Parameters, 
                            scen::Int64, 
                            lp::LPModel, 
                            clusters::Vector{Int64}, 
                            inserted::Set{CandType}, 
                            removed::Set{CandType}, 
                            best_viol::Float64)
    # Pega os clusters cujas linhas incidentes estão envolvidas em violações
    involved_clusters = Set{Int64}()
    for k in keys(inst.J)
        if isg(JuMP.value(lp.s[k]), 0.0)
            c1 = clusters[inst.bus_to_index[k[2]]]
            c2 = clusters[inst.bus_to_index[k[3]]]
            push!(involved_clusters, c1)
            push!(involved_clusters, c2)
        end
    end
    
    reinsert = Set{CandType}()
    for k in removed
        c1 = clusters[inst.bus_to_index[k[1][2]]]
        c2 = clusters[inst.bus_to_index[k[1][3]]]
        if c1 in involved_clusters || c2 in involved_clusters
            push!(reinsert, k)
        end
    end

    new_inserted = union(inserted, reinsert)
    update_lp!(inst, params, lp, new_inserted)
    viol = comp_viol(lp)

    num_applied = 1
    num_successes = 0
    if isl(viol, best_viol)
        best_viol = viol
        num_successes = 1

        l = length(removed)
        setdiff!(removed, reinsert)
        # best_viol = repair!(inst, params, scen, lp, removed, best_viol)
        @info "reinforce impr removed:$l -> $(length(removed))"
    else
        @info "reinforce did not succeed"
    end

    num_applied = 1
    num_successes = iseq(best_viol, 0.0) ? 1 : 0

    return best_viol, (num_applied, num_successes, 0)
end

function build_clusters_cycles(inst::Instance)
    elist = []
    for k in keys(inst.J)
        b1 = inst.bus_to_index[k[2]]
        b2 = inst.bus_to_index[k[3]]
        push!(elist, (b1, b2))
    end

    # Project into 2 dimensions
    start = time()
    g = SimpleGraph(Graphs.SimpleEdge.(elist));
    cycles = cycle_basis(g)
    # cycles = simplecycles_limited_length(g, 20)
    # Map each vertex → cycle index(es)
    println("Elapsed time: $(time() - start)")
    println("Num cycles:$(length(cycles))")

    readline()

    println("Cycles: $cycles")

    return nothing
end

function build_cycle_memberships(cycles)
    cycle_memberships = Dict(v => Vector{Int}() for v in 1:nv(g))

    for (i, cyc) in enumerate(cycles)
        for v in cyc
            push!(cycle_memberships[v], i)
        end
    end

    cycle_memberships
end

function build_graph(inst::Instance)
    elist = []
    for k in keys(inst.J)
        b1 = inst.bus_to_index[k[2]]
        b2 = inst.bus_to_index[k[3]]
        w = inst.J[k].gamma
        push!(elist, (b1, b2, w))
    end

    start = time()
    g = SimpleGraph(Graphs.SimpleEdge.(elist));

    # Project into 2 dimensions
    dists = Graphs.shortest_paths(g, 1)
    println(Graphs.distance(dists, 6))

    println("Elapsed time: $(time() - start)")

    return nothing
end