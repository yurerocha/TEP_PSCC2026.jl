"""
    detect_cycles(inst, model, is_drawing_en, filename="graph")

Detect cycles in the graph and return the free buses (those outside cycles).
"""
function detect_cycles(inst::Instance, 
                       md::MIPModel, 
                       is_drawing_en::Bool = true, 
                       filename::String = "graph")
    # TODO: Add "demand - generation" to the vertices.
    @info "Detect cycles"
    elist = []
    for c in inst.J
        push!(elist, (c.fr, c.to))
    end
    unique!(elist)
    # @info length(elist), elist

    g = SimpleGraph(Graphs.SimpleEdge.(elist))
    # g = SimpleDiGraph(Graphs.SimpleEdge.(elist))

    cycles = cycle_basis(g)
    # @info cycles
    # c = simplecycles(g)
    # c = simplecycles_limited_length(g, 10, 10)

    if is_drawing_en
        # vertices = collect(inst.I)
        # sort!(vertices)
        draw_cycles(dt, md, elist, cycles, filename)
        @info "Done drawing graph"
    end

    buses_per_cycle = Vector{Set{Int}}(undef, length(cycles))
    for (i, c) in enumerate(cycles)
        buses_per_cycle[i] = Set{Int}()
        for v in c
            push!(buses_per_cycle[i], v)
        end
    end

    # free_buses = setdiff(inst.I, busy_buses)
    # @info "Free buses: ", free_buses
    return cycles, buses_per_cycle
end

"""
    draw_cycles(inst, model, graph, cycles, filename)
"""
function draw_cycles(inst::Instance, 
                     md::MIPModel, 
                     elist::Vector{Tuple{Int64, Int64}}, 
                     cycles::Vector{Vector{Int64}}, 
                     filename::String)
    g = SimpleDiGraph(Graphs.SimpleEdge.(elist))
    flows = Dict{Circuit, Float64}()
    for l in 1:inst.num_J+inst.num_K
        c = get_circuit(dt, l)
        if c in keys(flows)
            flows[c] += value(md.f[l])
        else
            flows[c] = value(md.f[l])
        end
    end
    vertices = [round(Int, value(md.g[i]) - 
                (i in keys(inst.D) ? inst.D[i] : 0.0)) for i in 1:inst.num_I]
    edgeweights = [round(Int, flows[Circuit(src(e), dst(e))]) for e in edges(g)]
    # Generate n maximally distinguishable colors in LCHab space.
    vertexfillc = distinguishable_colors(nv(g), colorant"blue")
    @svg begin
        background("grey10")
        fontsize(8)
        sethue("white")
        println("Drawing first layer")
        @layer begin
            drawgraph(g,
                      layout = stress, 
                      vertexlabels = vertices,
                      edgegaps = 12,
                      edgestrokeweights = 0.5,
                      edgelabels = edgeweights,
                      vertexfillcolors = vertexfillc,
                      vertexshapesizes = 12
                      # vertexfillcolors = 
                      #     [RGB(rand(3)/2...) 
                      #        for i in 1:nv(g)]
            )
        end
        for (n, c) in enumerate(cycles)
            # @printf "\rLayer: %d of %d" n length(cycles)
            cycleedges = [Edge(c[i], c[mod1(i + 1, end)]) for i in eachindex(c)]
            @layer begin
                sethue(HSB(rescale(n, 1, length(cycles) + 1, 0, 360), 0.8, 0.6))
                setopacity(0.1)
                drawgraph(g, 
                          layout = stress,
                          # vertexlabels = (v) -> v in c && string(v),
                          vertexlabels = vertices,
                          edgegaps = 12,
                          edgelist = cycleedges,
                          edgestrokeweights = 3
                )
            end
        end
    end 2000 2000 filename * ".svg"
end

"""
    draw_solution(inst::Instance, 
                  params::Parameters, 
                  md::T, 
                  f::Vector{Float64}, 
                  viols::Vector{Tuple{Float64, Int64}}, 
                  filename::String = "solution") where T <: TepModel

Draw the graph of a solution.

Each edge is labeled with the flow value and each vertex is labeled with
"generation - demand".
"""
function draw_solution(inst::Instance, 
                       params::Parameters, 
                       md::T, 
                       f::Vector{Float64}, 
                       viols::Vector{Tuple{Float64, Int64}}, 
                       filename::String = "solution") where T <: TepModel
    flows = Dict{Circuit, Float64}()
    circuits = Dict{Circuit, Int64}()
    for l in 1:inst.num_J + inst.num_K
        c = get_circuit(inst, l)
        if c in keys(flows)
            flows[c] += value(f[l])
        else
            flows[c] = value(f[l])
            # Shift to the candidates
            k = inst.num_J + 1 + params.num_candidates * (l - 1)
            circuits[c] = k
        end
    end
    elist = []
    for l in 1:inst.num_J + inst.num_K
        c = get_circuit(inst, l)
        push!(elist, (c.fr, c.to))
    end
    # @show flows
    unique!(elist)

    g = SimpleDiGraph(Graphs.SimpleEdge.(elist))
    edgestrokecolors = [colorant"grey" for _ in elist]
    for (i, e) in enumerate(edges(g))
        for v in viols
            c = get_circuit(inst, v[2])
            if (src(e), dst(e)) == (c.fr, c.to)
                edgestrokecolors[i] = colorant"blue"
            end
        end
    end

    delta = [round(Int, (i in keys(md.g) ? value(md.g[i]) : 0.0) - 
             (i in keys(inst.D) ? inst.D[i] : 0.0)) for i in 1:inst.num_I]
    # edgelabels = ["$(round(Int, flows[Circuit(src(e), dst(e))])):$" 
    #               for e in edges(g)]
    edgelabels = []
    for e in edges(g)
        c = Circuit(src(e), dst(e))
        push!(edgelabels, "$(round(Int, flows[c])):$(circuits[c])")
    end
    # Generate n maximally distinguishable colors in LCHab space.
    # vertexfillc = distinguishable_colors(nv(g), colorant"blue")

    vertexfillc = Array{Colorant}(undef, nv(g))
    vertices = Array{String}(undef, nv(g))
    I = sort(collect(inst.I))
    @show length(I), length(delta)
    vertexshapesizes = []
    for i in eachindex(delta)
        d = delta[i]
        theta = value(md.theta[i])
        t = round(theta, digits=2)
        # vertices[i] = "$(I[i]) $d $t"
        vertices[i] = "$(I[i]) $t"
        # vertices[i] = "$(I[i]) $d"
        if d < 0 
            vertexfillc[i] = colorant"red"
        elseif d > 0
            vertexfillc[i] = colorant"green"
        else
            vertexfillc[i] = colorant"yellow"
        end
        push!(vertexshapesizes, 10 * log10(abs(d)))
    end
    @svg begin
        # background("grey")
        background("white")
        fontsize(6)
        sethue("black")
        text(filename, boxbottomcenter() + (0, -50), halign=:center)
        println("Drawing first layer")
        @layer begin
            drawgraph(g,
                    #   edgecurvature = 5,
                      layout = stress, 
                      vertexlabels = vertices,
                      edgegaps = 15,
                      edgestrokeweights = 1,
                      edgelabels = edgelabels,
                      edgestrokecolors = edgestrokecolors,
                      vertexfillcolors = vertexfillc,
                      vertexshapesizes = vertexshapesizes,
                      vertexlabeltextcolors = colorant"black"
                      # vertexfillcolors = 
                      #     [RGB(rand(3)/2...) 
                      #        for i in 1:nv(g)]
            )
        end
    end 1000 1000 filename * ".svg"
end