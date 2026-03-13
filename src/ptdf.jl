"""
    build_ptdf(inst::Instance, 
               params::Parameters, 
               scen::Int64 = 1; 
               T::Type{X} = Float16) where {X <: AbstractFloat}

Build the compact model for the TEP problem. 
    
In the original version of the model, there are no variables, since the 
generation is a parameter. In our version, we have added the generation as a 
variable.
"""
function build_ptdf(inst::Instance, 
                    params::Parameters, 
                    scen::Int64 = 1; 
                    T::Type{X} = Float16) where {X <: AbstractFloat}
    md = Model(Gurobi.Optimizer)

    bus_to_idx = build_bus_to_idx(inst)
    line_to_idx = build_line_to_idx(inst)

    d = build_d_injections(inst, bus_to_idx, params, scen)
    g, g_bus = build_g_injections(inst, bus_to_idx, params, scen, md)

    Gamma = build_susceptance_matrix(inst)

    S = build_incidence_matrix(inst, bus_to_idx; T = T)

    # I = SparseArrays.sparse(Matrix{T}(1.0 * LinearAlgebra.I(length(inst.I))))
    I = Matrix{T}(1.0 * LinearAlgebra.I(length(inst.I)))

    B, B_inv, beta = comp_beta(inst, bus_to_idx, params, Gamma, S, T = T)

    bus_inj = comp_bus_inj(d, g_bus)

    # TODO: Add auxiliar variable g_bus_sum[bus] == g_bus[bus] to avoid 
    # bottlenecks in the following operation
    f = beta * bus_inj

    # f = Vector{JuMP.AffExpr}()
    # s = Vector{JuMP.VariableRef}()
    # f_neg_cons = Vector{JuMP.ConstraintRef}()
    # f_pos_cons = Vector{JuMP.ConstraintRef}()

    s, f_neg_cons, f_pos_cons = add_thermal_limits_cons!(inst, params, md, f)

    add_load_balance_cons!(inst, md, d, g_bus)

    # The objective is to minimize slacks in the line flows
    # By constrolling the bounds on the slacks, we can simulate a model with or
    # without them
    obj = set_obj!(inst, params, scen, md, g, s)

    ptdf = PTDFModel{T}(md, obj, 
                        bus_to_idx, line_to_idx, 
                        Gamma, S, 
                        d, g, g_bus, 
                        B, B_inv, 
                        I, beta, bus_inj, f, 
                        s, f_neg_cons, f_pos_cons)

    # K = keys(inst.K)
    # K = collect(K)

    # Profile.clear()
    # @profile update_beta!(inst, bus_to_idx, params, ptdf, 1, K[1], T = T)
    # pprof()

    # update_beta!(inst, bus_to_idx, params, ptdf, 1, K[1], T = T)
    # old_beta = deepcopy(ptdf.beta)
    # rm_line!(inst, bus_to_idx, params, ptdf, 1, K[1], T = T)
    # add_line!(inst, bus_to_idx, params, ptdf, 1, K[1], T = T)

    # @warn norm(old_beta - ptdf.beta)

    config!(params, ptdf)

    return ptdf
end

function build_ptdf_system(inst::Instance, 
                           params::Parameters, 
                           # g_bus::Vector{Float64}, 
                           K, 
                           scen::Int64 = 1; 
                           T::Type{X} = Float16) where {X <: AbstractFloat}
    bus_to_idx = build_bus_to_idx(inst)
    line_to_idx = build_line_to_idx(inst)

    d = build_d_injections(inst, bus_to_idx, params, scen)

    Gamma = build_susceptance_matrix(inst, K)

    S = build_incidence_matrix(inst, bus_to_idx, T = T)

    I = Matrix{T}(1.0 * LinearAlgebra.I(length(inst.I)))

    B, B_inv, beta = comp_beta(inst, bus_to_idx, params, Gamma, S, T = T)

    # bus_inj = comp_bus_inj(d, g_bus)

    return PTDFSystem{T}(bus_to_idx, line_to_idx, 
                         Gamma, S, 
                         d, [], 
                         B, B_inv, 
                         I, beta, [])
end

