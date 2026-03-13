"""
    comp_gamma(params::Parameters, x::Float64, r::Float64 = 0.0)

Compute the susceptance of a circuit.
"""
function comp_gamma(params::Parameters, x::Float64, r::Float64 = 0.0)
    if params.model.is_dcp_power_model_en
        return -x / (r^2 + x^2)
    else
        return 1.0 / x
    end
    # return 1.0 / x
    # return -x / (r^2 + x^2)
end

"""
    get_circuit(inst, l)

Return the circuit at position l.
"""
function get_circuit(inst::Instance, l::Int64)
    if l <= inst.num_J
        return inst.J[l]
    else
        return inst.K[l - inst.num_J]
    end
end

"""
    get_existing_line(k::CandType)

Return k's corresponding existing line.
"""
function get_existing_line(k::CandType)
    return k[1]
end

"""
    get_inst_name(input::String)

Return the filename without the path and the extension.
"""
function get_inst_name(input::String)
    e = split(input, "/")[end]
    return String(split(e, ".")[1])
end

"""
    build_buses(data::Dict{String, Any})

Associate ids with indices in a vector.
"""
function build_buses(mpc::Dict{String, Any})
    I = Set{Int64}()
    # Get all ids
    for l in mpc["load"]
        push!(I, l[2]["load_bus"])
    end
    for g in mpc["gen"]
        # If generator is in service
        if g[2]["gen_status"] > 0
            push!(I, g[2]["gen_bus"])
        end
    end
    for b in mpc["branch"]
        # If branch is in service
        if b[2]["br_status"] > 0
            push!(I, b[2]["f_bus"])
            push!(I, b[2]["t_bus"])
        end
    end
    return I
end

"""
    build_loads(params::Parameters, 
                load::Dict{String, Any}, 
                shunt::Dict{String, Any})

Build Vector of loads (Pd) from the MATPOWER file.
"""
function build_loads(params::Parameters, 
                     load::Dict{String, Any}, 
                     shunt::Dict{String, Any})
    D = Dict{Int64, Float64}()
    for l in load
        D[l[2]["load_bus"]] = params.instance.load_gen_mult * l[2]["pd"]
    end
    if params.model.is_dcp_power_model_en
        for s in shunt
            bus = s[2]["shunt_bus"]
            if bus in keys(D)
                D[bus] += s[2]["gs"] * 1.0 ^ 2
            else
                D[bus] = s[2]["gs"] * 1.0 ^ 2
            end
        end
    end

    return D
end

"""
    build_gens(params::Parameters, 
               mpc::Dict{String, Any}, 
               cost_data::CostData, 
               inst_name::String)

Build generation data from MATPOWER file.
"""
function build_gens(params::Parameters, 
                    gen::Dict{String, Any}, 
                    cost_data::CostData, 
                    inst_name::String)
    # Get the inflation adjustment for the instance if it exists
    c_mult = 1.0
    if haskey(cost_data.gen_costs_mult, inst_name)
        c_mult = cost_data.gen_costs_mult[inst_name]
    else
        @warn "no inflation adjustment found for instance $inst_name\n" * 
                "\tusing 1.0 as multiplier."
    end
    
    G = Dict{Int64, GeneratorInfo}()
    for g in gen
        dt = g[2]
        # Machine out of service
        if dt["gen_status"] <= 0
            continue
        end
        lb = params.instance.load_gen_mult * dt["pmin"]
        ub = params.instance.load_gen_mult * dt["pmax"]
        costs = c_mult * abs.(reverse(dt["cost"]))
        G[dt["index"]] = GeneratorInfo(dt["gen_bus"], lb, ub, costs)
    end

    return G
end

"""
    build_existing_circuits(params::Parameters, 
                            baseMVA::Union{Int64, Float64}, 
                            cost_data::CostData)

Build existing lines, gamma values and capacities of circuits.
"""
function build_existing_circuits(params::Parameters, 
                                 mpc::Dict{String, Any}, 
                                 cost_data::CostData)
    J = Dict{Tuple3I, BranchInfo}()
    # min_gamma = 1e15
    # max_gamma = 0.0
    for b in mpc["branch"]
        dt = b[2]
        # Branch out of service
        if dt["br_status"] <= 0
            continue
        end
        
        x = dt["br_x"]
        r = dt["br_r"]

        if iseq(x, 0.0)
            fr = dt["f_bus"]
            to = dt["t_bus"]
            log(params, "x = 0.0 in circuit ($fr, $to).")
            # throw(ArgumentError("Error: x = 0.0 in circuit ($fr, $to)."))
            continue
        end

        j = (dt["index"], dt["f_bus"], dt["t_bus"])
        gamma = comp_gamma(params, x, r)
        cost = comp_existing_cost(mpc, cost_data, dt)
        # min_gamma = min(min_gamma, gamma)
        # max_gamma = max(max_gamma, gamma)
        J[j] = BranchInfo(dt["rate_a"], 
                          x, gamma, cost, 
                          (dt["angmin"], dt["angmax"]))
    end
    # @warn min_gamma, max_gamma
    # readline()

    return J
end

"""
    build_candidate_circuits!(params::Parameters, 
                              J::Dict{Tuple3I, BranchInfo})

Build candidate circuits based on exsting lines.
"""
function build_candidate_circuits(params::Parameters, 
                                  J::Dict{Tuple3I, BranchInfo})
    # TODO: K and J with the same key format
    K = Dict{CandType, BranchInfo}()
    rng = Random.MersenneTwister(params.instance.seed)
    # Candidate circuits are copies of the existing ones
    for (j, v) in J, l in 1:params.instance.num_candidates
        K[(j, l)] = deepcopy(v)
        # Compute the new costs based on the gamma values
        K[(j, l)].cost = comp_candidate_cost(params, v.cost, rng)
    end

    return K
end

"""
    function rm_g_quadratic_coeffs!(mpc)

Remove nonlinear coefficients of the generation terms to be used in the 
objective function.
"""
function rm_g_nonlinear_coeffs!(mpc::Dict{String, Any})
    for (i, g) in mpc["gen"]
        if length(g["cost"]) > 0
            new_g = reverse(g["cost"])[1:2]
            mpc["gen"][i]["cost"] = reverse(new_g)
        end
    end

    return nothing
end

"""
    read_reference_bus(params::Parameters, mpc::Dict{String, Any})

Read reference bus to instance.

Defaults to 1 if none is found.
"""
function read_reference_bus(params::Parameters, mpc::Dict{String, Any})
    ref_bus = params.instance.ref_bus
    for b in mpc["bus"]
        if b[2]["bus_type"] == mp_type_ref_bus
            ref_bus = b[2]["bus_i"]
            break
        end
    end

    return ref_bus
end

function read_cost_data(params::Parameters, costs_path::String)
    f = readlines(costs_path)

    # Parse circuit data
    i = findfirst(x -> contains(x, "# circuit data"), f) + 3
    @assert i != nothing "error: # circuit data section not found"

    voltage_classes = String[]
    reactances_km = Dict{String, Float64}()
    costs_km = Dict{String, Float64}()
    exp_lifetime = params.instance.expected_lifetime
    for line in f[i:end]
        if contains(line, "END")
            break
        end
        s = split(line)

        vclass = s[1]
        rkm = parse(Float64, s[2])
        ckm = parse(Float64, s[3])
        
        push!(voltage_classes, vclass)
        reactances_km[vclass] = rkm
        # Convert from M$/km-yr to $/km-hr
        costs_km[vclass] = ckm * 1e6 / (exp_lifetime * 365 * 24.0)
    end

    sort!(voltage_classes, by = x -> parse(Float64, x))

    # Parse transformer data
    i = findfirst(x -> contains(x, "# transformer data"), f) + 3
    @assert i != nothing "error: # transformer data section not found"

    transformers = Dict{Tuple{String, String}, Float64}()
    for line in f[i:end]
        if contains(line, "END")
            break
        end
        s = split(line)
        vclass1 = s[1]
        for j in 2:length(s)
            vclass2 = voltage_classes[j - 1]
            # Convert from $/MVA-yr to $/MVA-hr
            transformers[(vclass1, vclass2)] = 
                            parse(Float64, s[j]) / (exp_lifetime * 365 * 24.0)
        end
    end

    # Generation costs inflation adjustments
    i = findfirst(x -> contains(x, "# inflation adjustments"), f) + 3
    @assert i != nothing "error: # inflation adjustments section not found"

    gen_costs_mult = Dict{String, Float64}()
    for line in f[i:end]
        if contains(line, "END")
            break
        end
        s = split(line)
        gen_costs_mult[s[1]] = parse(Float64, s[3])
    end

    return CostData(voltage_classes, 
                    reactances_km, 
                    costs_km, 
                    transformers, 
                    gen_costs_mult)
end

"""
    vclass(cost_data::CostData, voltage_class::Int64)

Select the smallest voltage class that is greater than or equal to the given 
voltage class.
"""
function vclass(cost_data::CostData, voltage_class::String)
    i = findfirst(x -> parse(Float64, x) >= parse(Float64, voltage_class), 
                  cost_data.voltage_classes)
    @assert i != nothing "error: voltage class $voltage_class not assigned"

    return cost_data.voltage_classes[i]
end

function comp_length_km(baseMVA::Union{Int64, Float64}, 
                        cost_data::CostData, 
                        vclass::String, 
                        x_pu::Float64)
    x_ohms = x_pu * parse(Float64, vclass)^2 / baseMVA
    
    return x_ohms / cost_data.reactances_km[vclass]
end

function comp_circuit_cost(baseMVA::Union{Int64, Float64}, 
                           cost_data::CostData, 
                           voltage_class::String, 
                           x_pu::Float64)
    vcls = vclass(cost_data, voltage_class)
    length = comp_length_km(baseMVA, cost_data, vcls, x_pu)

    return cost_data.costs_km[vcls] * length
end

function comp_transformer_cost(baseMVA::Union{Int64, Float64}, 
                               cost_data::CostData, 
                               voltage_class_f::String, 
                               voltage_class_t::String)
    vcls_f = vclass(cost_data, voltage_class_f)
    vcls_t = vclass(cost_data, voltage_class_t)
    
    return cost_data.transformers[(vcls_f, vcls_t)] * baseMVA
end

function comp_existing_cost(mpc::Dict{String, Any},
                            cost_data::CostData, 
                            dt::Dict{String, Any})
    vcls_f = string(mpc["bus"]["$(dt["f_bus"])"]["base_kv"])
    vcls_t = string(mpc["bus"]["$(dt["t_bus"])"]["base_kv"])

    if vcls_f == vcls_t
        # Transmission line
        return comp_circuit_cost(mpc["baseMVA"], cost_data, vcls_f, dt["br_x"])
    else
        # Transformer
        return comp_transformer_cost(mpc["baseMVA"], cost_data, vcls_f, vcls_t)
    end
end

function comp_candidate_cost(params::Parameters, 
                             cost_existing_circuit::Float64, 
                             rng)
    # c = params.instance.cost_mult * abs(v.x)
    # K[(j, l)].cost = c / (params.instance.num_candidates + 1)

    # c = params.instance.cost_mult * 
    #         abs(v.x) / (params.instance.num_candidates + 1)
    # m = params.instance.cost_delta_mult * rand(rng, 1:10)
    # return c * (1 + m)

    m = params.instance.cost_delta_mult * rand(rng, 0:10)
    # cost_existing_circuit /= params.instance.num_candidates
    return  cost_existing_circuit * (1.0 + m)
end