function comp_load_sum(D::Dict{String, Any})
    return sum(l["pd"] for l in values(D))
end

function comp_gen_sum(G::Dict{String, Any}, g_key::String)
    return sum(g[g_key] for g in values(G) if g["gen_status"] > 0)
end

function comp_load_gen_rat(D::Dict{String, Any}, G::Dict{String, Any})
    return comp_load_sum(D) / comp_gen_sum(G, "pmax")
end

"""
    select_ren!(pglib_mpc::Dict{String, Any}, 
                cats_ren_avg::Float64, 
                cats_ren_ratio::Float64, 
                candidates::Set{String})

Randomly select generators from the PGLib-OPF system to act as renewable 
generators.

The generation profile from the associated CATS scenario is preserved.
"""
function select_ren!(pglib_mpc::Dict{String, Any}, 
                     cats_ren_avg::Float64, 
                     cats_ren_ratio::Float64, 
                     candidates::Set{String})
    gen_cap = comp_gen_sum(pglib_mpc["gen"], "pmax")
    ren_cap = 0.0
    ren_index = Set{String}()
    while !isempty(candidates) && isl(ren_cap / gen_cap, cats_ren_ratio)
        k = select_cand(pglib_mpc, cats_ren_avg, candidates)
        pop!(candidates, k)
        # Machine in service
        if pglib_mpc["gen"][k]["gen_status"] > 0
            ren_cap += pglib_mpc["gen"][k]["pmax"]
            push!(ren_index, k)
        end
    end

    return ren_index
end

function select_cand(pglib_mpc::Dict{String, Any}, 
                     cats_ren_avg::Float64, 
                     candidates::Set{String})
    sel_k = 0
    min_delta = const_infinite
    for k in candidates
        delta = abs(pglib_mpc["gen"][k]["pmax"] - cats_ren_avg) 
        if isl(delta, min_delta)
            sel_k = k
            min_delta = delta
        end
    end
    return sel_k
end

"""
    scale_ren_gen!(params::Parameters, 
                   scen::Int64, 
                   G::Dict{String, Any}, 
                   pglib_mpc::Dict{String, Any}, 
                   pglib_ren_gen_indices::Set{String}, 
                   cats_ren_gen::Dict{Int64, Float64}, 
                   cats_ren_cap::Float64)

Scale the capacity of renewable generators in a PGLIB-OPF scenario according 
to a constant involving renewable generation computed for the associated CATS
scenario.
"""
function scale_ren_gen!(params::Parameters, 
                        scen::Int64, 
                        G::Dict{String, Any}, 
                        pglib_mpc::Dict{String, Any}, 
                        pglib_ren_gen_indices::Set{String}, 
                        cats_ren_gen::Dict{Int64, Float64}, 
                        cats_ren_cap::Float64)
    for k in pglib_ren_gen_indices
        G[k]["pmin"] = 0.0
        g = comp_gen(cats_ren_gen[scen], 
                     pglib_mpc["baseMVA"], 
                     pglib_mpc["gen"][k]["pmax"], 
                     cats_ren_cap)
        G[k]["pmax"] = 
                round(g, digits = params.stochastic_instance.rounding_digits)
    end
    return nothing
end

"""
    comp_gen(ren_gen::Float64, 
             baseMVA::Float64, 
             gen::Float64, 
             cap::Float64)

Compute the generator capacity.
"""
function comp_gen(ren_gen::Float64, 
                  baseMVA::Float64, 
                  gen::Float64, 
                  cap::Float64)
    return ren_gen/baseMVA * gen/cap
end

"""
    scale_gen_lb!(params::Parameters, 
                  cats_mpc::Dict{String, Any}, 
                  pglib_mpc::Dict{String, Any}, 
                  G::Dict{String, Any})

Scale the lower bounds of generation in the PGLib-OPF base instance for the 
current scenario to match with the proportion of generation lower bound in the 
associated CATS scenario.
"""
function scale_gen_lb!(params::Parameters, 
                       cats_mpc::Dict{String, Any}, 
                       pglib_mpc::Dict{String, Any}, 
                       G::Dict{String, Any})
    gen_lb = comp_gen_sum(cats_mpc["gen"], "pmin")
    gen_ub = comp_gen_sum(cats_mpc["gen"], "pmax")
    cats_ratio = gen_lb / gen_ub

    gen_lb = comp_gen_sum(G, "pmin")
    gen_ub = comp_gen_sum(G, "pmax")
    new_ratio = gen_lb / gen_ub

    if iseq(gen_ub, 0.0) || iseq(new_ratio, 0.0)
        return nothing
    end

    m = cats_ratio / new_ratio
    for (k, g) in G
        if g["gen_status"] > 0
            G[k]["pmin"] = round(m * g["pmin"], 
                            digits = params.stochastic_instance.rounding_digits)
        end
    end

    if params.debugging_level == 1
        gen_lb = comp_gen_sum(G, "pmin")
        r1 = cats_ratio
        r2 = gen_lb / gen_ub
        @assert iseq(r1, r2, 1e-3) "diff gen ratio [cats and new scen $r1 $r2"
    end

    return nothing
end

"""
    scale_loads(params::Parameters, 
                cats_mpc::Dict{String, Any}, 
                pglib_mpc::Dict{String, Any}, 
                G::Dict{String, Any})

Scale the loads in the PGLib-OPF base instance for the current scenario to match 
with the proportion of load in the associated CATS scenario.
"""
function scale_loads(params::Parameters, 
                     cats_mpc::Dict{String, Any}, 
                     pglib_mpc::Dict{String, Any}, 
                     G::Dict{String, Any})
    cats_ratio = comp_load_gen_rat(cats_mpc["load"], cats_mpc["gen"])
    new_ratio = comp_load_gen_rat(pglib_mpc["load"], G)
    gen_cap = comp_gen_sum(G, "pmax")

    D = deepcopy(pglib_mpc["load"])
    if iseq(gen_cap, 0.0) || iseq(new_ratio, 0.0)
        return D
    end

    m = cats_ratio / new_ratio
    for (k, l) in pglib_mpc["load"]
        D[k]["pd"] = round(m * l["pd"], 
                            digits = params.stochastic_instance.rounding_digits)
    end

    if params.debugging_level == 1
        sum_load = comp_load_sum(D)
        r1 = cats_ratio
        r2 = sum_load / gen_cap
        @assert iseq(r1, r2, 1e-3) "diff load ratio [cats and new scen $r1 $r2"
    end

    return D
end

"""
    comp_ren_and_avg_cap(mpc::Dict{String, Any}, 
                         gen_data::DataFrame, 
                         ren_type::String)

Compute the renewable generation capacity.
"""
function comp_ren_and_avg_cap(mpc::Dict{String, Any}, 
                              gen_data::DataFrames.DataFrame, 
                              ren_type::String)
    gen_indices = [g for g in 1:size(gen_data)[1] if occursin(ren_type, 
                                               lowercase(gen_data.FuelType[g]))]
    
    # Renewable generation capacity
    rcap = sum(g["pmax"] for (i, g) in mpc["gen"] if g["index"] in gen_indices)
    cap = comp_gen_sum(mpc["gen"], "pmax")

    return rcap, rcap / length(gen_indices)
end

"""
    dict_hourly_data(hourly_data::DataFrames.DataFrame, 
                     scenarios::Vector{Int64}, 
                     ren_type::String)

Get the hourly data for scenarios and renewable type as a Dict.
"""
function dict_hourly_data(hourly_data::DataFrames.DataFrame, 
                          scenarios::Vector{Int64}, 
                          ren_type::String)
    return Dict([k => hourly_data[k, ren_type] for k in scenarios])
end

function scale_loads_gens!(inst::Instance, 
                           params::Parameters, 
                           pglib_mpc::Dict{String, Any})
    # Select the instance in a given percentile regarding their loads
    scenarios = deepcopy(inst.scenarios)
    sort!(scenarios, by = x -> sum(values(x.D)))
    scen = params.stochastic_instance.selection_percentile * length(scenarios)
    scen = round(Int64, scen)

    # Compute the multiplier for the selected instance according to parameter 
    # load_gen_mult
    pglib_sum = comp_load_sum(pglib_mpc["load"])
    sel_sum = sum(values(scenarios[scen].D))
    m = params.instance.load_gen_mult * pglib_sum / sel_sum
    @info "scen#$scen m:$m pglib_sum:$(round(pglib_sum, digits = 2))"

    # Update the instances
    for scen in eachindex(inst.scenarios)
        # d1 = sum(values(scenarios[scen].D))
        for k in keys(inst.scenarios[scen].D)
            inst.scenarios[scen].D[k] *= m
        end
        d = sum(values(inst.scenarios[scen].D))
        # msg = "scen$scen d:$(round(d1, digits = 2)) -> " * 
        #         "$(round(d2, digits = 2)) " * 
        #         "r:$(round(d2 / pglib_sum, digits = 2))"
        # println(msg)

        g_lb = 0.0
        g_ub = 0.0
        rd = params.stochastic_instance.rounding_digits
        for k in keys(inst.scenarios[scen].G)
            g = inst.scenarios[scen].G[k]
            inst.scenarios[scen].G[k].lower_bound = 
                                        round(m * g.lower_bound, digits = rd)
            inst.scenarios[scen].G[k].upper_bound = 
                                        round(m * g.upper_bound, digits = rd)
            g_lb += g.lower_bound
            g_ub += g.upper_bound
        end
        @info round.([d, g_lb, g_ub, d / g_ub, d / pglib_sum], digits = 2)
    end

    return nothing
end

function assert_feas_build_all(inst::Instance, params::Parameters)
    for scen in eachindex(inst.scenarios)
        @info "run feasibility test for scen#$scen"
        mip = build_mip(inst, params, scen)
        is_feas = 
                fix_start!(inst, params, scen, mip, Set{CandType}(keys(inst.K)))
        @assert is_feas "scen#$scen infeasible build all"
    end

    return nothing
end