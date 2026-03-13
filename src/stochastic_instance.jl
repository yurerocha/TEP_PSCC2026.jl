
"""
    build_stochastic_instance(params::Parameters, 
                              filepath::String, 
                              costs_path::String="MISO_TCEG_MTEP24.txt", 
                              days::Set{Int64} = Set([79, 171, 265, 355]))

Build stochastic instance considering a set of days
- Outono: 20 de março (equinócio de outono)
- Inverno: 20 de junho (solstício de inverno)
- Primavera: 22 de setembro (equinócio de primavera)
- Verão: 21 de dezembro (solstício de verão)
"""
function build_stochastic_instance(params::Parameters, 
                                   filepath::String, 
                                   costs_path::String="MISO_TCEG_MTEP24.txt", 
                                   days::Set{Int64} = Set([79, 171, 265, 355]))
    log(params, "Build stochastic instance", true)

    cost_data = read_cost_data(params, costs_path)
    inst_name = get_inst_name(filepath)

    # -------------- Set the scenarios based on the selected days --------------
    scenarios = Vector{Int64}()
    num_days = 24
    for d in days
        s = (d - 1) * num_days + 1
        e = d * num_days
        scenarios = vcat(scenarios, collect(range(s, e)))
    end
    # start_scen = 8649
    # N = 8760
    
    # ----------------------------- Read CATS data -----------------------------
    dir = "external/CATS-CaliforniaTestSystem"
    
    load_scens = CSV.read("$dir/Data/Load_Agg_Post_Assignment_v3_latest.csv", 
                          header = false, DataFrame)
    load_scens = [Dict([k => load_scens[i, k] for k in scenarios]) 
                                                for i in 1:size(load_scens)[1]]
    
    mpc = PowerModels.parse_file("$dir/MATPOWER/CaliforniaTestSystem.m")
    
    gen_data = CSV.read("$dir/GIS/CATS_gens.csv", DataFrame)

    load_mapping = map_buses_to_loads(mpc)
    
    hourly_data = CSV.read("$dir/Data/HourlyProduction2019.csv", DataFrame)
    solar_gen = dict_hourly_data(hourly_data, scenarios, "Solar")
    wind_gen = dict_hourly_data(hourly_data, scenarios, "Wind")

    # ----- Build base PGLib-OPF scen generators with renewable penetration ----
    solar_cap, solar_avg = comp_ren_and_avg_cap(mpc, gen_data, "solar")
    wind_cap, wind_avg = comp_ren_and_avg_cap(mpc, gen_data, "wind")
    gen_cap = sum(g["pmax"] for (i, g) in mpc["gen"])

    pglib_mpc = PowerModels.parse_file(filepath)
    candidates = Set(keys(pglib_mpc["gen"]))
    solar_gen_indices = select_ren!(pglib_mpc, solar_avg, 
                                    solar_cap / gen_cap, candidates)
    wind_gen_indices = select_ren!(pglib_mpc, wind_avg, 
                                   wind_cap / gen_cap, candidates)

    # ------------------------ Build multiple scenarios ------------------------
    inst = build_instance(params, filepath, costs_path)
    inst.scenarios = []
    inst.num_scenarios = length(scenarios)
    prob = 1.0 / inst.num_scenarios
    for (i, scen) in enumerate(scenarios)
        G = deepcopy(pglib_mpc["gen"])
        scale_ren_gen!(params, scen, G, pglib_mpc, 
                       solar_gen_indices, 
                       solar_gen, solar_cap)
        scale_ren_gen!(params, scen, G, pglib_mpc, 
                       wind_gen_indices, 
                       wind_gen, wind_cap)
        scale_gen_lb!(params, mpc, pglib_mpc, G)

        # Update the loads in the CATS mpc according to the current scenario
        update_loads!(scen, load_scens, mpc, load_mapping)

        D = scale_loads(params, mpc, pglib_mpc, G)

        # Parse to the Instance format
        D = build_loads(params, D, pglib_mpc["shunt"])
        G = build_gens(params, G, cost_data, inst_name)

        sumD = sum(d for d in values(D))
        sum_lb = sum(g.lower_bound for g in values(G))
        sum_ub = sum(g.upper_bound for g in values(G))

        # log(params, 
        #     "Scen#$i: $sumD, $sum_lb, $sum_ub, $(sumD / sum_ub)", 
        #     true)
        if params.debugging_level == 1
            @assert isl(sum_lb, sumD) "sum lb gen $sum_lb > sum load $sumD"
            @assert isl(sumD, sum_ub) "sum load $sumD > sum ub gen $sum_ub"
        end

        push!(inst.scenarios, Scenario(prob, D, G))
    end

    scale_loads_gens!(inst, params, pglib_mpc)

    if params.debugging_level == 1
        assert_feas_build_all(inst, params)
    end
    
    return inst
end

function build_cats_stochastic_instance(params::Parameters, 
                                    filepath::String, 
                                    costs_path::String="MISO_TCEG_MTEP24.txt", 
                                    days::Set{Int64} = Set([79, 171, 265, 355]))
    log(params, "Build CATS stochastic instance", true)

    cost_data = read_cost_data(params, costs_path)
    inst_name = get_inst_name(filepath)

    # -------------- Set the scenarios based on the selected days --------------
    scenarios = Vector{Int64}()
    num_days = 24
    for d in days
        s = (d - 1) * num_days + 1
        e = d * num_days
        scenarios = vcat(scenarios, collect(range(s, e)))
    end

    # ----------------------------- Read CATS data -----------------------------
    dir = "external/CATS-CaliforniaTestSystem"
    
    load_scens = CSV.read("$dir/Data/Load_Agg_Post_Assignment_v3_latest.csv", 
                          header = false, DataFrame)
    load_scens = [Dict([k => load_scens[i, k] for k in scenarios]) 
                                                for i in 1:size(load_scens)[1]]

    mpc = PowerModels.parse_file("$dir/MATPOWER/CaliforniaTestSystem.m")

    gen_data = CSV.read("$dir/GIS/CATS_gens.csv", DataFrame)

    load_mapping = map_buses_to_loads(mpc)

    hourly_data = CSV.read("$dir/Data/HourlyProduction2019.csv", DataFrame)
    solar_gen = dict_hourly_data(hourly_data, scenarios, "Solar")
    wind_gen = dict_hourly_data(hourly_data, scenarios, "Wind")

    pmax_og = [mpc["gen"][string(i)]["pmax"] for i in 1:size(gen_data)[1]]

    solar_cap, _ = comp_ren_and_avg_cap(mpc, gen_data, "solar")
    wind_cap, _ = comp_ren_and_avg_cap(mpc, gen_data, "wind")


    # ------------------------ Build multiple scenarios ------------------------
    inst = build_instance(params, filepath, costs_path)
    inst.scenarios = []
    inst.num_scenarios = length(scenarios)
    prob = 1.0 / inst.num_scenarios

    for k in scenarios
        # Change renewable generators' pg for the current scenario
        update_rgen!(k, mpc, gen_data, solar_gen, 
                        wind_gen, pmax_og, solar_cap, wind_cap)

        # Change load buses' Pd and Qd for the current scenario
        update_loads!(k, load_scens, mpc, load_mapping)

        D = build_loads(params, mpc["load"], mpc["shunt"])
        G = build_gens(params, mpc["gen"], cost_data, inst_name)

        sumD = sum(d for d in values(D))
        sum_lb = sum(g.lower_bound for g in values(G))
        sum_ub = sum(g.upper_bound for g in values(G))

        log(params, 
            "scen#$k: $sumD, $sum_lb, $sum_ub, $(sumD / sum_ub)", 
            true)
        if params.debugging_level == 1
            @assert isl(sum_lb, sumD) "sum lb gen $sum_lb > sum load $sumD"
            @assert isl(sumD, sum_ub) "sum load $sumD > sum ub gen $sum_ub"
        end

        push!(inst.scenarios, Scenario(prob, D, G))
    end

    if params.debugging_level == 1
        assert_feas_build_all(inst, params)
    end

    return inst
end
