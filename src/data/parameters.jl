const const_eps = 1e-5
const MAXINT = 2e9 # Gurobi MAXINT value
const const_infinite = 2.0e9
const mp_type_ref_bus = 3

abstract type AbstractSolutionStrategy end

struct Serial <: AbstractSolutionStrategy end

struct Deterministic <: AbstractSolutionStrategy end

struct Parallel <: AbstractSolutionStrategy end

Base.@kwdef mutable struct InstanceParameters
    seed::Int64 = 123 # Used for generating random number for the instance
    load_gen_mult::Float64 = 2.0 # Multiplier for the load and generation
    g_slack::Float64 = 0.15 # Generation slack with respect to the load
    cost_delta_mult::Float64 = 1e-2 # Value multiplied by the delta in cost
    num_candidates::Int64 = 2 # Number of candidates available per existing line
    # cost_mult::Float64 = 1e2 # Value multiplied by x to build the costs
    ref_bus::Int64 = 1 # Default reference bus used when none is found
    expected_lifetime::Int64 = 50 # Expected lifetime of the new lines in years
end

Base.@kwdef mutable struct StochasticInstanceParameters
    selection_percentile::Float64 = 0.8
    rounding_digits::Int64 = 2
end

Base.@kwdef mutable struct ModelParameters
    is_mip_en::Bool = true
    # penalty::Float64 = 100.0 # 1, 2, 3
    penalty::Float64 = 1.0
    # is_lp_model_s_var_set_req = true
    is_symmetry_en::Bool = false
    is_dcp_power_model_en::Bool = false # Build DCPPowerModel
    optimizer = Gurobi.Optimizer
end

Base.@kwdef mutable struct BinarySearchParameters
    is_en::Bool = true
    time_limit::Float64 = 600.0
    max_it::Int64 = 15 # 1.Calibrar 5, 10, 15: 10
    num_max_it_wo_impr::Int64 = 15 # 1.Calibrar 1, 3, 5: 1
end

Base.@kwdef mutable struct BeamSearchParameters
    time_limit::Float64 = 600.0
    num_children_per_parent::Int64 = 2 # w
    num_children_per_level::Int64 = 3 # N
    num_children_per_level_mult::Float64 = 0.5 # gamma
    candidates_per_batch_mult::Float64 = 5e-3 # 1.Cal. 0.25e-3, 5e-3, 1e-2: 5e-3
    num_max_it::Int64 = 2
    num_max_it_wo_impr::Int64 = 15 # 2.Calibrar 5, 10, 15
    is_shuffle_en::Bool = true
    restricted_list_ratio::Float64 = 1.0 # 1.Calibrar 0.5, 0.75, 1.0: 1.0
end

Base.@kwdef mutable struct ProgressiveHedgingParameters
    is_en::Bool = false
    time_limit::Float64 = 7200.0
    num_threads::Int64 = 49
    rho::Float64 = 1.0
    is_sep_rho_en::Bool = false # Calibrar true, false
    max_it::Int64 = 1000000
    convergence_eps::Float64 = 1e-3
    lb_threshold::Float64 = 0.25
    penalty_mult::Float64 = 1e6
end

Base.@kwdef mutable struct SolverParameters
    time_limit::Float64 = 7200.0
    num_threads::Int64 = 1
    log_level::Int64 = 0
end

Base.@kwdef mutable struct Parameters
    log_level::Int64 = 2
    log_dir::String = "logs"
    log_file::String = "log.txt"
    debugging_level::Int64 = 0
    solver::SolverParameters = SolverParameters()
    instance::InstanceParameters = InstanceParameters()
    stochastic_instance::StochasticInstanceParameters = 
                                                  StochasticInstanceParameters()
    model::ModelParameters = ModelParameters()
    binary_search::BinarySearchParameters = BinarySearchParameters()
    beam_search::BeamSearchParameters = BeamSearchParameters()
    progressive_hedging::ProgressiveHedgingParameters = 
                                                  ProgressiveHedgingParameters()
    # rng::MersenneTwister = Random.MersenneTwister(123)
    solution_strategy::AbstractSolutionStrategy = Serial()
end