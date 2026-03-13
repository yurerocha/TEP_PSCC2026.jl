# ------------------------------ Type declarations -----------------------------
const FloatVarRef = Union{Float64, JuMP.VariableRef}
# const VectorSet{T <: Number} = Union{Vector{T}, Set{T}}
const VectorSet = Union{Vector{Int64}, Set{Int64}}
const AffQuadExpr = Union{JuMP.AffExpr, JuMP.QuadExpr}
const Tuple3I = Tuple{Int64, Int64, Int64}
const CandType = Tuple{Tuple3I, Int64}

# -------------------------- Heuristic data structures -------------------------
struct Start
    inserted
    f::Dict{Any, Float64}
    g::Dict{Int, Float64}
    # theta::Vector{Float64}
end

struct Status
    name::String
    rm_ratio::Float64
    impr_ratio::Float64
    time::Float64

    Status(name::String) = new(name, 0.0, 0.0, 0.0)
    Status(name::String, 
           rm_ratio::Float64, 
           impr_ratio::Float64, 
           time::Float64) = new(name, rm_ratio, impr_ratio, time)

    function Status(name::String, 
                    num_inserted::Int64, num_start::Int64, 
                    cost::Float64, start_cost::Float64, 
                    start_time::Float64)
        rm_rat = num_inserted / num_start
        impr_rat = 1.0 - cost / start_cost

        return Status(name, rm_rat, impr_rat, time() - start_time)
    end
end

# -------------------------- Instance data structures --------------------------
struct CostData
    voltage_classes::Vector{String}
    reactances_km::Dict{String, Float64}
    costs_km::Dict{String, Float64} # Costs in $/km-hr
    transformers::Dict{Tuple{String, String}, Float64} # Costs in $/MVA-hr
    gen_costs_mult::Dict{String, Float64} # Inflation adjustments for gen costs
end

mutable struct GeneratorInfo
    bus::Int64
    lower_bound::Float64
    upper_bound::Float64
    costs::Vector{Float64}
end

# struct Circuit
#     fr::Int64 # "from" bus
#     to::Int64 # "to" bus
# end

mutable struct Scenario
    # id::Int64
    p::Float64 # Probability
    D::Dict{Int64, Float64} # Load
    G::Dict{Int64, GeneratorInfo}
end

mutable struct BranchInfo
    f_bar::Float64 # Capacity of circuits
    x::Float64 # Reactance of circuits
    gamma::Float64 # Susceptance of circuits
    cost::Float64 # Cost of candidate circuits
    Dtheta_bounds::Tuple{Float64, Float64}
end

mutable struct Instance
    name::String
    I::Set{Int64} # Buses
    J::Dict{Tuple3I, BranchInfo} # Existing lines
    K::Dict{CandType, BranchInfo} # Candidates
    key_to_index::Dict{CandType, Int64}
    costs::Vector{Float64}
    # f_bar::Vector{Float64} # Capacity of circuits
    # gamma::Vector{Float64} # Susceptance of circuits
    # Dtheta_limits::Vector{Tuple{Float64, Float64}} 
    # costs::Vector{Float64} # Cost of candidate circuits
    num_I::Int64 # Number of buses
    num_J::Int64 # Number of existing circuits
    num_K::Int64 # Number of candidate circuits
    ref_bus::Int64
    scenarios::Vector{Scenario}
    num_scenarios::Int64
end

# ---------------------------- Model data structures ---------------------------
# abstract type AbstractTEPModel end
abstract type TEPModel end

struct MIPModel <: TEPModel
    jump_model::JuMP.Model
    obj::AffQuadExpr
    x::Dict{CandType, JuMP.VariableRef}
    f::Dict{Any, JuMP.VariableRef}
    g::Dict{Int64, JuMP.VariableRef}
    g_bus::Dict{Int64, JuMP.AffExpr} # Sum of g for the same bus
    theta::Dict{Int64, JuMP.VariableRef}
    Dtheta::Dict{Tuple{Int64, Int64}, JuMP.VariableRef}
    fkl_cons::Dict{Int64, JuMP.ConstraintRef}

    MIPModel(params::Parameters) = 
                new(# direct_model(Gurobi.Optimizer()), 
                    JuMP.Model(params.model.optimizer), 
                    params.model.is_dcp_power_model_en ? QuadExpr() : AffExpr(), 
                    Dict{CandType, JuMP.VariableRef}(), 
                    Dict{Any, JuMP.VariableRef}(), 
                    Dict{Int64, JuMP.VariableRef}(), 
                    Dict{Int64, JuMP.AffExpr}(), 
                    Dict{Int64, JuMP.VariableRef}(),
                    Dict{Int64, JuMP.VariableRef}(), 
                    Dict{Int64, JuMP.ConstraintRef}())
end

mutable struct LPModel <: TEPModel
    jump_model::JuMP.Model
    obj::AffExpr
    has_fixed_s_vars::Bool
    s::Dict{Any, JuMP.VariableRef}
    f::Dict{Any, JuMP.VariableRef}
    g::Dict{Int64, JuMP.VariableRef}
    g_bus::Dict{Int64, JuMP.AffExpr} # Sum of g for the same bus
    theta::Dict{Int64, JuMP.VariableRef}
    Dtheta::Dict{Tuple{Int64, Int64}, JuMP.VariableRef}
    fkl_cons::Dict{Int64, JuMP.ConstraintRef}
    f_cons::Dict{Any, JuMP.ConstraintRef}

    LPModel(params::Parameters) = new(JuMP.Model(params.model.optimizer), 
                                      AffExpr(), 
                                      false, 
                                      Dict{Any, JuMP.VariableRef}(), 
                                      Dict{Any, JuMP.VariableRef}(), 
                                      Dict{Int64, JuMP.VariableRef}(), 
                                      Dict{Int64, JuMP.AffExpr}(), 
                                      Dict{Int64, JuMP.VariableRef}(), 
                                      Dict{Int64, JuMP.VariableRef}(), 
                                      Dict{Int64, JuMP.ConstraintRef}(), 
                                      Dict{Any, JuMP.ConstraintRef}())
end

const TepModel = Union{MIPModel, LPModel}

# ------------------------- PDDF model data structures -------------------------
# abstract type AbstractPTDFModel <: AbstractTEPModel end

"""
    PTDFModel{T <: AbstractFloat} <: TEPModel

- T: used to represent floating-point data types involved in computations 
that are required to be efficient.
"""
mutable struct PTDFModel{T <: AbstractFloat} <: TEPModel
    jump_model::JuMP.Model
    obj::AffExpr
    bus_to_idx::Dict{Any, Int64} # Map buses' ids to indices
    line_to_idx::Dict{Any, Int64} 
    Gamma::SparseArrays.SparseMatrixCSC{Float64, Int64} # m x m susceptances
    S::SparseArrays.SparseMatrixCSC{T, Int64} # m x n adjacency matrix
    d::Vector{Float64} # n vector of demands
    g::Dict{Int64, JuMP.VariableRef} 
    # g_bus::SparseArrays.SparseVector{JuMP.AffExpr, Int64} 
    g_bus::Vector{JuMP.AffExpr} 
    B::SparseArrays.SparseMatrixCSC{Float64, Int64} # n x n mat, where B = S'ΓS
    B_inv::Matrix{Float64} # n x n inverse of B
    # buffer::Matrix{T} # n x n
    I::Matrix{T} # n x n identity matrix
    beta::Matrix{T} # m x n, where β = ΓSB⁻¹
    bus_inj::Vector{JuMP.AffExpr}
    f::Vector{JuMP.AffExpr} # m x 1 vec of line flows
    s::Vector{JuMP.VariableRef} # m vector of slack variables
     # m x 1 vector of line flow conss
    f_neg_cons::Vector{JuMP.ConstraintRef}
    f_pos_cons::Vector{JuMP.ConstraintRef}
end

mutable struct PTDFSystem{T <: AbstractFloat}
    bus_to_idx::Dict{Any, Int64} # Map buses' ids to indices
    line_to_idx::Dict{Any, Int64}
    Gamma::SparseArrays.SparseMatrixCSC{Float64, Int64} # m x m susceptances
    S::SparseArrays.SparseMatrixCSC{T, Int64} # m x n adjacency matrix
    d::Vector{Float64} # n vector of demands
    g_bus::Vector{Float64}
    B::SparseArrays.SparseMatrixCSC{Float64, Int64} # n x n mat, where B = S'ΓS
    B_inv::Matrix{Float64} # n x n inverse of B
    # buffer::Matrix{T} # n x n
    I::Matrix{T} # n x n identity matrix
    beta::Matrix{T} # m x n, where β = ΓSB⁻¹
    bus_inj::Vector{Float64}
end

# ----------------------------- PH data structures -----------------------------

mutable struct ScenarioCache
    omega::Vector{Float64}
    state::Vector{Float64}
end

mutable struct Solution
    insert::Set{CandType}
    feas_insert::Set{CandType}
    count_use::Int64
    costs::Vector{Float64}
    g_costs::Vector{Float64}
    viols::Vector{Float64}

    Solution(inst::Instance) = new(Set{CandType}(keys(inst.K)), 
                                   Set{CandType}(keys(inst.K)), 
                                   0, 
                                   Vector{Float64}(undef, inst.num_scenarios), 
                                   Vector{Float64}(undef, inst.num_scenarios), 
                                   Vector{Float64}(undef, inst.num_scenarios))
end

# mutable struct ExecutionInfo
#     num_trials::Int64
#     num_successes::Int64
# end

mutable struct ScenarioStatus
    beam_search_runtime::Float64
    solver_runtime::Float64
    bin_search_rm_ratio::Float64
    beam_search_rm_ratio::Float64
    repair::Tuple{Int64, Int64} # applied, succeeded
    reinsert::Tuple{Int64, Int64, Int64} # applied, succeeded, iterations
end

mutable struct Cache
    it::Int64
    scenarios::Vector{ScenarioCache}
    x_hat::Vector{Float64}
    x_average::Vector{Float64}
    sol_lb::Solution
    sol_ub::Solution
    start_sol::Set{CandType}
    start_costs::Vector{Float64}
    # SEP rho heuristic by https://doi.org/10.1007/s10287-010-0125-4
    rho::Vector{Float64}
    sep_rho_x_min::Vector{Float64}
    sep_rho_x_max::Vector{Float64}
    deltas::Vector{Float64}
    best_convergence_delta::Float64
    t_deviation::Float64
    quality_deviation::Float64
    best_it::Int64
    best_sol::Set{CandType}
    best_cost::Float64
    hash_weights::Vector{Int64}
    hash_values::Vector{Float64}
    fixed_x_variables::Set{CandType}
    count_cycle_it::Vector{Int64}
    status::Vector{ScenarioStatus}

    Cache(inst::Instance, params::Parameters) = 
        new(0, 
            [ScenarioCache(zeros(inst.num_K), 
                           Vector{Float64}()) for _ in 1:inst.num_scenarios],
            # [zeros(num_vars) for _ in 1:num_scenarios], 
            # Vector{Vector{Float64}}(undef, num_vars), 
            # Vector{Vector{Float64}}(undef, num_vars), 
            Vector{Float64}(), 
            Vector{Float64}(), 
            # Set{CandType}(), 
            # Set{CandType}(),
            # 0, 0, 
            # Vector{Float64}(undef, inst.num_scenarios), 
            # Vector{Float64}(undef, inst.num_scenarios), 
            Solution(inst), Solution(inst), 
            Set{CandType}(keys(inst.K)), 
            [const_infinite for _ in eachindex(inst.scenarios)], 
            [params.progressive_hedging.rho / 2.0 for _ in eachindex(inst.K)], 
            Vector{Float64}(), 
            Vector{Float64}(), 
            Vector{Float64}(undef, inst.num_scenarios), 
            const_infinite, const_infinite, 100.0, 0, 
            Set{CandType}(), 
            const_infinite, 
            # Generate n numbers from 1:100 * n, randomly
            Distributions.sample(1:100 * inst.num_scenarios, 
                                    inst.num_scenarios, replace = false), 
            Vector{Float64}(undef, inst.num_K), 
            Set{CandType}(), 
            zeros(Int64, inst.num_K), 
            [ScenarioStatus(0.0, 0.0, 0.0, 0.0, (false, false), 
                (false, false, 0)) for _ in eachindex(inst.scenarios)]) 
end

@enum WorkerOption repair_sols comp_g_costs run_method

mutable struct WorkerCache
    option::WorkerOption
    scenarios::Vector{ScenarioCache}
    x_hat::Vector{Float64}
    sol_lb::Set{CandType}
    sol_ub::Set{CandType}
    inserted::Set{CandType}
    costs::Vector{Float64}
    rho::Vector{Float64}

    # Construct the struct according to the worker option
    function WorkerCache(cache::Cache, option::WorkerOption)
        if option == repair_sols
            return new(option, cache.scenarios, cache.x_hat, 
                        cache.sol_lb.insert, cache.sol_ub.insert, 
                        cache.start_sol, Vector{Float64}(), cache.rho) 
        elseif option == comp_g_costs
            return new(option, cache.scenarios, cache.x_hat, 
                        cache.sol_lb.insert, cache.sol_ub.insert, 
                        Set{CandType}(), Vector{Float64}(), cache.rho) 
        elseif option == run_method
            return new(option, cache.scenarios, cache.x_hat, 
                        Set{CandType}(), Set{CandType}(), 
                        cache.start_sol, cache.start_costs, cache.rho) 
        end
    end
end

mutable struct ControllerMessage
    cache::WorkerCache
    it::Int64
    scen::Int64
    time_limit::Float64
end

mutable struct SolutionInfo
    cost::Float64
    g_cost::Float64
    viol::Float64
    reinsert::Set{CandType}

    SolutionInfo(cost, g_cost, viol, reinsert) = 
                                            new(cost, g_cost, viol, reinsert)
end

mutable struct WorkerMessage
    state_values::Vector{Float64}
    it::Int64
    scen::Int64
    sol_info_lb::SolutionInfo
    sol_info_ub::SolutionInfo
    status::ScenarioStatus
end

# ------------------------- Beam Search data structures ------------------------
mutable struct Node # {T <: AbstractFloat}
    cost::Float64 # Build + generation
    viol::Float64
    inserted::Set{CandType}
    removed::Set{CandType}
    ignore::Set{CandType} # Cand lines whose removal leads to infeasibls sols
end

mutable struct BSControllerMessage
    it::Int64
    node_idx::Int64
    node::Node
    best_cost::Float64 # Build + generation
    k::Vector{Any}
    num_threads::Int64
end

mutable struct BSWorkerMessage
    node_idx::Int64
    lines::Set{CandType}
    is_feas::Bool
    cost::Float64 # Build + generation
    viol::Float64
end