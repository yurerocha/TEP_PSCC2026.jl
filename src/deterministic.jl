"""
    build_deterministic(inst::Instance, params::Parameters)

Build the deterministic equivalent model and return both problem and 
subproblems.
"""
function build_deterministic(inst::Instance, params::Parameters)
    mip = MIPModel(params)
    subproblems = Vector{MIPModel}(undef, inst.num_scenarios)

    for scen in 1:inst.num_scenarios
        # TODO: Change LP objective as well
        # TODO: Run heuristic in every it
        subproblem = build_mip(inst, params, scen, false)
        set_state!(inst, subproblem)
        add_subproblem!(inst, mip.jump_model, scen, subproblem.jump_model)
        subproblems[scen] = subproblem
    end
    add_non_anticipativity_cons!(inst, mip.jump_model, subproblems)
    add_obj_build_costs!(inst, mip.jump_model, subproblems)

    return mip, subproblems
end