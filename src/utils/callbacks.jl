"""
    barrier_callback(cb_data, cb_where::Cint)

Barrier callback to terminate Gurobi if the violation is worse than the best
violation found.

It is required to set jump_model.ext[:best_viol].
"""
function barrier_callback(cb_data, cb_where::Cint)
    if cb_where == GRB_CB_BARRIER
        dual_obj = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_BARRIER_DUALOBJ, dual_obj)
        if isg(dual_obj[], ext[:best_viol])
            @warn dual_obj, ext[:best_viol], "Terminate Gurobi"
            # readline()
            GRBterminate(backend(lp.jump_model).optimizer.model.inner)
        end
    end
end

"""
    empty_callback(cb_data, cb_where::Cint)

It can be used to remove the barrier callback.
"""
function empty_callback(cb_data, cb_where::Cint)
    return nothing
end

"""
    root_node_callback(cb_data, cb_where::Cint)

Get the root node best bound, the root node runtime, and the time to 
find the first incumbent solution.
"""
function root_node_callback(cb_data, cb_where::Cint)
    if cb_where == GRB_CB_MIPNODE
        node_count = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, node_count)
        if iseq(node_count[], 0.0)
            runtime = Ref{Cdouble}()
            root_bound = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, runtime)
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBND, root_bound)

            ext[:rt_runtime] = runtime[]
            ext[:rt_best_bound] = root_bound[]
        end
    elseif cb_where == GRB_CB_MIPSOL && 
           iseq(ext[:incumbent_time], params.solver.time_limit)
        # Prevents ext[:incumbent_time] from being updated after the first 
        # incumbent solution is found
        runtime = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, runtime)
        ext[:incumbent_time] = runtime[]
    end
end