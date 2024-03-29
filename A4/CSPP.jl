module Tools

export CSPP_INSTANCE

struct CSPP_INSTANCE

data::Array{Number}
Tmax::Number
origin::Int
destination::Int

end

end

module Solver

using JuMP, HiGHS

export IPSolver

function IPSolver(data, Tmax, origin, destination)
    start_node = Int.(data[:, 1])
    end_node = Int.(data[:, 2])
    arc_cost = data[:, 3]
    arc_time = data[:, 4]
    n_nodes = maximum(vcat(start_node, end_node))
    n_arcs = length(start_node)

    model = Model(HiGHS.Optimizer)

    @variable(model, x[i = 1:n_arcs], Bin)

    @objective(model, Min, sum(arc_cost .* x))

    @constraint(model, flow_origin, sum(x[start_node .== origin]) - sum(x[end_node .== origin]) == 1)
    @constraint(model, flow_balance[i = 1:n_nodes; (i != origin) && (i != destination)], sum(x[start_node .== i]) - sum(x[end_node .== i]) == 0)
    @constraint(model, flow_destination, sum(x[start_node .== destination]) - sum(x[end_node .== destination]) == -1)

    @constraint(model, time, sum(arc_time .* x) <= Tmax)

    optimize!(model)

    if primal_status(model) == FEASIBLE_POINT
        return objective_value(model), value.(x)
    else
        println("Infeasible instance!")
        return nothing
    end
end


function MyCGSolver(data, Tmax, origin, destination)
    start_node = Int.(data[:, 1])
    end_node = Int.(data[:, 2])
    arc_cost = data[:, 3]
    arc_time = data[:, 4]
    n_nodes = maximum(vcat(start_node, end_node))
    n_arcs = length(start_node)

    paths = []
    costs = []
    times = []

    m = Model(HiGHS.Optimizer)

    @variable(m, λ[paths] >= 0)
    @objective(m, Min, sum(costs .* λ))
    @constraint(m, resource, sum(times .* λ) <= Tmax)
    @constraint(m, convexity, sum(λ) == 1)
    optimize!(m)
end

end