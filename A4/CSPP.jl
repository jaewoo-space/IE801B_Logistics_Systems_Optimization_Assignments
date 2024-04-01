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
using Graphs, SimpleWeightedGraphs

export IPSolver, cal_cost_time, solve_rmp_and_pricing_problems

const EPS = 1e-10

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

function cal_cost_time(data, n_arcs, start_node, end_node, path)
    path_arcs = [findfirst(x->(start_node[x]==path[i]) && end_node[x]==path[i+1], 1:1:n_arcs) for i in 1:1:(length(path)-1)]
    tpv = sum(data[path_arcs, 3:4], dims=1)

    return path_arcs, tpv[1], tpv[2] # cost and time    
end

function solve_rmp_and_pricing_problems(path_ind, costs, times, paths, data, Tmax, start_node, end_node, arc_cost, arc_time, origin, destination, n_nodes, n_arcs)
    M = 100.0
    
    while true
        model = Model(HiGHS.Optimizer)

        @variable(model, y0 >= 0)                                   # artificial variable facilitating the solution procedure
        @variable(model, λ[path_ind] >= 0)                          # path variables
        @objective(model, Min, M*y0 + sum(costs .* λ))            # cost minimization
        @constraint(model, resource, sum(times .* λ) <= Tmax)       # resource constriaint, π1
        @constraint(model, convexity, y0 + sum(λ) == 1)             # convexity constraint, π0
        optimize!(model)                                            # solve the relaxaed RMP

        # results
        y0_val = value(y0)                                      
        λ_val = value.(λ)
        obj_val = objective_value(model)
        π0 = dual(convexity)
        π1 = dual(resource)

        if y0_val > EPS
            π0 = M
        end
        while true
            g = SimpleWeightedDiGraph(n_nodes)
            for i in 1:1:n_arcs
                add_edge!(g, start_node[i], end_node[i], arc_cost[i] - arc_time[i] * π1)
            end
            
            bf_state = bellman_ford_shortest_paths(g, origin)
            min_reduced_cost = bf_state.dists[destination] - π0
            println("min_reduced_cost:$(min_reduced_cost)")
            if min_reduced_cost < -EPS          # find a new column might improve the solution
                new_column = enumerate_paths(bf_state)[destination]
                path_arcs, cost, time = cal_cost_time(data, n_arcs, start_node, end_node, new_column)
                push!(path_ind, length(path_ind)+1)
                push!(costs, cost)
                push!(times, time)
                push!(paths, path_arcs)
                break
            elseif y0_val > EPS                # cannot find a new column, but the solution is infeasible
                π0 *= 10.0
            else
                return path_ind, costs, times, paths, λ_val, obj_val
            end
        end
    end
end

# TODO
function MyCGSolver(data, Tmax, origin, destination)
    start_node = Int.(data[:, 1])
    end_node = Int.(data[:, 2])
    arc_cost = data[:, 3]
    arc_time = data[:, 4]
    n_nodes = maximum(vcat(start_node, end_node))
    n_arcs = length(start_node)

    ## initialize
    path_ind = Int[]
    costs = Number[]
    times = Number[]
    paths = Vector{Int}[]

    # get an initial feasible solution
    g = SimpleWeightedDiGraph(n_nodes)
    for i in 1:1:n_arcs
        add_edge!(g, start_node[i], end_node[i], arc_cost[i])
    end
    bf_state = bellman_ford_shortest_paths(g, origin)
    new_column = enumerate_paths(bf_state)[destination]
    path_arcs, cost, time = cal_cost_time(data, n_arcs, start_node, end_node, new_column)
    push!(path_ind, length(path_ind)+1)
    push!(costs, cost)
    push!(times, time)
    push!(paths, path_arcs)

    path_ind, costs, times, paths, λ_val, obj_val = solve_rmp_and_pricing_problems(path_ind, costs, times, paths, data, Tmax, start_node, end_node, arc_cost, arc_time, origin, destination, n_nodes, n_arcs)

    return nothing
end

end