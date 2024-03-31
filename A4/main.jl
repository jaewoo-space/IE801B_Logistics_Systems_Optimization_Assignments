using HiGHS, JuMP, DelimitedFiles
include("./CSPP.jl")
using .Tools: CSPP_INSTANCE
using .Solver: IPSolver
using JLD
using Graphs, SimpleWeightedGraphs

function solve_pricing_problem(start_node, end_node, arc_cost, arc_time, π0, π1, origin, destination)
    g = SimpleWeightedDiGraph(n_nodes)
    for i in 1:1:n_arcs
        add_edge!(g, start_node[i], end_node[i], arc_cost[i] - arc_time[i] * π1)
    end
    
    bf_state = bellman_ford_shortest_paths(g, origin)
    min_reduced_cost = bf_state.dists[destination] - π0
    println("min_reduced_cost:$(min_reduced_cost)")
    if min_reduced_cost < -(1e-10)
        new_column = enumerate_paths(bf_state)[destination]
    else
        return nothing
    end 
end

function cal_cost_time(data, n_arcs, start_node, end_node, path)
    path_arcs = [findfirst(x->(start_node[x]==path[i]) && end_node[x]==path[i+1], 1:1:n_arcs) for i in 1:1:(length(path)-1)]
    tpv = sum(data[path_arcs, 3:4], dims=1)

    return path_arcs, tpv[1], tpv[2] # cost and time    
end

function solve_rmp(paths, costs, times, path_info)
    m = Model(HiGHS.Optimizer)

    @variable(m, y0 >= 0)
    @variable(m, λ[paths] >= 0)
    @objective(m, Min, 100*y0 + sum(costs .* λ))
    @constraint(m, resource, sum(times .* λ) <= Tmax)
    @constraint(m, convexity, y0 + sum(λ) == 1)
    optimize!(m)
    if primal_status(m) == FEASIBLE_POINT
        return value(y0), value.(λ), objective_value(m), dual(resource), dual(convexity)
    else
        return nothing
    end
end

prob = load("prob1.jld")["prob"]

data = prob.data
Tmax = prob.Tmax

start_node = Int.(data[:, 1])
end_node = Int.(data[:, 2])
arc_cost = data[:, 3]
arc_time = data[:, 4]
origin, destination = prob.origin, prob.destination
n_nodes = maximum(vcat(start_node, end_node))
n_arcs = length(start_node)


## initialize
paths = Int[1]
costs = Number[]
times = Number[]
path_info = Vector{Int}[]
M = 100.0

y0 = 0.0
λ = Float64[]
π0 = M
π1 = 0.0
z = Inf

new_column = solve_pricing_problem(start_node, end_node, arc_cost, arc_time, π0, π1, origin, destination)
path_arcs, cost, time = cal_cost_time(data, n_arcs, start_node, end_node, new_column)

push!(costs, cost)
push!(times, time)
push!(path_info, path_arcs)

function solve_node!(paths, costs, times, path_info)
    y0 = 0.0
    λ = Float64[]
    π1 = 100.0
    π0 = 0.0
    z = Inf

    while true
        master_sol = solve_rmp(paths, costs, times, path_info)
        
        if isnothing(master_sol)
            break
        end
        
        y0 = master_sol[1]
        λ = master_sol[2]
        z = master_sol[3]
        if y0 > 1e-5
            π0 = 100.0
        else
            π0 = master_sol[5]
        end
        π1 = master_sol[4]
        new_column = solve_pricing_problem(start_node, end_node, arc_cost, arc_time, π0, π1, origin, destination)
        println("new column: $(new_column)")
        
        if isnothing(new_column)
            break
        end

        path_arcs, cost, time = cal_cost_time(data, n_arcs, start_node, end_node, new_column)
        push!(paths, paths[end]+1)
        push!(costs, cost)
        push!(times, time)
        push!(path_info, path_arcs)
    end
    return y0, λ, π1, π0, z
end

y0, λ, π1, π0, z = solve_node!(paths, costs, times, path_info)

# TODO
# if λ is fractional, then start bnb
# if not, end of the branching