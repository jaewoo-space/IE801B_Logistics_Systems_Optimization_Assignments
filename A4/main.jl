using HiGHS, JuMP, DelimitedFiles
include("./CSPP.jl")
using .Tools: CSPP_INSTANCE
using .Solver: IPSolver
using JLD
using Graphs, SimpleWeightedGraphs

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
eps = 1e-10
paths = Int[1]
costs = Number[]
times = Number[]
M = 100.0
g = SimpleWeightedDiGraph(n_nodes)

π0 = M
π1 = 0.0

for i in 1:1:n_arcs
    add_edge!(g, start_node[i], end_node[i], arc_cost[i] - arc_time[i] * π1)
end

bf_state = bellman_ford_shortest_paths(g, origin)
min_reduced_cost = bf_state.dists[destination] - π0
new_column = enumerate_paths(bf_state)[destination]

function cal_cost(data, path)
    return 0
end

function cal_time(data, path)
    return 0
end

function solve_rmp!(paths, costs, times)

m = Model(HiGHS.Optimizer)

@variable(m, M*y0 >= 0)
@objective(m, Min, sum(costs .* λ))
@constraint(m, resource, sum(times .* λ) <= Tmax)
@constraint(m, convexity, sum(λ) == 1)
optimize!(m)

return 0

end
# paths, costs, times 


push!(costs, cal_cost(data, new_column))
push!(times, cal_cost(data, new_column))


## loop starts
# while true

# solve the RMP... if the solution is fractional, then BB
m = solve_rmp!()

# solve the pricing sub-problem
π0 = dual(convexity)
π1 = dual(resource)

for i in 1:1:n_arcs
    add_edge!(g, start_node[i], end_node[i], arc_cost[i] - arc_time[i] * π1)
end

bf_state = bellman_ford_shortest_paths(g, origin)
min_reduced_cost = bf_state.dists[destination] - π0
# terminal condition
# if min_reduced_cost > -eps
#     break
# end  
new_column = enumerate_paths(bf_state)[destination]

# paths, costs, times update
push!(paths, paths[end]+1)
push!(costs, cal_cost(data, new_column))
push!(times, cal_time(data, new_column))

# back to the RMP

# end



