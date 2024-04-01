using HiGHS, JuMP, DelimitedFiles
include("./CSPP.jl")
using .Tools: CSPP_INSTANCE
using .Solver: IPSolver, solve_rmp_and_pricing_problems, cal_cost_time
using JLD
using Graphs, SimpleWeightedGraphs

prob = load("prob1.jld")["prob"]

const data = prob.data
const Tmax = prob.Tmax

const start_node = Int.(data[:, 1])
const end_node = Int.(data[:, 2])
const arc_cost = data[:, 3]
const arc_time = data[:, 4]
const origin, destination = prob.origin, prob.destination
const n_nodes = maximum(vcat(start_node, end_node))
const n_arcs = length(start_node)

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

# TODO
# if λ is fractional, branch
# if not, end of the branch