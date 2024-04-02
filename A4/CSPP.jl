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

export IPSolver, MyCGSolver

const EPS = 1e-10

#* name: IPSolver
# inputs:
#   - data: the matrix containing the information of the network we are interested in
#   - Tmax: the maximum allowable time consumption (or any resource)
#   - origin: the node index where every path starts
#   - destination: the node index where every path ends
# outputs:
#   - cost_opt: optimal cost
#   - x_opt: optimal solution
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
    cost_opt, x_opt = objective_value(model), value.(x)
else
    println("Infeasible instance!")
    cost_opt, x_opt = nothing, nothing
end

path_opt = Set(findall(i->x_opt[i] > 1-EPS, 1:1:n_arcs))

return cost_opt, path_opt

end

#* name: solve_node!
function solve_node!(path_ind, costs, times, paths, one_arcs, zero_arcs, Tmax, start_node, end_node, arc_cost, arc_time, origin, destination, n_nodes, n_arcs)
    
## select paths satisfying the node conditions
path_ind_rmp = Int[]

if isempty(one_arcs)
    path_ind_rmp = copy(path_ind)
else
    for pInd in path_ind, arc in paths[pInd]
        if (arc in one_arcs) && !(pInd in path_ind_rmp)
            push!(path_ind_rmp, pInd)
        end
    end
end

# delete zero arcs
if !isempty(zero_arcs)
    for pInd in path_ind_rmp, arc in paths[pInd]
        if (arc in zero_arcs) && (pInd in path_ind_rmp)
            deleteat!(path_ind_rmp, findfirst(x->x==pInd, path_ind_rmp))
        end
    end
end

M = 100.0       # Big-M
while true
    # restricted master problem
    rmp = Model(HiGHS.Optimizer)

    @variable(rmp, y0 >= 0)                                   # artificial variable facilitating the solution procedure
    @variable(rmp, λ[path_ind_rmp] >= 0)                          # path variables
    @objective(rmp, Min, M*y0 + sum(costs[path_ind_rmp] .* λ))            # cost minimization
    @constraint(rmp, resource, sum(times[path_ind_rmp] .* λ) <= Tmax)       # resource constriaint, π1
    @constraint(rmp, convexity, y0 + sum(λ) == 1)             # convexity constraint, π0
    optimize!(rmp)                                            # solve the relaxaed RMP

    # results
    y0_val = value(y0)
    λ_val = value.(λ)
    obj_val = objective_value(rmp)
    π0 = dual(convexity)
    π1 = dual(resource)

    # if infeasible
    if y0_val > EPS
        π0 = M
    end

    # pricing sub-problem
    psp = Model(HiGHS.Optimizer)
    @variable(psp, x[i = 1:n_arcs], Bin)
    @objective(psp, Min, sum((arc_cost.-arc_time.*π1) .* x))
    @constraint(psp, flow_origin, sum(x[start_node .== origin]) - sum(x[end_node .== origin]) == 1)
    @constraint(psp, flow_balance[i = 1:n_nodes; (i != origin) && (i != destination)], sum(x[start_node .== i]) - sum(x[end_node .== i]) == 0)
    @constraint(psp, flow_destination, sum(x[start_node .== destination]) - sum(x[end_node .== destination]) == -1)
    if !isempty(one_arcs)
        @constraint(psp, node_ones[i = one_arcs], x[i] == 1)
    end
    if !isempty(zero_arcs)
        @constraint(psp, node_zeros[i = zero_arcs], x[i] == 0)
    end
    @constraint(psp, time, sum(arc_time .* x) <= Tmax)
    optimize!(psp)
    
    # no feasible path within node
    if primal_status(psp) != FEASIBLE_POINT
        λ_val_full = zeros(Float64, length(path_ind))
        λ_val_full[path_ind_rmp] = λ_val.data
        return λ_val_full, obj_val, y0_val, false
    end

    reduced_cost = objective_value(psp) - π0

    if y0_val > EPS     # RMP is infeasible
        x_opt = value.(x)
        new_path = Set(findall(i->x_opt[i]>(1-EPS), collect(1:1:n_arcs)))
        if new_path in paths[path_ind_rmp]    # no column to generate
            λ_val_full = zeros(Float64, length(path_ind))
            λ_val_full[path_ind_rmp] = λ_val.data
            return λ_val_full, obj_val, y0_val, true
        elseif reduced_cost < -EPS
            x_opt = value.(x)
            new_path = Set(findall(i->x_opt[i]>(1-EPS), collect(1:1:n_arcs)))
            new_cost = sum(arc_cost[collect(new_path)])
            new_time = sum(arc_time[collect(new_path)])
            # add column to the RMP
            push!(path_ind, length(path_ind)+1)
            push!(costs, new_cost)
            push!(times, new_time)
            push!(paths, new_path)
            push!(path_ind_rmp, length(path_ind))
        else
            M *= 10.0
        end
    else    # RMP is feasible
        if reduced_cost < -EPS      # find a new column
            x_opt = value.(x)
            new_path = Set(findall(i->x_opt[i]>(1-EPS), collect(1:1:n_arcs)))
            new_cost = sum(arc_cost[collect(new_path)])
            new_time = sum(arc_time[collect(new_path)])
            # add column to the RMP
            push!(path_ind, length(path_ind)+1)
            push!(costs, new_cost)
            push!(times, new_time)
            push!(paths, new_path)
            push!(path_ind_rmp, length(path_ind))
        else        # no column to generate
            λ_val_full = zeros(Float64, length(path_ind))
            λ_val_full[path_ind_rmp] = λ_val.data
            return λ_val_full, obj_val, y0_val, true
        end
    end
end

end

#* name: MyCGSolver
# inputs:
#   - data: the matrix containing the information of the network we are interested in
#   - Tmax: the maximum allowable time consumption (or any resource)
#   - origin: the node index where every path starts
#   - destination: the node index where every path ends
# outputs:
#   - cost_opt: optimal cost
#   - x_opt: optimal solution
function MyCGSolver(data, Tmax, origin, destination)

start_node = Int.(data[:, 1])
end_node = Int.(data[:, 2])
arc_cost = data[:, 3]
arc_time = data[:, 4]
n_nodes = maximum(vcat(start_node, end_node))
n_arcs = length(start_node)

## initialize
# containers for algorithms
pred_node = [0]
zero_arcs_node::Vector{Vector{Int}} = [[]]
one_arcs_node::Vector{Vector{Int}} = [[]]
state_node = [0] # 0: unsolved and not pruned out, 1: integer, 2: fractional, 3: infeasible, 4: pruned out
obj_val_node = [0.0]
path_ind_node = [Int[]]
costs_node = [Number[]]
times_node = [Number[]]
paths_node = [Set{Int}[]]

curr_node_ind = 1
best_node_ind = 0
best_path = Int[]
best_obj_val = Inf

# get an initial feasible solution
g = SimpleWeightedDiGraph(n_nodes)
for i in 1:1:n_arcs
    add_edge!(g, start_node[i], end_node[i], arc_cost[i])
end
bf_state = bellman_ford_shortest_paths(g, origin)
new_path = enumerate_paths(bf_state)[destination]
path_arcs = [findfirst(x->(start_node[x]==new_path[i]) && end_node[x]==new_path[i+1], 1:1:n_arcs) for i in 1:1:(length(new_path)-1)]
tpv = sum(data[path_arcs, 3:4], dims=1)
cost, time = tpv[1], tpv[2]

push!(path_ind_node[curr_node_ind], length(path_ind_node[curr_node_ind])+1)
push!(costs_node[curr_node_ind], cost)
push!(times_node[curr_node_ind], time)
push!(paths_node[curr_node_ind], Set(path_arcs))

# algorithm starts
while true
    # solve node by column generation
    λ_val, obj_val, y0_val, feasible = solve_node!(path_ind_node[curr_node_ind], costs_node[curr_node_ind], times_node[curr_node_ind], paths_node[curr_node_ind], one_arcs_node[curr_node_ind], zero_arcs_node[curr_node_ind], Tmax, start_node, end_node, arc_cost, arc_time, origin, destination, n_nodes, n_arcs)
    
    # update the current node status & best node
    if y0_val > EPS     # 3: infeasible
        state_node[curr_node_ind] = 3
        obj_val_node[curr_node_ind] = Inf
    elseif sum(λ_val .> 1-EPS) == 1     # 1: integer
        state_node[curr_node_ind] = 1
        obj_val_node[curr_node_ind] = obj_val
        if obj_val < best_obj_val # best solution update
            best_node_ind = curr_node_ind
            tpv = findfirst(i->λ_val[i]>1-EPS, path_ind_node[best_node_ind])
            best_path = paths_node[best_node_ind][tpv]
            best_obj_val = obj_val
        end
    elseif sum(λ_val .> 1-EPS) != 1     # fractional
        if best_obj_val - obj_val > EPS  # 2: not pruned out
            state_node[curr_node_ind] = 2
            obj_val_node[curr_node_ind] = obj_val
            ## branch, strategy: the first fractional arc
            if feasible
                # find the first fractional arc
                x_val = zeros(Float64, n_arcs)
                for pInd in path_ind_node[curr_node_ind], aInd in paths_node[curr_node_ind][pInd]
                    x_val[aInd] += λ_val[pInd]
                end
                arcs_for_nodeing = findall(i->(x_val[i] > EPS) && (x_val[i] - 1 < EPS), 1:1:n_arcs)
                arc_for_nodeing = arcs_for_nodeing[findfirst(a->!(a in one_arcs_node[curr_node_ind]) && !(a in zero_arcs_node[curr_node_ind]), arcs_for_nodeing)]
                # branch
                push!(pred_node, curr_node_ind); push!(pred_node, curr_node_ind);
                push!(zero_arcs_node, zero_arcs_node[curr_node_ind]); push!(zero_arcs_node, vcat(zero_arcs_node[curr_node_ind], arc_for_nodeing));
                push!(one_arcs_node, vcat(one_arcs_node[curr_node_ind], arc_for_nodeing)); push!(one_arcs_node, one_arcs_node[curr_node_ind]);
                push!(state_node, 0); push!(state_node, 0);
                push!(obj_val_node, 0.0); push!(obj_val_node, 0.0);
                push!(path_ind_node, path_ind_node[curr_node_ind]); push!(path_ind_node, path_ind_node[curr_node_ind]);
                push!(costs_node, costs_node[curr_node_ind]); push!(costs_node, costs_node[curr_node_ind]);
                push!(times_node, times_node[curr_node_ind]); push!(times_node, times_node[curr_node_ind]);
                push!(paths_node, paths_node[curr_node_ind]); push!(paths_node, paths_node[curr_node_ind]);
            end
        else    # 4: pruned out
            state_node[curr_node_ind] = 4
            obj_val_node[curr_node_ind] = obj_val
        end
    end

    # next node to solve
    curr_node_ind = 0
    for bInd in eachindex(state_node)
        if state_node[bInd] == 0
            curr_node_ind = bInd
            break
        end
    end

    if curr_node_ind == 0     # if no node to solve, terminate the algorithm
        break
    end
end

# solution history
solution_history = (pred_node, one_arcs_node, zero_arcs_node, state_node, obj_val_node)
best_solution = (best_node_ind, one_arcs_node[best_node_ind], zero_arcs_node[best_node_ind], best_obj_val, best_path)

return solution_history, best_solution

end

end