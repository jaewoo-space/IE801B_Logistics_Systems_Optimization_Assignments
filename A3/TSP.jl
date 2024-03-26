module Tools

using StatsBase: sample
using LinearAlgebra: norm
using Plots

export RandomProblemGenerator, VisualizeTour, VisualizeIntermediateSolution

#* name: RandomProblemGenerator
# inputs:
#   - sitesN: the number of sites
# outputs:
#   - sites: coordinates of the sites
function RandomProblemGenerator(sitesN::Int64)
    sites = rand(sitesN, 2)
    return sites
end

#* name: VisualizeTour
# inputs:
#   - sites: the coordinates of sites in 2D space
#   - tour: the tour the user want to visualize, a permutation of all the sites' indices
#   - title: title of the visual
# outputs:
#   - fig: visual plot
function VisualizeTour(sites, tour, title = "")
    sitesN = size(sites)[1]
    
    tour_plot = vcat(tour, tour[1])

    fig = plot(dpi=300)
    for i in 1:1:sitesN
        plot!([sites[tour_plot[i], 1], sites[tour_plot[i+1], 1]], [sites[tour_plot[i], 2], sites[tour_plot[i+1], 2]]; legend = false)
    end
    xaxis!([-0.15, 1.15])
    xlabel!("X")
    yaxis!([-0.15, 1.15])
    ylabel!("Y")
    xticks!([0.0, 0.5, 1.0])
    yticks!([0.0, 0.5, 1.0])
    scatter!(sites[:, 1], sites[:, 2], color=:black, markersize=3)
    title!(title)
    return fig
end

#* name: VisualizeLazyConstraints
# inputs:
#   - sites: the coordinates of sites in 2D space
#   - x_hist: a sequence of IP solutions obtained over the lazy constraints solution procedures
#   - title: title of the visual
# outputs:
#   - fig_list: list of visual plots
function VisualizeLazyConstraints(sites, x_hist, title = "")
    sitesN = size(sites)[1]
    
    fig_list = []

    for hInd in eachindex(x_hist)
        x_val = x_hist[hInd]
        selected_edges = [(i, findmax(x_val[i, :])[2]) for i in 1:1:sitesN]
        fig = plot(dpi=300, legend=false)
        for edge in selected_edges
            plot!([sites[edge[1], 1], sites[edge[2], 1]], [sites[edge[1], 2], sites[edge[2], 2]]; legend = false)
        end
        xaxis!([-0.15, 1.15])
        xlabel!("X")
        yaxis!([-0.15, 1.15])
        ylabel!("Y")
        xticks!([0.0, 0.5, 1.0])
        yticks!([0.0, 0.5, 1.0])
        scatter!(sites[:, 1], sites[:, 2], color=:black, markersize=3)
        title!(title)

        push!(fig_list, fig)
    end

    return fig_list
end

end

module ExactSolvers

using Gurobi, JuMP

export GurobiSolver, GurobiSolverWithLazyConstraints, GurobiSolverWithLazyConstraintsForVisuals

#* name: GurobiSolver
# inputs:
#   - C: {C_ij} is the length between the sites i and j
# outputs:
#   - opt_tour: optimal tour
#   - opt_cost: optimal tour length
#   - comp_time: solution time

function GurobiSolver(C)
    sitesN = size(C)[1]
    x_opt = Matrix{Bool}(undef, sitesN, sitesN)
    
    model = JuMP.Model(Gurobi.Optimizer)
    set_attribute(model, "OutputFlag", 0)

    @variables(model, begin
        x[1:1:sitesN, 1:1:sitesN], Bin
        t[i=2:1:sitesN] >= 0
    end)

    @constraints(model, begin
        outflow[i=1:1:sitesN], sum(x[i, :]) == 1
        inflow[j=1:1:sitesN], sum(x[:, j]) == 1
        selfloop[i=1:1:sitesN], sum(x[i, i]) == 0
        subtour[(i, j) = [(i, j) for i in 2:1:sitesN for j in 2:1:sitesN if i!=j]], t[j] - t[i] - sitesN*x[i, j] + (sitesN-1) >= 0
    end
    )

    @objective(model, Min, sum(C .* x))

    optimize!(model)

    x_opt = value.([x[i, j] for i in 1:1:sitesN, j in 1:1:sitesN])
    t_opt = value.([t[i] for i in 2:1:sitesN])
    opt_tour = sortperm(t_opt).+1; pushfirst!(opt_tour, 1)
    opt_cost = objective_value(model)

    comp_time = solve_time(model)

    return opt_tour, opt_cost, comp_time
end

#* name: GurobiSolverWithLazyConstraints
# inputs:
#   - C: {C_ij} is the length between the sites i and j
# outputs:
#   - opt_tour: optimal tour
#   - opt_cost: optimal tour length
#   - comp_time: solution time
function GurobiSolverWithLazyConstraints(C)
    sitesN = size(C)[1]
    
    model = JuMP.Model(Gurobi.Optimizer)
    set_attribute(model, "OutputFlag", 0)
    
    @variable(model, x[1:1:sitesN, 1:1:sitesN], Bin)

    @constraints(model, begin
        outflow[i=1:1:sitesN], sum(x[i, :]) == 1
        inflow[j=1:1:sitesN], sum(x[:, j]) == 1
        selfloop[i=1:1:sitesN], sum(x[i, i]) == 0
    end
    )

    @objective(model, Min, sum(C .* x))

    function find_subtour(x_val, s=1)
        subtour = [s]
        while true
            _, next_site = findmax(x_val[subtour[end], :])
            if next_site == s
                break
            else
                push!(subtour, next_site)
            end
        end
        return subtour
    end

    function call_back_function(cb_data)
        status = callback_node_status(cb_data, model)
        
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end
        
        x_val = (y->callback_value(cb_data, y)).(x).data
        
        subtour = find_subtour(x_val)
        S = length(subtour)
        if S > sitesN-1
            return
        end
        
        ex = AffExpr()
        for i in subtour, j in subtour
            add_to_expression!(ex, 1, x[i, j])
        end

        con = @build_constraint(ex <= S-1)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        
        return
    end

    MOI.set(model, MOI.LazyConstraintCallback(), call_back_function)

    optimize!(model)

    x_opt = value.([x[i, j] for i in 1:1:sitesN, j in 1:1:sitesN])
    opt_tour = find_subtour(x_opt)
    opt_cost = objective_value(model)
    comp_time = solve_time(model)
    
    return opt_tour, opt_cost, comp_time
end

#* name: GurobiSolverWithLazyConstraintsForGIF
# inputs:
#   - C: {C_ij} is the length between the sites i and j
# outputs:
#   - x_hist: sequence of IP solutions obtained over the solution procedures
function GurobiSolverWithLazyConstraintsForGIF(C)
    sitesN = size(C)[1]
    
    model = JuMP.Model(Gurobi.Optimizer)
    set_attribute(model, "OutputFlag", 0)
    
    x_hist = []
    
    @variable(model, x[1:1:sitesN, 1:1:sitesN], Bin)

    @constraints(model, begin
        outflow[i=1:1:sitesN], sum(x[i, :]) == 1
        inflow[j=1:1:sitesN], sum(x[:, j]) == 1
        selfloop[i=1:1:sitesN], sum(x[i, i]) == 0
    end
    )

    @objective(model, Min, sum(C .* x))

    function find_subtour(x_val, s=1)
        subtour = [s]
        while true
            _, next_site = findmax(x_val[subtour[end], :])
            if next_site == s
                break
            else
                push!(subtour, next_site)
            end
        end
        return subtour
    end

    function call_back_function(cb_data)
        status = callback_node_status(cb_data, model)
        
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end
        
        x_val = (y->callback_value(cb_data, y)).(x).data
        
        subtour = find_subtour(x_val)
        S = length(subtour)
        if S > sitesN-1
            return
        end
        
        push!(x_hist, x_val)
        ex = AffExpr()
        for i in subtour, j in subtour
            add_to_expression!(ex, 1, x[i, j])
        end
        
        push!(x_hist, x_val)
        
        con = @build_constraint(ex <= S-1)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        
        return
    end

    MOI.set(model, MOI.LazyConstraintCallback(), call_back_function)

    optimize!(model)

    push!(x_hist, value.([x[i, j] for i in 1:1:sitesN, j in 1:1:sitesN]))

    return x_hist
end

end

module Heuristics

const eps = 1e-5

using LinearAlgebra

export NearestNeighbor, Greedy, FarthestInsertion, Christofides, TwoOptSwapSlow, TwoOptSwap

#* name: NearestNeighbor
# inputs:
#   - C: {C_ij} is the length between the sites i and j
# outputs:
#   - opt_tour: optimal tour
#   - opt_cost: optimal tour length
function NearestNeighbor(C)
    sitesN = size(C)[1]
    C_tpv = copy(C); C_tpv[diagind(C_tpv)] .= Inf

    opt_tour = [rand(collect(1:1:sitesN))]
    remaining_sites = collect(setdiff(Set(collect(1:1:sitesN)), Set(opt_tour)))
    for _ in 1:1:(sitesN-1)
        next_site = remaining_sites[argmin(C_tpv[opt_tour[end], remaining_sites])]
        push!(opt_tour, next_site)
        remaining_sites = collect(setdiff(Set(remaining_sites), Set([next_site])))
    end
    opt_cost = sum([C[opt_tour[i], opt_tour[i+1]] for i in 1:1:(sitesN-1)]) + C[opt_tour[end], opt_tour[1]]
    return opt_tour, opt_cost
end

#* name: Greedy
# inputs:
#   - C: {C_ij} is the length between the sites i and j
# outputs:
#   - opt_tour: optimal tour
#   - opt_cost: optimal tour length
function Greedy(C)

    function is_subtour(connected_sites, site_cnt, new_edge)
        i = new_edge[1]
        j = new_edge[2]
        if (site_cnt[i] != 1) || (site_cnt[j] != 1)
            return false
        else
            subtour = [i, connected_sites[i][1]]
            while true
                tpv1 = connected_sites[subtour[end]]
                if tpv1[2] == 0
                    break
                else
                    if tpv1[1] == subtour[end-1]
                        push!(subtour, tpv1[2])
                    else
                        push!(subtour, tpv1[1])
                    end
                end
            end

            if subtour[end] == j
                return true
            else
                return false
            end
        end
    end

    sitesN = size(C)[1]

    edges = Tuple{Int, Int}[(i, j) for i in 1:1:sitesN for j in 1:1:sitesN if i<j]
    edges_cost = [C[edge[1], edge[2]] for edge in edges]
    sorted_arg = sortperm(edges_cost)

    connected_sites = [zeros(Int, 2) for sInd in 1:1:sitesN]
    site_cnt = zeros(Int, sitesN)

    opt_cost = 0.0
    for _ in 1:1:(sitesN-1)
        while true
            eInd = sorted_arg[1]
            if (site_cnt[edges[eInd][1]] > 1) || (site_cnt[edges[eInd][2]] > 1) || is_subtour(connected_sites, site_cnt, edges[eInd])
                popfirst!(sorted_arg)
            else
                opt_cost += edges_cost[sorted_arg[1]]
                i = edges[eInd][1]; j = edges[eInd][2]
                site_cnt[i] += 1; site_cnt[j] += 1
                connected_sites[i][site_cnt[i]] = j; connected_sites[j][site_cnt[j]] = i
                popfirst!(sorted_arg)
                break
            end
        end                
    end

    tpv = findall(x->x==1, site_cnt)
    opt_cost += C[tpv[1], tpv[2]]

    connected_sites[tpv[1]][2] = tpv[2]
    connected_sites[tpv[2]][2] = tpv[1]
    

    opt_tour = [1, connected_sites[1][1]]
    for sInd in 1:1:(sitesN-1)
        prev_site = opt_tour[end-1]
        curr_site = opt_tour[end]
        if connected_sites[curr_site][1] == prev_site
            push!(opt_tour, connected_sites[curr_site][2])
        else
            push!(opt_tour, connected_sites[curr_site][1])
        end
    end

    return opt_tour, opt_cost
end

#* name: FarthestInsertion
# inputs:
#   - C: {C_ij} is the length between the sites i and j
# outputs:
#   - opt_tour: optimal tour
#   - opt_cost: optimal tour length
function FarthestInsertion(C)
    sitesN = size(C)[1]
    
    opt_tour = Int[rand(1:1:sitesN)]
    push!(opt_tour, argmax(C[opt_tour[end], :]))
    opt_cost = 2 * C[opt_tour[1], opt_tour[end]]

    remaining_sites = collect(setdiff(Set(collect(1:1:sitesN)), Set(opt_tour)))
    for _ in 1:1:(sitesN-2)
        sInd = argmin(C[opt_tour, remaining_sites])[2]
        new_site = remaining_sites[sInd]
        
        min_dl = Inf
        insert_at = 0
        for sInd in 1:1:(length(opt_tour)-1)
            dl = C[opt_tour[sInd], new_site] + C[new_site, opt_tour[sInd+1]] - C[opt_tour[sInd], opt_tour[sInd+1]]
            if dl < min_dl
                min_dl = dl
                insert_at = sInd
            end
        end
        dl = C[opt_tour[end], new_site] + C[new_site, opt_tour[1]] - C[opt_tour[end], opt_tour[1]]
        if dl < min_dl
            min_dl = dl
            insert_at = length(opt_tour)
        end
        
        opt_tour = vcat(opt_tour[1:insert_at], [new_site], opt_tour[(insert_at+1):end])
        opt_cost += min_dl
        remaining_sites = collect(setdiff(Set(remaining_sites), [new_site]))
    end

    return opt_tour, opt_cost
end

# TODO
function Christofides(C)
    opt_tour = 1
    opt_cost = 1
    return opt_tour, opt_cost
end

#* name: TwoOptSwapSlow
# inputs:
#   - C: {C_ij} is the length between the sites i and j
#   - init_tour: the initial tour where the algorithm starts
# outputs:
#   - opt_tour: optimal tour
#   - opt_cost: optimal tour length
function TwoOptSwapSlow(C, init_tour)
    sitesN = size(C)[1]
    opt_tour = copy(init_tour)
    opt_cost = sum([C[opt_tour[sInd], opt_tour[sInd+1]] for sInd in 1:1:(sitesN-1)]) + C[opt_tour[end], opt_tour[1]]
    flag = false
    while !flag
        flag = true
        for sInd in 1:1:(sitesN-1)
            for ssInd in (sInd+1):1:sitesN
                new_tour = vcat(opt_tour[1:sInd], opt_tour[ssInd:-1:(sInd+1)], opt_tour[(ssInd+1):end])
                new_cost = sum([C[new_tour[sInd], new_tour[sInd+1]] for sInd in 1:1:(sitesN-1)]) + C[new_tour[end], new_tour[1]]
                if new_cost < opt_cost
                    opt_tour = new_tour
                    opt_cost = new_cost
                    flag = false
                    break
                end
            end
        end
    end

    return opt_tour, opt_cost
end

#* name: TwoOptSwap
# inputs:
#   - C: {C_ij} is the length between the sites i and j
#   - init_tour: the initial tour where the algorithm starts
#   - mode: 1 - update the solution whenever meet an improvement, otherwise - update the solution with the best swap
# outputs:
#   - opt_tour: optimal tour
#   - opt_cost: optimal tour length
function TwoOptSwap(C, init_tour, mode=1)
    global eps
    sitesN = size(C)[1]
    
    opt_tour = copy(init_tour)
    opt_cost = sum([C[opt_tour[sInd], opt_tour[sInd+1]] for sInd in 1:1:(sitesN-1)]) + C[opt_tour[end], opt_tour[1]]
    flag = false
    if mode==1
        while !flag
            flag = true
            for sInd in 1:1:(sitesN-1)
                for ssInd in (sInd+1):1:sitesN
                    dl = C[opt_tour[sInd], opt_tour[ssInd]] + C[opt_tour[sInd+1], opt_tour[ssInd+1 - (ssInd==sitesN)*sitesN]] - C[opt_tour[ssInd], opt_tour[ssInd+1 - (ssInd==sitesN)*sitesN]] - C[opt_tour[sInd], opt_tour[sInd+1]]
                    if dl < -eps
                        opt_tour = vcat(opt_tour[1:sInd], opt_tour[ssInd:-1:(sInd+1)], opt_tour[(ssInd+1):end])
                        opt_cost += dl
                        flag = false
                        break
                    end
                end
            end
        end
    else
        while true
            flag = true
            dl_min = -eps
            candidates = Vector{Int}(undef, 2)
            for sInd in 1:1:(sitesN-1)
                for ssInd in (sInd+1):1:sitesN
                    dl = C[opt_tour[sInd], opt_tour[ssInd]] + C[opt_tour[sInd+1], opt_tour[ssInd+1 - (ssInd==sitesN)*sitesN]] - C[opt_tour[ssInd], opt_tour[ssInd+1 - (ssInd==sitesN)*sitesN]] - C[opt_tour[sInd], opt_tour[sInd+1]]
                    if dl < dl_min
                        dl_min = dl
                        candidates[1] = sInd
                        candidates[2] = ssInd
                        flag = false
                    end
                end
            end
            if flag
                break
            end
            opt_tour = vcat(opt_tour[1:candidates[1]], opt_tour[candidates[2]:-1:(candidates[1]+1)], opt_tour[(candidates[2]+1):end])
            opt_cost += dl_min
        end
    end

    return opt_tour, opt_cost
end

end

module MetaHeuristics

end