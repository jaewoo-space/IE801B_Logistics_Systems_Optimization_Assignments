module Tools

using StatsBase: sample
using LinearAlgebra: norm
using Plots

export RandomProblemGenerator, VisualizeTour, VisualizeIntermediateSolution

#* name: RandomProblemGenerator
# inputs:
#   - sitesN: the number of sites
#   - sqrL: the lenght of the one side of the squre where the sites are distributed
# outputs:
#   - sites: coordinates of the sites
#   - C: {C_ij} is the length between the sites i and j
function RandomProblemGenerator(sitesN::Int64)
    sites = rand(sitesN, 2)
    return sites
end

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
#   - x_opt: optimal solution of the TSP formulated as IP
#   - opt_tour: 
#   - opt_cost: 
#   - comp_time: 

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

using LinearAlgebra

export NearestNeighbor, Greedy, FarthestInsertion, Christofides

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

# incomplete
function Greedy(C)
    sitesN = size(C)[1]

    site_cnt = zeros(Int, sitesN)
    edges = [(i, j) for i in 1:1:sitesN for j in 1:1:sitesN if i<j]
    edges_cost = [C[edge[1], edge[2]] for edge in edges]
    selected_edges = []
    for _ in 1:1:sitesN

        # TODO: subtour check part
        while true
            eInd = argmin(edges_cost)
            break
        end

        push!(selected_edges, edges[eInd])
        i, j = edges[eInd]

        popat!(edges, eInd)
        popat!(edges_cost, eInd)
        if site_cnt[i] == 1
            tpv1 = Set(findall(x->(x[1] == i) || (x[2] == i), edges))
            tpv2 = collect(setdiff(Set(eachindex(edges)), tpv1))
            edges = edges[tpv2]
            edges_cost = edges_cost[tpv2]
        end
        if site_cnt[j] == 1
            tpv1 = Set(findall(x->(x[1] == j) || (x[2] == j), edges))
            tpv2 = collect(setdiff(Set(eachindex(edges)), tpv1))
            edges = edges[tpv2]
            edges_cost = edges_cost[tpv2]
        end
        
        site_cnt[i] += 1; site_cnt[j] += 1
    end
    
    opt_cost = sum([C[edge[1], edge[2]] for edge in selected_edges])

    opt_tour = [selected_edges[1][1], selected_edges[1][2]]
    popat!(selected_edges, 1)
    for _ in 1:1:(sitesN-2)
        tpv = findfirst(x->(x[1]==opt_tour[end]) || (x[2]==opt_tour[end]), selected_edges)
        if selected_edges[tpv][1] == opt_tour[end]
            push!(opt_tour, selected_edges[tpv][2])
            popat!(selected_edges, tpv)
        else
            push!(opt_tour, selected_edges[tpv][1])
            popat!(selected_edges, tpv)
        end
    end
    
    return opt_tour, opt_cost
end

function FarthestInsertion(C)
    sitesN = size(C)[1]
    
    opt_tour = [rand(1:1:sitesN)]
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

function Christofides(C)
    opt_tour = 1
    opt_cost = 1
    return opt_tour, opt_cost
end

end

module MetaHeuristics

end