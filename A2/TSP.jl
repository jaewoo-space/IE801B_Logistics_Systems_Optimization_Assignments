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

module Solvers

using Gurobi, JuMP

export GurobiSolver, GurobiSolverWithLazyConstraints, GurobiSolverWithLazyConstraintsForVisuals

#* name: GurobiSolver
# inputs:
#   - C: {C_ij} is the length between the sites i and j
# outputs:
#   - x_opt: optimal solution of the TSP formulated as IP
#   - 

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