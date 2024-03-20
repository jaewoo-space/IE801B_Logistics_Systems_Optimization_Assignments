module Tools

using StatsBase: sample
using LinearAlgebra: norm
using Plots

export SymmetricProblemGenerator, VisualizeTour, VisualizeIntermediateSolution

#* name: SymmetricProblemGenerator
# inputs:
#   - sitesN: the number of sites
#   - sqrL: the lenght of the one side of the squre where the sites are distributed
# outputs:
#   - sites: coordinates of the sites
#   - C: {C_ij} is the length between the sites i and j
function SymmetricProblemGenerator(sitesN::Int64, cu=1)
    sites = cu * rand(sitesN, 2)
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

function VisualizeIntermediateSolution(sites, x_val, title = "")
    sitesN = size(sites)[1]
    
    tour_plot = vcat(tour, tour[1])

    fig = plot(dpi=300)
    for i in 1:1:sitesN
        plot!([sites[tour_plot[i], 1], sites[tour_plot[i+1], 1]], [sites[tour_plot[i], 2], sites[tour_plot[i+1], 2]]; legend = false)
    end
    xaxis!([0.0, 1.0])
    xlabel!("X")
    yaxis!([0.0, 1.0])
    xlabel!("Y")
    xticks!([0.0, 0.5, 1.0])
    yticks!([0.5, 1.0])
    scatter!(sites[:, 1], sites[:, 2])
    title!(title)
    return fig
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

    solution_summary(model, verbose=true)

    x_opt = Matrix{Bool}(undef, sitesN, sitesN)
    
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
        x_val = (y->callback_value(cb_data, y)).(x)
        subtour = find_subtour(x_val)
        if length(subtour) < sitesN
            con = @build_constraint(sum(x[i, j] for i in subtour, j in subtour) <= length(subtour)-1)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
        
        return
    end

    MOI.set(model, MOI.LazyConstraintCallback(), call_back_function)

    optimize!(model)

    return x_hist
end
end