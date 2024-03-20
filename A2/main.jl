using JuMP, Gurobi, LinearAlgebra
include("TSP.jl")
using .Solvers: GurobiSolver, GurobiSolverWithLazyConstraints, GurobiSolverWithLazyConstraintsForGIF
using .Tools: SymmetricProblemGenerator, VisualizeTour, VisualizeIntermediateSolution
using Plots
using Random: seed!
using Concorde: solve_tsp
using DelimitedFiles: writedlm, readdlm
using ProgressBars

seed!(3141592)

sitesN = collect(10:5:100)
# println("========== Generate Instances ==========")
# mkpath("instances")

# for n in ProgressBar(sitesN)
#     sites = SymmetricProblemGenerator(n)
#     writedlm("./instances/tsp$(n)_inst.csv", sites)
# end

println("========== Solve Problems ==========")
mkpath("results")

for iInd in ProgressBar(eachindex(sitesN))
    sites = readdlm("./instances/tsp$(sitesN[iInd])_inst.csv")
    C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN[iInd], j in 1:1:sitesN[iInd]]
    if iInd == 1
        GurobiSolver(C)
        GurobiSolverWithLazyConstraints(C)
        solve_tsp(sites[:, 1].*10000, sites[:, 2].*10000; dist="EUC_2D")
    end

    results = Matrix{Any}(undef, 3, sitesN[iInd]+2)
    

    ## All subtour constraints included
    opt_tour, opt_cost, comp_time = GurobiSolver(C)
    results[1, 1] = round(comp_time, digits=3)
    results[1, 2] = round(opt_cost, digits=4)
    results[1, 3:end] = opt_tour

    ## Lazy constraints
    opt_tour, opt_cost, comp_time = GurobiSolverWithLazyConstraints(C)
    results[2, 1] = round(comp_time, digits=3)
    results[2, 2] = round(opt_cost, digits=4)
    results[2, 3:end] = opt_tour

    ## Concorde
    ti = time()
    opt_tour, opt_cost = solve_tsp(sites[:, 1].*10000, sites[:, 2].*10000; dist="EUC_2D")
    te = time()
    results[3, 1] = round(te - ti, digits=3)
    results[3, 2] = round(opt_cost/10000, digits=4)
    results[3, 3:end] = opt_tour

    writedlm("./results/tsp$(sitesN[iInd])_results.csv", results)
end

# println("========== Generate Visuals ==========")
# for iInd in ProgressBar(eachindex(sitesN))

#     fig1 = VisualizeTour(sites, opt_tour, "time = $(comp_times[iInd, 1]) secs")
#     savefig(fig1, "./results/tsp$(sitesN[iInd])_gurobi.png")

#     fig2 = VisualizeTour(sites, opt_tour, "time = $(comp_times[iInd, 2]) secs")
#     savefig(fig2, "./results/tsp$(sitesN[iInd])_lazy_constraints.png")

#     fig3 = VisualizeTour(sites, opt_tour, "time = $(comp_times[iInd, 3]) secs")
#     savefig(fig3, "./results/tsp$(sitesN[iInd])_concorde.png")

# end

# println("========== Generate .gif ==========")
# mkpath("./results/gif")

# sites = readdlm("./instances/tsp20_inst.txt")
# C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN[iInd], j in 1:1:sitesN[iInd]]

# x_hist = GurobiSolverWithLazyConstraintsForGIF(C, true)
    
#     anim = @animate for i in eachindex(x_hist)
        
# end