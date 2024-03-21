using JuMP, Gurobi, LinearAlgebra
include("TSP.jl")
using .Solvers: GurobiSolver, GurobiSolverWithLazyConstraints, GurobiSolverWithLazyConstraintsForGIF
using .Tools: RandomProblemGenerator, VisualizeTour, VisualizeLazyConstraints
using Plots
using Random: seed!
using Concorde: solve_tsp
using DelimitedFiles: writedlm, readdlm
using Statistics: mean, std
using ProgressBars

seed!(31415)

## Exp 1
println("========== Experiment #1 ==========")
mkpath("./instances/exp1")
mkpath("./results/exp1")

sitesN = collect(5:1:40)

println("----- Generate Instances -----")
for n in ProgressBar(sitesN)
    sites = RandomProblemGenerator(n)
    writedlm("./instances/exp1/tsp$(n)_inst.csv", sites)
end

println("----- Solve Instances -----")
for iInd in ProgressBar(eachindex(sitesN))
    sites = readdlm("./instances/exp1/tsp$(sitesN[iInd])_inst.csv")
    C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN[iInd], j in 1:1:sitesN[iInd]]

    if iInd == 1
        GurobiSolver(C)
        GurobiSolverWithLazyConstraints(C)
        solve_tsp(sites[:, 1].*10000, sites[:, 2].*10000; dist="EUC_2D")
    end

    results = Matrix{Any}(undef, 3, sitesN[iInd]+2)
    
    ## Concorde
    ti = time()
    opt_tour, opt_cost = solve_tsp(sites[:, 1].*10000, sites[:, 2].*10000; dist="EUC_2D")
    te = time()
    results[1, 1] = round(te - ti, digits=3)
    results[1, 2] = round(opt_cost/10000, digits=4)
    results[1, 3:end] = opt_tour

    ## All at once
    opt_tour, opt_cost, comp_time = GurobiSolver(C)
    results[2, 1] = round(comp_time, digits=3)
    results[2, 2] = round(opt_cost, digits=4)
    results[2, 3:end] = opt_tour

    ## Lazy constraints
    opt_tour, opt_cost, comp_time = GurobiSolverWithLazyConstraints(C)
    results[3, 1] = round(comp_time, digits=3)
    results[3, 2] = round(opt_cost, digits=4)
    results[3, 3:end] = opt_tour

    writedlm("./results/exp1/tsp$(sitesN[iInd])_results.csv", results)
end


println("----- Generate Results Plot -----")
time_plot = reduce(hcat, [readdlm("./results/exp1/tsp$(sitesN[iInd])_results.csv")[:, 1] for iInd in eachindex(sitesN)])
writedlm("./results/exp1/results_time.csv", time_plot')

fig = plot(dpi=300)
plot!(sitesN, time_plot[1, :], label="Concorde", linewidth=2)
plot!(sitesN, time_plot[2, :], label="All at once", linewidth=2)
plot!(sitesN, time_plot[3, :], label="Lazy constraints", linewidth=2)
xlabel!("# Nodes")
ylabel!("Computation Time (sec)")
savefig(fig, "./results_plot.png")


println("----- Visualize Lazy Constraints -----")
sitesN_gif = 10
sites_gif = readdlm("./instances/exp1/tsp$(sitesN_gif)_inst.csv")
C_gif = [norm(sites_gif[i, :] .- sites_gif[j, :]) for i in 1:1:sitesN_gif, j in 1:1:sitesN_gif]
x_hist_gif = GurobiSolverWithLazyConstraintsForGIF(C_gif)

fig_list = VisualizeLazyConstraints(sites_gif, x_hist_gif, "")
[push!(fig_list, fig_list[end]) for i in 1:1:4]

anim = @animate for fInd in eachindex(fig_list)
    plot(fig_list[fInd])
end

Plots.buildanimation(anim, "./animation_tsp$(sitesN_gif).gif", fps=4)

## Exp 2
println("========== Experiment #2 ==========")
mkpath("./instances/exp2")
mkpath("./results/exp2")

instN = 100
sitesN = 30

for iInd in ProgressBar(1:1:instN)
    sites = RandomProblemGenerator(sitesN)
    writedlm("./instances/exp2/tsp_inst$(iInd).csv", sites)
end

println("----- Solve Random Problems -----")
results = zeros(instN+2, 3)

for iInd in ProgressBar(1:1:instN)
    sites = readdlm("./instances/exp2/tsp_inst$(iInd).csv")
    C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN, j in 1:1:sitesN]

    if iInd == 1
        GurobiSolver(C)
        GurobiSolverWithLazyConstraints(C)
        solve_tsp(sites[:, 1].*10000, sites[:, 2].*10000; dist="EUC_2D")
    end

    ## Concorde
    ti = time()
    opt_tour, opt_cost = solve_tsp(sites[:, 1].*10000, sites[:, 2].*10000; dist="EUC_2D")
    te = time()
    results[iInd, 1] = round(te - ti, digits=3)
    
    ## All at once
    opt_tour, opt_cost, comp_time = GurobiSolver(C)
    results[iInd, 2] = round(comp_time, digits=3)
    
    ## Lazy constraints
    opt_tour, opt_cost, comp_time = GurobiSolverWithLazyConstraints(C)
    results[iInd, 3] = round(comp_time, digits=3)
end

results[instN+1, :] = mean(results[1:instN, :], dims=1)
results[instN+2, :] = std(results[1:instN, :], dims=1)

writedlm("./results/exp2/results_time.csv", results)