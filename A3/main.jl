using JuMP, Gurobi, LinearAlgebra
include("TSP.jl")
using .Heuristics: TwoOptSwapSlow, TwoOptSwap, Greedy
using .Tools: RandomProblemGenerator, VisualizeTour, VisualizeLazyConstraints
using Plots
using Random: seed!
using Concorde: solve_tsp
using DelimitedFiles: writedlm, readdlm
using Statistics: mean, std
using ProgressBars
using Base.Threads

seed!(31415)

println("=========== Instances Generation ==========")

sitesN =  [20, 50, 100, 1000, 10000]
instN = 1000

mkpath("./instances/optimal")
opt_cost = zeros(instN, length(sitesN))

for nInd in eachindex(sitesN)
    println("---------- N = $(sitesN[nInd]) ----------")
    mkpath("./instances/N$(sitesN[nInd])")
    
    for iInd in ProgressBar(1:1:instN)
        sites = RandomProblemGenerator(sitesN[nInd])
        writedlm("./instances/N$(sitesN[nInd])/inst$iInd.csv", sites)

        _, tpv = solve_tsp(sites[:, 1].*10000, sites[:, 2].*10000; dist="EUC_2D")
        opt_cost[iInd, nInd] = tpv/10000
    end
end
writedlm("./instances/opt_sol.csv", opt_cost)

## optimal solution for each instance


## Exp 1

# for n in sitesN
#     for iInd in 1:1:instN
#         sites = readdlm("./instances/N$n/inst$Iind.csv")
#         C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:n, j in 1:1:n]
#         init_tour = collect(1:1:n)
        
#         ti = time()
#         _, opt_cost = TwoOptSwapSlow(C, init_tour)
#         te = time()


#     end
# end

# C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN, j in 1:1:sitesN]

# opt_tour, opt_cost = Greedy(C)

# fig = VisualizeTour(sites, opt_tour, "opt_cost=$(round(opt_cost, digits=2))")

# init_tour = opt_tour

# opt_tour, opt_cost = TwoOptSwap(C, init_tour, 1)

# fig = VisualizeTour(sites, opt_tour, "opt_cost=$(round(opt_cost, digits=2))")

# opt_tour, opt_cost = TwoOptSwap(C, init_tour, 2)

# fig = VisualizeTour(sites, opt_tour, "opt_cost=$(round(opt_cost, digits=2))")

# opt_tour, opt_cost = TwoOptSwap(C, init_tour, 3)

# fig = VisualizeTour(sites, opt_tour, "opt_cost=$(round(opt_cost, digits=2))")

