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
using Base.Threads: @threads

seed!(3141592)

# println("=========== Instances Generation ==========")

sitesN =  [20, 50, 100, 1000]
instN = 100

# opt_cost = zeros(instN, length(sitesN))

# for nInd in eachindex(sitesN)
#     println("---------- N = $(sitesN[nInd]) ----------")
#     mkpath("./instances/N$(sitesN[nInd])")
    
#     @threads for iInd in ProgressBar(1:1:instN)
#         sites = RandomProblemGenerator(sitesN[nInd])
#         writedlm("./instances/N$(sitesN[nInd])/inst$iInd.csv", sites)

#         _, tpv = solve_tsp(sites[:, 1].*10000, sites[:, 2].*10000; dist="EUC_2D")
#         opt_cost[iInd, nInd] = tpv/10000
#     end
# end
# writedlm("./instances/opt_sol.csv", opt_cost)

## Exp

seed!(2718281)
mkpath("./results")
opt_costs = [zeros(instN, length(sitesN)) for _ in 1:4]
comp_times = [zeros(instN, length(sitesN)) for _ in 1:4]
for nInd in eachindex(sitesN)
    println("---------- N = $(sitesN[nInd]) ----------")
    for iInd in ProgressBar(1:1:instN)
        sites = readdlm("./instances/N$(sitesN[nInd])/inst$iInd.csv")
        C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN[nInd], j in 1:1:sitesN[nInd]]
        init_tour = collect(1:1:sitesN[nInd])

        ti = time()
        _, opt_cost = TwoOptSwapSlow(C, init_tour)
        te = time()
        opt_costs[1][iInd, nInd] = round(opt_cost, digits=4)
        comp_times[1][iInd, nInd] = te-ti

        ti = time()
        _, opt_cost = TwoOptSwap(C, init_tour, 1)
        te = time()
        opt_costs[2][iInd, nInd] = round(opt_cost, digits=4)
        comp_times[2][iInd, nInd] = te-ti

        ti = time()
        _, opt_cost = TwoOptSwap(C, init_tour, 2)
        te = time()
        opt_costs[3][iInd, nInd] = round(opt_cost, digits=4)
        comp_times[3][iInd, nInd] = te-ti

        ti = time()
        _, opt_cost = TwoOptSwap(C, init_tour, 3)
        te = time()
        opt_costs[4][iInd, nInd] = round(opt_cost, digits=4)
        comp_times[4][iInd, nInd] = te-ti
    end
    writedlm("./results/1a_cost.csv", opt_costs[1])
    writedlm("./results/1a_time.csv", comp_times[1])

    writedlm("./results/1b_cost.csv", opt_costs[2])
    writedlm("./results/1b_time.csv", comp_times[2])

    writedlm("./results/2b_cost.csv", opt_costs[3])
    writedlm("./results/2b_time.csv", comp_times[3])

    writedlm("./results/3_cost.csv", opt_costs[4])
    writedlm("./results/3_time.csv", comp_times[4])
end

## Visualize

