using JuMP, LinearAlgebra
include("TSP.jl")
using .Heuristics: TwoOptSwapSlow, TwoOptSwap, MyTwoOptSwap, TwoOptSwapForGIF, MyTwoOptSwapForGIF
using .Tools: RandomProblemGenerator, VisualizeTourHist
using Plots
using Random: seed!, shuffle
using Concorde: solve_tsp
using DelimitedFiles: writedlm, readdlm
using Statistics: mean, std
using ProgressBars
using Base.Threads: @threads


## Instance generation
# seed!(3141592)

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

#
sites_test = readdlm("./instances/N20/inst20.csv");
C_test = [norm(sites_test[i, :] .- sites_test[j, :]) for i in 1:1:20, j in 1:1:20];
init_tour_test = collect(1:1:20);
TwoOptSwapSlow(C_test, init_tour_test);
TwoOptSwap(C_test, init_tour_test, 1);
MyTwoOptSwap(C_test, init_tour_test);

## Experiments
mkpath("./results")

# for 1a
seed!(2718281)
opt_costs = zeros(instN, 3)
comp_times = zeros(instN, 3)
println("=========== Exp. for 1a ===========")
for nInd in 1:1:3, iInd in ProgressBar(1:1:instN)
    sites = readdlm("./instances/N$(sitesN[nInd])/inst$iInd.csv")
    C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN[nInd], j in 1:1:sitesN[nInd]]
    init_tour = collect(1:1:sitesN[nInd])

    ti = time()
    _, opt_cost = TwoOptSwapSlow(C, init_tour)
    te = time()
    opt_costs[iInd, nInd] = round(opt_cost, digits=4)
    comp_times[iInd, nInd] = te-ti
end
writedlm("./results/1a_cost.csv", opt_costs)
writedlm("./results/1a_time.csv", comp_times)

# for the others
seed!(2718281)
opt_costs = [zeros(instN, length(sitesN)) for _ in 1:3]
comp_times = [zeros(instN, length(sitesN)) for _ in 1:3]
println("========== Exp. for other questions ==========")
for nInd in eachindex(sitesN)
    println("---------- N = $(sitesN[nInd]) ----------")
    for iInd in ProgressBar(1:1:instN)
        sites = readdlm("./instances/N$(sitesN[nInd])/inst$iInd.csv")
        C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN[nInd], j in 1:1:sitesN[nInd]]
        init_tour = collect(1:1:sitesN[nInd])

        ti = time()
        _, opt_cost = TwoOptSwap(C, init_tour, 1)
        te = time()
        opt_costs[1][iInd, nInd] = round(opt_cost, digits=4)
        comp_times[1][iInd, nInd] = te-ti

        ti = time()
        _, opt_cost = TwoOptSwap(C, init_tour, 2)
        te = time()
        opt_costs[2][iInd, nInd] = round(opt_cost, digits=4)
        comp_times[2][iInd, nInd] = te-ti

        ti = time()
        _, opt_cost = MyTwoOptSwap(C, init_tour)
        te = time()
        opt_costs[3][iInd, nInd] = round(opt_cost, digits=4)
        comp_times[3][iInd, nInd] = te-ti
    end
    writedlm("./results/1b_cost.csv", opt_costs[1])
    writedlm("./results/1b_time.csv", comp_times[1])

    writedlm("./results/2b_cost.csv", opt_costs[2])
    writedlm("./results/2b_time.csv", comp_times[2])

    writedlm("./results/3_cost.csv", opt_costs[3])
    writedlm("./results/3_time.csv", comp_times[3])
end

## Visuals
sites_gif = readdlm("./instances/N20/inst20.csv")
C_gif = [norm(sites_gif[i, :] .- sites_gif[j, :]) for i in 1:1:20, j in 1:1:20]
init_tour_gif = collect(1:1:20)

opt_tour_hist_gif = Vector{Vector{Vector{Int}}}(undef, 3)
opt_tour_hist_gif[1], _ = TwoOptSwapForGIF(C_gif, init_tour_gif, 1)
opt_tour_hist_gif[2], _ = TwoOptSwapForGIF(C_gif, init_tour_gif, 2)
opt_tour_hist_gif[3], _ = MyTwoOptSwapForGIF(C_gif, init_tour_gif)

mkpath("./results/gif")
for eInd in 1:1:3
    if eInd == 3
        fig_list = VisualizeTourHist(sites_gif, opt_tour_hist_gif[eInd]); [push!(fig_list, fig_list[end]) for i in 1:1:80]
    else
        fig_list = VisualizeTourHist(sites_gif, opt_tour_hist_gif[eInd]); [push!(fig_list, fig_list[end]) for i in 1:1:10]
    end
    for fInd in eachindex(fig_list)
        savefig(fig_list[fInd], "./results/gif/fig$eInd-$fInd.png")
    end
    anim = @animate for fig in fig_list
        plot(fig)
    end
    if eInd == 3
        Plots.buildanimation(anim, "./results/animation$eInd.gif", fps=32)
    else
        Plots.buildanimation(anim, "./results/animation$eInd.gif", fps=4)
    end
end