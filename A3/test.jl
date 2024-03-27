using JuMP, Gurobi, LinearAlgebra
include("TSP.jl")
using .Heuristics: TwoOptSwapSlow, TwoOptSwap, MyTwoOptSwap, TwoOptSwapForGIF, MyTwoOptSwapForGIF
using .Tools: RandomProblemGenerator, VisualizeTourHist
using Plots
using Random: seed!, shuffle
using Concorde: solve_tsp
using DelimitedFiles: writedlm, readdlm
using Statistics: mean, std
using ProgressBars

seed!(2718281)

N = 1000

sites = readdlm("./instances/N$N/inst1.csv")
C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:N, j in 1:1:N]
init_tour = collect(1:1:N)

@time MyTwoOptSwap(C, init_tour)
@time TwoOptSwap(C, init_tour, 2)