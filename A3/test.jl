using JuMP, Gurobi, LinearAlgebra
include("TSP.jl")
using .Heuristics: TwoOptSwapSlow, TwoOptSwap, MyTwoOptSwap
using .Tools: RandomProblemGenerator, VisualizeTour, VisualizeLazyConstraints
using Plots
using Random: seed!
using Concorde: solve_tsp
using DelimitedFiles: writedlm, readdlm
using Statistics: mean, std
using ProgressBars
using Base.Threads: @threads

N = 1000
sites = readdlm("./instances/N$N/inst1.csv")
C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:N, j in 1:1:N]

seed!(2718281)
# init_tour = collect(1:1:1000)
# init_tour, _ = FarthestInsertion(C)
using StatsBase: sample
init_tour = sample(1:1:N, N, replace=false)

# @time opt_tour, opt_cost = TwoOptSwap(C, init_tour, 1)
# @time opt_tour, opt_cost = TwoOptSwap(C, init_tour, 2)
@time opt_tour, opt_cost = MyTwoOptSwap(C, init_tour)