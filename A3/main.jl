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

seed!(31415)

sitesN =  30

sites = RandomProblemGenerator(sitesN)
C = [norm(sites[i, :] .- sites[j, :]) for i in 1:1:sitesN, j in 1:1:sitesN]

opt_tour, opt_cost = Greedy(C)

fig = VisualizeTour(sites, opt_tour, "opt_cost=$(round(opt_cost, digits=2))")

