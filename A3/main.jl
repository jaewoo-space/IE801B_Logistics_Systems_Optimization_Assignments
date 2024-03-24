using JuMP, Gurobi, LinearAlgebra
include("TSP.jl")
using .Heuristics: NearestNeighbor, Greedy, FarthestInsertion
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

opt_tour, opt_cost = NearestNeighbor(C)

fig = VisualizeTour(sites, opt_tour, "opt_cost=$(opt_cost)")

opt_tour, opt_cost = FarthestInsertion(C)

fig = VisualizeTour(sites, opt_tour, "opt_cost=$(opt_cost)")
