include("ESPPRC.jl")
using .Solvers: Feillet
using CVRPLIB, DelimitedFiles

prob_name_list = ["P-n16-k8", "A-n32-k5", "B-n64-k9"]
pInd = 1

cvrp, _, _ = readCVRPLIB(prob_name_list[pInd])
dual = readdlm("./data/dual_var_"*prob_name_list[pInd]*".csv")

capacity = cvrp.capacity
weights = cvrp.weights
demand = cvrp.demand
depot = cvrp.depot
customers = cvrp.customers

Feillet(capacity, demand, depot, customers, weights .- dual[1:(end-1)])