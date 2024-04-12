include("ESPPRC.jl")
using .Solvers: Feillet
using CVRPLIB, DelimitedFiles

cvrp1, _, _ = readCVRPLIB("P-n16-k8")
dual1 = readdlm("./data/dual_var_P-n16-k8.csv")

cvrp2, _, _ = readCVRPLIB("A-n32-k5")
dual2 = readdlm("./data/dual_var_A-n32-k5.csv")

cvrp3, _, _ = readCVRPLIB("B-n64-k9")
dual3 = readdlm("./data/dual_var_B-n64-k9.csv")


capacity = cvrp1.capacity
weights = cvrp1.weights
demand = cvrp1.demand
depot = cvrp1.depot
customers = cvrp1.customers

nodes = vcat(depot, customers)
N = length(nodes)
weights
demand
customers
label = zeros(Float64, length(nodes), 2 + 1 + N)