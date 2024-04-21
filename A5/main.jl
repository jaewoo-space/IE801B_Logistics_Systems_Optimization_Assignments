include("ESPPRC.jl")
using .Solvers: Feillet
using CVRPLIB, DelimitedFiles

mkpath("./results")

prob_name_list = ["P-n16-k8", "A-n32-k5", "B-n64-k9"]

exm_cnt = zeros(Int64, 3)
for pInd in eachindex(prob_name_list)
    println("****************************************************************")
    println("**************************  $(prob_name_list[pInd])  **************************")
    println("****************************************************************")

    ## 1) load the instance and duals
    cvrp, _, _ = readCVRPLIB(prob_name_list[pInd])
    dual = readdlm("./data/dual_var_"*prob_name_list[pInd]*".csv")
    
    ## 2) preprocessing
    # capacity, demand, depot, and customers
    capacity = cvrp.capacity
    demand = vcat(cvrp.demand, cvrp.demand[1])
    depot = cvrp.depot
    customers = cvrp.customers

    # weights
    N = length(demand)
    weights = zeros(Float64, N-1, N)
    # a. add the destination column
    weights[:, 1:(N-1)] = cvrp.weights
    weights[:, N] = weights[:, 1]
    # b. screening out the depot -> depot arc
    weights[1, N] = Inf

    ## 3) solve ESPPRC
    result, exm_cnt[pInd] = Feillet(capacity, demand, depot, customers, weights .- dual[1:end-1])

    ## 4) save the result
    writedlm("./results/result_$(prob_name_list[pInd]).csv", result, '\t')
end

writedlm("./results/examinations_count.csv", exm_cnt, '\t')