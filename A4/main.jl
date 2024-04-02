include("./CSPP.jl")
using .Tools: CSPP_INSTANCE
using .Solver: IPSolver, MyCGSolver
using JLD

println("========== Solution Procedure ==========")
## problem 1
println("------ Solve Problem 1 ------")
prob1 = load("prob1.jld")["prob"]

data1 = prob1.data
Tmax1 = prob1.Tmax
origin1, destination1 = prob1.origin, prob1.destination
solution_history_cg1, opt_sol_cg1 = MyCGSolver(data1, Tmax1, origin1, destination1)
opt_cost_ip1, opt_path_ip1 = IPSolver(data1, Tmax1, origin1, destination1)

## problem 2
println("------ Solve Problem 2 ------")
prob2 = load("prob2.jld")["prob"]

data2 = prob2.data
Tmax2 = prob2.Tmax
origin2, destination2 = prob2.origin, prob2.destination
solution_history_cg2, opt_sol_cg2 = MyCGSolver(data2, Tmax2, origin2, destination2)
opt_cost_ip2, opt_path_ip2 = IPSolver(data2, Tmax2, origin2, destination2)


## display the solution
println("==================== Results Summary ====================")
println("---------- Problem 1 ----------")
println("link/node-based formulation: selected arcs - $([Int.(data1[i, 1:2]) for i in collect(opt_path_ip1)]), opt_cost - $opt_cost_ip1")
println("column generation: selected arcs - $([Int.(data1[i, 1:2]) for i in collect(opt_sol_cg1[5])]), opt_cost - $(opt_sol_cg1[4])")

println("---------- Problem 2 ----------")
println("link/node-based formulation: selected arcs - $([Int.(data2[i, 1:2]) for i in collect(opt_path_ip2)]), opt_cost - $opt_cost_ip2")
println("column generation: selected arcs - $([Int.(data2[i, 1:2]) for i in collect(opt_sol_cg2[5])]), opt_cost - $(opt_sol_cg2[4])")