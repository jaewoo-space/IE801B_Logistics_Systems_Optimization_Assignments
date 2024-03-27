using DelimitedFiles
using Statistics: mean, std

opt_cost = readdlm("./instances/opt_sol.csv")

# 1a
cost1 = readdlm("./results/1a_cost.csv")
time1 = readdlm("./results/1a_time.csv")

mean(time1, dims=1)
std(time1, dims=1)

cost2 = readdlm("./results/1b_cost.csv")
time2 = readdlm("./results/1b_time.csv")
opt_gap2 = (cost2.-opt_cost)./opt_cost

mean(time2, dims=1)
std(time2, dims=1)
mean(opt_gap2, dims=1) .* 100
std(opt_gap2, dims=1) .* 100

# 1b

cost3 = readdlm("./results/2b_cost.csv")
time3 = readdlm("./results/2b_time.csv")
opt_gap3 = (cost3.-opt_cost)./opt_cost

mean(time3, dims=1)
std(time3, dims=1)
mean(opt_gap3, dims=1) .* 100
std(opt_gap3, dims=1) .* 100


cost4 = readdlm("./results/3_cost.csv")
time4 = readdlm("./results/3_time.csv")
opt_gap4 = (cost4.-opt_cost)./opt_cost

mean(time4, dims=1)
std(time4, dims=1)
mean(opt_gap4, dims=1) .* 100
std(opt_gap4, dims=1) .* 100