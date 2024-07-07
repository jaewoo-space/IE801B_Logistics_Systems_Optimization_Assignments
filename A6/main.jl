include("./Functions.jl")

using DelimitedFiles, StatsBase, Random, Plots, JLD
using .Tools
using .Solvers

Random.seed!(2718)


#! Instance 1

rf = 0.2    # the number of facilities ratio to the number of locations

# load data and generate instance
data = readdlm("sortcap.csv")
n = size(data)[1]
facilities = sample(1:n, round(Int, n*rf), replace=false)
customers = sort(collect(setdiff(Set(collect(1:n)), Set(facilities))))

# instance coefficients
a = [data[j, 5]/1e02 for j in customers]
b = [data[i, 4]/1e02 for i in facilities]
c = [haversine((data[i, 3], data[i, 2]), (data[j, 3], data[j, 2]))/1e05 for i in facilities, j in customers]
F = [data[i, 6]/1e03 for i in facilities]
prob1 = Dict(
    "a" => a,
    "b" => b,
    "c" => c,
    "F" => F
)

# exact solver
obj_exact, x_exact, y_exact = Exact(a, b, c, F)

# heuristic solver
obj_heuristic, x_heuristic, y_heuristic, Zl_list, Zu_list = Klincewicz(a, b, c, F, 1e-4)

# visuals generation
fig1_exact, fig2_exact = VisualizeResults(data[facilities, [2,3]], data[customers, [2,3]], x_exact, y_exact)
fig1_heuristic, fig2_heuristic = VisualizeResults(data[facilities, [2,3]], data[customers, [2,3]], x_heuristic, y_heuristic)
fig1_exact = plot(fig1_exact, xlabel="Longitude", ylabel="Latitude")
fig2_exact = plot(fig2_exact, title="obj=$(round(obj_exact, digits=2))", xlabel="Longitude", ylabel="Latitude")
fig2_heuristic = plot(fig2_heuristic, title="obj=$(round(obj_heuristic, digits=2))", xlabel="Longitude", ylabel="Latitude")

fig3_heuristic = plot()
plot!(Zu_list, color=:red, label="Upper bound", legend=:bottomright, ylimits=(400.0, 550.0), xlabel="Iterations")
plot!(Zl_list, color=:blue, label="Lower bound")
xlabel!("Iterations")

# save the results
savefig(fig1_exact, "prob1_instance.png")
savefig(fig2_exact, "prob1_exact.png")
savefig(fig2_heuristic, "prob1_heuristic.png")
savefig(fig3_heuristic, "prob1_bounds.png")

#! Instance 2

n = 100     # the number of random locations

rf = 0.75   # the number of facilities ratio to the number of locations

# generate random locations and allocate their roles
locations = rand(n, 2)
facilities = sample(1:n, round(Int, n*rf), replace=false)
customers = sort(collect(setdiff(Set(collect(1:n)), Set(facilities))))
N = length(facilities)
M = length(customers)

# instance coefficients
a = rand(1:10, M)
b = rand(20:10:50, N)
c = round.(Int, [sqrt(sum((locations[i, :] .- locations[j, :]).^2))*10 for i in facilities, j in customers])
F = rand(50:50:150, N)
prob2 = Dict(
    "a" => a,
    "b" => b,
    "c" => c,
    "F" => F
)

# exact solver
obj_exact, x_exact, y_exact = Exact(a, b, c, F)

# heuristic solver
obj_heuristic, x_heuristic, y_heuristic, Zl_list, Zu_list = Klincewicz(a, b, c, F, 1e-4)

# visuals generation
fig1_exact, fig2_exact = VisualizeResults(locations[facilities, :], locations[customers, :], x_exact, y_exact)
fig1_heuristic, fig2_heuristic = VisualizeResults(locations[facilities, :], locations[customers, :], x_heuristic, y_heuristic)
fig1_exact = plot(fig1_exact, xlabel="longitude", ylabel="latitude")
fig2_exact = plot(fig2_exact, title="obj=$(round(obj_exact, digits=2))", xlabel="X", ylabel="Y")
fig2_heuristic = plot(fig2_heuristic, title="obj=$(round(obj_heuristic, digits=2))", xlabel="X", ylabel="Y")
fig3_heuristic = plot()
plot!(Zu_list, color=:red, label="Upper bound", legend=:bottomright, ylimits=(100.0, 300.0), xlabel="Iterations")
plot!(Zl_list, color=:blue, label="Lower bound")

# save the results
savefig(fig1_exact, "prob2_instance.png")
savefig(fig2_exact, "prob2_exact.png")
savefig(fig2_heuristic, "prob2_heuristic.png")
savefig(fig3_heuristic, "prob2_bounds.png")

# save the instances
save("instances.jld", "prob1", prob1, "prob2", prob2)