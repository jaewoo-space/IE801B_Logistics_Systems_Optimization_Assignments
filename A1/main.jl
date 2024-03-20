include("ShortestPathProblems.jl")
using .Tools: ProblemGenerator
using .LabelCorrecting: ModifiedBellmanFordAlgorithm
using Statistics: mean
using DelimitedFiles: readdlm, writedlm
using ProgressBars
using Random, DelimitedFiles
using Base.Threads: @threads
using ProgressBars

## Instances generation
println("========== Instances Generation ==========")
Random.seed!(0)

instN = 1000
groupN = 19

density = [0.2, 0.5, 0.8, 0.2, 0.5, 0.8, 0.2, 0.5, 0.8, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5]
calV = [10, 10, 10, 50, 50, 50, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
cl = [0, 0, 0, 0, 0, 0, 0, 0, 0, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5]
cu = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50]
n_ratio = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.02, 0.05, 0.08, 0.1, 0.01, 0.02, 0.05, 0.08, 0.1]

mkpath("./instances")

inst_info_tot = zeros(Int64, groupN, 3)
for gInd in 1:1:groupN
    println("Generate Instances of Group$(gInd)!")
    calE = round(Int64, calV[gInd]*(calV[gInd]-1) * density[gInd])
    neg_c = round(Int64, calE*n_ratio[gInd])
    inst_info_tot[gInd, :] = [calV[gInd], calE, neg_c]

    mkpath("./instances/group$gInd")

    writedlm("./instances/group$gInd/Info.csv", [calV[gInd], calE, cl[gInd], cu[gInd], neg_c])

    # @threads for eInd in ProgressBar(1:1:instN)
    for eInd in ProgressBar(1:1:instN)
        A = ProblemGenerator(calV[gInd], calE, cl[gInd], cu[gInd], neg_c)
        writedlm("./instances/group$gInd/test$eInd.csv", A)
    end
end
writedlm("./instances/inst_info_tot.csv", inst_info_tot)

## Solve problems
println("========== Solve Problems ==========")
results_tot = zeros(groupN, 6)

mkpath("./results")

for gInd in 1:1:groupN
    println("Solve Group$(gInd) Instances!")
    n_ex = zeros(Float64, instN, 2)
    time = zeros(Float64, instN, 2)
    for i in ProgressBar(1:1:instN)
        A = convert.(Int64, readdlm("./instances/group$gInd/test$i.csv"))
        s = 1
        d1, pred1, n_ex[i, 1], time[i, 1] = ModifiedBellmanFordAlgorithm(A, s, true)
        d2, pred2, n_ex[i, 2], time[i, 2] = ModifiedBellmanFordAlgorithm(A, s, false)
    end

    results_tot[gInd, 1:2] = round.(mean(n_ex, dims=1), digits=1)
    results_tot[gInd, 4:5] = mean(time[2:end, :], dims=1)
    results_tot[gInd, 3] = round((results_tot[gInd, 1] - results_tot[gInd, 2]) / results_tot[gInd, 2] * 100, digits=1)
    results_tot[gInd, 6] = round((results_tot[gInd, 4] - results_tot[gInd, 5]) / results_tot[gInd, 5] * 100, digits=1)

    mkpath("./results/group$gInd")

    writedlm("./results/group$(gInd)/label_examination_count.csv", n_ex)
    writedlm("./results/group$(gInd)/computation_time.csv", time)
end

writedlm("./results/results_tot.csv", results_tot)