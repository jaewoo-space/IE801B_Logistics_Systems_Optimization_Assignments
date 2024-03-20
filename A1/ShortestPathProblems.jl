module Tools

using StatsBase: sample
using Random: shuffle, shuffle!
using DataStructures

export ProblemGenerator

#* name: ProblemGenerator
# inputs:
#   - calV: # of vertices
#   - calE: # of edges
#   - cl: lower bound of the length
#   - cu: upper bound of the length
#   - neg_c: # of negative lengths
# outputs:
#   - A: shortest path problem instance
function ProblemGenerator(calV, calE, cl::Int64 = 0, cu::Int64 = 10, neg_c::Int64 = 0)
    all_edge_candidates = [(i, j) for i in 1:1:calV for j in 1:1:calV if i!=j]
    # generate a tree
    tpv1 = vcat([1], shuffle(collect(2:1:calV)))
    tpv2 = [rand(tpv1[1:i]) for i in 1:1:(calV-1)]
    edges = [[tpv2[i], tpv1[i+1]] for i in 1:1:(calV-1)]
    # generate a graph DataStructures
    edge_candidates = collect(setdiff(Set(all_edge_candidates), Set(edges)))
    edges = vcat(edges, collect.(sample(edge_candidates, calE-calV+1, replace = false)))
    edge_candidates = collect(setdiff(Set(all_edge_candidates), Set(edges)))
    edges = reduce(hcat, edges)'

    if neg_c == 0
        # length allocation
        A = hcat(edges, rand(cl:1:cu, calE))
        A = A[sortperm(A[:, 1]), :]
        for i in 1:1:calV
            idx = findall(x->x==i, A[:, 1])
            A[idx, :] = A[idx, :][sortperm(A[idx, 2]), :]
        end
    else
        c = vcat(rand(0:1:cu, calE-neg_c), rand(cl:1:-1, neg_c))
        A = hcat(edges, c)
        A = A[sortperm(A[:, 1]), :]
        for i in 1:1:calV
            idx = findall(x->x==i, A[:, 1])
            A[idx, :] = A[idx, :][sortperm(A[idx, 2]), :]
        end
    end
    return A
end
end

module LabelSetting

export DijkstraAlgorithm

#* name: DijkstraAlgorithm
# inputs:
#   - A: forward star representation of a graph, (tail, head, cost)
#   - s: source node index, default = 1
# outputs:
#   - d: distance from the source of each node
#   - pred: predecssor of each node

function DijkstraAlgorithm(A, s::Int64 = 1)
    # check the nonnegative length assumption
    if sum(A[:, 3] .< 0) > 0
        println("Invalid instance!")
        return nothing
    end

    # preprocessing
    calV = maximum(A[:, 1:2])                                   # cardinality of the vertices set V (or, the number of the vertices)
    calE = size(A)[1]                                           # cardinality of the edges set E (or, the number of the edges)
    pnt = [findfirst(x->x==v, A[:, 1]) for v in 1:1:(calV+1)]   # point of each node
    pnt[calV+1] = calE+1
    tpv1 = findall(x->isnothing(x), pnt)
    tpv2 = deleteat!(collect(1:1:(calV+1)), tpv1)
    for v in tpv1
        pnt[v] = pnt[tpv2[findfirst(x->x>v, tpv2)]]
    end

    #* initialization
    inf = sum(A[:, 3]) + 1                                          # define a large number, playing a role as infinity in efficient
    i = s
    Sbar = copy(V); popat!(Sbar, s)                                 # unlabeled set
    S = [s]                                                         # labeled set
    d = fill(Inf, calV); d[s] = 0                                   # labels
    pred = Vector{Int64}(undef, calV); pred[s] = 0                  # predecessor
    
    #* algorithm
    while true
        # distance update
        for ind in pnt[i]:1:(pnt[i+1]-1)
            j = A[ind, 2]
            c_ij = A[ind, 3]
            if d[j] > d[i] + c_ij
                d[j] = d[i] + c_ij
                pred[j] = i
            end
        end
        # terminal condition
        if length(S) == calV
            break
        end
        # node selection & update S, Sbar
        i = argmin(d[Sbar]); i = popat!(Sbar, i)
        push!(S, i)
    end

    return d, pred
end

end

module LabelCorrecting

using DataStructures

export BellmanFordAlgorithm, ModifiedBellmanFordAlgorithm

#* name: BellmanFordAlgorithm
# inputs:
#   - A: forward star representation of a graph, (tail, head, cost)
#   - s: source node index, default = 1
# outputs:
#   - d: distance from the source of each node
#   - pred: predecssor of each node
function BellmanFordAlgorithm(A, s::Int64 = 1)
    # preprocessing
    calV = maximum(A[:, 1:2])                                   # cardinality of the vertices set V (or, the number of the vertices)
    calE = size(A)[1]                                           # cardinality of the edges set E (or, the number of the edges)
    pnt = [findfirst(x->x==v, A[:, 1]) for v in 1:1:(calV+1)]   # point of each node
    pnt[calV+1] = calE+1
    tpv1 = findall(x->isnothing(x), pnt)
    tpv2 = deleteat!(collect(1:1:(calV+1)), tpv1)
    for v in tpv1
        pnt[v] = pnt[tpv2[findfirst(x->x>v, tpv2)]]
    end

    #* initialization
    d = fill(Inf, calV); d[s] = 0                               # labels
    pred = Vector{Int64}(undef, calV); pred[s] = 0              # predecessor
    flag_nc = false
    #* algorithm
    while true
        # check optimality condition 
        tpv3 = d[A[:, 1]] .+ A[:, 3] .- d[A[:, 2]]
        tpv4 = findall(x->x<0, tpv3)
        if isempty(tpv4)
            break
        end
        # update label
        for ind in tpv4
            i = A[ind, 1]
            j = A[ind, 2]
            c_ij = A[ind, 3]
            d[j] = d[i] + c_ij
            # negative cycle detection
            if d[j] < neg_thr
                flag_nc = true
                break
            end
            pred[j] = i
        end
    end

    return d, pred, flag_nc
end

#* name: ModifiedBellmanFordAlgorithm
# inputs:
#   - A: forward star representation of a graph, (tail, head, cost)
#   - s: source node index, default = 1
# outputs:
#   - d: distance from the source of each node
#   - pred: predecssor of each node
function ModifiedBellmanFordAlgorithm(A, s::Int64 = 1, use_dequeue::Bool = true)
    #* preprocessing
    calV = maximum(A[:, 1:2])                                   # cardinality of the vertices set V (or, the number of the vertices)
    calE = size(A)[1]                                           # cardinality of the edges set E (or, the number of the edges)
    pnt = [findfirst(x->x==v, A[:, 1]) for v in 1:1:(calV+1)]   # point of each node
    pnt[calV+1] = calE+1
    tpv1 = findall(x->isnothing(x), pnt)
    tpv2 = deleteat!(collect(1:1:(calV+1)), tpv1)
    for v in tpv1
        pnt[v] = pnt[tpv2[findfirst(x->x>v, tpv2)]]
    end

    #* initialization
    neg_thr = sum(A[:, 3] .* (A[:, 3] .< 0))                    # threshold for negative cycle detection
    d = fill(Inf, calV); d[s] = 0                               # labels
    pred = Vector{Int64}(undef, calV); pred[s] = 0              # predecessor
    if use_dequeue
        LIST = Deque{Int64}(); push!(LIST, s)
        LIST_hist = Queue{Int64}(); enqueue!(LIST_hist, s)
    else
        LIST = Queue{Int64}(); enqueue!(LIST, s) 
    end
    flag_nc = false

    #* algorithm
    n_ex = 0    # the number of node examinations
    time::Float64 = 0
    if use_dequeue
        time = @elapsed while !isempty(LIST)
            i = popfirst!(LIST)
            for ind in pnt[i]:1:(pnt[i+1]-1)
                j = A[ind, 2]
                c_ij = A[ind, 3]
                n_ex += 1
                # check optimality condition
                if d[j] > d[i] + c_ij
                    # update label
                    d[j] = d[i] + c_ij
                    # negative cycle detection
                    if d[j] < neg_thr
                        flag_nc = true
                        break
                    end
                    pred[j] = i
                    if !(j in LIST)
                        if j in LIST_hist
                            pushfirst!(LIST, j)
                        else
                            push!(LIST, j)
                            enqueue!(LIST_hist, j)
                        end
                    end
                end
            end
        end
    else
        time = @elapsed while !isempty(LIST)
            i = dequeue!(LIST)
            for ind in pnt[i]:1:(pnt[i+1]-1)
                j = A[ind, 2]
                c_ij = A[ind, 3]
                n_ex += 1
                # check optimality condition
                if d[j] > d[i] + c_ij 
                    # update label
                    d[j] = d[i] + c_ij
                    # negative cycle detection
                    if d[j] < neg_thr
                        flag_nc = true
                        break
                    end
                    pred[j] = i
                    if !(j in LIST)
                        enqueue!(LIST, j)
                    end

                end
            end
        end
    end

    return d, pred, n_ex, time, flag_nc
end

end