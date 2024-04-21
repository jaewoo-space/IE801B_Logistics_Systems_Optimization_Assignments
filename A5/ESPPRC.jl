module Solvers

using LinearAlgebra: I

export Feillet

const EPS = 1e-6

#* name: is_dominated
# description: check whether label_exm is dominated by label_cmp
# inputs:
#   - label_exm: label to examin
#   - label_cmp: label to compare
#   - label_length: the length of a label
# outputs:
#   - whether label_exm is dominated by label_cmp or not
function is_dominated(label_exm, label_cmp, label_length)

if sum((label_cmp .- label_exm) .< EPS) == label_length
    return true
else
    return false
end

end

#* name: nondominated_labels
# description: find dominated labels among labels
# inputs:
#   - labels: the list of the labels to check
#   - label_length: the length of a label
# outputs:
#   - nondominated labels
function nondominated_labels(labels, label_length)

n_labels = length(labels)

labels_nd::Array{Bool} = fill(true, n_labels)

for lInd1 in 1:1:(n_labels-1)
    if labels_nd[lInd1]
        for lInd2 in (lInd1+1):1:n_labels
            if labels_nd[lInd2]
                if is_dominated(labels[lInd1], labels[lInd2], label_length)
                    labels_nd[lInd1] = false
                    break
                end
                if is_dominated(labels[lInd2], labels[lInd1], label_length)
                    labels_nd[lInd2] = false
                end
            end
        end
    end
end

return labels[labels_nd]

end

#* name: Feillet
# description: algorithm for solving ESPPRC, proposed by Feillet et al. (2004)
# inputs:
#   - capacity: vehicle capacity
#   - demand: demand at each node
#   - depot: depot node index
#   - customers: list of customer indices
#   - weights: arc weight matrix
# outputs:
#   - non-dominated labels list at the destination
#   - the number of node examinations
function Feillet(capacity, demand, depot, customers, weights)

N = length(customers)

## 1) initialization
# define containers and supporting variables
labels_list::Array{Array{Tuple}} = [[] for _ in 1:1:(N+1)]
push!(labels_list[depot], (0.0, 0, zeros(Bool, N)..., 0.0))
E = [depot]
F_ij::Array{Tuple} = []
eye = I(N)
e = [eye[:, i] for i in 1:1:N]

## 2) algorithm
exm_cnt = 1
while !isempty(E)
    # a. choose a node to examine
    i = popfirst!(E)
    println("============ # Examinations: $exm_cnt, node = $i ==========")
    # b. extend labels of each successor node
    for j in customers
        # if i -> j is an available arc
        if (j != i) && (weights[i, j] < Inf)
            F_ij = []
            for label in labels_list[i]
                if !label[1+j]
                    new_label = (label[1]+demand[j], label[2]+1, Bool.(label[3:(end-1)].+e[j-1])..., round(label[end]+weights[i, j], digits=6))
                    # capacity constraint
                    if new_label[1] <= capacity
                        push!(F_ij, new_label)
                    end
                end
            end
            # find non-dominated labels
            labels_nd = nondominated_labels(vcat(labels_list[j], F_ij), N+3)
            
            # if the non-dominated labels set of a node is updated, add the node to the examination waiting list
            if Set(labels_list[j]) != Set(labels_nd)
                labels_list[j] = labels_nd
                if !(j in E) && (j != customers[end])
                    push!(E, j)
                end
            end
        end
    end

    println("Nodes to examin: $E")
    tpv = length(labels_list[customers[end]])
    println("# labels at the destination: $tpv")
    exm_cnt += 1
end

return labels_list[customers[end]], exm_cnt

end

end