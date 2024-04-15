module Solvers

using LinearAlgebra: I

export Feillet

const EPS = 1e-6

function dominated_labels(labels, label_length)

labels_dom::Array{Tuple} = []

for label_check in labels
    for label_comp in labels
        if (sum((label_comp .- label_check) .<= EPS) == label_length) && (sum((label_comp .- label_check) .< EPS) > 0)
            push!(labels_dom, label_check)
        end
    end
end

return labels_dom

end

function Feillet(capacity, demand, depot, customers, weights)

nodes = vcat(depot, customers)
N = length(nodes)

# initialization
labels_list::Array{Array{Tuple}} = [[] for _ in 1:1:N]
push!(labels_list[depot], (0.0, 0, zeros(Bool, N)..., 0.0))
E = [depot]

# supporting variables
F::Array{Tuple} = []
eye = I(N)
e = [eye[:, i] for i in 1:1:N]

while !isempty(E)
    i = pop!(E)
    for j in customers
        if j != i
            F = []
            for label in labels_list[i]
                if !label[2+j]
                    new_label = (label[1]+demand[j], label[2]+1, (label[3:(end-1)].+e[j])..., label[end]+weights[i, j])
                    push!(F, new_label)
                end
            end
            # EFF
            labels_dom = dominated_labels(vcat(labels_list[j], F), N+3)
            new_labels = setdiff(Set(F), Set(labels_dom))
            if !isempty(new_labels)
                labels_list[j] = vcat(labels_list[j], collect(new_labels))
                push!(E, j)
            end
        end
    end
end

return lables_list

end

end