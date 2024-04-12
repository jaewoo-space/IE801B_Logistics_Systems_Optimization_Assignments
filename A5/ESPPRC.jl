module Solvers

using LinearAlgebra: I

function nondominate(labels)

return 0

end

function Feillet(capacity, demand, depot, customers, weights)

nodes = vcat(depot, customers)
N = length(nodes)

labels_list::Array{Array{Tuple}} = [[] for _ in 1:1:N]
push!(labels_list[depot], (0.0, 0, zeros(Int, N)..., 0.0))

E = [depot]
F::Array{Array{Tuple}} = []
eye = I(N)
e = [I[:, i] for i in 1:1:N]

while !isempty(E)
    i = pop!(E)
    for j in customers if j != i
        F = []
        for label in labels_list[i]
            if label[2+j] == 0
                new_label = (label[1]+demand[j], label[2]+1, (label[3:(end-1)].+e[j])..., label[end]+weights[i, j])
                push!(F, new_label)
            end
        end
        # EFF
        # new_labels = nondominate()
        # labels_list[j] = 
    end
end

return 0

end

end