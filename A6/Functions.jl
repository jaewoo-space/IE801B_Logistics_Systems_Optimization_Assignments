module Tools

using Plots

export haversine, VisualizeResults

const RE = 6371e03
const D2R = pi/180

function haversine(coordinate1, coordinate2)

lat1, lon1 = D2R .* coordinate1
lat2, lon2 = D2R .* coordinate2

a = sin((lat2-lat1)/2)^2 + cos(lat1) * cos(lat2) * sin((lon2-lon1)/2)^2
c = 2 * atan(sqrt(a), sqrt(1-a))
d = RE * c

return d

end

function VisualizeResults(facilities, customers, x, y)

N = size(facilities)[1]
M = size(customers)[1]


fig1 = plot(dpi=300)
scatter!(facilities[:, 1], facilities[:, 2], color =:red, markersize = 3, label="Facilities")
scatter!(customers[:, 1], customers[:, 2], color =:blue, markersize = 3, label="Customers", legend=:outerbottomright)
    
fig2 = plot(dpi=300)
scatter!(facilities[:, 1], facilities[:, 2], color =:red, markersize = 3, label="Closed facilities")
scatter!(customers[:, 1], customers[:, 2], color =:blue, markersize = 3, label="Customers")
scatter!(facilities[findall(x->x>0, y), 1], facilities[findall(x->x>0, y), 2], color =:yellow, marker =:star, markersize = 6, label="Open facilities", legend=:outerbottomright)
for i in 1:1:N
    for j in 1:1:M
        if x[i, j] > 0
            plot!([customers[j, 1], facilities[i, 1]], [customers[j, 2], facilities[i, 2]], color =:black, linestyle =:dash, label=false)
        end
    end
end

return fig1, fig2

end

end


module Solvers

export Klincewicz, Exact

using JuMP, HiGHS

const EPS = 1e-5

#!
function add_heuristic(a, b, c, F, N, M)

facilities = Set(1:1:N)
K = Set([])
Kc = Set(facilities)
Zu = Inf
x_best::Array{Int} = zeros(Int, N, M)
y_best::Array{Int} = zeros(Int, N)
x_curr::Array{Int} = zeros(Int, N, M)
y_curr::Array{Int} = zeros(Int, N)

while true
    ## Step 1
    w = zeros(Float64, N, M)
    R = fill(-Inf, N)
    for i in Kc
        for j in 1:1:M
            w[i, j] = max(minimum(c[collect(K), j].- c[i, j], init=0), 0)
        end
        Omega = sum(w[i, :])
        R[i] = Omega * min(b[i]/sum(a[w[i, :].>0], init=0), 1) - F[i]
    end

    ## Step 2
    i_add = argmax(R)
    push!(K, i_add)
    Kc = setdiff(facilities, K)
    K_list = collect(K)

    ## Step 3
    if (sum(b[K_list]) - sum(a)) < EPS
        continue
    end
    
    ## Step 4
    if length(K_list) == 1
        order = sortperm([minimum(c[K_list, j]) for j in 1:1:M], rev=true)
    else
        cost_first_second = [sort(c[K_list, j])[1:2] for j in 1:1:M]
        cost_diff = [cost[2] for cost in cost_first_second] .- [cost[1] for cost in cost_first_second]
        order = sortperm(cost_diff, rev=true)
    end

    ## Step 5
    remaining_capacity = b[K_list]
    TC = 0.0
    flag = 0
    x_curr = zeros(Int, N, M)
    for j in order
        tpv = sortperm(c[K_list, j])
        for i in tpv
            if a[j] - remaining_capacity[i] < EPS
                x_curr[K_list[i], j] = 1
                remaining_capacity[i] -= a[j]
                flag = 1
                break
            end
        end

        if flag == 0
            flag = 2
            break
        end

        flag = 0

    end

    if flag == 2
        continue
    end

    ## Step 6
    y_curr = [maximum(x_curr[i, :]) for i in 1:1:N]
    TC = sum(c.*x_curr) + sum(F.*y_curr)
    if TC - Zu < -EPS
        x_best = x_curr
        y_best = y_curr
        Zu = TC
    else
        break
    end
    
    K_ban = K_list[b[K_list] .- remaining_capacity .< EPS]

    if !isempty(K_ban)
        setdiff!(facilities, Set(K_ban))
    end
end

return Zu, x_best, y_best

end

#!
function final_heuristic(a, b, c, N, M, x, y)

K = findall(x -> x > EPS, y)
x_mod::Array{Int} = copy(x)

flag = true
while flag
    flag = false

    cost_first_second = [sort(c[K, j])[1:2] for j in 1:1:M]
    cost_diff = [cost[2] for cost in cost_first_second] .- [cost[1] for cost in cost_first_second]
    remaining_capacity = b[K] .- [sum(a.*x[i, :]) for i in K]
    order = sortperm(cost_diff, rev=true)
    for j in order
        can_move = K[findall(x -> x < EPS, a[j] .- remaining_capacity)]
        min_move = can_move[argmin(c[can_move, j])]
        if c[min_move, j] - sum(c[:, j] .* x[:, j]) < -EPS
            x_mod[:, j] = zeros(Int, N)
            x_mod[min_move, j] = 1
            flag = true
        end
    end
end

return x_mod

end

#!
function solve_dual(a, b, c, F, N, M, lambda)

model = JuMP.Model(HiGHS.Optimizer)

model = Model(HiGHS.Optimizer)
set_silent(model)

@variable(model, x[i = 1:N, j=1:M], Bin)
@variable(model, y[i = 1:N], Bin)

@objective(model, Min, sum(F.*y) + sum((c .+ lambda .* a') .* x) - sum(lambda.*b))

@constraint(model, single_source[j=1:M], sum(x[:, j]) == 1)
@constraint(model, open_facility[i=1:N, j=1:M], x[i, j] <= y[i])

optimize!(model)

return objective_value(model), round.(Int, value.(x)), round.(Int, value.(y))

end

#!
function is_feasible(a, b, N, x)

if isempty(findall(x->x>EPS, sum(a' .* x, dims=2) .- b))
    return true
else
    return false
end

end

#!
function Klincewicz(a, b, c, F, w = 0.25, epsilon = 1e-3, max_solve = 200)

N::Int = length(b)
M::Int = length(a)
x_best::Array{Int} = zeros(Int, N, M)
y_best::Array{Int} = zeros(Int, N)
x_dual::Array{Int} = zeros(Int, N, M)
y_dual::Array{Int} = zeros(Int, N)
final_heuristic_flag::Bool = false
flag0::Bool = true
dual_solve_cnt::Int = 0
U::Array{Int} = collect(1:1:N)
capacity_constraint::Array{Float64} = zeros(Float64, N)
Zl_list::Array{Float64} = []
Zu_list::Array{Float64} = []

## Step 0: initialize
# calculate upper bound 
Zu, x_best, y_best = add_heuristic(a, b, c, F, N, M)
# calculate lower bound
lambda = zeros(Float64, N)
Zl, x_dual, y_dual = solve_dual(a, b, c, F, N, M, lambda)
push!(Zl_list, Zl)
push!(Zu_list, Zu)
println("Zl: $Zl, Zu: $Zu")
# if relaxed problem is feasible, stop
if is_feasible(a, b, N, x_dual)
    final_heuristic_flag = true
    flag0 = false
    x_best = x_dual
    y_best = y_dual
    obj = sum(x_best.*c) + sum(y_best.*F)
    return obj, x_dual, y_dual
# if not, update lambda
else
    capacity_constraint = reduce(vcat, sum(a'.*x_dual, dims=2)) .- b
    lambda[U] = [max(lambda[i]+w*(Zu-Zl)*(capacity_constraint[i])/sqrt(sum(capacity_constraint.^2)), 0) for i in U]
end

while flag0
    ## Step 1
    #  solve relaxed problem
    Zl_new, x_dual, y_dual = solve_dual(a, b, c, F, N, M, lambda)
    dual_solve_cnt += 1
    # if necessary, update Zl
    if Zl - Zl_new < -EPS
        Zl = Zl_new
        println("Zl: $Zl, Zu: $Zu")
    end
    push!(Zl_list, Zl)
        
    ## Step 2: if Step 1 obtains an infeasible solution
    if !is_feasible(a, b, N, x_dual)
        if dual_solve_cnt == max_solve
            break
        end
    else
        # initialize counter
        dual_solve_cnt = 0
        # update Zu if necessary
        Zu_new = sum(x_dual.*c) + sum(y_dual.*F)
        if Zu_new - Zu < -EPS
            Zu = Zu_new
            x_best = x_dual
            y_best = y_dual
            final_heuristic_flag = true
            println("Zl: $Zl, Zu: $Zu")
        ## Step 4: if Step 1 obtains a feasible solution
        # if a better soluton already obtained
        else
            push!(Zu_list, Zu)
            break
        end
        # if the gap is tight enough,
        if Zu/Zl <= 1 + epsilon
            push!(Zu_list, Zu)
            break
        end
        # otherwise,unmark multipliers
        U = collect(1:1:N)
    end
    push!(Zu_list, Zu)

    ## Step 3
    # constraint violation
    capacity_constraint = reduce(vcat, sum(a'.*x_dual, dims=2)) .- b
    # not marked (U) or constraint violated... multipliers to update
    i_to_update = collect(union(Set(U), Set(findall(x -> x > EPS, capacity_constraint))))
    # update multipliers
    lambda_new = copy(lambda)
    lambda_new[i_to_update] = [max(lambda_new[i]+w*(Zu-Zl)*(capacity_constraint[i])/sqrt(sum(capacity_constraint[i_to_update].^2)), 0) for i in i_to_update]
    lambda = lambda_new
    # mark decreased multipliers
    U = collect(setdiff(Set(U), Set(findall(x -> x < -EPS, capacity_constraint))))
end

## Step 5
if final_heuristic_flag
    x_best = final_heuristic(a, b, c, N, M, x_best, y_best)
    y_best = [maximum(x_best[i, :]) for i in 1:1:N]
end

obj = sum(x_best.*c) + sum(y_best.*F)

return obj, x_best, y_best, Zl_list, Zu_list

end

#!

function Exact(a, b, c, F)

N::Int = length(b)
M::Int = length(a)

model = JuMP.Model(HiGHS.Optimizer)

@variable(model, x[i = 1:N, j=1:M], Bin)
@variable(model, y[i = 1:N], Bin)

@objective(model, Min, sum(F.*y) + sum(c .* x))

@constraint(model, capacity[i=1:N], sum(a.*x[i, :]) <= b[i])
@constraint(model, single_source[j=1:M], sum(x[:, j]) == 1)
@constraint(model, open_facility[i=1:N, j=1:M], x[i, j] <= y[i])

optimize!(model)

return objective_value(model), round.(Int, value.(x)), round.(Int, value.(y))

end

end