using JuMP, HiGHS

function f()
    model = Model(HiGHS.Optimizer)
    @variable(model, y >= 0)
    @constraint(model, y <= 1)
    @objective(model, Max, y)
    optimize!(model)
    return model
end


tpv = f()

value(all_variables(tpv)[1])