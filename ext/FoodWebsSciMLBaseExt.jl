module FoodWebsSciMLBaseExt

using FoodWebs: FoodWebs, TypeIIModel, community_matrix, lotka_volterra, type2_dynamics
using LinearAlgebra: Diagonal, mul!, norm
using SciMLBase: SciMLBase, ODEFunction, ODEProblem, VectorContinuousCallback

# Absent species (`p.present` false) take no part in the dynamics: their rate is 0 and their
# abundance enters no other species' rate. Implicit solvers can leave them at rounding errors
# instead of exactly 0.
function lotka_volterra_rhs!(dx, x, p, t)
    mul!(dx, p.A, p.present .* x)
    dx .= p.present .* x .* (p.r .+ dx)
    return nothing
end

# diag(r + A y) + diag(y) A diag(present), with y = present .* x
function lotka_volterra_jacobian!(J, x, p, t)
    y = p.present .* x
    J .= y .* p.A .* p.present'
    g = p.r .+ p.A * y
    for i in eachindex(x)
        p.present[i] && (J[i, i] += g[i])
    end
    return nothing
end

function FoodWebs.lotka_volterra_problem(A::AbstractMatrix, r::AbstractVector, x0::AbstractVector, tspan)
    S = length(x0)
    size(A) == (S, S) && length(r) == S ||
        throw(DimensionMismatch("A must be $S × $S and r of length $S for $S species, got $(size(A)) and $(length(r))"))
    f = ODEFunction(lotka_volterra_rhs!; jac = lotka_volterra_jacobian!)
    p = (A = float.(A), r = float.(collect(r)), present = collect(x0 .> 0))
    return ODEProblem(f, float.(collect(x0)), float.(tspan), p)
end
function FoodWebs.lotka_volterra_problem(community::NamedTuple, x0::AbstractVector, tspan)
    return FoodWebs.lotka_volterra_problem(community.A, community.r, x0, tspan)
end

# Integrate `problem(start)`, the ODEProblem from the abundances above their thresholds, setting a
# species to 0 when it falls below its threshold while declining. `per_capita(x, p, i)` is species
# i's per capita growth rate at x.
function integrate_with_extinctions(alg, problem, per_capita, x0, tspan, extinction_threshold; kwargs...)
    S = length(x0)
    threshold = extinction_threshold isa Real ? fill(float(extinction_threshold), S) : float.(collect(extinction_threshold))
    length(threshold) == S ||
        throw(DimensionMismatch("expected $S extinction thresholds, got $(length(threshold))"))
    any(<(0), x0) && throw(ArgumentError("abundances must be nonnegative"))
    start = [x0[i] > threshold[i] ? float(x0[i]) : 0.0 for i in 1:S]
    extinction_times = [start[i] > 0 ? Inf : float(first(tspan)) for i in 1:S]
    # a species crossing its threshold while declining is set to 0, where it stays; one crossing
    # it while growing is left alone. DiffEqBase passes the index of the event, or in newer
    # versions the crossing direction of every event, 0 for those that did not occur.
    condition(out, x, t, integrator) = (out .= x .- threshold)
    function extinguish!(integrator, events)
        x = integrator.u
        for i in (events isa Integer ? (events,) : findall(!=(0), events))
            per_capita(x, integrator.p, i) < 0 || continue
            x[i] = 0
            integrator.p.present[i] = false
            extinction_times[i] = integrator.t
        end
        return nothing
    end
    callback = VectorContinuousCallback(condition, extinguish!, S; save_positions = (false, false))
    solution = SciMLBase.solve(problem(start), alg; callback, kwargs...)
    x = reduce(hcat, solution.u)
    for i in 1:S
        x[i, solution.t .>= extinction_times[i]] .= 0
    end
    return (
        t = solution.t, x, survivors = isinf.(extinction_times), extinction_times,
        solution.retcode, solution,
    )
end

function FoodWebs.lotka_volterra(
        alg, A::AbstractMatrix, r::AbstractVector, x0::AbstractVector, tspan;
        extinction_threshold = 1.0e-9, kwargs...
    )
    S = length(x0)
    per_capita(x, p, i) = r[i] + sum(A[i, j] * x[j] for j in 1:S if p.present[j])
    problem(start) = FoodWebs.lotka_volterra_problem(A, r, start, tspan)
    return integrate_with_extinctions(alg, problem, per_capita, x0, tspan, extinction_threshold; kwargs...)
end
function FoodWebs.lotka_volterra(alg, community::NamedTuple, x0::AbstractVector, tspan; kwargs...)
    return lotka_volterra(alg, community.A, community.r, x0, tspan; kwargs...)
end

function FoodWebs.pulse_response(
        alg, A::AbstractMatrix, r::AbstractVector, x̂::AbstractVector, δ::AbstractVector, times; kwargs...
    )
    length(δ) == length(x̂) || throw(DimensionMismatch("δ has length $(length(δ)), x̂ has $(length(x̂))"))
    # each species' growth rate at x̂ is zero, relative to the size of its terms
    scale = abs.(r) .+ abs.(A) * abs.(x̂)
    all(abs.(r .+ A * x̂) .<= sqrt(eps()) .* scale) ||
        throw(ArgumentError("x̂ is not an equilibrium: r + A x̂ is not zero"))
    integrate(x0, tmax) = lotka_volterra(alg, A, r, x0, (0.0, tmax); saveat = times, kwargs...)
    return response(integrate, Diagonal(x̂) * A, x̂, δ, times)
end

# The deviation from x̂ after the pulse δ, integrated by `integrate(x0, tmax)`, against the
# linearization with Jacobian J
function response(integrate, J, x̂, δ, times)
    x0 = x̂ .+ δ
    any(<(0), x0) && throw(ArgumentError("the pulse makes an abundance negative"))
    trajectory = integrate(x0, float(maximum(times)))
    magnitude = norm(δ)
    deviation = [norm(x .- x̂) for x in eachcol(trajectory.x)] ./ magnitude
    linear = [norm(exp(t * J) * δ) for t in trajectory.t] ./ magnitude
    k = argmax(deviation)
    return (; trajectory.t, deviation, linear, amplification = deviation[k], time = trajectory.t[k], trajectory)
end
function FoodWebs.pulse_response(alg, community::NamedTuple, x̂::AbstractVector, δ::AbstractVector, times; kwargs...)
    return FoodWebs.pulse_response(alg, community.A, community.r, x̂, δ, times; kwargs...)
end

# The type II model with absent species (`p.present` false) taking no part, as for Lotka–Volterra
function type2_rhs!(dx, x, p, t)
    y = p.present .* x
    FoodWebs.type2_per_capita!(dx, y, p.model, p.present)
    dx .= y .* dx
    return nothing
end

# diag(present g) + diag(y) ∂g/∂x diag(present), with y = present .* x
function type2_jacobian!(J, x, p, t)
    y = p.present .* x
    FoodWebs.type2_per_capita_jacobian!(J, y, p.model, p.present)
    J .= y .* J .* p.present'
    g = FoodWebs.type2_per_capita!(similar(y), y, p.model, p.present)
    for i in eachindex(x)
        p.present[i] && (J[i, i] += g[i])
    end
    return nothing
end

function FoodWebs.type2_problem(m::TypeIIModel, x0::AbstractVector, tspan)
    length(x0) == length(m.x) || throw(DimensionMismatch("x0 must have length $(length(m.x))"))
    f = ODEFunction(type2_rhs!; jac = type2_jacobian!)
    p = (model = m, present = collect(x0 .> 0))
    return ODEProblem(f, float.(collect(x0)), float.(tspan), p)
end

function FoodWebs.type2_dynamics(alg, m::TypeIIModel, x0::AbstractVector, tspan; extinction_threshold = 1.0e-9, kwargs...)
    g = zeros(length(m.x))
    per_capita(x, p, i) = FoodWebs.type2_per_capita!(g, p.present .* x, m, p.present)[i]
    problem(start) = FoodWebs.type2_problem(m, start, tspan)
    return integrate_with_extinctions(alg, problem, per_capita, x0, tspan, extinction_threshold; kwargs...)
end

function FoodWebs.pulse_response(alg, m::TypeIIModel, δ::AbstractVector, times; kwargs...)
    length(δ) == length(m.x) || throw(DimensionMismatch("δ has length $(length(δ)), the model has $(length(m.x)) species"))
    integrate(x0, tmax) = type2_dynamics(alg, m, x0, (0.0, tmax); saveat = times, kwargs...)
    return response(integrate, community_matrix(m), m.x, δ, times)
end

end
