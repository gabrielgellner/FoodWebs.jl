module FoodWebsSciMLBaseExt

using FoodWebs: FoodWebs, lotka_volterra
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

function FoodWebs.lotka_volterra(
        alg, A::AbstractMatrix, r::AbstractVector, x0::AbstractVector, tspan;
        extinction_threshold = 1.0e-9, kwargs...
    )
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
            r[i] + sum(A[i, j] * x[j] for j in 1:S if integrator.p.present[j]) < 0 || continue
            x[i] = 0
            integrator.p.present[i] = false
            extinction_times[i] = integrator.t
        end
        return nothing
    end
    callback = VectorContinuousCallback(condition, extinguish!, S; save_positions = (false, false))
    prob = FoodWebs.lotka_volterra_problem(A, r, start, tspan)
    solution = SciMLBase.solve(prob, alg; callback, kwargs...)
    x = reduce(hcat, solution.u)
    for i in 1:S
        x[i, solution.t .>= extinction_times[i]] .= 0
    end
    return (
        t = solution.t, x, survivors = isinf.(extinction_times), extinction_times,
        solution.retcode, solution,
    )
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
    x0 = x̂ .+ δ
    any(<(0), x0) && throw(ArgumentError("the pulse makes an abundance negative"))
    trajectory = lotka_volterra(alg, A, r, x0, (0.0, float(maximum(times))); saveat = times, kwargs...)
    magnitude = norm(δ)
    deviation = [norm(x .- x̂) for x in eachcol(trajectory.x)] ./ magnitude
    J = Diagonal(x̂) * A
    linear = [norm(exp(t * J) * δ) for t in trajectory.t] ./ magnitude
    k = argmax(deviation)
    return (; trajectory.t, deviation, linear, amplification = deviation[k], time = trajectory.t[k], trajectory)
end
function FoodWebs.pulse_response(alg, community::NamedTuple, x̂::AbstractVector, δ::AbstractVector, times; kwargs...)
    return FoodWebs.pulse_response(alg, community.A, community.r, x̂, δ, times; kwargs...)
end

end
