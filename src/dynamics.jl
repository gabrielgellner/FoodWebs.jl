"""
    return_time(J)

The characteristic return time `-1 / λ` of an equilibrium with Jacobian `J`, where `λ` is the
[`spectral_abscissa`](@ref): the time over which a perturbation shrinks by a factor `e` in the
long run (Neubert & Caswell 1997). `Inf` when the equilibrium is not stable.
"""
function return_time(J::AbstractMatrix)
    λ = spectral_abscissa(J)
    return λ < 0 ? -1 / λ : Inf
end

"""
    amplification_envelope(J, t)

The largest growth `‖exp(tJ)‖₂` at time `t` of a perturbation of an equilibrium with Jacobian `J`,
over all perturbations (the amplification envelope of Neubert & Caswell 1997). Each
perturbation's own growth `‖exp(tJ) δ‖ / ‖δ‖` lies below it. Its slope at `t = 0` is the
[`numerical_abscissa`](@ref), and its decay rate as `t → ∞` the [`spectral_abscissa`](@ref).
"""
amplification_envelope(J::AbstractMatrix, t::Real) = opnorm(exp(t * Matrix(J)))

"""
    max_amplification(J)

The maximum of the [`amplification_envelope`](@ref) of a stable equilibrium with Jacobian `J`, and
the time it is reached: a named tuple `(amplification, time)` (Neubert & Caswell 1997). The
maximum is 1 at time 0 when the equilibrium is not reactive (`numerical_abscissa(J) ≤ 0`).
Otherwise it is found on a grid of times growing by 5% from `10⁻³` of the shorter of the return
time and `1 / numerical_abscissa(J)`, up to at least 20 return times and until the envelope is
below 1, and refined by golden-section search. Throws an `ArgumentError` if `J` is not stable.
"""
function max_amplification(J::AbstractMatrix)
    T = return_time(J)
    isfinite(T) || throw(ArgumentError("the equilibrium is not stable, so the amplification is unbounded"))
    ν = numerical_abscissa(J)
    ν <= 0 && return (amplification = 1.0, time = 0.0)
    M = Matrix(J)
    envelope(t) = opnorm(exp(t * M))
    times, values = [0.0], [1.0]
    t = 1.0e-3 * min(T, 1 / ν)
    while true
        push!(times, t)
        push!(values, envelope(t))
        (t >= 20T && values[end] < 1) && break
        t *= 1.05
    end
    k = argmax(values)
    a, b = times[max(k - 1, 1)], times[min(k + 1, end)]
    φ = (sqrt(5) - 1) / 2
    for _ in 1:60
        c, d = b - φ * (b - a), a + φ * (b - a)
        envelope(c) >= envelope(d) ? (b = d) : (a = c)
    end
    tmax = (a + b) / 2
    amplification = envelope(tmax)
    return amplification >= values[k] ? (; amplification, time = tmax) : (amplification = values[k], time = times[k])
end

"""
    lotka_volterra_problem(A, r, x0, tspan)
    lotka_volterra_problem(community, x0, tspan)

The `ODEProblem` of the Lotka–Volterra dynamics `dx/dt = x (r + A x)` from `x0` over `tspan`,
with its analytic Jacobian `diag(r + A x) + diag(x) A` for stiff solvers (load SciMLBase.jl to
enable this method). A `community` is a named tuple with fields `A` and `r`, such as one from
[`sample_communities`](@ref).
"""
function lotka_volterra_problem end

"""
    lotka_volterra(alg, A, r, x0, tspan; extinction_threshold = 1e-9, kwargs...)
    lotka_volterra(alg, community, x0, tspan; kwargs...)

Integrate the Lotka–Volterra dynamics `dx/dt = x (r + A x)` from `x0` over `tspan` with the
SciML solver `alg`, such as `Tsit5()` from OrdinaryDiffEqTsit5.jl or `Rodas5P()` from
OrdinaryDiffEqRosenbrock.jl (load SciMLBase.jl to enable this method). A species whose abundance
falls below `extinction_threshold` (a number, or one per species) goes extinct: it is set to 0
and stays there. Species that start at or below the threshold are absent from the start. Other
keywords, such as `saveat`, `abstol` and `reltol`, go to `solve`.

Returns a named tuple: the saved times `t`; the abundances `x`, one column per time;
`survivors`, whether each species is present at the end; `extinction_times`, the time each
species went extinct (`Inf` for survivors); `retcode`, the solver's return code; and the
`solution` itself.
"""
function lotka_volterra end

"""
    pulse_response(alg, A, r, x̂, δ, times; kwargs...)
    pulse_response(alg, community, x̂, δ, times; kwargs...)
    pulse_response(alg, m::TypeIIModel, δ, times; kwargs...)

The response of the Lotka–Volterra dynamics to a pulse `δ` away from the equilibrium `x̂`
(`r + A x̂ = 0`), integrated by [`lotka_volterra`](@ref) with solver `alg` and saved at `times`
(load SciMLBase.jl to enable this method). Keywords go to [`lotka_volterra`](@ref). For a
[`TypeIIModel`](@ref) the pulse is away from its equilibrium `m.x`, integrated by
[`type2_dynamics`](@ref), and the linearization uses its [`community_matrix`](@ref).

Returns a named tuple: `t`; `deviation`, the distance `‖x(t) - x̂‖ / ‖δ‖`; `linear`, the same
distance in the linearization, `‖exp(tJ) δ‖ / ‖δ‖` with `J = diag(x̂) A`; `amplification`, the
largest deviation; `time`, when it is reached; and the `trajectory` from [`lotka_volterra`](@ref).
For small pulses `deviation` follows `linear`, which lies below the
[`amplification_envelope`](@ref) of `J`. Throws an `ArgumentError` if `x̂` is not an equilibrium
or `x̂ + δ` has a negative abundance.
"""
function pulse_response end
