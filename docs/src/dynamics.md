# Dynamics

Each point of a [`FeasibleCommunities`](@ref) polytope is a Lotka–Volterra model
`dx/dt = x (r + A x)` with equilibrium `x`, and its community matrix `J = diag(x) A` describes how
small perturbations of that equilibrium evolve. This page measures the transient response from
`J` alone, then integrates the dynamics to follow large perturbations and unstable communities.

## Transients from the community matrix

A perturbation `δ` of a stable equilibrium evolves as `exp(tJ) δ` while it is small. The
[`return_time`](@ref) `-1 / λ`, from the [`spectral_abscissa`](@ref) `λ`, sets its decay in the
long run. Before that, perturbations can grow: the [`amplification_envelope`](@ref)
`‖exp(tJ)‖` bounds the growth of every perturbation at time `t`, starts with slope equal to the
[`numerical_abscissa`](@ref), and reaches its [`max_amplification`](@ref) at some time before it
decays (Neubert & Caswell 1997).

```@example dynamics
using FoodWebs, PolytopeSamplers, Random, CairoMakie

rng = Xoshiro(3)
web = niche_network(rng, 15, 0.15; basal_paths = true)
x = exp10.(2 .* rand(rng, 15) .- 1)
fc = FeasibleCommunities(web, x; efficiency = (0.0, 0.3))
samples = sample_communities(rng, fc; target_ess = 200)
communities = samples.communities
stable = filter(c -> spectral_abscissa(c.J) < 0, communities[1:10:end])
length(stable), length(communities[1:10:end])
```

```@example dynamics
peaks = [max_amplification(c.J) for c in stable]
fig = Figure(size = (800, 350))
ax1 = Axis(fig[1, 1]; xlabel = "numerical abscissa", ylabel = "maximum amplification", xscale = log10, yscale = log10)
scatter!(ax1, [numerical_abscissa(c.J) for c in stable], [p.amplification for p in peaks]; markersize = 6)
ax2 = Axis(fig[1, 2]; xlabel = "return time", ylabel = "time of the maximum", xscale = log10, yscale = log10)
scatter!(ax2, [return_time(c.J) for c in stable], [p.time for p in peaks]; markersize = 6)
fig
```

Every stable community here is reactive, and its envelope peaks within a third of its return
time, often several orders of magnitude sooner. The numerical abscissa sets the early growth and
the spectral abscissa the late decay of the same response.

## Integrating the dynamics

Loading SciMLBase.jl enables [`lotka_volterra_problem`](@ref), [`lotka_volterra`](@ref) and
[`pulse_response`](@ref). They take any SciML solver: `Tsit5()` from OrdinaryDiffEqTsit5.jl for
nonstiff problems, or `Rodas5P()` from OrdinaryDiffEqRosenbrock.jl, which uses the analytic
Jacobian and suits webs whose biomasses span orders of magnitude.

[`pulse_response`](@ref) follows a pulse away from the equilibrium and compares it with the
linearization. A small pulse follows `‖exp(tJ) δ‖`, below the envelope:

```@example dynamics
using OrdinaryDiffEqRosenbrock: Rodas5P

c = stable[1]
T = return_time(c.J)
times = exp10.(range(log10(T) - 4, log10(T) + 1; length = 200))
δ = 1.0e-3 .* x .* randn(rng, 15)
p = pulse_response(Rodas5P(), c, x, δ, times)
fig = Figure(size = (600, 350))
ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "‖x(t) − x̂‖ / ‖δ‖", xscale = log10, yscale = log10)
lines!(ax, times, [amplification_envelope(c.J, t) for t in times]; label = "envelope", color = :gray)
lines!(ax, p.t, p.linear; label = "linearization", linewidth = 3)
lines!(ax, p.t, p.deviation; label = "Lotka–Volterra", linestyle = :dash, color = :black)
axislegend(ax; position = :lb)
fig
```

[`lotka_volterra`](@ref) integrates from any start and removes a species when its abundance
falls below `extinction_threshold` while declining. From a start near the equilibrium of an
unstable community, the dynamics move to another state, here with fewer species:

```@example dynamics
unstable = communities[findfirst(c -> spectral_abscissa(c.J) > 0, communities)]
s = lotka_volterra(Rodas5P(), unstable, x .* (1 .+ 1.0e-3 .* randn(rng, 15)), (0.0, 2.0e3))
fig = Figure(size = (600, 350))
ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "abundance / equilibrium", yscale = log10)
for i in 1:15
    alive = s.x[i, :] .> 0
    lines!(ax, s.t[alive], s.x[i, alive] ./ x[i])
end
fig
```

```@example dynamics
count(s.survivors), filter(isfinite, s.extinction_times)
```

## Saturating functional responses

[`type2_model`](@ref) turns a point of the polytope into a multi-prey Holling type II model
([`TypeIIModel`](@ref)) with the same equilibrium and link fluxes. Each predator gets a handling
time; its attack rates rise to `λ / (1 - τ)`, where `τ` is the fraction of its time spent
handling prey at the equilibrium ([`handling_fractions`](@ref)), so the fluxes stay the same, as
do the growth rates, self-regulation and efficiencies. Saturation changes the Jacobian
([`community_matrix`](@ref community_matrix(::TypeIIModel))): a satiated predator takes a smaller
share of a growing prey, which gives each prey positive density dependence and a positive effect
on the other prey of its predators, and the predator responds less to its prey.

`match_jacobian = true` raises the self-regulation and efficiencies to cancel all but the effects
between prey of a shared predator. Across handling fractions, saturation alone destabilizes the
stable communities sampled above within `τ = 0.1`. With the compensation their spectral abscissa
barely moves, except in one community that the remaining effects between prey destabilize beyond
`τ = 0.3`:

```@example dynamics
points = samples.points[:, 1:10:end]
θs = [θ for (θ, c) in zip(eachcol(points), communities[1:10:end]) if spectral_abscissa(c.J) < 0][1:5]
τs = 0:0.05:0.6
abscissa(θ, τ; kwargs...) = spectral_abscissa(community_matrix(type2_model(fc, θ; handling_fraction = τ, kwargs...)))
fig = Figure(size = (800, 350))
ax1 = Axis(fig[1, 1]; xlabel = "handling fraction τ", ylabel = "spectral abscissa", title = "saturation")
ax2 = Axis(fig[1, 2]; xlabel = "handling fraction τ", title = "saturation, compensated")
for θ in θs
    lines!(ax1, τs, [abscissa(θ, τ) for τ in τs])
    lines!(ax2, τs, [abscissa(θ, τ; match_jacobian = true) for τ in τs])
end
fig
```

[`type2_dynamics`](@ref) integrates a type II model with the extinctions of
[`lotka_volterra`](@ref), and [`pulse_response`](@ref) follows a pulse away from its equilibrium.

```@docs
return_time
amplification_envelope
max_amplification
lotka_volterra_problem
lotka_volterra
pulse_response
```

```@docs
TypeIIModel
type2_model
community_matrix(::TypeIIModel)
handling_fractions
type2_problem
type2_dynamics
```
