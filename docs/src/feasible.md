# Feasible communities

Lotka–Volterra dynamics `dxᵢ/dt = xᵢ (rᵢ + Σⱼ Aᵢⱼ xⱼ)` have the equilibrium `x > 0` when
`r + A x = 0`. At a fixed `x` and web this is linear in the interaction strengths and in `r`, so
the parameters consistent with the equilibrium, with linear ecological constraints, form a
polytope. [`FeasibleCommunities`](@ref) describes it.

The coordinates are a loss and a gain magnitude for each feeding link, a self-regulation
magnitude for each self-regulated species, and the growth rates when they are given as intervals
(the default). Species that eat each other get a loss and a gain at each of their entries of
`A`, and cannibalism merges into self-regulation, so every web can be used.

```@example feasible
using FoodWebs, Random

rng = Xoshiro(2026)
web = generalized_cascade_network(rng, 15, 0.15; basal_paths = true)
x = steady_state_biomass(web; efficiency = 0.2)
fc = FeasibleCommunities(web, x; efficiency = (0.0, 0.3), bounds = (0.0, 5.0))
```

```@example feasible
dimension(fc), isbounded(fc), structurally_singular(fc), diagnose(fc)
```

## Constraints

- `rates`: fixed growth rates, or intervals by role (basal species in `[0, 5]` and consumers in
  `[-1, 0]` by default), which makes them coordinates.
- `self_regulation`: which species have a self-regulation term (`:all` by default, as a magnitude
  of 0 is allowed; `:all_but_top` reproduces Gellner, McCann & Hastings).
- `efficiency`: a band `e_min λ ≤ γ ≤ e_max λ` on each link's gain against its loss.
- `self_share`: self-regulation as a share of each species' turnover.
- `bounds`: a box on every magnitude.

[`halfspaces`](@ref) gives the polytope as plain arrays, `lower ≤ θ ≤ upper`, `A θ ≤ b`,
`C θ = d`, and [`interaction_matrix`](@ref), [`growth_rates`](@ref),
[`community_matrix`](@ref) and [`fluxes`](@ref) turn a point into a model.

## Structure without a solver

[`dimension`](@ref) counts the free directions, [`isbounded`](@ref) decides boundedness from the
web and the constraints (exactly, as the tests check against linear programs),
[`structurally_singular`](@ref) finds webs whose interaction matrix is singular at every point,
and [`diagnose`](@ref) lists violated necessary conditions, such as an energy budget that cannot
balance:

```@example feasible
chain = food_chain_network(2)
diagnose(FeasibleCommunities(chain, [1.0, 1.0]; rates = [0.1, -1.0], efficiency = 0.5))
```

Testing whether the polytope has any point, and sampling it, need a linear program and a sampler,
which the PolytopeSamplers.jl extension provides (see [Sampling](@ref)).

## Priors

A uniform distribution on the polytope is one prior. [`CommunityPrior`](@ref) puts a
distribution on each kind of coordinate, and [`logdensity`](@ref) gives the density restricted
to the polytope, which a sampler for non-uniform targets uses. A density on the growth rates
times a uniform on the interaction strengths conditions a joint prior on the equilibrium holding,
weighting each rate vector by the volume of interaction strengths it allows; drawing the rates
first and sampling the interactions on their slice does not, so the two priors differ. A uniform
distribution stays uniform under linear changes of coordinates, so uniform in the magnitudes is
uniform in the fluxes or in the entries of the community matrix.

## Sampling

Loading PolytopeSamplers.jl enables `PolytopeSamplers.Polytope(fc)`, `isempty(fc)` (a linear
program), [`sample_communities`](@ref) and [`feasible_growth_rates`](@ref).
[`sample_communities`](@ref) runs chains until the effective sample size reaches `target_ess`,
with the billiard walk for the uniform distribution and coordinate hit-and-run for a prior,
which was 2 to 3.5 times faster than CRHMC on webs of 63 and 210 dimensions.

```@example feasible
using PolytopeSamplers
using LinearAlgebra: eigvals

isempty(fc)
```

```@example feasible
res = sample_communities(rng, fc; target_ess = 200)
stable = count(c -> maximum(real, eigvals(c.J)) < 0, res.communities) / length(res.communities)
res.converged, stable
```

The published constrained food web studies fix the growth rates instead, drawing them by role
with [`draw_growth_rates`](@ref) and redrawing until the polytope has a point;
[`feasible_growth_rates`](@ref) runs this loop. The biomasses must be on a scale where the rates
divided by them fit within `bounds`:

```@example feasible
found = feasible_growth_rates(rng, web, 20 .* x; efficiency = (0.0, 0.3), bounds = (0.0, 5.0))
found.draws, dimension(found.communities)
```

## Three ensembles

The same web gives three ensembles of community matrices, each assuming more than the last:

1. [`random_predprey`](@ref) draws the link magnitudes independently, with a fixed diagonal and
   no equilibrium. Its stability depends on the web's sign structure and the magnitudes' moments
   (see [`elliptic_law`](@ref)).
2. Feasible communities, `J = diag(x) A` for points of `FeasibleCommunities(web, x)`, require
   the equilibrium `x`, which ties each row of `A` to the growth rate and biomass of its species.
3. Feasible and energetic communities add an efficiency band, so that no consumer gains more
   than a share of what its prey lose.

```@example feasible
using Distributions: Uniform

stability(Js) = count(J -> maximum(real, eigvals(J)) < 0, Js) / length(Js)
unconstrained = stability([random_predprey(rng, web, Uniform(0, 1), 0.3) for _ in 1:1000])
feasible = sample_communities(rng, FeasibleCommunities(web, x; bounds = (0.0, 5.0)); target_ess = 200)
energetic = sample_communities(rng, FeasibleCommunities(web, x; efficiency = (0.0, 0.3), bounds = (0.0, 5.0)); target_ess = 200)
(; unconstrained, feasible = stability([c.J for c in feasible.communities]),
    energetic = stability([c.J for c in energetic.communities]))
```

The gallery note `gallery/equilibria.qmd` compares the feasible ensembles across biomass
pyramids, from Eltonian to inverted, and on the observed biomasses of the Tuesday Lake webs.

```@docs
sample_communities
feasible_growth_rates
draw_growth_rates
```

## Reference

```@docs
FeasibleCommunities
coordinates
halfspaces
interaction_matrix
growth_rates
community_matrix
fluxes
dimension
isbounded
structurally_singular
diagnose
CommunityPrior
logdensity
```
