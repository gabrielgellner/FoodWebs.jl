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
which the PolytopeSamplers.jl extension provides.

## Priors

A uniform distribution on the polytope is one prior. [`CommunityPrior`](@ref) puts a
distribution on each kind of coordinate, and [`logdensity`](@ref) gives the density restricted
to the polytope, which a sampler for non-uniform targets uses. A density on the growth rates
times a uniform on the interaction strengths conditions a joint prior on the equilibrium holding,
weighting each rate vector by the volume of interaction strengths it allows; drawing the rates
first and sampling the interactions on their slice does not, so the two priors differ. A uniform
distribution stays uniform under linear changes of coordinates, so uniform in the magnitudes is
uniform in the fluxes or in the entries of the community matrix.

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
