# Food web models

Each generator returns a [`PredationMatrix`](@ref) for `S` species at a target connectance `C`.
Webs are resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` (default 3%) of
`C S^2`, widened for small webs to the whole number of links nearest `C S^2`. With
`basal_paths = true` they are also resampled until every species has a feeding path up from a
basal species, which dynamical models usually need. A target that cannot be met throws an
`ArgumentError` instead of resampling forever.

| Model | Generator | Prey selection |
|---|---|---|
| Food chain | [`food_chain_network`](@ref) | each species eats the one below |
| May | [`may_network`](@ref) | random pairs, random direction |
| Cascade | [`cascade_network`](@ref) | lower-ranked species, fixed probability |
| Generalized cascade | [`generalized_cascade_network`](@ref) | lower-ranked species, a Beta-distributed probability per predator |
| Niche | [`niche_network`](@ref) | an interval of niche values |
| Generalized niche | [`generalized_niche_network`](@ref) | part interval, part cascade |
| Nested hierarchy | [`nested_hierarchy_network`](@ref) | prey of the consumers sharing the first prey |
| Minimum potential niche | [`minimum_potential_niche_network`](@ref) | an interval with forbidden links inside |
| Probabilistic niche | [`probabilistic_niche_network`](@ref) | a kernel around the feeding center |

Mean connectance with acceptance turned off (`tol = Inf`) matches the target for every model,
which checks each model's parameterization separately from the acceptance step.

```@example models
using FoodWebs, Random

rng = Xoshiro(1)
webs = Dict(
    "cascade" => cascade_network(rng, 40, 0.1),
    "niche" => niche_network(rng, 40, 0.1),
    "nested hierarchy" => nested_hierarchy_network(rng, 40, 0.1),
)
Dict(name => round(connectance(w), digits = 3) for (name, w) in webs)
```

```@docs
PredationMatrix
food_chain_network
may_network
cascade_network
generalized_cascade_network
niche_network
generalized_niche_network
nested_hierarchy_network
minimum_potential_niche_network
probabilistic_niche_network
probabilistic_niche_probabilities
```
