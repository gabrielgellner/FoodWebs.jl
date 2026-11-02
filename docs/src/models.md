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

## Allometric diet breadth model

The allometric diet breadth model of Petchey et al. (2008) predicts who eats whom in a given set
of species from their body masses, where the models above generate random webs. Each consumer
ranks its resources by profitability, energy over handling time, and eats the most profitable
ones up to the diet that maximizes its rate of energy intake. Encounter rates, energy contents
and handling times are power laws of the consumer and resource masses. With
[`RatioHandling`](@ref), a consumer cannot eat resources heavier than a fraction `b` of its own
mass, so the web has no loops and no cannibalism, and each diet is an interval of resource
masses centered near `b / 2` of the consumer's mass. [`PowerHandling`](@ref) lets consumers eat
heavier resources but gives every consumer the same ranking of resources.

[`fit_diet_breadth`](@ref) fits the model to an observed web with body masses, holding the
number of links at the observed number and maximizing the fraction of them in the right place.
Petchey et al. predict 57% of the Benguela links with ratio handling times and 42% with power
handling times; the fit below finds parameters that do better on the cheddar version of the
web.

```@example models
fw = foodweb("cheddar/benguela")
ratio = fit_diet_breadth(Xoshiro(1), fw; handling = :ratio)
power = fit_diet_breadth(Xoshiro(1), fw; handling = :power)
(ratio = round(ratio.fraction, digits = 2), power = round(power.fraction, digits = 2))
```

The fitted model comes close to the observed omnivory and spread of generality. Its web has no
loops, since every consumer eats lighter species, while 41% of the Benguela species are in
loops:

```@example models
using Statistics: mean

S, C = size(fw, 1), connectance(fw)
niche = [structural_properties(niche_network(rng, S, C; basal_paths = true)) for _ in 1:100]
props = (:T, :B, :GenSD, :VulSD, :Loop, :Omniv, :TL)
Dict(
    "observed" => round.(getindex.(Ref(structural_properties(fw)), props), digits = 2),
    "diet breadth (ratio)" => round.(getindex.(Ref(structural_properties(ratio.web)), props), digits = 2),
    "niche (mean)" => round.(Tuple(mean(getindex.(niche, p)) for p in props), digits = 2),
)
```

```@docs
diet_breadth_network
RatioHandling
PowerHandling
fit_diet_breadth
```

Petchey, O. L., Beckerman, A. P., Riede, J. O. & Warren, P. H. Size, foraging, and food web
structure. PNAS 105, 4191–4196 (2008).
