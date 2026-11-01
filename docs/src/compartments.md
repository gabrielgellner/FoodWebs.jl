# Compartments

Compartments are groups of species that interact more among themselves than with others. In real
food webs they are mostly energy channels and habitat splits, and their modularity is close to
what the generalized niche model predicts (Guimerà et al. 2010). They limit how far extinctions
and perturbations spread (Stouffer & Bascompte 2011), while the effect of block structure on the
eigenvalues of a community matrix depends on the sign and strength of the links between the
compartments (Grilli, Rogers & Allesina 2016).

[`web_modularity`](@ref) scores a partition of a web, or of a weight matrix such as energy fluxes,
by how many more links fall inside the groups than expected from the species' degrees. Its
directed form (Leicht & Newman 2008) is the default; the undirected form of Newman & Girvan
(2004) is the one most published food web values use.

```@example compartments
using FoodWebs

# two food chains, 1 → 2 → 3 and 4 → 5 → 6, joined by a top predator 7 that eats 3 and 6
web = falses(7, 7)
for (prey, predator) in [(1, 2), (2, 3), (4, 5), (5, 6), (3, 7), (6, 7)]
    web[prey, predator] = true
end
channels = [1, 1, 1, 2, 2, 2, 3]
layers = [1, 2, 3, 1, 2, 3, 3]
web_modularity(web, channels), web_modularity(web, layers)
```

The two channels score 1/3. Grouping by trophic layer scores no better than chance here, as
every link runs between the groups. Grouping the Tuesday Lake web by the category of each species
(producers, invertebrates and vertebrate ectotherms), which follows its trophic layers, scores
below zero:

```@example compartments
fw = foodweb("cheddar/tl84")
web_modularity(fw, fw.traits["category"]), web_modularity(fw, fw.traits["category"]; directed = false)
```

The resolution `γ` (Reichardt & Bornholdt 2006) sets the size of the groups that modularity
favours; values above 1 favour smaller groups. The modularity of a web's best partition is never
negative and is high even for random webs (Guimerà, Sales-Pardo & Amaral 2004), so its
significance has to be judged against null webs.

## Finding compartments

[`compartments`](@ref) searches for a partition of high modularity, by the Louvain method or by
simulated annealing, and combines repeated searches by consensus. In the two-chain web the best
partition puts the top predator with the top of one chain, which scores slightly more than the
two channels:

```@example compartments
using Random: Xoshiro

compartments(Xoshiro(1), web)
```

Many partitions of a real web have nearly the highest modularity (Good, de Montjoye & Clauset
2010), so a single search is one of many answers. The modularity of each search and the fraction
of searches that put each pair of species together show how well defined the compartments are:

```@example compartments
r = compartments(Xoshiro(2026), fw; runs = 20)
maximum(r.groups), r.modularity, extrema(r.run_modularity)
```

```@example compartments
# the fraction of species pairs whose grouping differs between searches
count(x -> 0 < x < 1, r.coassignment) / length(r.coassignment)
```

## Three views of a compartment

Guimerà et al. (2010) define a compartment in three ways: species densely linked to each other
(`view = :density`, the default), species that share many prey (`:predator`), and species that
share many predators (`:prey`), the last two by the bipartite modularity of Guimerà, Sales-Pardo &
Amaral (2007). The density view finds energy channels, which follow the sink webs; the predator
and prey views find guilds of species at similar trophic levels, which cut across the sink webs,
as Guimerà et al. found. Species with no prey (or no predators) share them with no one and stay
alone in the predator (prey) view:

```@example compartments
binary = first(remove_cannibalism(fw)) .!= 0
map((:density, :predator, :prey)) do view
    r = compartments(Xoshiro(2026), binary; view, runs = 20)
    h = compartment_homogeneity(binary, r.best)
    (; view, Q = round(maximum(r.run_modularity); digits = 3), groups = maximum(r.best),
        trophic = round(h.trophic; digits = 1), sink = round(h.sink; digits = 1))
end
```

## Significance against null webs

Random webs reach high modularity from fluctuations alone (Guimerà, Sales-Pardo & Amaral 2004),
so [`modularity_zscore`](@ref) runs the same search on every null web and compares the best
modularity found. Against webs rewired to keep each species' links, Tuesday Lake is strongly
compartmented; against the niche model it is less modular, and against the generalized niche
model, whose diets need not be contiguous, it is typical, as Guimerà et al. (2010) found for
other food webs:

```@example compartments
rewired = modularity_zscore(Xoshiro(1), fw; samples = 50)
niche = modularity_zscore(Xoshiro(1), fw; samples = 50, null = niche_network)
generalized = modularity_zscore(Xoshiro(1), fw; samples = 50,
    null = (rng, S, C) -> generalized_niche_network(rng, S, C, 0.6))
[(null, round(r.mean; digits = 3), round(r.z; digits = 1)) for (null, r) in
    (("rewired", rewired), ("niche", niche), ("generalized niche", generalized))]
```

```@example compartments
rewired.observed
```

## Species roles and homogeneity

[`species_roles`](@ref) places each species by its within-module degree `z` and its
participation coefficient `P`, the spread of its links over the groups (Guimerà & Amaral 2005).
Connectors, with `P > 0.62`, link compartments, and in Tuesday Lake they sit higher in the web
than the peripheral species:

```@example compartments
groups = compartments(Xoshiro(2026), fw; runs = 20).best
roles = species_roles(fw, groups)
levels = trophic_levels(fw)
[(role, count(==(role), roles.role), round(sum(levels[roles.role .== role]) / count(==(role), roles.role); digits = 2))
    for role in unique(roles.role)]
```

[`compartment_homogeneity`](@ref) measures whether compartments follow trophic levels and sink
webs (Guimerà et al. 2010). Tuesday Lake's compartments follow its sink webs, the energy
channels from the basal species up; its species categories follow trophic levels instead and cut
across the sink webs:

```@example compartments
compartment_homogeneity(fw, groups), compartment_homogeneity(fw, fw.traits["category"])
```

## Links between compartments

The effect of compartments on the eigenvalues of a community matrix depends on the links between
them: the mean of the between-group entries sets the outlier eigenvalues, and their variance and
pair correlation the bulk (Grilli, Rogers & Allesina 2016). [`block_structure`](@ref) splits a
community or flux matrix into these within- and between-group statistics, with the kinds of pairs
in each block and, for each species, the share of its interaction magnitude with other groups.
For a random community matrix on Tuesday Lake, 113 of the 263 predator–prey pairs cross
compartments:

```@example compartments
using Distributions: Uniform

web = PredationMatrix(first(remove_cannibalism(fw)) .!= 0)
J = random_predprey(Xoshiro(1), web, Uniform(0, 1), 0.3)
blocks = block_structure(J, groups)
k = length(blocks.sizes)
within = sum(blocks.predator_prey[r, r] for r in 1:k)
between = sum(blocks.predator_prey[r, s] for r in 1:k for s in (r + 1):k)
within, between
```

The 33 coupling species span the web, from basal species eaten in several compartments to
consumers at trophic level 3.7; the species with no links outside their compartment are all
basal:

```@example compartments
couplers = blocks.coupling .> 0
count(couplers), extrema(levels[couplers]), extrema(levels[.!couplers])
```

## Groups beyond compartments

Modularity finds only compartments, groups linked mostly among themselves. The group model of
Allesina & Pascual (2009) is a stochastic block model: a species of group `r` is eaten by one of
group `s` with probability `p[r, s]`. Dense blocks on the diagonal of `p` are compartments; dense
blocks off it are guilds, species that eat the same prey or are eaten by the same predators.
[`fit_group_model`](@ref) searches for the partition of highest likelihood for each number of
groups and selects the number by the integrated classification likelihood, which, unlike the AIC
of the paper, picks one group for a random web. Tuesday Lake divides into trophic guilds: two
basal groups, two of herbivores and two of predators, with the dense blocks above the diagonal:

```@example compartments
model = fit_group_model(Xoshiro(2026), binary; k = 1:8)
model.k, [round(sum(levels[model.groups .== r]) / count(==(r), model.groups); digits = 1) for r in 1:model.k]
```

```@example compartments
round.(model.p; digits = 2)
```

## Energy channels

Energy channels are an ecological definition of compartments: the pathways of energy from each
basal source up the web (Moore & Hunt 1988). [`source_reliance`](@ref) gives the fraction of each
species' energy that comes from each source, following its diet (Rooney et al. 2008), and
[`energy_channels`](@ref) puts each species in the channel of the source it relies on most, then
merges channels while modularity rises (Zhao et al. 2018). In two chains joined by a top predator
that takes three times as much from the first, the predator gets three quarters of its energy
from the first source and joins its channel:

```@example compartments
F = zeros(7, 7)
F[1, 2] = F[2, 3] = F[4, 5] = F[5, 6] = 1
F[3, 7], F[6, 7] = 3, 1
channels = energy_channels(F)
channels.groups, channels.reliance[7, :]
```

In Tuesday Lake, the channels from its 31 basal species merge into eight, about half as modular
as the partition of highest modularity, which need not follow the flow of energy:

```@example compartments
channels = energy_channels(binary)
(channels = (maximum(channels.groups), round(channels.modularity; digits = 3)),
    highest = (maximum(groups), round(web_modularity(binary, groups); digits = 3)))
```

## Compartments of interaction strengths

Modularity needs non-negative weights, and a community matrix is signed: a predator–prey pair has
one negative and one positive entry. [`pair_weights`](@ref) gives each pair a symmetric weight,
`√|J_ij J_ji|` by default, which does not depend on the abundance units of each species, or
`|J_ij| + |J_ji|`. Two nulls then separate the sources of compartments: rewiring tests the
topology, and [`shuffle_weights`](@ref) (`null = :shuffle`) permutes the strengths among the
interacting pairs, testing whether they are arranged in compartments beyond what the topology
gives. In the random community matrix of Tuesday Lake the topology is strongly compartmented, but
the strengths, drawn independently, add nothing:

```@example compartments
W = pair_weights(J)
strengths = modularity_zscore(Xoshiro(2), W; null = :shuffle, samples = 30, runs = 3)
topology = modularity_zscore(Xoshiro(2), web; samples = 30, runs = 3)
(strengths = round(strengths.z; digits = 2), topology = round(topology.z; digits = 1))
```

## Signed modularity

`pair_weights(J; method = :signed)` keeps the sign of each pair, `sign(J_ij J_ji) √|J_ij J_ji|`:
negative for predator–prey pairs, positive for competition and mutualism. The signed modularity
of Gómez, Jensen & Arenas (2009) (`signed = true`) rewards positive pairs within groups and
negative pairs between them. In a purely trophic community matrix every pair is negative, so it
finds the anti-modular structure of Grilli, Rogers & Allesina (2016): in Tuesday Lake, two groups
that alternate up the trophic levels, basal species with the third level and herbivores with
the fourth, with nearly every predator–prey pair between them:

```@example compartments
signed_weights = pair_weights(J; method = :signed)
s = compartments(Xoshiro(2), signed_weights; signed = true, method = :annealing, runs = 3)
rounded = round.(levels)
[(group = g, levels = [count(i -> s.groups[i] == g && rounded[i] == l, eachindex(levels)) for l in 1:4]) for g in 1:maximum(s.groups)]
```

## Block-structured random matrices

[`block_random_matrix`](@ref) draws community matrices whose species fall into groups, with the
connectance, means, spreads and pair correlation of each pair of groups set separately; the
means may differ by direction, as when one group eats another. [`modular_random_matrix`](@ref)
is the two-group case of Grilli, Rogers & Allesina (2016), set by the overall connectance and the
modularity `Q` through [`modular_connectances`](@ref). [`block_spectrum`](@ref) predicts the
eigenvalues from the block moments: a bulk from the variances and pair correlations, and outliers
from the means.

The ratio of the predicted spectral abscissa with block structure to that without shows the
regimes of Grilli et al.: with a negative mean interaction, modularity stabilizes slightly and
anti-modularity (`Q < 0`) destabilizes strongly; with a positive mean, block structure is neutral
for two groups of equal size, since the leading outlier `(S/2)(μ_w + |μ_b|)` is then `SCμ` for
every `Q`, and destabilizing for unequal groups:

```@example compartments
function ratio(Q, μ, α; S = 1000, C = 0.2)
    n = round(Int, α * S)
    predicted(q) = let (Cw, Cb) = modular_connectances(C, q, α)
        block_spectrum([n, S - n]; connectance = [Cw Cb; Cb Cw], μ).abscissa
    end
    return predicted(Q) / predicted(0.0)
end
[(μ, α, [round(ratio(Q, μ, α); digits = 2) for Q in (-0.3, -0.1, 0.1, 0.3)])
    for (μ, α) in ((-1.0, 0.5), (-1.0, 0.25), (1.0, 0.5), (1.0, 0.25))]
```

The predictions also apply to a sampled matrix and its groups, from the moments of
[`block_structure`](@ref):

```@example compartments
using LinearAlgebra: eigvals

M = block_random_matrix(Xoshiro(3), [300, 500]; connectance = [0.15 0.1; 0.1 0.2],
    μ = [-0.3 -1.5; 1.0 -0.3], ρ = -0.4, diagonal = -1)
p = block_spectrum(M, [fill(1, 300); fill(2, 500)])
p.outliers, p.abscissa, maximum(real, eigvals(M))
```

## Drawing compartments

With Makie loaded, [`foodwebplot`](@ref) colours species by group, keeps each group together within
its trophic band, draws the links between groups in red, and marks species roles by shape;
`show_within = false` keeps only the coupling links:

```@example compartments
using CairoMakie

fig, ax = foodwebplot(binary; groups, roles = roles.role, show_within = false, labels = nothing,
    size = (750, 480))
fig
```

[`blockplot`](@ref) orders a matrix by group and draws the block boundaries, here the random
community matrix of Tuesday Lake from above and the co-assignment of 20 searches:

```@example compartments
r = compartments(Xoshiro(2026), fw; runs = 20)
fig = Figure(size = (900, 400))
blockplot!(Axis(fig[1, 1]; title = "community matrix", aspect = DataAspect(), yreversed = true), J, groups)
blockplot!(Axis(fig[1, 2]; title = "co-assignment", aspect = DataAspect(), yreversed = true), r.coassignment, r.groups)
fig
```

[`eigenplot`](@ref) with `groups` draws the bulk and outliers of [`block_spectrum`](@ref) in place
of the elliptic law, here for an anti-modular matrix whose outlier to the right sets its
stability:

```@example compartments
M = modular_random_matrix(Xoshiro(4), 600; α = 0.25, C = 0.2, Q = -0.3, μ = -1.0)
fig, ax = eigenplot(M; groups = [fill(1, 150); fill(2, 450)])
fig
```

```@docs
blockplot
blockplot!
web_modularity
compartments
modularity_zscore
species_roles
compartment_homogeneity
block_structure
group_model
fit_group_model
source_reliance
energy_channels
pair_weights
shuffle_weights
block_random_matrix
modular_random_matrix
modular_connectances
block_spectrum
```
