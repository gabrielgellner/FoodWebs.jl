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

```@docs
web_modularity
compartments
modularity_zscore
species_roles
compartment_homogeneity
block_structure
block_random_matrix
modular_random_matrix
modular_connectances
block_spectrum
```
