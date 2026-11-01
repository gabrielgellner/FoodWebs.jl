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

```@docs
web_modularity
compartments
```
