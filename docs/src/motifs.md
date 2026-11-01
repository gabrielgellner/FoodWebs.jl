# Motifs

A motif is a connected subgraph of a few species, counted up to relabeling. [`motif_census`](@ref)
counts the 13 three-species motifs or the 199 four-species classes, labeled by the mfinder IDs of
Milo et al. (2002). The three-species motifs carry the names of Stouffer et al. (2007): S1 food
chain, S2 omnivory, S3 loop, S4 exploitative competition, S5 apparent competition, and D1–D8
with two-way links. The four-species diamond, two consumers sharing a prey and a predator, is
`M904`.

```@example motifs
using FoodWebs

web = foodweb("williams_martinez_2000/bridge_brook_lake")
Dict(motif_name(id) => n for (id, n) in zip(motif_ids(3), motif_census(web)) if n > 0)
```

The census matches igraph and NetworkX on every bundled web.

## Over- and under-representation

[`motif_zscores`](@ref) compares the census with webs randomized by [`rewire_links`](@ref), which
keeps each species' numbers of prey, predators, one-way and two-way links, and its cannibalism.
All eight bundled webs studied by Stouffer et al. (2007) fall in their published groups by the
signs of the S2, S4 and S5 z-scores.

```@example motifs
using Random

p = motif_zscores(Xoshiro(1), web; samples = 100)
Dict(motif_name(id) => round(z, digits = 1) for (id, z) in zip(p.ids, p.z) if isfinite(z))
```

## Species roles

[`motif_roles`](@ref) counts each species' occurrences in the 30 positions of the
three-species motifs (Stouffer et al. 2012), such as the top, middle and bottom of a food chain.

```@docs
motif_census
motif_ids
motif_name
motif_roles
rewire_links
motif_zscores
```
