# Feeding loops

A feeding loop of a community matrix `J` is a directed cycle of nonzero off-diagonal entries,
`i₁ → i₂ → ⋯ → iₖ → i₁`, and its weight is the geometric mean of the absolute entries around it.
Neutel, Heesterbeek & de Ruiter (2002) found that the heaviest loop of three or more species sets
how much self-regulation a soil food web needs to be stable, and that long omnivory loops are
light because they contain weak top-down effects. [`feeding_loops`](@ref) lists the loops of `J`
up to a given length, heaviest first, and [`heaviest_loop`](@ref) finds the heaviest one with a
pruned search that reaches longer loops.

The de Ruiter, Neutel & Moore (1995) community matrix follows from the link fluxes `F` and the
biomasses `x`: predator `j` affects prey `i` by `-F[i, j] / x[j]` and prey `i` affects predator
`j` by `e[i] F[i, j] / x[i]`. For Tuesday Lake with fluxes from metabolism:

```@example loops
using FoodWebs, CairoMakie, Random, Statistics

fw = foodweb("cheddar/tl84")
F = fluxes_from_metabolism(fw)
x = biomass(fw)
e = assimilation_efficiencies(fw)
J = -F ./ x' .+ (e .* F)' ./ x'
loops = feeding_loops(J; min_length = 3, max_length = 4)
length(loops), heaviest_loop(J)
```

Each loop is a named tuple with its `species` in loop order, the `links` of `J` around it, its
`length`, `sign` and `weight`. Loops of three or more species appear in both directions, which
run through different entries; `both_directions = false` keeps the heavier one.
[`loopplot`](@ref) draws loops over the web, each as a chain of arrows in loop order:

```@example loops
top = feeding_loops(J; min_length = 3, max_length = 4, both_directions = false)[1:5]
fig, ax = loopplot(fw, top; size = (900, 500))
fig
```

## Self-regulation and loop weight

[`min_self_regulation`](@ref) is Neutel's stability measure: the smallest `s` for which the
off-diagonal part of `J` with diagonal `-s d` is stable. Here `d` is each species' metabolic rate
`0.71 M^-0.25`. Shuffling the top-down entries among the links, and the bottom-up entries among
the links, keeps the strengths and the web but breaks their pattern. The shuffled matrices have
heavier loops and need from 6 to over 100 times as much self-regulation:

```@example loops
d = 0.71 .* fw.traits["M"] .^ -0.25
links = [(i, j) for i in axes(F, 1), j in axes(F, 2) if i != j && F[i, j] > 0]
function shuffled(rng, J, links)
    K = zero(J)
    down, up = shuffle(rng, [J[i, j] for (i, j) in links]), shuffle(rng, [J[j, i] for (i, j) in links])
    for (k, (i, j)) in enumerate(links)
        K[i, j] += down[k]
        K[j, i] += up[k]
    end
    return K
end
rng = Xoshiro(1)
shuffles = [shuffled(rng, J, links) for _ in 1:30]
fig = Figure(size = (500, 400))
ax = Axis(fig[1, 1]; xlabel = "heaviest loop weight", ylabel = "minimum self-regulation s", xscale = log10, yscale = log10)
scatter!(ax, [heaviest_loop(K; max_length = 6).weight for K in shuffles], [min_self_regulation(K; scale = d) for K in shuffles]; label = "shuffled")
scatter!(ax, [heaviest_loop(J; max_length = 6).weight], [min_self_regulation(J; scale = d)]; label = "observed", markersize = 14)
axislegend(ax; position = :rb)
fig
```

## Species and entries of loops

[`loop_membership`](@ref) marks the species in each loop and [`loop_mask`](@ref) the entries of
`J` on any of a set of loops, so an analysis can change strengths inside the loops and outside
them separately. Here the entries of the 20 heaviest loops are shuffled among themselves, and
separately the other off-diagonal entries among themselves, each time keeping top-down and
bottom-up entries apart. The median self-regulation needed over 30 shuffles of each, against the
observed `s`:

```@example loops
S = size(J, 1)
inside = loop_mask(feeding_loops(J; min_length = 3, max_length = 4)[1:20], S)
outside = [i != j && J[i, j] != 0 && !inside[i, j] for i in 1:S, j in 1:S]
function scramble(rng, J, mask)
    K = copy(J)
    for sgn in (-1, 1)
        entries = findall(mask .& (sign.(J) .== sgn))
        K[entries] = shuffle(rng, J[entries])
    end
    return K
end
rng = Xoshiro(2)
(
    entries = (inside = count(inside), outside = count(outside)),
    observed = min_self_regulation(J; scale = d),
    inside = median(min_self_regulation(scramble(rng, J, inside); scale = d) for _ in 1:30),
    outside = median(min_self_regulation(scramble(rng, J, outside); scale = d) for _ in 1:30),
)
```

Shuffling the 38 entries of the heaviest loops raises the self-regulation needed fivefold, and
shuffling the 488 other entries twentyfold.

The heaviest loops sit low in the web, as Neutel et al. (2002) report. Across the loops of three
species, weight falls with the mean trophic level of the species in the loop:

```@example loops
levels = trophic_levels(fw)
triangles = feeding_loops(J; min_length = 3, max_length = 3, both_directions = false)
fig = Figure(size = (500, 400))
ax = Axis(fig[1, 1]; xlabel = "mean trophic level", ylabel = "loop weight", yscale = log10)
scatter!(ax, [mean(levels[l.species]) for l in triangles], [l.weight for l in triangles])
fig
```

```@docs
feeding_loops
heaviest_loop
loop_membership
loop_mask
min_self_regulation
loopplot
loopplot!
```

## References

- de Ruiter, P. C., Neutel, A.-M. & Moore, J. C. Energetics, patterns of interaction strengths,
  and stability in real ecosystems. Science 269, 1257–1260 (1995).
- Neutel, A.-M., Heesterbeek, J. A. P. & de Ruiter, P. C. Stability in real food webs: weak links
  in long loops. Science 296, 1120–1123 (2002).
- Neutel, A.-M. et al. Reconciling complexity with stability in naturally assembling food webs.
  Nature 449, 599–602 (2007).
