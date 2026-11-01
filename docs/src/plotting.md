# Plotting

[`trophic_layout`](@ref) places each species at the height of its trophic level and orders the
species across each band of levels to reduce edge crossings. With a Makie backend loaded,
[`foodwebplot`](@ref) draws a web from that layout, with arrows from prey to predator, and
[`eigenplot`](@ref) draws the eigenvalues of a matrix with its elliptic-law prediction.

```@example plotting
using FoodWebs, CairoMakie

fw = foodweb("williams_martinez_2000/chesapeake_bay")
fig, ax = foodwebplot(fw; size = (700, 500), title = fw.name)
fig
```

Any trophic levels can set the heights. The MacKay levels show where a loopy, omnivorous web
mixes its levels:

```@example plotting
fw = foodweb("williams_martinez_2000/coachella_valley")
fig = Figure(size = (1100, 450))
foodwebplot!(Axis(fig[1, 1]; title = "prey-averaged levels"), fw)
foodwebplot!(Axis(fig[1, 2]; title = "MacKay levels"), fw; levels = trophic_levels(fw; method = :mackay))
fig
```

Node size, edge opacity and arrowhead size scale with the size of the web, and labels are left
out above 60 species; each can be set by keyword. The repository's `gallery/` holds Quarto notes
drawing every bundled web and the random-matrix reference cases (`just gallery`).

```@docs
trophic_layout
FoodWebs.edge_crossings
foodwebplot
foodwebplot!
eigenplot
eigenplot!
```
