# Equilibrium biomass

Equilibrium-based community models need the biomass of each species. FoodWebs.jl reads it from
data where a web has it, or computes it from steady-state energy flow through the web.

## Observed biomass

Webs read from data keep the extra columns of their `nodes.csv` as `fw.traits`. For the cheddar
webs with body mass `M` and numerical abundance `N`, [`biomass`](@ref) is their product; species
without measurements, such as detritus, are `missing`.

```@example biomass
using FoodWebs

fw = foodweb("cheddar/tl84")
x_observed = biomass(fw)
extrema(x_observed)
```

## Steady-state energy flow

[`steady_state_biomass`](@ref) passes production up the web: basal species produce at a fixed
rate, each consumer produces a transfer efficiency times the production it takes from its prey,
and biomass is production divided by turnover. It needs no trophic levels, so omnivory and
feeding loops are handled by the web itself. With turnover scaling as `M^(-1/4)`, the biomass
ratio per level is `TE · R^(1/4)` for a predator–prey mass ratio `R`, so pyramids invert when
`TE · R^(1/4) > 1` (McCauley et al. 2018).

```@example biomass
x_flow = steady_state_biomass(fw; efficiency = 0.1, turnover = :allometric)
pyramid_slope(fw, x_flow), pyramid_slope(fw, x_observed)
```

The slope of log biomass against trophic level is negative for the energy-flow prediction and
positive for the observed biomass: Tuesday Lake is top-heavy, with more biomass at higher levels
than steady-state transfer at 10% efficiency and allometric turnover gives.

## Pyramid measures

[`biomass_ratios`](@ref) gives the consumer-to-resource ratio on every link,
[`pyramid_slope`](@ref) the mean change per trophic level, and [`is_pyramid`](@ref) and
[`is_inverted_pyramid`](@ref) whether biomass never rises (or never falls) up a link.

```@docs
biomass
steady_state_biomass
biomass_ratios
pyramid_slope
is_pyramid
is_inverted_pyramid
```
