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

## Body mass

Body mass `M` and abundance `N` test allometric theory along a web's links. Use
[`mass_ratios`](@ref) for the predator-to-prey mass ratio on every link and
[`mass_abundance_slope`](@ref) for the regression of `log10(N)` on `log10(M)` across species.
[`mass_abundance_links`](@ref) treats each link as an arrow in the plane of log mass and log
abundance (Cohen, Schittler, Raffaelli & Reuman 2009), and [`mass_level_slope`](@ref) gives
the trophic levels gained per tenfold increase in mass.

In Tuesday Lake in 1984 abundance falls as `M^-0.83` (Cohen, Jonsson & Carpenter 2003). On all
but 2 of its 264 links the predator outweighs its prey, and the median ratio is about 1,900.

```@example biomass
using Statistics: median

fw = foodweb("cheddar/tl84")
ratios = [r.ratio for r in mass_ratios(fw)]
(
    slope = round(mass_abundance_slope(fw).slope, digits = 3),
    heavier_prey = count(<(1), ratios),
    median_ratio = round(median(ratios)),
    mean_link_length = round(sum(l.length for l in mass_abundance_links(fw)) / length(ratios), digits = 2),
    levels_per_decade = round(mass_level_slope(fw).slope, digits = 3),
)
```

```@docs
mass_ratios
mass_abundance_slope
mass_abundance_links
mass_level_slope
```
