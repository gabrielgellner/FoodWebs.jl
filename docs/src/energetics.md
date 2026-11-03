# Energy fluxes

[`fluxes_from_metabolism`](@ref) estimates the biomass flux on every link at equilibrium from
biomasses and body masses, following fluxweb (Gauzens et al. 2019). Each consumer loses
`0.71 M^-0.25` per unit biomass to metabolism, and its assimilated intake pays for those losses
and for what its predators take. Its intake splits among its prey in proportion to their biomass,
so the balance is one linear equation per consumer. [`steady_state_biomass`](@ref) works the
other way, passing production up from the basal species.

[`assimilation_efficiencies`](@ref) gives the efficiency with which each species is assimilated
when eaten: 0.545 for plants, 0.906 for animals and 0.158 for detritus (Lang et al. 2017). For a
cheddar web it reads the types from the `"category"` trait.

```@example energetics
using FoodWebs

fw = foodweb("cheddar/tl84")
F = fluxes_from_metabolism(fw)
producer = fw.traits["category"] .== "producer"
(herbivory = sum(F[producer, :]), carnivory = sum(F[.!producer, :]))
```

Tuesday Lake records body mass in kg. The metabolic intercept depends on the units of mass and
time, so pass `intercept` to match other units.

## Fluxes in feasible communities

Fluxes and biomasses fix the losses of a [`FeasibleCommunities`](@ref) through `fixed_fluxes`,
and the same efficiencies fix the gains. Each species keeps its growth rate and self-regulation,
tied by its equilibrium row, so the polytope has one dimension per self-regulated species.

```@example energetics
x = biomass(fw)
fc = FeasibleCommunities(
    fw, x; fixed_fluxes = F, efficiency = assimilation_efficiencies(fw),
    rates = fill((-10.0, 10.0), size(fw, 1))
)
dimension(fc), isbounded(fc)
```

With self-regulation only from cannibalism, the community matrix is the Jacobian of fluxweb's
`create.jacob`, and its [`spectral_abscissa`](@ref) is fluxweb's `stability.value`. fluxweb gives
detritus a positive diagonal, which this model has no counterpart for.

```@docs
fluxes_from_metabolism
assimilation_efficiencies
```

## Flow networks

A [`FlowNetwork`](@ref) holds the flows between compartments at steady state, with each
compartment's input, export and respiration. Its indices follow enaR (Borrett & Lau 2014):
[`throughflow`](@ref), [`total_system_throughflow`](@ref), [`total_system_throughput`](@ref),
[`average_path_length`](@ref), [`finn_cycling_index`](@ref) (cycled throughflow over total
system throughflow, Finn 1980), [`ascendency`](@ref) with development capacity and overhead in
bits, and [`mixed_trophic_impact`](@ref) (Ulanowicz & Puccia 1990). The flow-weighted
[`trophic_levels`](@ref trophic_levels(::FlowNetwork)) and [`omnivory_index`](@ref) follow
NetIndices (Kones et al. 2009), without cannibalism.

The fluxes of Tuesday Lake give a budget at steady state: producers take in what their consumers
eat, and each consumer respires its metabolic losses and exports what it does not assimilate.

```@example energetics
e = assimilation_efficiencies(fw)
M = fw.traits["M"]
intake = vec(sum(F; dims = 1))
consumer = intake .> 0
net = FlowNetwork(
    F; inputs = .!consumer .* vec(sum(F; dims = 2)),
    exports = intake .- vec(sum(e .* F; dims = 1)),
    respiration = [consumer[i] ? 0.71 * M[i]^-0.25 * x[i] : 0.0 for i in eachindex(x)],
    names = fw.species,
)
a = ascendency(net)
(
    path_length = average_path_length(net),
    cycling = finn_cycling_index(net),
    relative_ascendency = a.ascendency / a.capacity,
    max_omnivory = maximum(omnivory_index(net)),
)
```

```@docs
FlowNetwork
throughflow
total_system_throughflow
total_system_throughput
average_path_length
finn_cycling_index
ascendency
trophic_levels(::FlowNetwork)
omnivory_index
mixed_trophic_impact
```

## References

- Borrett, S. R. & Lau, M. K. enaR: An R package for Ecosystem Network Analysis. Methods Ecol.
  Evol. 5, 1206–1213 (2014).
- Finn, J. T. Flow analysis of models of the Hubbard Brook ecosystem. Ecology 61, 562–571 (1980).
- Kones, J. K., Soetaert, K., van Oevelen, D. & Owino, J. O. Are network indices robust indicators
  of food web functioning? A Monte Carlo approach. Ecol. Model. 220, 370–382 (2009).
- Ulanowicz, R. E. & Puccia, C. J. Mixed trophic impacts in ecosystems. Coenoses 5, 7–16 (1990).
- Gauzens, B. et al. fluxweb: An R package to easily estimate energy fluxes in food webs.
  Methods Ecol. Evol. 10, 270–279 (2019).
- Lang, B., Ehnes, R. B., Brose, U. & Rall, B. C. Temperature and consumer type dependencies of
  energy flows in natural communities. Oikos 126, 1717–1725 (2017).
