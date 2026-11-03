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

## References

- Gauzens, B. et al. fluxweb: An R package to easily estimate energy fluxes in food webs.
  Methods Ecol. Evol. 10, 270–279 (2019).
- Lang, B., Ehnes, R. B., Brose, U. & Rall, B. C. Temperature and consumer type dependencies of
  energy flows in natural communities. Oikos 126, 1717–1725 (2017).
