# Trophic levels and structure

## Trophic levels

[`trophic_levels`](@ref) offers four definitions:

- `:prey_averaged` (default): basal species at 1, every other species 1 above the mean level of
  its prey (Levine 1980).
- `:shortest`: 1 plus the number of links in the shortest food chain down to a basal species.
- `:short_weighted`: the mean of the shortest and prey-averaged levels (Williams & Martinez
  2004).
- `:mackay`: the least-squares fit of "each link spans one level" over all links (MacKay,
  Johnson & Sansom 2020), defined for every web, including webs with no basal species.

The first three need a feeding path from a basal species to every species
([`disconnected_from_basal`](@ref)). The MacKay levels place a basal species one level below the
mean level of its predators, so a resource eaten directly by top predators sits above level 1.
[`trophic_incoherence`](@ref) is the mean squared amount by which links miss spanning one
level: 0 for a perfectly layered web, 1 for a web with no net direction.

```@example levels
using FoodWebs

# B1 and B2 eat nothing, H eats B1, P eats H and B2
web = Bool[0 0 1 0; 0 0 0 1; 0 0 0 1; 0 0 0 0]
trophic_levels(web), trophic_levels(web; method = :mackay), trophic_incoherence(web)
```

The prey-averaged levels call P an omnivore at 2.5. The MacKay levels lift B2 to level 2 so that
every link spans one level, and the web is perfectly coherent.

## Cannibalism and two-way links

Many community-matrix studies simplify a web before using it: [`remove_cannibalism`](@ref)
drops self-links and [`orient_two_way_links`](@ref) keeps one link of each pair of species that
eat each other, by a rule the caller chooses. No rule is established, so comparing results
under several rules, and against keeping both links, shows how much the choice matters.

## Structural properties

[`structural_properties`](@ref) returns the properties of Williams & Martinez (2000) and the
characteristic path length of Dunne, Williams & Martinez (2002). On the bundled Williams &
Martinez webs they reproduce 45 of the 53 comparable values of their Table 2 to its printed
precision, and the path lengths of Dunne et al. Table 1 for all seven webs; the tests list the
exceptions and their causes.

```@example levels
structural_properties(foodweb("williams_martinez_2000/bridge_brook_lake"))
```

```@docs
connectance
basal_species
top_species
intermediate_species
cannibals
disconnected_from_basal
remove_cannibalism
orient_two_way_links
trophic_levels
trophic_incoherence
species_in_loops
omnivores
generality
vulnerability
trophic_similarity
mean_max_similarity
characteristic_path_length
structural_properties
```
