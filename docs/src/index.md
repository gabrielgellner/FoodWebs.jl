# FoodWebs.jl

FoodWebs.jl generates structural food web models, reads empirical food webs, measures their
structure and builds random community matrices with their structure, for studying the stability
of large ecological communities.

- [Food web models](models.md): the cascade, niche and related generators, with
  Williams & Martinez style acceptance.
- [Empirical food webs](data.md): 26 bundled webs with citations, and scripts that download
  the Web of Life and Mangal webs.
- [Trophic levels and structure](properties.md): prey-averaged, short-weighted and MacKay
  trophic levels, trophic incoherence, and the structural properties of Williams & Martinez
  (2000), tested against their published tables.
- [Motifs](motifs.md): three- and four-species motif censuses, species roles and z-scores
  against a rewiring null model.
- [Community matrices and stability](matrices.md): random predator–prey Jacobians, spectral
  and numerical abscissas, and elliptic-law predictions with outliers.
- [Plotting](plotting.md): trophic-level layouts and Makie plots of webs and spectra.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/gabrielgellner/FoodWebs.jl")
```

FoodWebs.jl requires Julia 1.13. Plotting needs a Makie backend such as CairoMakie.

## Conventions

- **Food webs** are square binary matrices where `web[i, j] == 1` means species `j` eats
  species `i`: prey are rows, predators are columns. [`PredationMatrix`](@ref) stores a
  generated web and [`FoodWeb`](@ref) an empirical one with species names.
- **Connectance** is `C = L / S^2`, with `L` links among `S` species (Williams & Martinez 2000).
- **Community matrices** use the Jacobian convention: `J[i, j]` is the effect of species `j` on
  the growth of species `i`, so for predator `j` and prey `i`, `J[i, j] < 0 < J[j, i]`.
- **Random generators** take an optional leading `rng::AbstractRNG`, so results are
  reproducible from a seeded generator.

## A first example

```@example index
using FoodWebs, Distributions, Random

web = niche_network(Xoshiro(1), 30, 0.15)
connectance(web), length(basal_species(web)), length(top_species(web))
```

```@example index
J = random_predprey(Xoshiro(2), web, Normal(), 0.5)
spectral_abscissa(J), numerical_abscissa(J)
```

The package is licensed under the GNU General Public License, version 3 or later.
