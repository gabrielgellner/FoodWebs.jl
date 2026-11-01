# Changelog

All notable changes to this project will be documented in this file.

## [0.7.0] - 2026-11-01

### Features

- **feasible**: Add the `per_capita_bounds` keyword to `FeasibleCommunities`, which bounds each loss, gain and self-regulation term by its per capita rate, `lo / x[j] ≤ θ ≤ hi / x[j]`, so that rescaling the biomasses `x` leaves the spectra unchanged; magnitude `bounds` cap the per capita rates at `hi x[j]` and so depend on the scale and shape of `x`. The docs gain a section on how the biomasses enter the model, and the docs examples and `gallery/equilibria.qmd` now use per capita bounds: under an efficiency cap of 0.3 the stable fraction falls from 0.82 to about 0.72 between Eltonian (ρ = 0.2) and inverted (ρ = 5) pyramids, against 0.89 to 0.71 in v0.6.0, which overstated the effect of the pyramid (#67)

## [0.6.0] - 2026-11-01

### Bug Fixes

- **randcmat**: Give species that eat each other both links' effects in `random_predprey`: each entry of a two-way pair is `f |b'| - |a|`, of either sign. Before, the lower index was always treated as the predator, which reversed the links of niche-family webs where a prey has a higher niche value than its predator, so results from those webs change (#56)

### Documentation

- **feasible**: Add a Three ensembles section to the feasible communities page, comparing `random_predprey`, feasible and energetically constrained communities, and the `gallery/equilibria.qmd` note on the stable fraction against the shape of the biomass pyramid, for generated pyramids and the observed Tuesday Lake biomasses (#54)

### Features

- **topology**: Add `remove_cannibalism` and `orient_two_way_links`, which keeps one link of each two-way pair by trophic level (`:up`), index, at random, or drops both; both return the web and the removed links and keep the input type (#57)
- **randcmat**: Add `expected_community_matrix`, the mean of `random_predprey` given the web, and the `expected` keyword of `eigenplot` to mark its eigenvalues, which predict the largest imaginary part of cascade and niche spectra to within about 1% (#59)
- **biomass**: Keep node traits from `nodes.csv` in the new `FoodWeb.traits` field; add `biomass` (body mass times abundance), `steady_state_biomass` from energy flow through the web with fixed, per-species or allometric turnover, and the pyramid measures `biomass_ratios`, `pyramid_slope`, `is_pyramid` and `is_inverted_pyramid`, with a new docs page (#51)
- **feasible**: Add `FeasibleCommunities`, the polytope of Lotka–Volterra interaction strengths, self-regulation and growth rates with a given equilibrium, for every web including cannibalism and mutual feeding; `halfspaces`, `coordinates`, `interaction_matrix`, `growth_rates`, `community_matrix` and `fluxes` convert between the polytope and models; `dimension`, `isbounded`, `structurally_singular` and `diagnose` describe it without a solver; `CommunityPrior` and `logdensity` put densities on it (#52)
- **feasible**: Add the PolytopeSamplers.jl extension, with `Polytope(fc)`, `isempty(fc)` by linear program, `sample_communities` for uniform and prior targets, and `feasible_growth_rates`, the redraw loop of the published studies; add `draw_growth_rates` (#53)

## [0.5.0] - 2026-11-01

### Documentation

- Add a Documenter site deployed to GitHub Pages at https://gabrielgellner.github.io/FoodWebs.jl, with pages on the models, empirical webs, trophic levels and structure, motifs, community matrices and plotting, live examples and figures, and every exported docstring; `just docs` builds it locally (#49)

### Features

- **elliptic law**: Predict the outlier that a nonzero mean adds to the spectrum, `EllipticLaw.outlier` at `θ + ρSσ²/θ` from the center (O'Rourke & Renfrew 2014); `EllipticLaw` gains this field, so its positional constructor takes one more argument; add `eigenplot` and `eigenplot!` to draw eigenvalues with the elliptic-law ellipse, outlier and stability boundary; test the random-matrix reference cases of Allesina's notes and the imaginary outliers of cascade-structured webs; add the `gallery/spectra.qmd` note (#20, #41)
- **motifs**: Add `motif_census` for the 13 three-species and 199 four-species motifs by mfinder ID, `motif_name`, `motif_roles` for the 30 species positions, `rewire_links` for the null model keeping each species' one-way and two-way links, and `motif_zscores`; the census matches igraph and NetworkX on every bundled web (#42, #43)
- **motifs**: Return the normalized profile from `motif_zscores`; the eight bundled webs studied by Stouffer et al. (2007) fall in their published groups by the signs of the omnivory and competition motifs (#44, #45)

## [0.4.0] - 2026-10-31

### Breaking Changes

- License the package under GPL-3.0-or-later, so that code adapted from GPL-licensed ecology software can be brought in; versions up to 0.3.0 remain available under MIT (#38)

### Bug Fixes

- **topology**: Accept a web when its number of links is within `tol` of `C S^2` or is the whole number of links nearest it, so that small webs reach every connectance; throw an `ArgumentError` up front when the target allows too few links to connect `S` species (#26)

### Features

- **utils**: Add MacKay–Johnson–Sansom trophic levels, `trophic_levels(web; method = :mackay)`, defined for every web including those without basal species, and `trophic_incoherence` (`F₀`); both accept weighted webs (#14)
- **topology**: Add the nested-hierarchy model, `nested_hierarchy_network`, following Cattin et al. (2004) as described in Stouffer et al. (2005) (#16)
- **topology**: Add the minimum potential niche model, `minimum_potential_niche_network`, with niche intervals drawn at a potential connectance that keeps the expected connectance at `C` after forbidding links (#17)
- **topology**: Add the probabilistic niche model: `probabilistic_niche_probabilities` gives the link probabilities, with a kernel exponent `p` (the Gaussian is `p = 2`), and `probabilistic_niche_network` draws webs from given parameters or for a target `S` and `C` (#18)
- **data**: Bundle 26 empirical food webs, the 19 webs of the R package cheddar and the 7 webs of Williams & Martinez (2000), each with its citation, source and license in a catalog and a Frictionless data package; add the `FoodWeb` type, `foodweb`, `foodweb_ids` and `read_foodweb`, and scripts that download the Web of Life and Mangal food webs (#19, #35)
- **properties**: Add the structural properties of Williams & Martinez (2000), `structural_properties` and its parts (`intermediate_species`, `cannibals`, `species_in_loops`, `omnivores`, `generality`, `vulnerability`, `trophic_similarity`, `mean_max_similarity`), `characteristic_path_length`, and shortest and short-weighted trophic levels; tested against their Table 2 and the path lengths of Dunne et al. (2002) (#19, #37)
- **plotting**: End arrows at node edges in screen space, draw all arrows in one call, scale node size, edge opacity and arrowhead size with the size of the web, and colour labels by node luminance; `trophic_layout` measures bands from their lowest species, so layouts can differ from 0.3.0; add a Quarto gallery of the bundled webs, rendered by `just gallery` (#36, #39)

## [0.3.0] - 2026-10-31

### Features

- **plotting**: Lay out food webs by trophic level with `trophic_layout`, ordering species within each level to reduce edge crossings, and draw them with `foodwebplot` and `foodwebplot!` in a Makie package extension (#23, #27)

## [0.2.0] - 2026-10-30

### Breaking Changes

- Require Julia 1.13 (#15)
- **topology**: Measure connectance as `L / S^2` everywhere; generators take `C` in this convention, with niche `β = 1/(2C) - 1` and cascade link probability `2CS/(S - 1)` (#13)
- **topology**: Store food webs as `Matrix{Bool}`; `PredationMatrix` is an `AbstractMatrix{Bool}` that checks its input is square and binary (#13, #21)
- **topology**: Link each pair in `may_network` with probability `2CS/(S - 1)` in a random direction, with no self links or mutual predation (#13)
- **randcmat**: Return a Jacobian from `random_predprey`, with `J[i, j]` the effect of species `j` on species `i`, so `J[prey, pred] < 0 < J[pred, prey]`; eigenvalues and reactivity are unchanged, eigenvectors are transposed (#13)
- **elliptic law**: Rename `ellipsis_sample` to `elliptic_law` and `EllipsisLawSolution` to `EllipticLaw`, declared `public`, and drop its unused `μ2` field (#21)
- **utils**: Throw an `ArgumentError` from `trophic_levels` for species with no feeding path up from a basal species, whose prey-averaged levels are undefined (#13)

### Features

- **topology**: Accept generated webs Williams & Martinez style: resample until the web is connected and its connectance is within `tol` (default 3%) of the target, up to `maxiter` attempts (#13)
- **topology**: Take an optional leading `rng::AbstractRNG` in every generator, `random_predprey` and `shuffle_pairs` (#13)
- **topology**: Reject webs with species cut off from basal species with `basal_paths = true` (#13)
- **topology**: Check `S`, `C` and `diet` arguments in the generators (#13)
- **utils**: Add `disconnected_from_basal`, the species with no feeding path up from a basal species (#13)
- **stability**: Add `spectral_abscissa` and `numerical_abscissa`; `λ_stability` and `ν_stability` remain as aliases (#21)
- Export `PredationMatrix` and `shuffle_pairs` (#13)

### Bug Fixes

- **topology**: Draw one feeding center per species in `niche_network`; a center was drawn per candidate prey, so diets were not intervals (#13)
- **topology**: Draw one `Beta(1, β)` probability per predator in `generalized_cascade_network`; a draw per pair reduced it to the plain cascade model (#13)
- **topology**: Make `generalized_niche_network` run (`round(Int, vector)` threw), use the `web[prey, pred]` orientation, sample cascade links without replacement, exclude self links and retry until connected (#13)
- **randcmat**: Keep links below the diagonal in `random_predprey`; links where the prey ranks above its predator, common in niche webs, were dropped (#13)
- **utils**: Solve for trophic levels directly; `pinv` returned meaningless values for species without a path from a basal species (#13)

### Performance

- **stability**: Compute `numerical_abscissa` with `eigmax(Symmetric(...))` (#21)

### Miscellaneous

- Replace LightGraphs with Graphs.jl and drop unused dependencies (DifferentialEquations, Plots, PyPlot, NLsolve, Parameters, ForwardDiff, RecursiveArrayTools) (#13)
- Test on GitHub Actions (Linux and macOS), with Aqua and ExplicitImports checks, and test model properties with per-testset seeded RNGs (#13, #15)
- Format with Runic and lint with JETLS in CI; add a justfile (#15, #22)
- Add a git-cliff changelog and release recipes (#24)
