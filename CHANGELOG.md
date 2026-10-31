# Changelog

All notable changes to this project will be documented in this file.

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
