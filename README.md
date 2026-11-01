# FoodWebs.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://gabrielgellner.github.io/FoodWebs.jl/dev/)

Models and utilities for ecological food webs: network topologies, random community
matrices, stability metrics and elliptic-law predictions.

```julia
using Pkg
Pkg.add(url = "https://github.com/gabrielgellner/FoodWebs.jl")
```

## Conventions

- **Food webs** are square binary matrices where `web[i, j] == 1` means species `j` eats
  species `i` (prey are rows, predators are columns).
- **Connectance** is `C = L / S^2` (Williams & Martinez 2000).
- **Community matrices** use the Jacobian convention: `J[i, j]` is the effect of species `j`
  on the growth of species `i`, so for predator `j` and prey `i`, `J[i, j] < 0 < J[j, i]`.
- **Random generators** take an optional leading `rng::AbstractRNG` and resample until the web
  is connected with connectance within a relative `tol` (default 3%) of the target. Pass
  `basal_paths = true` to also reject webs with species that have no feeding path up from a
  basal species (about 1–2% of niche webs), which is useful for dynamical models.

```julia
using FoodWebs, Distributions, Random
web = niche_network(Xoshiro(1), 30, 0.15)
J = random_predprey(Xoshiro(2), web, Normal(), 0.5)
spectral_abscissa(J), numerical_abscissa(J)
```

## Plotting

`trophic_layout(web)` gives node positions with each species at the height of its trophic level,
ordered across each level to reduce edge crossings. With a Makie backend loaded, `foodwebplot`
draws the web from those positions, with arrows from prey to predator.

```julia
using FoodWebs, CairoMakie, Random
fig, ax = foodwebplot(niche_network(Xoshiro(1), 20, 0.12; basal_paths = true))
```

## Development

Tasks are run with [just](https://github.com/casey/just):

```sh
just setup      # install Runic and JETLS as Julia apps, instantiate environments, enable hooks
just fmt        # format with Runic
just ci         # formatting check, JETLS lint, and tests (run by the pre-push hook)
just gallery    # render the Quarto gallery of the bundled webs to gallery/webs.html
just docs       # build the Documenter site to docs/build
```

Formatting uses [Runic](https://github.com/fredrikekre/Runic.jl). To have `git blame` skip the
bulk formatting commit, run `git config blame.ignoreRevsFile .git-blame-ignore-revs`.

## License

Copyright (C) 2015–2026 Gabriel Gellner.

FoodWebs.jl is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version. It is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See [LICENSE](LICENSE) for the full text.

Versions up to 0.3.0 were released under the MIT License. The bundled data in `data/` keep their
own terms; see [data/README.md](data/README.md).
