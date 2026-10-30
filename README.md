# FoodWebs.jl

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
  is connected with connectance within a relative `tol` (default 3%) of the target.

```julia
using FoodWebs, Distributions, Random
web = niche_network(Xoshiro(1), 30, 0.15)
J = random_predprey(Xoshiro(2), web, Normal(), 0.5)
λ_stability(J), ν_stability(J)
```
