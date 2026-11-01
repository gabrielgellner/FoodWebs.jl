# Community matrices and stability

[`random_predprey`](@ref) builds a random Jacobian with a web's structure: for each link where
predator `j` eats prey `i`, `J[i, j] = -|a|` and `J[j, i] = f |b|`, with `a`, `b` drawn from a
distribution and `f` the conversion efficiency, and `-1` on the diagonal.

The [`spectral_abscissa`](@ref) decides local stability and the asymptotic return rate; the
[`numerical_abscissa`](@ref) decides reactivity, whether some perturbations grow at first.

## The elliptic law

For a large random matrix whose off-diagonal entries have variance `σ²` and pairs `(M[i, j],
M[j, i])` with correlation `ρ`, the eigenvalues fill an ellipse centered on the mean diagonal with
semi-axes `√(Sσ²)(1 ± ρ)` (Sommers et al. 1988; Allesina & Tang 2012). [`elliptic_law`](@ref)
estimates this from a matrix's moments, including the zero pairs of a sparse matrix, so the
connectance enters through `σ²`. A nonzero mean adds an outlier at `θ + ρSσ²/θ` from the center,
`θ = (S - 1)μ` (O'Rourke & Renfrew 2014).

```@example matrices
using FoodWebs, Distributions, Random, CairoMakie

rng = Xoshiro(1)
J = random_predprey(rng, may_network(rng, 300, 0.1), Normal(), 0.3)
fig, ax = eigenplot(J; title = "predator–prey, f = 0.3")
fig
```

Predator–prey pairs have negative correlation, so the ellipse is narrow along the real axis, and
with `f < 1` the negative mean adds an outlier on the stable side.

The elliptic law assumes the pairs are exchangeable. Structured food webs break this: in a
cascade web every predator's effect lies above the diagonal and every prey's effect below, and
the spectrum gains large imaginary outliers.

```@example matrices
web = cascade_network(rng, 300, 0.1)
fig = Figure(size = (900, 450))
eigenplot!(Axis(fig[1, 1]; title = "cascade web"), random_predprey(rng, web, Normal(), 1.0))
eigenplot!(Axis(fig[1, 2]; title = "predators shuffled"), random_predprey(rng, web, Normal(), 1.0; shuffle_prey = true))
fig
```

```@docs
random_predprey
spectral_abscissa
numerical_abscissa
λ_stability
ν_stability
elliptic_law
FoodWebs.EllipticLaw
offdiag_pairs
offdiag_pairs_nz
shuffle_pairs
```
