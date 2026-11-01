# # Block-structured random matrices
# Random community matrices whose species fall into groups, with the pair statistics of each pair
# of groups set separately, after Grilli, Rogers & Allesina (2016), and predictions of their
# eigenvalues: a bulk from the variances and pair correlations, and outliers from the means.

# a scalar or a k × k matrix of block parameters, as a k × k matrix
block_parameter(x::Real, k) = fill(float(x), k, k)
function block_parameter(x::AbstractMatrix, k)
    size(x) == (k, k) || throw(DimensionMismatch("block parameters must be $k×$k, got $(size(x))"))
    return Matrix{Float64}(x)
end

"""
    block_random_matrix([rng], sizes; connectance, μ, σ = 1, ρ = 0, diagonal = 0, cascade = 0)

A random community matrix of `sum(sizes)` species in groups of `sizes`, the first `sizes[1]`
species in group 1 and so on. Each pair of species `i` in group `r` and `j` in group `s` interacts
with probability `connectance[r, s]`, and then `(M[i, j], M[j, i])` is bivariate normal with means
`(μ[r, s], μ[s, r])`, standard deviations `(σ[r, s], σ[s, r])` and correlation `ρ[r, s]`;
otherwise both entries are 0. Each parameter is a number or a `k × k` matrix; `connectance` and
`ρ` describe pairs and must be symmetric, while `μ` and `σ` can differ by direction, as when the
species of one group eat those of another (`μ[r, s] < 0 < μ[s, r]`). The diagonal is `diagonal`.

Interactions come in pairs: either both entries of a pair are nonzero or neither is, as for the
feeding links of a community matrix, where a predator depresses its prey and the prey benefits
the predator. One-way effects, such as donor control, are not generated.

With `cascade = ξ ≠ 0`, each species gets a random size, and for a pair in which `i` is smaller
than `j`, the means are multiplied by `1 + ξ` for the effect of `j` on `i` and by `1 − ξ` for the
effect of `i` on `j`, so that with `μ < 0` and `ξ > 1` the larger species depresses the smaller and
benefits from it, as in a cascade food web (the food-web variant of Grilli, Rogers & Allesina
2016, with `ξ = 3`); for symmetric `μ` the mean over the pair stays `μ`.

Grilli, J., Rogers, T. & Allesina, S. Modularity and stability in ecological communities. Nat.
Commun. 7, 12031 (2016).
"""
function block_random_matrix(
        rng::AbstractRNG, sizes::AbstractVector{<:Integer}; connectance, μ, σ = 1, ρ = 0, diagonal::Real = 0,
        cascade::Real = 0
    )
    all(>(0), sizes) || throw(ArgumentError("group sizes must be positive"))
    k = length(sizes)
    Cm, μm, σm, ρm = (block_parameter(x, k) for x in (connectance, μ, σ, ρ))
    for (name, x) in (("connectance", Cm), ("ρ", ρm))
        x == x' || throw(ArgumentError("$name must be symmetric, as it describes pairs"))
    end
    all(c -> 0 <= c <= 1, Cm) || throw(ArgumentError("connectances must be in [0, 1]"))
    all(>=(0), σm) || throw(ArgumentError("σ must be non-negative"))
    all(r -> -1 <= r <= 1, ρm) || throw(ArgumentError("ρ must be in [-1, 1]"))
    S = sum(sizes)
    groups = reduce(vcat, (fill(r, n) for (r, n) in enumerate(sizes)))
    body = cascade == 0 ? zeros(S) : rand(rng, S)
    M = zeros(S, S)
    for i in 1:S, j in (i + 1):S
        r, s = groups[i], groups[j]
        rand(rng) < Cm[r, s] || continue
        z1, z2 = randn(rng), randn(rng)
        # correlated standard normal deviations for M[i, j] and M[j, i]
        x, y = z1, ρm[r, s] * z1 + sqrt(1 - ρm[r, s]^2) * z2
        # with a cascade, (1 + ξ) on the effect of the larger species on the smaller
        fij, fji = body[i] <= body[j] ? (1 + cascade, 1 - cascade) : (1 - cascade, 1 + cascade)
        M[i, j] = fij * μm[r, s] + σm[r, s] * x
        M[j, i] = fji * μm[s, r] + σm[s, r] * y
    end
    for i in 1:S
        M[i, i] = diagonal
    end
    return M
end
block_random_matrix(sizes::AbstractVector{<:Integer}; kwargs...) = block_random_matrix(Random.default_rng(), sizes; kwargs...)

"""
    modular_connectances(C, Q, α)

The within- and between-group connectances `(C_w, C_b)` of two groups holding the fractions `α`
and `1 − α` of the species, for overall connectance `C` and modularity `Q` against an
Erdős–Rényi reference (Grilli, Rogers & Allesina 2016, equations 3 and 4):
`C_w = C (1 + Q / (α² + (1 − α)²))` and `C_b = C (1 − Q / (2α(1 − α)))`. Throws an `ArgumentError`
when either falls outside `[0, 1]`, which bounds the attainable `Q`.
"""
function modular_connectances(C::Real, Q::Real, α::Real)
    0 < α < 1 || throw(ArgumentError("α must be in (0, 1), got $α"))
    0 <= C <= 1 || throw(ArgumentError("C must be in [0, 1], got $C"))
    Cw = C * (1 + Q / (α^2 + (1 - α)^2))
    Cb = C * (1 - Q / (2α * (1 - α)))
    (0 <= Cw <= 1 && 0 <= Cb <= 1) ||
        throw(ArgumentError("Q = $Q is not attainable with C = $C and α = $α: it gives C_w = $Cw and C_b = $Cb"))
    return (Cw, Cb)
end

"""
    modular_random_matrix([rng], S; α = 1/2, C, Q, μ, σ = 1, ρ = 0, diagonal = 0, cascade = 0)

The two-group random community matrix of Grilli, Rogers & Allesina (2016): `S` species, the first
`round(αS)` in group 1, overall connectance `C` and modularity `Q`, with the connectances of
[`modular_connectances`](@ref) and the same pair distribution (mean `μ`, standard deviation `σ`,
correlation `ρ`) within and between groups. `Q > 0` is modular, `Q < 0` anti-modular (bipartite-
like), and `Q = 0` unstructured. See [`block_random_matrix`](@ref) for `diagonal` and `cascade`.
"""
function modular_random_matrix(
        rng::AbstractRNG, S::Integer; α::Real = 1 / 2, C::Real, Q::Real, μ::Real, σ::Real = 1, ρ::Real = 0,
        diagonal::Real = 0, cascade::Real = 0
    )
    Cw, Cb = modular_connectances(C, Q, α)
    n = round(Int, α * S)
    0 < n < S || throw(ArgumentError("α = $α leaves a group empty with S = $S"))
    return block_random_matrix(rng, [n, S - n]; connectance = [Cw Cb; Cb Cw], μ, σ, ρ, diagonal, cascade)
end
modular_random_matrix(S::Integer; kwargs...) = modular_random_matrix(Random.default_rng(), S; kwargs...)

"""
    block_spectrum(M, groups)
    block_spectrum(sizes; connectance, μ, σ = 1, ρ = 0, diagonal = 0)

Predicted eigenvalues of a block-structured community matrix, from the moments of its blocks:
those of [`block_structure`](@ref) for a matrix `M` and its partition `groups`, or those of the
parameters of [`block_random_matrix`](@ref) without a cascade (`μ_x = C μ`,
`σ_x² = C(σ² + (1 − C)μ²)`, and for equal directions `ρ_x = (ρσ² + (1 − C)μ²) / (σ² + (1 − C)μ²)`,
Grilli, Rogers & Allesina 2016). Returns a named tuple:

- `center`, `horizontal`, `vertical`: the bulk, an ellipse centred on the real axis. Writing the
  matrix as its block means plus a remainder of mean zero, the remainder gives the bulk, centred
  at the size-weighted mean of `d_r − μ_rr`, where `d_r` is the mean diagonal of group `r`. With
  every pair correlation 0, the bulk is the disc of Aljadeff, Renfrew & Stern (2015), whose
  squared radius is the largest eigenvalue of `n_s σ_rs²`. Otherwise it is the ellipse of the
  elliptic law for the mean block variance and the variance-weighted mean correlation, exact for
  one group, for two groups of equal size (Grilli et al. 2016) and for blocks with the same
  variance and correlation, an approximation otherwise;
- `exact`: whether the bulk is one of these exact cases;
- `means`: the eigenvalues of the block means, `R[r, s] = n_s μ_rs` off the diagonal and
  `R[r, r] = (n_r − 1) μ_rr + d_r`, the eigenvalues of the expected matrix on vectors that are
  constant within groups;
- `outliers`: the real parts of the `means` that lie outside the bulk, moved by the correction of
  O'Rourke & Renfrew (2014) for the pooled moments, `θ + ρ̄ S σ̄² / θ` with `θ` measured from the
  bulk centre;
- `abscissa`: the predicted spectral abscissa, the rightmost of the bulk edge and the outliers.

Aljadeff, J., Renfrew, D. & Stern, M. Eigenvalues of block structured asymmetric random
matrices. J. Math. Phys. 56, 103502 (2015). O'Rourke, S. & Renfrew, D. Low rank perturbations of
large elliptic random matrices. Electron. J. Probab. 19, 1–65 (2014).
"""
function block_spectrum(M::AbstractMatrix{<:Real}, groups::AbstractVector)
    b = block_structure(M, groups)
    variance = replace(b.variance, NaN => 0.0)
    correlation = replace(b.correlation, NaN => 0.0)
    return spectrum_from_moments(b.sizes, replace(b.mean, NaN => 0.0), variance, correlation, b.diagonal)
end
function block_spectrum(sizes::AbstractVector{<:Integer}; connectance, μ, σ = 1, ρ = 0, diagonal::Real = 0)
    k = length(sizes)
    Cm, μm, σm, ρm = (block_parameter(x, k) for x in (connectance, μ, σ, ρ))
    mean = Cm .* μm
    variance = Cm .* (σm .^ 2 .+ (1 .- Cm) .* μm .^ 2)
    # the covariance of a pair present with probability C, with means (μ_rs, μ_sr)
    covariance = [
        Cm[r, s] * (ρm[r, s] * σm[r, s] * σm[s, r] + μm[r, s] * μm[s, r]) - Cm[r, s]^2 * μm[r, s] * μm[s, r]
            for r in 1:k, s in 1:k
    ]
    correlation = [
        variance[r, s] > 0 && variance[s, r] > 0 ? covariance[r, s] / sqrt(variance[r, s] * variance[s, r]) : 0.0
            for r in 1:k, s in 1:k
    ]
    return spectrum_from_moments(collect(sizes), mean, variance, correlation, fill(float(diagonal), k))
end

function spectrum_from_moments(sizes, μ, v, ρ, d)
    k = length(sizes)
    S = sum(sizes)
    center = sum(sizes[r] * (d[r] - μ[r, r]) for r in 1:k) / S
    weights = [sizes[r] * sizes[s] for r in 1:k, s in 1:k]
    σ̄2 = sum(weights .* v) / S^2
    ρ̄ = σ̄2 > 0 ? sum(weights .* ρ .* v) / sum(weights .* v) : 0.0
    # one group, two groups of equal size, or the same moments in every block (the elliptic law)
    homogeneous = all(≈(v[1]), v) && all(≈(ρ[1]), ρ)
    exact = all(iszero, ρ) || k == 1 || (k == 2 && sizes[1] == sizes[2]) || homogeneous
    horizontal, vertical = if all(iszero, ρ)
        radius = sqrt(max(maximum(real, eigvals([sizes[s] * v[r, s] for r in 1:k, s in 1:k])), 0.0))
        (radius, radius)
    else
        (sqrt(S * σ̄2) * (1 + ρ̄), sqrt(S * σ̄2) * (1 - ρ̄))
    end
    R = [r == s ? (sizes[r] - 1) * μ[r, r] + d[r] : sizes[s] * μ[r, s] for r in 1:k, s in 1:k]
    means = eigvals(R)
    outliers = Float64[]
    for λ in means
        θ = real(λ) - center
        abs(imag(λ)) < 1.0e-9 * max(1, abs(λ)) && abs(θ) > horizontal || continue
        push!(outliers, center + θ + ρ̄ * S * σ̄2 / θ)
    end
    abscissa = max(center + horizontal, maximum(outliers; init = -Inf))
    return (; center, horizontal, vertical, exact, means, outliers, abscissa)
end

"""
    blockplot(M, groups; size = (520, 450), title = "", colormap, linecolor = :black, linewidth = 1)
    blockplot!(ax, M, groups; colormap, linecolor = :black, linewidth = 1)

A heatmap of the square matrix `M` (a web, a community matrix or the `coassignment` of
[`compartments`](@ref)) with the species ordered by group, groups in order of first appearance in
`groups`, and lines at the boundaries between groups, so that block structure shows as blocks
(Guimerà et al. 2010, Fig. 1D). Rows run down and columns across, as the matrix is written. A
matrix with negative entries gets a diverging colour scale symmetric about 0 (`:RdBu`), any other
a sequential one (`:Blues`). Requires Makie: load a Makie backend such as CairoMakie to enable it.
`blockplot` returns the figure and axis, with a colour bar; `blockplot!` draws into an existing
axis and returns the heatmap.
"""
function blockplot end

"""
    blockplot!(ax, M, groups; kwargs...)

Draw [`blockplot`](@ref) into the Makie axis `ax`.
"""
function blockplot! end
