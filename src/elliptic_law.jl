# # Elliptic law
# For a random `S`×`S` matrix whose off-diagonal entries have variance `V` and pairs
# `(M[i, j], M[j, i])` with correlation `ρ`, the eigenvalues are uniformly distributed in an
# ellipse centered at the mean diagonal `d` with semi-axes `√(SV)(1 + ρ)` along the real axis
# and `√(SV)(1 - ρ)` along the imaginary axis (Sommers et al. 1988; Allesina & Tang 2012). A
# nonzero mean `μ` of the off-diagonal entries adds the rank-one matrix `μ (11ᵀ - I)`, whose
# eigenvalue `θ = (S - 1)μ` appears as an outlier near `d + θ + ρSV/θ` when `|θ| > √(SV)`
# (O'Rourke & Renfrew 2014).

"""
    EllipticLaw

The elliptic law prediction for a matrix, returned by [`elliptic_law`](@ref).

# Fields
- `S`: the number of species
- `μ`: the mean of the off-diagonal entries
- `μ_diag`: the mean of the diagonal entries, the center of the ellipse
- `σ2`: the variance of the off-diagonal entries
- `ρ`: the correlation of the off-diagonal pairs `(M[i, j], M[j, i])`
- `x_semi_axis`: the semi-axis along the real axis, `√(Sσ2)(1 + ρ)`
- `y_semi_axis`: the semi-axis along the imaginary axis, `√(Sσ2)(1 - ρ)`
- `outlier`: the predicted real outlier eigenvalue `μ_diag + θ + ρSσ2/θ` from the mean, with
  `θ = (S - 1)μ`, or `NaN` when `|θ| ≤ √(Sσ2)` and the mean leaves no outlier
"""
struct EllipticLaw
    S::Int
    μ::Float64
    μ_diag::Float64
    σ2::Float64
    ρ::Float64
    x_semi_axis::Float64
    y_semi_axis::Float64
    outlier::Float64
end

"""
    offdiag_pairs(M)

An `S(S - 1)/2 × 2` matrix whose rows are the off-diagonal pairs `[M[i, j], M[j, i]]`, `i < j`,
of the square matrix `M`.
"""
function offdiag_pairs(M::AbstractMatrix)
    S = LinearAlgebra.checksquare(M)
    upper = [(i, j) for i in 1:(S - 1) for j in (i + 1):S]
    return [[M[i, j] for (i, j) in upper] [M[j, i] for (i, j) in upper]]
end

"""
    offdiag_pairs_nz(M)

A `2 × n` matrix whose columns are the off-diagonal pairs `[M[i, j], M[j, i]]` for the
upper-triangle entries with `M[i, j] != 0`.
"""
function offdiag_pairs_nz(M::AbstractMatrix)
    S = LinearAlgebra.checksquare(M)
    pairs = [(M[i, j], M[j, i]) for i in 1:(S - 1) for j in (i + 1):S if M[i, j] != 0]
    return permutedims([first.(pairs) last.(pairs)])
end

"""
    shuffle_pairs([rng], M)

Randomly swap each off-diagonal pair `(M[i, j], M[j, i])` with probability 1/2, leaving the
diagonal alone. This gives a matrix whose upper and lower triangles have the same marginal
distribution while keeping the pair structure.
"""
function shuffle_pairs(rng::AbstractRNG, M::AbstractMatrix)
    S = LinearAlgebra.checksquare(M)
    shuffled = copy(M)
    #TODO: it would be nice to have some kind of degree of shuffling in the future
    for i in 1:(S - 1), j in (i + 1):S
        if rand(rng, Bool)
            shuffled[i, j], shuffled[j, i] = M[j, i], M[i, j]
        end
    end
    return shuffled
end
shuffle_pairs(M::AbstractMatrix) = shuffle_pairs(Random.default_rng(), M)

"""
    elliptic_law(M)

The [`EllipticLaw`](@ref) prediction for the eigenvalue support of `M`, estimated from the
moments of its off-diagonal pairs, including the zero pairs of a sparse matrix, so that the
connectance enters the variance: a matrix whose pairs are nonzero with probability `C` and have
variance `σ²` when nonzero has `σ2 = Cσ²`.

The bulk is centered on the mean diagonal, and a nonzero mean of the off-diagonal entries adds
the outlier of O'Rourke & Renfrew (2014). The prediction assumes the pairs are exchangeable
across the matrix; structured food webs, such as cascade webs with predators' effects above the
diagonal and prey's effects below, have outliers that this does not capture.

O'Rourke, S. & Renfrew, D. Low rank perturbations of large elliptic random matrices. Electron.
J. Probab. 19, 1–65 (2014).

When the two marginal variances of the pairs differ, the variance used is their geometric
mean, which is also the denominator of the pair correlation.
"""
function elliptic_law(M::AbstractMatrix{<:Real})
    S = LinearAlgebra.checksquare(M)
    pairs = offdiag_pairs(M)
    μ = mean(pairs)
    μ_diag = mean(diag(M))
    # population covariance of the pairs
    C = cov(pairs; dims = 1, corrected = false)
    #TODO: is the geometric mean the right way to combine unequal marginal variances in the
    # √(SV) part?
    σ2 = sqrt(C[1, 1] * C[2, 2])
    ρ = C[1, 2] / σ2
    θ = (S - 1) * μ
    outlier = abs(θ) > sqrt(S * σ2) ? μ_diag + θ + ρ * S * σ2 / θ : NaN
    return EllipticLaw(S, μ, μ_diag, σ2, ρ, sqrt(S * σ2) * (1 + ρ), sqrt(S * σ2) * (1 - ρ), outlier)
end

"""
    eigenplot(M; size = (500, 500), title = "", kwargs...)
    eigenplot!(ax, M; kwargs...)

Plot the eigenvalues of the square matrix `M` in the complex plane, optionally with the
[`elliptic_law`](@ref) prediction. Requires Makie: load a Makie backend such as CairoMakie to
enable it. `eigenplot` returns the figure and axis; `eigenplot!` draws into an existing axis and
returns the axis.

Keywords:
- `ellipse`: draw the elliptic-law ellipse estimated from the off-diagonal pairs of `M`
  (default `true`);
- `outlier`: mark the predicted outlier from a nonzero mean, when there is one (default `true`);
- `stability`: draw the stability boundary `Re λ = 0` (default `true`);
- `law`: the [`EllipticLaw`](@ref) to draw, by default `elliptic_law(M)`;
- `expected`: a matrix whose eigenvalues to mark as crosses, such as the
  [`expected_community_matrix`](@ref) of a structured web, whose eigenvalues give the
  imaginary spread of the spectrum (default `nothing`), with `expectedcolor`;
- `color`, `markersize`: for the eigenvalues;
- `ellipsecolor`, `outliercolor`: for the prediction;
- `equal_aspect`: show a square window with equal scales on both axes, so the ellipse keeps its
  shape (default `true`); with `false` Makie fits each axis to the data.
"""
function eigenplot end

"""
    eigenplot!(ax, M; kwargs...)

Draw the eigenvalues of `M` into the Makie axis `ax`; see [`eigenplot`](@ref) for the keywords.
"""
function eigenplot! end
