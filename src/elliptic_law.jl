# # Elliptic law
# For a random `S`×`S` matrix whose off-diagonal entries have variance `V` and pairs
# `(M[i, j], M[j, i])` with correlation `ρ`, the eigenvalues are uniformly distributed in an
# ellipse centered at the mean diagonal `d` with semi-axes `√(SV)(1 + ρ)` along the real axis
# and `√(SV)(1 - ρ)` along the imaginary axis (Sommers et al. 1988; Allesina & Tang 2012).

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
"""
struct EllipticLaw
    S::Int
    μ::Float64
    μ_diag::Float64
    σ2::Float64
    ρ::Float64
    x_semi_axis::Float64
    y_semi_axis::Float64
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
moments of its off-diagonal pairs.

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
    return EllipticLaw(S, μ, μ_diag, σ2, ρ, sqrt(S * σ2) * (1 + ρ), sqrt(S * σ2) * (1 - ρ))
end
