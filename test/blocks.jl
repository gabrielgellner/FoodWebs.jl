using FoodWebs
using Test
using Random: Xoshiro
using LinearAlgebra: eigvals
using Statistics: mean

# the eigenvalue of `λ` with real part nearest to `x`
nearest(λ, x) = real(λ[argmin(abs.(real.(λ) .- x))])

@testset "block random matrices" begin
    rng = Xoshiro(87)
    # prey in group 1 and their predators in group 2: the predators depress the prey (μ[1, 2] < 0)
    # and gain from them (μ[2, 1] > 0), more variably in one direction than the other
    sizes = [150, 250]
    connectance = [0.2 0.1; 0.1 0.3]
    μ = [-0.5 -1.0; 0.4 -0.2]
    σ = [1.0 0.8; 0.5 1.0]
    ρ = [-0.3 -0.6; -0.6 0.2]
    M = block_random_matrix(rng, sizes; connectance, μ, σ, ρ, diagonal = -2)
    groups = [fill(1, 150); fill(2, 250)]
    @test size(M) == (400, 400) && all(M[i, i] == -2 for i in 1:400)
    # interactions come in pairs
    @test all((M[i, j] != 0) == (M[j, i] != 0) for i in 1:400, j in 1:400 if i != j)
    # the block moments of the sample are those predicted from the parameters
    b = block_structure(M, groups)
    C = connectance
    @test b.connectance ≈ C atol = 0.01
    @test b.mean ≈ C .* μ atol = 0.02
    @test b.variance ≈ C .* (σ .^ 2 .+ (1 .- C) .* μ .^ 2) atol = 0.03
    cov12 = C[1, 2] * (ρ[1, 2] * σ[1, 2] * σ[2, 1] + μ[1, 2] * μ[2, 1]) - C[1, 2]^2 * μ[1, 2] * μ[2, 1]
    @test b.correlation[1, 2] ≈ cov12 / sqrt(b.variance[1, 2] * b.variance[2, 1]) atol = 0.04
    # opposite signs with probability about 0.89 × 0.79 from the means and spreads
    @test b.predator_prey[1, 2] > 0.7 * count(!iszero, M[1:150, 151:400])

    @test M == block_random_matrix(Xoshiro(87), sizes; connectance, μ, σ, ρ, diagonal = -2)
    @test_throws ArgumentError block_random_matrix(sizes; connectance = [0.2 0.1; 0.3 0.2], μ = 0)
    @test_throws ArgumentError block_random_matrix(sizes; connectance = 0.2, μ = 0, ρ = [0 0.5; 0 0])
    @test_throws ArgumentError block_random_matrix(sizes; connectance = 1.5, μ = 0)
    @test_throws DimensionMismatch block_random_matrix(sizes; connectance = 0.2, μ = zeros(3, 3))
end

@testset "modular connectances" begin
    # equations 3 and 4 of Grilli et al. (2016); Q = 0 is unstructured
    @test modular_connectances(0.2, 0.0, 0.3) == (0.2, 0.2)
    Cw, Cb = modular_connectances(0.2, 0.25, 0.5)
    @test Cw ≈ 0.2 * 1.5 && Cb ≈ 0.2 * 0.5
    # the overall connectance is kept: the fraction of within-group pairs is α² + (1 − α)²
    α = 0.25
    Cw, Cb = modular_connectances(0.2, -0.3, α)
    @test (α^2 + (1 - α)^2) * Cw + 2α * (1 - α) * Cb ≈ 0.2
    M = modular_random_matrix(Xoshiro(1), 1000; α, C = 0.2, Q = -0.3, μ = -1)
    @test count(!iszero, M) / (1000 * 999) ≈ 0.2 atol = 0.005
    # Q beyond the attainable range leaves a connectance outside [0, 1]
    @test_throws ArgumentError modular_connectances(0.2, 0.6, 0.5)
    @test_throws ArgumentError modular_connectances(0.8, 0.4, 0.5)
end

@testset "block spectrum against simulation" begin
    rng = Xoshiro(88)
    S = 800
    for (α, μ, ρ, Q) in ((1 / 4, -1.0, 0.0, 0.3), (1 / 4, -1.0, 0.0, -0.3), (1 / 4, 1.0, -0.5, 0.2), (1 / 2, -0.5, -0.5, 0.25))
        M = modular_random_matrix(rng, S; α, C = 0.2, Q, μ, ρ)
        λ = eigvals(M)
        n = round(Int, α * S)
        Cw, Cb = modular_connectances(0.2, Q, α)
        p = block_spectrum([n, S - n]; connectance = [Cw Cb; Cb Cw], μ, ρ)
        # the prediction from the sampled block moments agrees with that from the parameters
        e = block_spectrum(M, [fill(1, n); fill(2, S - n)])
        @test e.abscissa ≈ p.abscissa rtol = 0.03
        # outliers within 2% or 3 (outliers near the bulk fluctuate by an amount of order 1); the
        # rightmost eigenvalue within 3% when an outlier sets it, and within 10% when the bulk edge
        # does, as finite matrices overshoot it and the pooled ellipse is an approximation here
        for o in p.outliers
            @test isapprox(nearest(λ, o), o; rtol = 0.02, atol = 3)
        end
        bulk = p.abscissa == p.center + p.horizontal
        @test maximum(real, λ) ≈ p.abscissa rtol = (bulk ? 0.1 : 0.03)
    end
    # the exact cases: equal group sizes and blocks with the same moments
    @test block_spectrum([300, 300]; connectance = [0.3 0.1; 0.1 0.3], μ = -1, ρ = -0.5).exact
    @test block_spectrum([200, 400]; connectance = 0.2, μ = 0, ρ = -0.5).exact
    @test !block_spectrum([200, 400]; connectance = [0.3 0.1; 0.1 0.3], μ = -1, ρ = -0.5).exact

    # food-web blocks with means that differ by direction
    sizes = [300, 500]
    connectance = [0.15 0.1; 0.1 0.2]
    μ = [-0.3 -1.5; 1.0 -0.3]
    M = block_random_matrix(rng, sizes; connectance, μ, ρ = -0.4, diagonal = -1)
    λ = eigvals(M)
    p = block_spectrum(sizes; connectance, μ, ρ = -0.4, diagonal = -1)
    for o in p.outliers
        @test isapprox(nearest(λ, o), o; rtol = 0.03, atol = 3)
    end
    @test maximum(real, λ) ≈ p.abscissa rtol = (p.abscissa == p.center + p.horizontal ? 0.1 : 0.03)

    # one group: the outlier of the elliptic law
    M = block_random_matrix(rng, [600]; connectance = 0.2, μ = 0.5, ρ = 0.3)
    p = block_spectrum(M, ones(Int, 600))
    @test only(p.outliers) ≈ elliptic_law(M).outlier rtol = 0.01
end

@testset "modularity and stability (Grilli et al. 2016)" begin
    # the spectral abscissa with block structure against that without. For two groups of equal
    # size and μ < 0, modularity stabilizes moderately and anti-modularity destabilizes strongly.
    # With μ > 0 the leading outlier sits at (S/2)(μ_w + |μ_b|), which is SCμ for every Q when
    # α = 1/2, so equal groups are neutral, while unequal groups (α = 1/4) destabilize: the
    # predicted outlier rises from 100 to 112 at Q = 0.3 and to 111 at Q = −0.5
    rng = Xoshiro(89)
    abscissa(Q, μ; α = 1 / 2) = mean(maximum(real, eigvals(modular_random_matrix(rng, 500; α, C = 0.2, Q, μ))) for _ in 1:3)
    base = abscissa(0.0, -1.0)
    @test abscissa(0.3, -1.0) < base
    @test abscissa(-0.3, -1.0) > 2 * base
    base = abscissa(0.0, 1.0)
    @test abscissa(0.3, 1.0) ≈ base rtol = 0.03
    @test abscissa(-0.3, 1.0) ≈ base rtol = 0.03
    base = abscissa(0.0, 1.0; α = 1 / 4)
    @test abscissa(0.3, 1.0; α = 1 / 4) > 1.05 * base && abscissa(-0.5, 1.0; α = 1 / 4) > 1.05 * base

    # the cascade variant: with μ < 0 and ξ = 3 each pair has one negative and one positive mean
    M = block_random_matrix(rng, [200]; connectance = 0.2, μ = -1, σ = 0.1, cascade = 3)
    pairs = [(M[i, j], M[j, i]) for i in 1:200 for j in (i + 1):200 if M[i, j] != 0]
    @test all(p -> p[1] * p[2] < 0, pairs)
    @test mean(p -> p[1] + p[2], pairs) / 2 ≈ -1 atol = 0.05
end
