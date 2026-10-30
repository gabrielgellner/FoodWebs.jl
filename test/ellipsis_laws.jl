using FoodWebs
using Test
using Random: Xoshiro
using LinearAlgebra: eigvals, diag, I

@testset "off-diagonal pairs" begin
    M = [1 2 3; 4 5 6; 7 8 9]
    @test offdiag_pairs(M) == [2 4; 3 7; 6 8]
    @test offdiag_pairs_nz([0 1 0; 2 0 3; 0 4 0]) == [1 3; 2 4]
    @test size(offdiag_pairs_nz(zeros(3, 3))) == (2, 0)
    @test size(offdiag_pairs(randn(Xoshiro(1), 20, 20))) == (20 * 19 ÷ 2, 2)
end

@testset "shuffle_pairs" begin
    rng = Xoshiro(11)
    M = randn(rng, 30, 30)
    Ms = shuffle_pairs(rng, M)
    @test diag(Ms) == diag(M)
    # every pair is kept, possibly swapped
    @test all(Set((Ms[i, j], Ms[j, i])) == Set((M[i, j], M[j, i])) for i in 1:30 for j in 1:30)
    @test Ms != M
    @test shuffle_pairs(Xoshiro(1), M) == shuffle_pairs(Xoshiro(1), M)
end

@testset "elliptic law" begin
    # entries with variance σ² and pair correlation ρ: eigenvalues fill an ellipse centered
    # at the diagonal d with semi-axes √(Sσ²)(1 ± ρ)
    rng = Xoshiro(12)
    S, σ, ρ, d = 600, 1 / sqrt(600), 0.5, -1.0
    M = d * Matrix(1.0I, S, S)
    for i in 1:(S - 1), j in (i + 1):S
        z1, z2 = randn(rng), randn(rng)
        M[i, j] = σ * z1
        M[j, i] = σ * (ρ * z1 + sqrt(1 - ρ^2) * z2)
    end
    sol = ellipsis_sample(M)
    @test sol.S == S
    @test sol.ρ ≈ ρ atol = 0.02
    @test sol.σ2 ≈ σ^2 rtol = 0.02
    @test sol.μ_diag ≈ d
    λ = eigvals(M)
    @test maximum(real, λ) - d ≈ sol.x_semi_axis rtol = 0.1
    @test maximum(imag, λ) ≈ sol.y_semi_axis rtol = 0.1
end
