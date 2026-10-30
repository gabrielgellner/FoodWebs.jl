using FoodWebs
using Test
using Random: Xoshiro

@testset "resilience" begin
    A1 = [-1.0 0; 0 -10]
    A2 = [-1.0 15; 0 -10]
    @test λ_stability(A1) == λ_stability(A2)
    @test ν_stability(A1) ≈ -1.0
    # symmetric part [-1 7.5; 7.5 -10] has largest eigenvalue (-11 + √306) / 2
    @test ν_stability(A2) ≈ (-11 + sqrt(306)) / 2

    # complex eigenvalues: rotation with damping
    @test λ_stability([-0.5 -2.0; 2.0 -0.5]) ≈ -0.5

    # the numerical abscissa bounds the spectral abscissa, with equality for normal matrices
    rng = Xoshiro(10)
    for _ in 1:50
        A = randn(rng, 10, 10)
        @test ν_stability(A) >= λ_stability(A) - 1.0e-10
        Asym = A + A'
        @test ν_stability(Asym) ≈ λ_stability(Asym)
    end
end
