using FoodWebs
using Test
using Random: Xoshiro

@testset "stability" begin
    @test λ_stability === spectral_abscissa
    @test ν_stability === numerical_abscissa

    A1 = [-1.0 0; 0 -10]
    A2 = [-1.0 15; 0 -10]
    @test spectral_abscissa(A1) == spectral_abscissa(A2)
    @test numerical_abscissa(A1) ≈ -1.0
    # symmetric part [-1 7.5; 7.5 -10] has largest eigenvalue (-11 + √306) / 2
    @test numerical_abscissa(A2) ≈ (-11 + sqrt(306)) / 2

    # complex eigenvalues: rotation with damping
    @test spectral_abscissa([-0.5 -2.0; 2.0 -0.5]) ≈ -0.5

    # the numerical abscissa bounds the spectral abscissa, with equality for normal matrices
    rng = Xoshiro(10)
    for _ in 1:50
        A = randn(rng, 10, 10)
        @test numerical_abscissa(A) >= spectral_abscissa(A) - 1.0e-10
        Asym = A + A'
        @test numerical_abscissa(Asym) ≈ spectral_abscissa(Asym)
    end
end
