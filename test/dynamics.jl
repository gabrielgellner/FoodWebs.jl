using FoodWebs
using Test
using LinearAlgebra: Diagonal, norm, triu
using Random: Xoshiro
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqTsit5: Tsit5
using PolytopeSamplers: PolytopeSamplers
using SciMLBase: SciMLBase

@testset "return time and amplification" begin
    @test return_time([-1.0 0; 0 -10]) ≈ 1
    @test return_time([-1.0 15; 0 -10]) ≈ 1
    @test return_time([0.5 0; 0 -1]) == Inf

    # a normal matrix has no transient growth: its envelope is exp(λ t)
    N = [-1.0 2; 2 -5]
    λ = spectral_abscissa(N)
    for t in (0.0, 0.3, 2.0)
        @test amplification_envelope(N, t) ≈ exp(λ * t)
    end
    @test max_amplification(N) == (amplification = 1.0, time = 0.0)
    @test_throws ArgumentError max_amplification([0.5 0; 0 -1])

    rng = Xoshiro(20)
    for _ in 1:10
        # eigenvalues -0.5, -1, ..., -4 on the diagonal, and large entries above it make it reactive
        J = Diagonal(-0.5 .* (1:8)) + 3 .* triu(randn(rng, 8, 8), 1)
        @test spectral_abscissa(J) ≈ -0.5 && numerical_abscissa(J) > 0
        # the envelope starts at 1 with slope ν and decays at rate λ
        h = 1.0e-6
        @test amplification_envelope(J, 0.0) ≈ 1
        @test (amplification_envelope(J, h) - 1) / h ≈ numerical_abscissa(J) rtol = 1.0e-3
        @test log(amplification_envelope(J, 61.0) / amplification_envelope(J, 60.0)) ≈ spectral_abscissa(J) rtol = 1.0e-3
        # the maximum lies above the envelope on a fine grid, and every perturbation lies below it
        m = max_amplification(J)
        @test m.amplification > 1 && m.time > 0
        @test m.amplification ≈ amplification_envelope(J, m.time)
        grid = range(0, 40; length = 4001)
        @test m.amplification >= maximum(t -> amplification_envelope(J, t), grid) - 1.0e-9
        δ = randn(rng, 8)
        @test all(norm(exp(t * J) * δ) / norm(δ) <= amplification_envelope(J, t) + 1.0e-9 for t in grid)
    end
end

# the bounded polytopes of the feasible tests, with biomasses over two orders of magnitude
dynamics_communities(seed) = FeasibleCommunities(
    niche_network(Xoshiro(seed + 500), 12, 0.15; basal_paths = true),
    exp10.(2 .* rand(Xoshiro(seed + 900), 12) .- 1); efficiency = (0.0, 0.3)
)

@testset "Lotka-Volterra dynamics" begin
    # logistic growth, x(t) = K / (1 + (K / x0 - 1) exp(-r t)) with K = r / a
    r, a, x0 = 0.8, 0.4, 0.1
    s = lotka_volterra(Tsit5(), fill(-a, 1, 1), [r], [x0], (0.0, 10.0); saveat = 0:0.5:10, abstol = 1.0e-12, reltol = 1.0e-10)
    @test s.retcode == SciMLBase.ReturnCode.Success
    @test s.t == 0:0.5:10
    @test vec(s.x) ≈ [(r / a) / (1 + (r / a / x0 - 1) * exp(-r * t)) for t in s.t] rtol = 1.0e-8
    @test s.survivors == [true] && s.extinction_times == [Inf]

    # the analytic Jacobian matches finite differences
    rng = Xoshiro(21)
    A, r = randn(rng, 5, 5), randn(rng, 5)
    prob = lotka_volterra_problem(A, r, rand(rng, 5), (0.0, 1.0))
    x = rand(rng, 5)
    f(x) = (dx = similar(x); prob.f(dx, x, prob.p, 0.0); dx)
    Jx = zeros(5, 5)
    prob.f.jac(Jx, x, prob.p, 0.0)
    h = 1.0e-7
    @test Jx ≈ reduce(hcat, (f(x .+ h .* (1:5 .== j)) .- f(x)) ./ h for j in 1:5) rtol = 1.0e-5
    @test Jx ≈ Diagonal(r .+ A * x) + Diagonal(x) * A

    # a consumer without prey declines as x0 exp(r t) and goes extinct when it reaches the
    # threshold; one starting below the threshold is absent from the start
    s = lotka_volterra(
        Tsit5(), [0.0 0; 0 -1], [-1.0, 1.0], [1.0, 0.5], (0.0, 30.0);
        extinction_threshold = [1.0e-6, 1.0e-6], abstol = 1.0e-12, reltol = 1.0e-10
    )
    @test s.survivors == [false, true]
    @test s.extinction_times[1] ≈ log(1.0e6) rtol = 1.0e-6
    @test s.x[1, end] == 0 && s.x[2, end] ≈ 1
    s = lotka_volterra(Tsit5(), [-1.0;;], [1.0], [1.0e-10], (0.0, 30.0))
    @test s.survivors == [false] && s.extinction_times == [0.0] && all(==(0), s.x)
    @test_throws ArgumentError lotka_volterra(Tsit5(), [-1.0;;], [1.0], [-1.0], (0.0, 1.0))
    @test_throws DimensionMismatch lotka_volterra(Tsit5(), [-1.0;;], [1.0], [1.0, 1.0], (0.0, 1.0))
    @test_throws DimensionMismatch lotka_volterra(Tsit5(), [-1.0;;], [1.0], [1.0], (0.0, 1.0); extinction_threshold = [1.0, 1.0])

    # with A + A' negative definite the equilibrium is globally stable (Goh 1977), and the
    # dynamics return to it from starts off by up to a factor of 10
    A = [-1.0 -0.5 0; 0.2 -1 -0.5; 0 0.2 -1]
    x̂ = [2.0, 1.0, 0.5]
    r = -A * x̂
    rng = Xoshiro(22)
    for _ in 1:5
        s = lotka_volterra(Rodas5P(), A, r, x̂ .* exp10.(2 .* rand(rng, 3) .- 1), (0.0, 200.0); abstol = 1.0e-12, reltol = 1.0e-10)
        @test all(s.survivors)
        @test s.x[:, end] ≈ x̂ rtol = 1.0e-6
    end

    for seed in 1:3
        fc = dynamics_communities(seed)
        samples = sample_communities(Xoshiro(seed), fc; target_ess = 200).communities
        stable = samples[findfirst(c -> spectral_abscissa(c.J) < 0, samples)]
        T = return_time(stable.J)
        # a stable community returns to its equilibrium from a small perturbation, with either
        # solver; from a large one it can reach another attractor
        start = fc.x .* (1 .+ 0.01 .* randn(Xoshiro(seed + 1), 12))
        for alg in (Tsit5(), Rodas5P())
            s = lotka_volterra(alg, stable, start, (0.0, 100T); abstol = 1.0e-12, reltol = 1.0e-10)
            @test all(s.survivors)
            @test s.x[:, end] ≈ fc.x rtol = 1.0e-6
        end
        @test lotka_volterra(Tsit5(), stable, start, (0.0, T)).x == lotka_volterra(Tsit5(), stable.A, stable.r, start, (0.0, T)).x

        # small pulses follow the linearization, with an error that shrinks with the pulse
        δ = fc.x .* randn(Xoshiro(seed + 2), 12)
        times = range(0, 3T; length = 61)
        err(ε) = (
            p = pulse_response(Rodas5P(), stable, fc.x, ε .* δ, times; abstol = 1.0e-14, reltol = 1.0e-12);
            maximum(abs.(p.deviation .- p.linear))
        )
        @test 5 < err(1.0e-3) / err(1.0e-4) < 20
        p = pulse_response(Rodas5P(), stable, fc.x, 1.0e-4 .* δ, times; abstol = 1.0e-14, reltol = 1.0e-12)
        @test p.deviation[1] ≈ 1 && p.linear[1] ≈ 1
        @test p.amplification == maximum(p.deviation) && p.time in times
        @test all(p.linear .<= [amplification_envelope(stable.J, t) for t in p.t] .+ 1.0e-9)

        # an unstable community leaves its equilibrium for another state, with abundances that
        # stay nonnegative and extinct species at zero; at the default tolerances Rodas5P damps a
        # weakly unstable oscillation (abscissa 0.02, frequency 3.4) back to the equilibrium
        k = findfirst(c -> spectral_abscissa(c.J) > 0, samples)
        k === nothing && continue
        s = lotka_volterra(
            Rodas5P(), samples[k], fc.x .* (1 .+ 1.0e-3 .* randn(Xoshiro(seed + 3), 12)), (0.0, 1.0e3);
            abstol = 1.0e-12, reltol = 1.0e-10
        )
        @test s.retcode == SciMLBase.ReturnCode.Success
        @test all(>=(0), s.x)
        @test s.survivors == (s.x[:, end] .> 0)
        @test norm(s.x[:, end] .- fc.x) > 1.0e-3 * norm(fc.x)
    end

    fc = dynamics_communities(1)
    c = sample_communities(Xoshiro(1), fc; target_ess = 100).communities[1]
    @test_throws ArgumentError pulse_response(Tsit5(), c, 2 .* fc.x, 0.1 .* fc.x, 0:1)
    @test_throws ArgumentError pulse_response(Tsit5(), c, fc.x, -2 .* fc.x, 0:1)
    @test_throws DimensionMismatch pulse_response(Tsit5(), c, fc.x, [0.1], 0:1)
end
