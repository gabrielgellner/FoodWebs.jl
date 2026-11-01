using FoodWebs
using Test
using Random: Xoshiro
using LinearAlgebra: eigvals, diag, I
using Statistics: mean, quantile
using Distributions: Normal
using CairoMakie: CairoMakie

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
    sol = elliptic_law(M)
    @test sol isa FoodWebs.EllipticLaw
    @test sol.S == S
    @test sol.ρ ≈ ρ atol = 0.02
    @test sol.σ2 ≈ σ^2 rtol = 0.02
    @test sol.μ_diag ≈ d
    λ = eigvals(M)
    @test maximum(real, λ) - d ≈ sol.x_semi_axis rtol = 0.1
    @test maximum(imag, λ) ≈ sol.y_semi_axis rtol = 0.1
end

# A sparse `S`×`S` matrix with diagonal `d` whose off-diagonal pairs are zero with probability
# `1 - C`, and otherwise have mean `μ`, variance `σ²` and correlation `ρ`, drawn from normal or
# uniform marginals (Allesina, Theoretical Community Ecology, ch. 4)
function sparse_pairs(rng, S, C, d; σ = 1.0, ρ = 0.0, μ = 0.0, uniform = false)
    M = zeros(S, S)
    draw() = uniform ? sqrt(12) * (rand(rng) - 0.5) : randn(rng)
    for i in 1:(S - 1), j in (i + 1):S
        rand(rng) < C || continue
        z1, z2 = draw(), draw()
        M[i, j] = μ + σ * z1
        M[j, i] = μ + σ * (ρ * z1 + sqrt(1 - ρ^2) * z2)
    end
    M[CartesianIndex.(1:S, 1:S)] .= d
    return M
end

# the largest real part, measured from the diagonal, and the largest imaginary part
function extents(M, d)
    λ = eigvals(M)
    return maximum(real, λ) - d, maximum(imag, λ)
end

# An observed edge of the spectrum agrees with a predicted semi-axis to within a tenth of the
# radius scale √(Sσ²). Edge eigenvalues fluctuate on that scale, so a narrow ellipse's short
# axis cannot be held to a relative tolerance of its own length.
near_edge(observed, predicted, sol) = abs(observed - predicted) <= 0.1 * sqrt(sol.S * sol.σ2)

@testset "May's sparse matrices" begin
    # the eigenvalues of a matrix with connectance C and entry variance σ² fill a disc of
    # radius √(SCσ²), whatever the distribution of the entries
    rng = Xoshiro(13)
    S, C, σ, d = 500, 0.2, 1.0, -10.0
    for uniform in (false, true)
        M = sparse_pairs(rng, S, C, d; σ, uniform)
        sol = elliptic_law(M)
        @test sol.x_semi_axis ≈ sqrt(S * C * σ^2) rtol = 0.03
        @test sol.y_semi_axis ≈ sqrt(S * C * σ^2) rtol = 0.03
        re, im = extents(M, d)
        @test near_edge(re, sqrt(S * C * σ^2), sol)
        @test near_edge(im, sqrt(S * C * σ^2), sol)
        @test isnan(sol.outlier)
    end
end

@testset "correlated pairs and interaction types" begin
    # with pair correlation ρ the disc stretches to semi-axes √(SCσ²)(1 ± ρ)
    rng = Xoshiro(14)
    S, C, d = 500, 0.2, -10.0
    abscissa = Dict{Float64, Float64}()
    for ρ in (-0.5, 0.0, 0.5)
        M = sparse_pairs(rng, S, C, d; ρ)
        sol = elliptic_law(M)
        @test sol.ρ ≈ ρ atol = 0.03
        @test sol.x_semi_axis ≈ sqrt(S * C) * (1 + ρ) rtol = 0.05
        @test sol.y_semi_axis ≈ sqrt(S * C) * (1 - ρ) rtol = 0.05
        re, im = extents(M, d)
        @test near_edge(re, sol.x_semi_axis, sol)
        @test near_edge(im, sol.y_semi_axis, sol)
        abscissa[ρ] = re + d
    end
    # Allesina & Tang (2012): at equal S, C and σ, predator-prey pairs (ρ < 0) are the easiest
    # to stabilize and competition or mutualism pairs (ρ > 0) the hardest
    @test abscissa[-0.5] < abscissa[0.0] < abscissa[0.5]

    # random_predprey pairs a negative and a positive entry of independent half-normal sizes,
    # so with f = 1 the pair correlation is -E|a| E|b| / E[a²] = -2/π
    web = may_network(rng, 400, 0.1)
    J = random_predprey(rng, web, Normal(), 1.0)
    sol = elliptic_law(J)
    @test sol.ρ ≈ -2 / π atol = 0.03
    @test abs(sol.μ) < 0.01
    re, im = extents(J, -1.0)
    @test near_edge(re, sol.x_semi_axis, sol)
    @test near_edge(im, sol.y_semi_axis, sol)
end

@testset "outliers from a nonzero mean" begin
    # a mean μ of the off-diagonal entries adds the outlier d + θ + ρSσ²/θ, θ = (S - 1)μ
    # (O'Rourke & Renfrew 2014); the correction ρSσ²/θ matters when the pairs are correlated
    rng = Xoshiro(15)
    S, C, d = 500, 0.2, -5.0
    for (ρ, μ) in ((0.0, 0.5), (0.5, 0.5), (-0.5, 0.5), (0.5, -0.6))
        M = sparse_pairs(rng, S, C, d; ρ, μ)
        sol = elliptic_law(M)
        λ = eigvals(M)
        observed = μ > 0 ? maximum(real, λ) : minimum(real, λ)
        @test abs(observed - sol.outlier) < 1.0
        # the outlier lies outside the bulk
        @test abs(sol.outlier - d) > sol.x_semi_axis
    end
    M = sparse_pairs(rng, S, C, d; ρ = 0.5, μ = 0.5)
    sol = elliptic_law(M)
    θ = (S - 1) * sol.μ
    @test abs(maximum(real, eigvals(M)) - sol.outlier) < abs(maximum(real, eigvals(M)) - (d + θ))
end

@testset "structured predator-prey matrices" begin
    # Structured webs break the exchangeable pairs the elliptic law assumes. The eigenvalues of
    # the expected community matrix given the web give the spread of imaginary parts, including
    # the largest, which the elliptic law misses; the real parts, and so stability, come from
    # the fluctuations around it.
    rng = Xoshiro(16)
    S = 300
    for (gen, f) in ((cascade_network, 1.0), (cascade_network, 0.3), (niche_network, 1.0), (niche_network, 0.3))
        web = gen(rng, S, 0.1)
        J = random_predprey(rng, web, Normal(), f)
        λ, μ = eigvals(J), eigvals(expected_community_matrix(web, Normal(), f))
        # the largest imaginary part to within about 1% (8 seeds), and the fluctuations widen
        # the 90% quantile of |Im λ| by 2-16%
        @test maximum(imag, λ) ≈ maximum(imag, μ) rtol = 0.03
        @test 0.98 < quantile(abs.(imag.(λ)), 0.9) / quantile(abs.(imag.(μ)), 0.9) < 1.25
        @test maximum(imag, λ) > 2 * elliptic_law(J).y_semi_axis
        # the expected matrix is stable; the fluctuations push the spectral abscissa above zero
        @test maximum(real, μ) < 0 < maximum(real, λ)
        # shuffling who is the predator makes the pairs exchangeable again; in niche webs the
        # uneven numbers of links per species still widen the edge beyond the ellipse
        if gen === cascade_network
            Js = random_predprey(rng, web, Normal(), f; shuffle_prey = true)
            sol = elliptic_law(Js)
            @test near_edge(extents(Js, -1.0)[2], sol.y_semi_axis, sol)
        end
    end
end

@testset "eigenplot" begin
    rng = Xoshiro(17)
    M = sparse_pairs(rng, 100, 0.2, -5.0; μ = 0.5)
    fig, ax = eigenplot(M; title = "spectrum")
    @test fig isa CairoMakie.Figure
    path = tempname() * ".png"
    CairoMakie.save(path, fig)
    @test filesize(path) > 0
    # eigenvalues, ellipse, outlier and stability line
    @test count(p -> p isa CairoMakie.Scatter, ax.scene.plots) == 2
    @test count(p -> p isa CairoMakie.Lines, ax.scene.plots) == 1
    # each part can be left out, and an existing axis or a given law can be used
    fig = CairoMakie.Figure()
    ax = CairoMakie.Axis(fig[1, 1])
    @test eigenplot!(ax, M; ellipse = false, outlier = false, stability = false, equal_aspect = false) === ax
    @test count(p -> p isa CairoMakie.Scatter, ax.scene.plots) == 1
    @test isempty(filter(p -> p isa CairoMakie.Lines, ax.scene.plots))
    @test eigenplot!(CairoMakie.Axis(fig[1, 2]), M; law = elliptic_law(M)) isa CairoMakie.Axis
    # the eigenvalues of an expected matrix as a further layer
    web = cascade_network(rng, 60, 0.1)
    J = random_predprey(rng, web, Normal(), 1.0)
    fig, ax = eigenplot(J; expected = expected_community_matrix(web, Normal(), 1.0))
    @test count(p -> p isa CairoMakie.Scatter, ax.scene.plots) == 2
end
