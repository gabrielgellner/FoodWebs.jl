using FoodWebs
using Test
using Random: Xoshiro
using Distributions: Exponential, Uniform, pdf
using PolytopeSamplers: PolytopeSamplers, CRHMC, Polytope
using Statistics: mean, std, var
using LinearAlgebra: Diagonal

# A three-species chain (2 eats 1, 3 eats 2) with fixed rates, a self-regulated resource and an
# efficiency cap of 0.5. The equilibrium rows give λ₁₂ = (1 - 3σ) / 2, γ₂₁ = (λ₂₃ + 0.2) / 3 and
# γ₃₂ = 0.05, so the polytope is the set of (σ, λ₂₃) with γ₂₁ ≤ λ₁₂ / 2, γ₃₂ ≤ λ₂₃ / 2 and every
# magnitude in [0, 5]. Drawing (σ, λ₂₃) uniformly and keeping the feasible draws samples it
# exactly, as the other coordinates are linear in them.
const CHAIN_X = [3.0, 2.0, 1.0]
const CHAIN_R = [1.0, -0.2, -0.1]
chain_communities() = FeasibleCommunities(
    food_chain_network(3), CHAIN_X; rates = CHAIN_R, self_regulation = :basal, efficiency = (0.0, 0.5),
    bounds = (0.0, 5.0)
)
function chain_rejection(rng, n)
    kept = Tuple{Float64, Float64, Float64, Float64}[]
    while length(kept) < n
        σ, λ23 = rand(rng) / 3, 5 * rand(rng)
        λ12, γ21, γ32 = (1 - 3σ) / 2, (λ23 + 0.2) / 3, 0.05
        γ21 <= λ12 / 2 && γ32 <= λ23 / 2 && push!(kept, (σ, λ12, γ21, λ23))
    end
    return kept
end

@testset "polytope and emptiness" begin
    fc = chain_communities()
    P, H = Polytope(fc), halfspaces(fc)
    @test P.lower == H.lower && P.upper == H.upper && P.A == H.A && P.C == H.C && P.d == H.d
    @test !isempty(fc)
    # rates that fail a necessary condition, and a budget that cannot balance
    @test isempty(FeasibleCommunities(food_chain_network(2), ones(2); rates = [-1.0, 0.5], self_regulation = :none))
    @test isempty(FeasibleCommunities(food_chain_network(2), ones(2); rates = [0.1, -1.0], efficiency = 0.5))
end

@testset "uniform sampling" begin
    fc = chain_communities()
    res = sample_communities(Xoshiro(31), fc; target_ess = 1000)
    @test res.converged
    θ = res.points
    k(kind, i, j) = findfirst(c -> c == (kind = kind, i = i, j = j), coordinates(fc))
    σ, λ23 = θ[k(:self, 1, 1), :], θ[k(:loss, 2, 3), :]
    exact = chain_rejection(Xoshiro(32), 200_000)
    ess = res.diagnostics.min_ess
    for (sampled, reference) in ((σ, first.(exact)), (λ23, last.(exact)))
        # means within 4 Monte Carlo standard errors, and matching spread
        @test abs(mean(sampled) - mean(reference)) < 4 * std(reference) / sqrt(ess)
        @test std(sampled) ≈ std(reference) rtol = 0.1
    end
    # every sample is an equilibrium with the chain's signs
    for c in res.communities[1:50:end]
        @test c.r .+ c.A * CHAIN_X ≈ zeros(3) atol = 1.0e-7
        @test c.A[1, 2] <= 0 <= c.A[2, 1] && c.A[2, 3] <= 0 <= c.A[3, 2]
        @test c.J ≈ CHAIN_X .* c.A
    end
end

@testset "two-species segment" begin
    # 2 eats 1, only the resource self-regulates: σ is uniform on [0, (r₁ - x₂ γ / e) / x₁]
    x, r, e = [2.0, 0.5], [3.0, -0.4], 0.25
    fc = FeasibleCommunities(Bool[0 1; 0 0], x; rates = r, self_regulation = :basal, efficiency = (0.0, e))
    top = (r[1] - x[2] * (-r[2] / x[1]) / e) / x[1]
    res = sample_communities(Xoshiro(33), fc; target_ess = 1000)
    σ = res.points[findfirst(c -> c.kind === :self, coordinates(fc)), :]
    @test minimum(σ) >= -1.0e-9 && maximum(σ) <= top + 1.0e-9
    @test abs(mean(σ) - top / 2) < 4 * top / sqrt(12 * res.diagnostics.min_ess)
    @test var(σ) ≈ top^2 / 12 rtol = 0.1
end

@testset "sampling a prior" begin
    fc = chain_communities()
    prior = CommunityPrior(; interaction = Exponential(1.0))
    exact = chain_rejection(Xoshiro(34), 200_000)
    # importance weights: the exponential density of the losses and gains that vary
    w = [pdf(Exponential(1.0), λ12) * pdf(Exponential(1.0), γ21) * pdf(Exponential(1.0), λ23) for (_, λ12, γ21, λ23) in exact]
    weighted(v) = sum(w .* v) / sum(w)
    k(kind, i, j) = findfirst(c -> c == (kind = kind, i = i, j = j), coordinates(fc))
    for sampler in (nothing, CRHMC())
        res = sample_communities(Xoshiro(35), fc, prior; target_ess = 1000, sampler)
        @test res.converged
        λ23 = res.points[k(:loss, 2, 3), :]
        reference = weighted(last.(exact))
        spread = sqrt(weighted((last.(exact) .- reference) .^ 2))
        @test abs(mean(λ23) - reference) < 4 * spread / sqrt(res.diagnostics.min_ess)
    end
    # the prior pulls the loss λ₂₃ down from its uniform mean
    @test weighted(last.(exact)) < mean(last.(exact))
end

@testset "per capita invariances" begin
    # J = diag(x) A is similar to B = A diag(x), and a coordinate at A[i, j] has the per capita
    # rate θ x[j]. Mapping each sampled point to these rates gives a point of the matching polytope
    # with the same spectrum.
    web = generalized_cascade_network(Xoshiro(37), 12, 0.15; basal_paths = true)
    x = 2.0 .^ (trophic_levels(web) .- 1)
    per_capita(fc, θ) = [c.kind === :growth ? θ[k] : θ[k] * fc.x[c.j] for (k, c) in enumerate(coordinates(fc))]
    inside(fc, θ) = isfinite(logdensity(CommunityPrior(), fc, θ; atol = 1.0e-7))

    # a pyramid is x = 1 with the per capita efficiency band [0, 0.3 x_prey / x_pred] on each link
    pyramid = FeasibleCommunities(web, x; efficiency = (0.0, 0.3), per_capita_bounds = (0.0, 5.0))
    hi = [web[i, j] && i != j ? 0.3 * x[i] / x[j] : 0.0 for i in 1:12, j in 1:12]
    flat = FeasibleCommunities(web, ones(12); efficiency = (zero(hi), hi), bounds = (0.0, 5.0))
    res = sample_communities(Xoshiro(38), pyramid; target_ess = 100)
    for θ in eachcol(res.points[:, 1:25:end])
        B = per_capita(pyramid, θ)
        @test inside(flat, B)
        J = community_matrix(pyramid, θ)
        @test community_matrix(flat, B) ≈ Diagonal(1 ./ x) * J * Diagonal(x)
    end

    # under per capita bounds, rescaling x rescales the coordinates and leaves J unchanged
    scaled = FeasibleCommunities(web, 7 .* x; efficiency = (0.0, 0.3), per_capita_bounds = (0.0, 5.0))
    for θ in eachcol(res.points[:, 1:25:end])
        θ7 = [c.kind === :growth ? θ[k] : θ[k] / 7 for (k, c) in enumerate(coordinates(pyramid))]
        @test inside(scaled, θ7)
        @test community_matrix(scaled, θ7) ≈ community_matrix(pyramid, θ)
    end
end

@testset "growth rate draws" begin
    web = cascade_network(Xoshiro(36), 20, 0.15)
    r = draw_growth_rates(Xoshiro(1), web)
    basal = basal_species(web)
    @test all(r[basal] .> 0) && all(r[setdiff(1:20, basal)] .< 0)
    @test r == draw_growth_rates(Xoshiro(1), web)
    @test all(draw_growth_rates(Xoshiro(2), web; basal = Uniform(10, 11))[basal] .>= 10)

    # redraw until feasible; per capita bounds fit the rates at any scale of x, here biomasses of
    # at most 0.001
    web = generalized_cascade_network(Xoshiro(4), 15, 0.15; basal_paths = true)
    x = steady_state_biomass(web; efficiency = 0.5) ./ 1000
    found = feasible_growth_rates(Xoshiro(3), web, x; efficiency = (0.0, 0.5), per_capita_bounds = (0.0, 5.0))
    @test !isempty(found.communities) && found.draws >= 1
    @test found.r == feasible_growth_rates(Xoshiro(3), web, x; efficiency = (0.0, 0.5), per_capita_bounds = (0.0, 5.0)).r
    # magnitude bounds (0, 5) cap every per capita rate at 0.005, far below the basal rates
    @test_throws ErrorException feasible_growth_rates(Xoshiro(3), web, x; efficiency = (0.0, 0.5), bounds = (0.0, 5.0), maxiter = 20)
    # basal species that must lose biomass can never balance
    @test_throws ErrorException feasible_growth_rates(Xoshiro(5), web, x; basal = Uniform(-2, -1), maxiter = 5)
end
