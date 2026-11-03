using FoodWebs
using Test
using LinearAlgebra: Diagonal, norm
using Random: Xoshiro
using OrdinaryDiffEqRosenbrock: Rodas5P
using PolytopeSamplers: PolytopeSamplers
using SciMLBase: SciMLBase

# A prey with self-regulation 0.5 at x = 2, eaten by a predator at x = 1 with loss rate 0.4 and
# efficiency 0.5, and the point θ of these values
function predator_prey()
    fc = FeasibleCommunities(Bool[0 1; 0 0], [2.0, 1.0]; rates = [1.4, -0.4], self_regulation = :basal)
    return fc, [0.5, 0.4, 0.2]
end

# the right-hand side x .* g of a type II model
type2_rates(m, x) = x .* FoodWebs.type2_per_capita!(zeros(length(x)), x, m, trues(length(x)))

@testset "type II models by hand" begin
    fc, θ = predator_prey()
    P = halfspaces(fc)
    @test P.C * θ ≈ P.d
    N, Pr = fc.x
    # τ = λ h N = 0.8 h
    h = 0.5
    m = type2_model(fc, θ; handling = [0 h; 0 0])
    @test handling_fractions(m) ≈ [0, 0.4]
    a = m.attack[1, 2]
    @test a ≈ 0.4 / 0.6
    @test m.efficiency[1, 2] ≈ 0.5 && m.self_regulation == [0.5, 0.0] && m.r == [1.4, -0.4]
    # the Rosenzweig–MacArthur Jacobian with prey self-regulation σ = 0.5
    D = 1 + a * h * N
    J = [
        N * (-0.5 + a^2 * h * Pr / D^2) -a * N / D
        0.5 * a * Pr / D^2 0
    ]
    @test community_matrix(m) ≈ J
    @test type2_rates(m, fc.x) ≈ zeros(2) atol = 1.0e-14
    # handling times of 0 give the Lotka–Volterra model
    m0 = type2_model(fc, θ; handling = zeros(2, 2))
    @test m0.attack[1, 2] == 0.4 && community_matrix(m0) ≈ community_matrix(fc, θ)
    # the same τ from a handling fraction or a maximum consumption rate 1 / h
    @test community_matrix(type2_model(fc, θ; handling_fraction = 0.4)) ≈ J
    @test community_matrix(type2_model(fc, θ; max_consumption = 1 / h)) ≈ J

    # with one prey the refit is exact: σ rises by λ² h P and e to e / (1 - τ)
    mm = type2_model(fc, θ; handling = [0 h; 0 0], match_jacobian = true)
    @test community_matrix(mm) ≈ community_matrix(fc, θ)
    @test mm.self_regulation[1] ≈ 0.5 + 0.4^2 * h * Pr
    @test mm.efficiency[1, 2] ≈ 0.5 / 0.6
    @test type2_rates(mm, fc.x) ≈ zeros(2) atol = 1.0e-14
    @test mm.attack == m.attack && mm.handling == m.handling
end

@testset "type II models of sampled communities" begin
    for seed in 1:3
        web = niche_network(Xoshiro(seed + 600), 12, 0.15; basal_paths = true)
        x = exp10.(2 .* rand(Xoshiro(seed + 700), 12) .- 1)
        fc = FeasibleCommunities(web, x; efficiency = (0.0, 0.3))
        points = sample_communities(Xoshiro(seed), fc; target_ess = 50).points
        eats = fc.eats
        for θ in eachcol(points[:, 1:3])
            τ = 0.2 .+ 0.5 .* rand(Xoshiro(seed + 800), 12)
            m = type2_model(fc, θ; handling_fraction = τ)
            consumers = [any(eats[:, j]) for j in 1:12]
            @test handling_fractions(m)[consumers] ≈ τ[consumers]
            # the same equilibrium and fluxes as the Lotka–Volterra point
            @test maximum(abs, type2_rates(m, x)) < 1.0e-10 * maximum(abs, x)
            D = 1 ./ (1 .- handling_fractions(m))
            @test m.attack .* x .* x' ./ D' ≈ fluxes(fc, θ)
            # the Jacobian against central differences of x .* g
            J = community_matrix(m)
            δ = 1.0e-6 .* x
            Jfd = reduce(hcat, (type2_rates(m, x .+ δ[k] .* (1:12 .== k)) .- type2_rates(m, x .- δ[k] .* (1:12 .== k))) ./ (2δ[k]) for k in 1:12)
            @test J ≈ Jfd rtol = 1.0e-6
            # the change from the Lotka–Volterra Jacobian: + F_ij λ_kj h_kj for prey i and k of a
            # shared predator j, and - x_j λ_kj h_kj G_j for predator j and its prey k
            λ = fluxes(fc, θ) ./ (x .* x')
            F = fluxes(fc, θ)
            G = [sum(m.efficiency[i, j] * λ[i, j] * x[i] for i in 1:12) for j in 1:12]
            Δprey = zeros(12, 12)
            shared = falses(12, 12)
            for j in 1:12, i in 1:12, k in 1:12
                eats[i, j] && eats[k, j] || continue
                Δprey[i, k] += F[i, j] * λ[k, j] * m.handling[k, j]
                shared[i, k] = true
            end
            Δ = copy(Δprey)
            for j in 1:12, k in 1:12
                eats[k, j] && (Δ[j, k] -= x[j] * λ[k, j] * m.handling[k, j] * G[j])
            end
            JLV = community_matrix(fc, θ)
            @test J - JLV ≈ Δ atol = 1.0e-10 * maximum(abs, JLV)

            # the refit leaves only the effects between distinct prey of a shared predator, and
            # only raises self-regulation and efficiencies
            mm = type2_model(fc, θ; handling_fraction = τ, match_jacobian = true)
            @test maximum(abs, type2_rates(mm, x)) < 1.0e-10 * maximum(abs, x)
            @test all(mm.self_regulation .>= m.self_regulation) && all(mm.efficiency .>= m.efficiency)
            R = community_matrix(mm) - JLV
            kept = [i != k && shared[i, k] for i in 1:12, k in 1:12]
            @test maximum(abs, R[.!kept]) < 1.0e-9 * maximum(abs, JLV)
            @test R[kept] ≈ Δprey[kept]
            @test norm(R) <= norm(J - JLV)
        end
    end
end

@testset "type II dynamics" begin
    fc, θ = predator_prey()
    m = type2_model(fc, θ; handling = [0 0.5; 0 0])
    # the equilibrium stays put, and small pulses follow the linearization
    s = type2_dynamics(Rodas5P(), m, fc.x, (0.0, 50.0); abstol = 1.0e-12, reltol = 1.0e-10)
    @test s.retcode == SciMLBase.ReturnCode.Success
    @test s.x[:, end] ≈ fc.x
    δ = [0.3, -0.2]
    err(ε) = (
        p = pulse_response(Rodas5P(), m, ε .* δ, range(0, 20; length = 41); abstol = 1.0e-14, reltol = 1.0e-12);
        maximum(abs.(p.deviation .- p.linear))
    )
    @test 5 < err(1.0e-3) / err(1.0e-4) < 20
    # a predator without its prey declines and goes extinct
    s = type2_dynamics(Rodas5P(), m, [0.0, 1.0], (0.0, 200.0); extinction_threshold = 1.0e-6)
    @test s.survivors == [false, false]
    @test s.extinction_times[2] ≈ log(1.0e6) / 0.4 rtol = 1.0e-3
    # with handling times of 0 the dynamics are Lotka–Volterra
    m0 = type2_model(fc, θ; handling = zeros(2, 2))
    c = (A = interaction_matrix(fc, θ), r = growth_rates(fc, θ))
    x0 = [1.0, 2.0]
    @test type2_dynamics(Rodas5P(), m0, x0, (0.0, 10.0); saveat = 0:1:10, abstol = 1.0e-12, reltol = 1.0e-10).x ≈
        lotka_volterra(Rodas5P(), c, x0, (0.0, 10.0); saveat = 0:1:10, abstol = 1.0e-12, reltol = 1.0e-10).x rtol = 1.0e-8
end

@testset "type II argument errors" begin
    fc, θ = predator_prey()
    @test_throws ArgumentError type2_model(fc, θ)
    @test_throws ArgumentError type2_model(fc, θ; handling = zeros(2, 2), handling_fraction = 0.1)
    @test_throws DimensionMismatch type2_model(fc, θ[1:2]; handling_fraction = 0.1)
    @test_throws DimensionMismatch type2_model(fc, θ; handling = zeros(3, 3))
    @test_throws ArgumentError type2_model(fc, θ; handling = [0 -1; 0 0])
    @test_throws ArgumentError type2_model(fc, θ; handling_fraction = 1.0)
    @test_throws DimensionMismatch type2_model(fc, θ; handling_fraction = [0.1])
    # consumption at the equilibrium, 0.8, above the maximum of 0.5
    @test_throws ArgumentError type2_model(fc, θ; max_consumption = 0.5)
    @test_throws ArgumentError type2_model(fc, θ; max_consumption = 0.0)
    # a gain with no loss carries no flux
    fc2 = FeasibleCommunities(Bool[0 1; 0 0], [2.0, 1.0]; rates = [1.0, -0.2], self_regulation = :basal)
    @test_throws ArgumentError type2_model(fc2, [0.5, 0.0, 0.2]; handling_fraction = 0.1)
    # at τ = 0.6 the refit efficiency would be 0.5 / 0.4 > 1
    @test_throws ArgumentError type2_model(fc, θ; handling_fraction = 0.6, match_jacobian = true)
end
