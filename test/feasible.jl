using FoodWebs
using Test
using Random: Xoshiro
using Distributions: Exponential, Normal, Uniform, logpdf
using HiGHS: HiGHS
using JuMP: @constraint, @objective, @variable, MOI, Model, objective_value, optimize!, set_silent,
    termination_status, value
using LinearAlgebra: Diagonal

# The point of the polytope that maximizes c ⋅ θ, or nothing if it is empty or unbounded that way
function lp_point(P, c)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    θ = @variable(model, [1:length(P.lower)])
    for k in eachindex(P.lower)
        @constraint(model, θ[k] >= P.lower[k])
        isfinite(P.upper[k]) && @constraint(model, θ[k] <= P.upper[k])
    end
    @constraint(model, P.A * θ .<= P.b)
    @constraint(model, P.C * θ .== P.d)
    @objective(model, Max, c' * θ)
    optimize!(model)
    termination_status(model) == MOI.OPTIMAL || return nothing
    return value.(θ)
end

# A point away from the faces: the mean of vertices maximizing random objectives
function lp_interior(rng, P; vertices = 20)
    points = [lp_point(P, randn(rng, length(P.lower))) for _ in 1:vertices]
    any(isnothing, points) && return nothing
    return sum(points) / vertices
end

# A nonempty polytope is unbounded exactly when it has a nonzero direction of recession: a d with
# C d = 0 and A d ≤ 0 that is zero on coordinates bounded on both sides. Here every coordinate
# is bounded below, so d ≥ 0, and it is nonzero when sum(d) can reach 1.
function has_recession_direction(P)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    n = length(P.lower)
    all(isfinite, P.lower) || error("expected every coordinate bounded below")
    d = @variable(model, [1:n], lower_bound = 0)
    for k in 1:n
        isfinite(P.upper[k]) && @constraint(model, d[k] == 0)
    end
    @constraint(model, P.A * d .<= 0)
    @constraint(model, P.C * d .== 0)
    @constraint(model, sum(d) <= 1)
    @objective(model, Max, sum(d))
    optimize!(model)
    return objective_value(model) > 1.0e-7
end

isnonempty(P) = lp_point(P, zeros(length(P.lower))) !== nothing

@testset "coordinates and construction" begin
    # 2 eats 1, 3 eats 2
    chain = food_chain_network(3)
    fc = FeasibleCommunities(chain, [3.0, 2.0, 1.0]; rates = [1.0, -0.2, -0.1], self_regulation = :basal)
    @test coordinates(fc) == [
        (kind = :self, i = 1, j = 1), (kind = :loss, i = 1, j = 2), (kind = :gain, i = 2, j = 1),
        (kind = :loss, i = 2, j = 3), (kind = :gain, i = 3, j = 2),
    ]
    # rates as intervals add a growth coordinate per species, by role
    fc = FeasibleCommunities(chain, [3.0, 2.0, 1.0])
    P = halfspaces(fc)
    growth = [k for (k, c) in enumerate(coordinates(fc)) if c.kind === :growth]
    @test P.lower[growth] == [0, -1, -1] && P.upper[growth] == [5, 0, 0]
    # cannibals always self-regulate; two-way pairs get a loss and a gain at each entry
    fc = FeasibleCommunities(Bool[1 1; 1 0], [1.0, 1.0]; self_regulation = :none)
    kinds = [(c.kind, c.i, c.j) for c in coordinates(fc)]
    @test (:self, 1, 1) in kinds && !((:self, 2, 2) in kinds)
    @test Set(k for k in kinds if k[1] !== :growth && k[1] !== :self) ==
        Set([(:loss, 1, 2), (:gain, 1, 2), (:loss, 2, 1), (:gain, 2, 1)])
    @test sprint(show, fc) == "FeasibleCommunities on 2 species with 7 coordinates"

    @test_throws ArgumentError FeasibleCommunities(chain, [1.0, 0.0, 1.0])
    @test_throws DimensionMismatch FeasibleCommunities(chain, [1.0, 1.0])
    @test_throws ArgumentError FeasibleCommunities(chain, ones(3); efficiency = 1.5)
    @test_throws ArgumentError FeasibleCommunities(chain, ones(3); self_share = (0.6, 0.4))
    @test_throws ArgumentError FeasibleCommunities(chain, ones(3); self_regulation = :some)
    @test_throws ArgumentError FeasibleCommunities(chain, ones(3); bounds = (-1.0, 1.0))
    @test_throws ArgumentError FeasibleCommunities(chain, ones(3); per_capita_bounds = (2.0, 1.0))
    @test_throws ArgumentError FeasibleCommunities(chain, ones(3); per_capita_bounds = (-1.0, 1.0))
end

@testset "per capita bounds" begin
    # a coordinate at A[i, j] has the per capita rate θ x[j], so (lo, hi) bounds it by lo / x[j]
    # and hi / x[j]; growth rates keep their intervals
    chain = food_chain_network(3)
    x = [4.0, 2.0, 0.5]
    fc = FeasibleCommunities(chain, x; per_capita_bounds = (0.5, 2.0))
    P, P0 = halfspaces(fc), halfspaces(FeasibleCommunities(chain, x))
    for (k, c) in enumerate(coordinates(fc))
        if c.kind === :growth
            @test (P.lower[k], P.upper[k]) == (P0.lower[k], P0.upper[k])
        else
            @test P.lower[k] ≈ 0.5 / x[c.j] && P.upper[k] ≈ 2.0 / x[c.j]
        end
    end
    # with magnitude bounds, the tighter of the two applies
    P = halfspaces(FeasibleCommunities(chain, x; bounds = (0.0, 1.0), per_capita_bounds = (0.0, 2.0)))
    for (k, c) in enumerate(coordinates(fc))
        c.kind === :growth || @test P.upper[k] ≈ min(1.0, 2.0 / x[c.j])
    end
    # rescaling x scales the box by its inverse and leaves the equilibrium rows of the per capita
    # rates unchanged
    P1, P5 = halfspaces(fc), halfspaces(FeasibleCommunities(chain, 5 .* x; per_capita_bounds = (0.5, 2.0)))
    interaction = [c.kind !== :growth for c in coordinates(fc)]
    @test P5.upper[interaction] ≈ P1.upper[interaction] ./ 5
    @test P5.C[:, interaction] ≈ 5 .* P1.C[:, interaction]
    # finite per capita bounds always bound the polytope; a mutual-feeding loop with no
    # efficiency band is unbounded without them
    loop = Bool[0 1; 1 0]
    @test !isbounded(FeasibleCommunities(loop, ones(2); self_regulation = :none))
    @test isbounded(FeasibleCommunities(loop, ones(2); self_regulation = :none, per_capita_bounds = (0.0, 3.0)))
    @test !isbounded(FeasibleCommunities(loop, ones(2); self_regulation = :none, per_capita_bounds = (0.0, Inf)))
end

@testset "the PolytopeSamplers.jl example" begin
    # a 15-species generalized cascade web with fixed rates and no self-regulation for top
    # species: 78 coordinates in a 63-dimensional slice
    rng = Xoshiro(2026)
    web = generalized_cascade_network(rng, 15, 0.15; basal_paths = true)
    r = [i in basal_species(web) ? 5 * rand(rng) : -rand(rng) for i in 1:15]
    fc = FeasibleCommunities(web, rand(rng, 15) .+ 1; rates = r, self_regulation = :all_but_top)
    @test length(coordinates(fc)) == 78
    @test dimension(fc) == 63
end

@testset "points are equilibria" begin
    rng = Xoshiro(26)
    # (keywords, whether the rates are fixed)
    configurations = [
        ((; self_regulation = :all, bounds = (0.0, 10.0)), false),
        ((; self_regulation = :all_but_top, efficiency = (0.05, 0.3), bounds = (0.0, 10.0)), false),
        ((; efficiency = 0.2, self_share = (0.1, 0.6), bounds = (0.0, 10.0)), false),
        ((; efficiency = (0.0, 0.5), bounds = (0.0, 10.0)), true),
    ]
    tested = 0
    for gen in (cascade_network, niche_network, nested_hierarchy_network), (config, fixed) in configurations
        web = gen(rng, 20, 0.15; basal_paths = true)
        x = steady_state_biomass(web; efficiency = 0.2) .* exp.(0.3 .* randn(rng, 20))
        kw = fixed ? (; rates = [i in basal_species(web) ? 2.0 : -0.1 for i in 1:20], config...) : config
        fc = FeasibleCommunities(web, x; kw...)
        P = halfspaces(fc)
        θ = lp_interior(rng, P)
        θ === nothing && continue
        tested += 1
        A, r = interaction_matrix(fc, θ), growth_rates(fc, θ)
        @test r .+ A * x ≈ zeros(20) atol = 1.0e-7
        @test community_matrix(fc, θ) ≈ Diagonal(x) * A
        # one-way links keep their signs; cannibals and self-regulated species have A_ii ≤ 0
        for i in 1:20, j in 1:20
            (i != j && web[i, j] && !web[j, i]) || continue
            @test A[i, j] <= 0 <= A[j, i]
        end
        @test all(A[i, i] <= 1.0e-12 for i in 1:20)
        # efficiency bands on every link
        if haskey(kw, :efficiency)
            lo, hi = kw.efficiency isa Tuple ? kw.efficiency : (kw.efficiency, kw.efficiency)
            for (k, c) in enumerate(coordinates(fc))
                c.kind === :loss || continue
                g = findfirst(d -> d.kind === :gain && d.i == c.j && d.j == c.i, coordinates(fc))
                @test lo * θ[k] - 1.0e-7 <= θ[g] <= hi * θ[k] + 1.0e-7
            end
            # the budget: Σ (λ - γ) x_i x_j + Σ σ x² = x ⋅ r
            budget = sum(
                c.kind === :loss ? θ[k] * x[c.i] * x[c.j] : c.kind === :gain ? -θ[k] * x[c.i] * x[c.j] :
                    c.kind === :self ? θ[k] * x[c.i]^2 : 0.0 for (k, c) in enumerate(coordinates(fc))
            )
            @test budget ≈ sum(x .* r) atol = 1.0e-6
        end
        # fluxes
        F = fluxes(fc, θ)
        @test all(F[c.i, c.j] ≈ θ[k] * x[c.i] * x[c.j] for (k, c) in enumerate(coordinates(fc)) if c.kind === :loss)
    end
    @test tested >= 8
end

@testset "two species eating each other and cannibalism" begin
    # 1 basal; 2 and 3 eat 1 and each other; 3 is a cannibal
    web = Bool[0 1 1; 0 0 1; 0 1 1]
    fc = FeasibleCommunities(web, [4.0, 1.0, 1.0]; self_regulation = :basal, efficiency = (0.0, 0.5), bounds = (0.0, 20.0))
    @test isempty(diagnose(fc))
    P = halfspaces(fc)
    rng = Xoshiro(27)
    signs = Set{Int}()
    for _ in 1:30
        θ = lp_point(P, randn(rng, length(P.lower)))
        θ === nothing && continue
        A = interaction_matrix(fc, θ)
        @test growth_rates(fc, θ) .+ A * [4.0, 1.0, 1.0] ≈ zeros(3) atol = 1.0e-7
        @test A[3, 3] <= 0
        push!(signs, sign(round(A[2, 3], digits = 9)))
    end
    # the net effect within the pair can take either sign
    @test 1 in signs && -1 in signs
end

@testset "two-species consumer-resource segment" begin
    # 2 eats 1; fixed rates; only the resource self-regulates. Row 2: x₁ γ + r₂ = 0, so
    # γ = -r₂ / x₁. Row 1: -x₁ σ - x₂ λ + r₁ = 0, so λ = (r₁ - x₁ σ) / x₂, and γ ≤ e λ bounds σ
    # by (r₁ - x₂ γ / e) / x₁.
    x, r, e = [2.0, 0.5], [3.0, -0.4], 0.25
    fc = FeasibleCommunities(Bool[0 1; 0 0], x; rates = r, self_regulation = :basal, efficiency = (0.0, e))
    @test dimension(fc) == 1
    P = halfspaces(fc)
    σ = findfirst(c -> c.kind === :self, coordinates(fc))
    γ = findfirst(c -> c.kind === :gain, coordinates(fc))
    top = lp_point(P, [k == σ ? 1.0 : 0.0 for k in eachindex(P.lower)])
    @test top[γ] ≈ -r[2] / x[1]
    @test top[σ] ≈ (r[1] - x[2] * (-r[2] / x[1]) / e) / x[1]
    @test lp_point(P, [k == σ ? -1.0 : 0.0 for k in eachindex(P.lower)])[σ] ≈ 0 atol = 1.0e-9
end

@testset "boundedness against linear programs" begin
    rng = Xoshiro(28)
    compared = 0
    outcomes = Bool[]
    for _ in 1:60
        gen = rand(rng, (cascade_network, niche_network, nested_hierarchy_network))
        web = gen(rng, 8, 0.2)
        efficiency = rand(rng, (nothing, 1.0, (0.1, 1.0), (0.1, 0.6)))
        self_share = rand(rng, (nothing, (0.2, 0.9), (0.0, 1.0)))
        self_regulation = rand(rng, (:all, :all_but_top, :basal))
        fc = FeasibleCommunities(web, rand(rng, 8) .+ 0.5; efficiency, self_share, self_regulation)
        P = halfspaces(fc)
        isnonempty(P) || continue
        compared += 1
        bounded = isbounded(fc)
        @test bounded == !has_recession_direction(P)
        push!(outcomes, bounded)
    end
    @test compared >= 30
    # both outcomes occur, so the comparison checks each rule
    @test true in outcomes && false in outcomes
    # finite bounds always bound it
    @test isbounded(FeasibleCommunities(niche_network(rng, 10, 0.2), ones(10); bounds = (0.0, 3.0)))
end

@testset "structural singularity and diagnosis" begin
    # 2 and 3 both eat only 1 and nothing eats them: their rows share one column
    web = Bool[0 1 1; 0 0 0; 0 0 0]
    @test structurally_singular(FeasibleCommunities(web, ones(3); self_regulation = :basal))
    @test !structurally_singular(FeasibleCommunities(web, ones(3); self_regulation = :all))

    # a basal species with a fixed negative rate cannot pay for its losses; a top species
    # without self-regulation cannot have a positive rate
    fc = FeasibleCommunities(food_chain_network(2), ones(2); rates = [-1.0, 0.5], self_regulation = :none)
    problems = diagnose(fc)
    @test length(problems) == 2
    @test !isnonempty(halfspaces(fc))
    # the budget under efficiency bands
    fc = FeasibleCommunities(food_chain_network(2), ones(2); rates = [0.1, -1.0], efficiency = 0.5)
    @test any(occursin("budget", p) for p in diagnose(fc))
    @test !isnonempty(halfspaces(fc))
    # dense May webs often have no basal species
    rng = Xoshiro(29)
    for gen in (cascade_network, niche_network, nested_hierarchy_network, minimum_potential_niche_network)
        web = gen === minimum_potential_niche_network ? gen(rng, 25, 0.12, 0.3) : gen(rng, 25, 0.12)
        fc = FeasibleCommunities(web, steady_state_biomass(remove_cannibalism(web)[1]; efficiency = 0.2))
        @test isempty(diagnose(fc))
    end
end

@testset "every bundled web" begin
    for id in filter(id -> startswith(id, "cheddar/") || startswith(id, "williams_martinez_2000/"), foodweb_ids())
        fw = foodweb(id)
        fc = FeasibleCommunities(fw, ones(size(fw, 1)))
        @test length(coordinates(fc)) == 2 * count(fw[i, j] for i in axes(fw, 1), j in axes(fw, 2) if i != j) +
            size(fw, 1) + size(fw, 1)
        @test dimension(fc) >= 0
    end
end

@testset "priors" begin
    fc = FeasibleCommunities(food_chain_network(3), [3.0, 2.0, 1.0]; efficiency = (0.0, 0.5), bounds = (0.0, 5.0))
    P = halfspaces(fc)
    θ = lp_interior(Xoshiro(30), P)
    @test logdensity(CommunityPrior(), fc, θ) == 0
    prior = CommunityPrior(; interaction = Exponential(1.0), growth = (basal = Uniform(0, 5), consumer = Normal(-0.5, 0.2)))
    expected = sum(
        c.kind in (:loss, :gain) ? logpdf(Exponential(1.0), θ[k]) :
            c.kind === :growth ? (c.i == 1 ? logpdf(Uniform(0, 5), θ[k]) : logpdf(Normal(-0.5, 0.2), θ[k])) : 0.0
            for (k, c) in enumerate(coordinates(fc))
    )
    @test logdensity(prior, fc, θ) ≈ expected
    # outside the polytope
    @test logdensity(prior, fc, θ .+ 1) == -Inf
    @test CommunityPrior(; interaction = Exponential(1.0), gain = Exponential(2.0)).gain == Exponential(2.0)
end
