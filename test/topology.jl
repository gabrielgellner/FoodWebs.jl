using FoodWebs
using Test
using Random: Xoshiro
using LinearAlgebra: diag, istriu
using Graphs: DiGraph, is_connected, is_cyclic
using Statistics: mean, std

# generators taking (rng, S, C; kwargs...)
const GENERATORS = Dict(
    "may" => may_network,
    "cascade" => cascade_network,
    "generalized cascade" => generalized_cascade_network,
    "niche" => niche_network,
    "generalized niche" => (rng, S, C; kw...) -> generalized_niche_network(rng, S, C, 0.5; kw...),
)

# species are sorted by niche value, so an interval diet is a contiguous run of row indices
is_interval(idx) = isempty(idx) || idx == first(idx):last(idx)

@testset "PredationMatrix" begin
    web = PredationMatrix([0 1; 0 0])
    @test size(web) == (2, 2)
    @test web[1, 2] == 1
    @test web isa AbstractMatrix{Int}
    @test_throws DimensionMismatch PredationMatrix([0 1 0; 0 0 1])
    @test_throws ArgumentError PredationMatrix([0 2; 0 0])
    @test_throws ArgumentError PredationMatrix([0 -1; 0 0])
end

@testset "food chain" begin
    web = food_chain_network(5)
    @test connectance(web) == 4 / 25
    @test basal_species(web) == [1]
    @test top_species(web) == [5]
    @test trophic_levels(web) ≈ 1:5
    @test food_chain_network(1) == fill(0, 1, 1)
    @test_throws ArgumentError food_chain_network(0)
end

@testset "common properties: $name" for (name, gen) in GENERATORS
    rng = Xoshiro(1)
    for (S, C) in ((10, 0.2), (40, 0.1), (100, 0.15))
        web = gen(rng, S, C)
        @test web isa PredationMatrix
        @test size(web) == (S, S)
        @test all(x -> x == 0 || x == 1, web)
        @test is_connected(DiGraph(web.links))
        @test abs(connectance(web) - C) <= 0.03C
        # May networks are random graphs, not food web models, so they need not have basal
        # species
        name == "may" || @test !isempty(basal_species(web))
        # cascade webs are acyclic, so every species has a defined trophic level. Niche webs can
        # contain loops of species that only eat each other, cut off from any basal species.
        occursin("cascade", name) && @test all(isfinite, trophic_levels(web))
    end
    # acceptance tolerance is respected
    @test abs(connectance(gen(rng, 60, 0.1; tol = 0.005)) - 0.1) <= 0.005 * 0.1
    # reproducible from a seeded rng
    @test gen(Xoshiro(42), 30, 0.15) == gen(Xoshiro(42), 30, 0.15)
    # impossible targets fail loudly instead of looping forever: 0.101 * 20^2 is not an integer
    @test_throws ErrorException gen(rng, 20, 0.101; tol = 0.0, maxiter = 5)
    @test_throws ArgumentError gen(rng, 1, 0.1)
    @test_throws ArgumentError gen(rng, 20, 0.0)
    @test_throws ArgumentError gen(rng, 20, -0.1)
    @test_throws ArgumentError gen(rng, 20, 0.6)
end

@testset "unbiased connectance: $name" for (name, gen) in GENERATORS
    # with acceptance off the mean connectance over many webs should hit the target, which
    # checks each model's parameterization rather than the rejection step. The niche models
    # lose links at large C as ranges are truncated at 1, so stay at small C.
    rng = Xoshiro(2)
    for C in (0.05, 0.1, 0.2)
        mC = mean(connectance(gen(rng, 100, C; tol = Inf)) for _ in 1:200)
        @test isapprox(mC, C, rtol = 0.05)
    end
end

@testset "may network" begin
    rng = Xoshiro(3)
    for _ in 1:20
        web = may_network(rng, 40, 0.2)
        @test all(iszero, diag(web))
        # no mutual predation
        @test all(web[i, j] + web[j, i] <= 1 for i in 1:40 for j in 1:40)
    end
    # link direction is unbiased, so roughly half the links are above the diagonal
    webs = [may_network(rng, 100, 0.1) for _ in 1:20]
    frac_upper = mean(sum(w[i, j] for i in 1:99 for j in (i + 1):100) / sum(w) for w in webs)
    @test isapprox(frac_upper, 0.5, atol = 0.03)
    @test_throws ArgumentError may_network(rng, 10, 0.46)
end

@testset "cascade networks" begin
    rng = Xoshiro(4)
    for gen in (cascade_network, generalized_cascade_network), _ in 1:20
        web = gen(rng, 40, 0.2)
        # species only eat lower ranked species, so the web is acyclic with no cannibals
        @test istriu(web.links, 1)
        @test !is_cyclic(DiGraph(web.links))
        @test 1 in basal_species(web)
    end
    # the maximum possible connectance, S(S-1)/2 links, is allowed
    @test connectance(cascade_network(rng, 10, 0.45)) == 0.45
    @test_throws ArgumentError cascade_network(rng, 10, 0.46)

    # one Beta(1, β) draw per predator gives diet fractions with far more spread than the
    # plain cascade's binomial ones
    S, C = 200, 0.2
    frac(web) = [sum(web[1:(j - 1), j]) / (j - 1) for j in (S ÷ 2):S]
    @test std(frac(generalized_cascade_network(rng, S, C))) > 3 * std(frac(cascade_network(rng, S, C)))
end

@testset "niche network" begin
    rng = Xoshiro(5)
    for _ in 1:50
        web = niche_network(rng, 50, 0.15)
        @test all(is_interval(findall(==(1), web[:, i])) for i in 1:50)
        # the species with the lowest niche value has range 0
        @test 1 in basal_species(web)
    end
end

@testset "basal paths" begin
    rng = Xoshiro(13)
    # roughly 1-2% of niche webs contain loops of species cut off from any basal species
    @test any(!isempty(disconnected_from_basal(niche_network(rng, 40, 0.25))) for _ in 1:500)
    # random link directions mean dense May networks almost never have a basal species
    @test_throws ErrorException may_network(rng, 40, 0.25; basal_paths = true, maxiter = 100)
    for (name, gen) in GENERATORS, _ in 1:100
        name == "may" && continue
        web = gen(rng, 40, 0.25; basal_paths = true)
        @test isempty(disconnected_from_basal(web))
        @test all(isfinite, trophic_levels(web))
    end
    @test niche_network(Xoshiro(1), 30, 0.15; basal_paths = true) ==
          niche_network(Xoshiro(1), 30, 0.15; basal_paths = true)
end

@testset "generalized niche network" begin
    rng = Xoshiro(6)
    for diet in (0.0, 0.25, 0.5, 0.75, 1.0), _ in 1:10
        web = generalized_niche_network(rng, 50, 0.15, diet)
        @test all(iszero, diag(web))
        @test 1 in basal_species(web)
    end
    # with diet = 1 the diets are niche intervals, except that excluding the cannibalistic
    # link can split the interval in two around the species itself
    for _ in 1:20
        web = generalized_niche_network(rng, 50, 0.15, 1.0)
        @test all(is_interval(sort(union(findall(==(1), web[:, i]), i))) ||
                  is_interval(findall(==(1), web[:, i])) for i in 1:50)
    end
    # lowering diet contiguity makes diets less interval-like
    n_interval(web) = count(i -> is_interval(findall(==(1), web[:, i])), 1:50)
    @test mean(n_interval(generalized_niche_network(rng, 50, 0.15, 0.0)) for _ in 1:20) <
          mean(n_interval(generalized_niche_network(rng, 50, 0.15, 1.0)) for _ in 1:20)
    @test_throws ArgumentError generalized_niche_network(rng, 20, 0.1, 1.5)
    @test_throws ArgumentError generalized_niche_network(rng, 20, 0.1, -0.1)
end
