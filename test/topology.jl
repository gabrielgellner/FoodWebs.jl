using FoodWebs
using Test
using Random: Xoshiro
using LinearAlgebra: diag, istriu
using Graphs: DiGraph, is_connected, is_cyclic
using Statistics: mean, quantile, std

# generators taking (rng, S, C; kwargs...)
const GENERATORS = Dict(
    "may" => may_network,
    "cascade" => cascade_network,
    "generalized cascade" => generalized_cascade_network,
    "niche" => niche_network,
    "generalized niche" => (rng, S, C; kw...) -> generalized_niche_network(rng, S, C, 0.5; kw...),
    "nested hierarchy" => nested_hierarchy_network,
    "minimum potential niche" => (rng, S, C; kw...) -> minimum_potential_niche_network(rng, S, C, 0.3; kw...),
    "probabilistic niche" => probabilistic_niche_network,
)

# species are sorted by niche value, so an interval diet is a contiguous run of row indices
is_interval(idx) = isempty(idx) || idx == first(idx):last(idx)

@testset "PredationMatrix" begin
    web = PredationMatrix([0 1; 0 0])
    @test size(web) == (2, 2)
    @test web[1, 2] == 1
    @test web isa AbstractMatrix{Bool}
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
        @test abs(count(web) - C * S^2) <= 0.03 * C * S^2
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
    # 0.1 * 8^2 = 6.4 links allows at most 6, too few to connect 8 species
    @test_throws ArgumentError gen(rng, 8, 0.1)
    @test_throws ArgumentError gen(rng, 1, 0.1)
    @test_throws ArgumentError gen(rng, 20, 0.0)
    @test_throws ArgumentError gen(rng, 20, -0.1)
    @test_throws ArgumentError gen(rng, 20, 0.6)
end

@testset "small webs: $name" for (name, gen) in GENERATORS
    # no whole number of links lies within 3% of C S^2 here, so the generators accept the
    # nearest one
    rng = Xoshiro(7)
    for (S, C) in ((6, 0.18), (7, 0.15), (9, 0.12)), _ in 1:10
        @test count(gen(rng, S, C)) == round(Int, C * S^2)
    end
    # with tol = 0 every web has the nearest whole number of links
    @test count(gen(rng, 20, 0.101; tol = 0.0)) == 40
end

@testset "acceptance unchanged for large webs: $name" for (name, gen) in GENERATORS
    # when tol C S^2 >= 1/2 the window is the relative one, so a seeded generator returns the
    # first draw that is connected with a link count within tol of C S^2
    S, C = 40, 0.1
    for seed in 1:5
        rng = Xoshiro(seed)
        expected = nothing
        while expected === nothing
            web = try
                gen(rng, S, C; tol = Inf, maxiter = 1)
            catch err
                err isa ErrorException || rethrow()
                continue
            end
            abs(count(web) - C * S^2) <= 0.03 * C * S^2 && (expected = web)
        end
        @test gen(Xoshiro(seed), S, C) == expected
    end
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
        @test all(
            is_interval(sort(union(findall(==(1), web[:, i]), i))) ||
                is_interval(findall(==(1), web[:, i])) for i in 1:50
        )
    end
    # lowering diet contiguity makes diets less interval-like
    n_interval(web) = count(i -> is_interval(findall(==(1), web[:, i])), 1:50)
    @test mean(n_interval(generalized_niche_network(rng, 50, 0.15, 0.0)) for _ in 1:20) <
        mean(n_interval(generalized_niche_network(rng, 50, 0.15, 1.0)) for _ in 1:20)
    @test_throws ArgumentError generalized_niche_network(rng, 20, 0.1, 1.5)
    @test_throws ArgumentError generalized_niche_network(rng, 20, 0.1, -0.1)
end

# mean Jaccard similarity of the prey sets of consumer pairs that share at least one prey
function prey_overlap(web)
    S = size(web, 1)
    J = Float64[]
    for a in 1:S, b in (a + 1):S
        shared = count(web[:, a] .& web[:, b])
        shared > 0 && push!(J, shared / count(web[:, a] .| web[:, b]))
    end
    return mean(J)
end

@testset "nested hierarchy network" begin
    rng = Xoshiro(10)
    # diet sizes are drawn as in the niche model
    diets(web) = [count(web[:, j]) for j in axes(web, 2)]
    nhm = reduce(vcat, diets(nested_hierarchy_network(rng, 100, 0.1; tol = Inf)) for _ in 1:300)
    niche = reduce(vcat, diets(niche_network(rng, 100, 0.1; tol = Inf)) for _ in 1:300)
    @test isapprox(mean(nhm), mean(niche), rtol = 0.05)
    @test isapprox(std(nhm), std(niche), rtol = 0.05)
    @test all(abs(quantile(nhm, q) - quantile(niche, q)) <= 2 for q in (0.25, 0.5, 0.75, 0.9))

    # every consumer's first prey has a lower niche value, and diets need not be intervals
    webs = [nested_hierarchy_network(rng, 100, 0.08) for _ in 1:50]
    @test all(any(web[1:(j - 1), j]) for web in webs for j in 2:100 if any(web[:, j]))
    @test any(!is_interval(findall(web[:, j])) for web in webs for j in 1:100)
    # Stouffer et al. (2005, Fig. 4) report 92% of links going to species with lower niche
    # values at S = 100, C = 0.08, and very little cannibalism
    L = sum(sum, webs)
    @test 0.88 <= sum(sum(web[i, j] for j in 1:100 for i in 1:(j - 1)) for web in webs) / L <= 0.95
    @test sum(sum(diag(web)) for web in webs) / L < 0.01

    # prey picked from the consumer group's diets make consumers that share prey more alike
    # than in the generalized cascade model, which has the same diet sizes and random prey
    @test mean(prey_overlap(nested_hierarchy_network(rng, 60, 0.1)) for _ in 1:30) >
        mean(prey_overlap(generalized_cascade_network(rng, 60, 0.1)) for _ in 1:30) + 0.02
end

@testset "minimum potential niche network" begin
    # with f = 0 this is the niche model, draw for draw
    for seed in 1:5
        @test minimum_potential_niche_network(Xoshiro(seed), 40, 0.1, 0.0) ==
            niche_network(Xoshiro(seed), 40, 0.1)
    end

    # with f > 0 the web keeps both ends of every niche interval of the same draw and a
    # fraction 1 - f of the species inside it. Compare first draws, before any resampling.
    rng = Xoshiro(11)
    S, C, f = 60, 0.1, 0.4
    Cp = FoodWebs.potential_connectance(S, C, f)
    subset, ends_kept, kept, inside = true, true, 0, 0
    for seed in 1:300
        mpn = try
            minimum_potential_niche_network(Xoshiro(seed), S, C, f; tol = Inf, maxiter = 1)
        catch err
            err isa ErrorException || rethrow()
            continue
        end
        niche = niche_network(Xoshiro(seed), S, Cp; tol = Inf, maxiter = 1)
        subset &= all(mpn .<= niche)
        for j in 1:S
            prey = findall(niche[:, j])
            isempty(prey) && continue
            ends_kept &= mpn[first(prey), j] && mpn[last(prey), j]
            kept += count(mpn[prey[2:(end - 1)], j])
            inside += max(length(prey) - 2, 0)
        end
    end
    @test subset
    @test ends_kept
    @test isapprox(kept / inside, 1 - f, atol = 0.02)

    # forbidden links break up diet intervals
    webs = [minimum_potential_niche_network(rng, 60, 0.1, 0.5) for _ in 1:20]
    @test any(!is_interval(findall(web[:, j])) for web in webs for j in 1:60)

    # the potential connectance accounts for the kept ends of each interval, so mean
    # connectance stays on target for several f. Drawing intervals at C / (1 - f) would overshoot
    # by 2-28% here, and at C = 0.1, f = 0.8 it would need connectance 0.5.
    for f in (0.2, 0.5, 0.8), C in (0.05, 0.1)
        mC = mean(connectance(minimum_potential_niche_network(rng, 100, C, f; tol = Inf)) for _ in 1:200)
        @test isapprox(mC, C, rtol = 0.05)
    end

    @test_throws ArgumentError minimum_potential_niche_network(rng, 20, 0.1, 1.0)
    @test_throws ArgumentError minimum_potential_niche_network(rng, 20, 0.1, -0.1)
    # 0.3 links per species pair cannot survive forbidding 80% of the interior links
    @test_throws ArgumentError minimum_potential_niche_network(rng, 50, 0.3, 0.8)
end

@testset "probabilistic niche probabilities" begin
    η = [0.1, 0.3, 0.5, 0.7, 0.9]
    c = [0.0, 0.3, 0.6, 0.5, 0.7]
    r = [0.0, 0.2, 0.2, 0.4, 0.8]
    P = probabilistic_niche_probabilities(η, c, r; a = 0.8)
    @test size(P) == (5, 5)
    @test all(0 .<= P .<= 0.8)
    # a species with no range eats nothing
    @test all(iszero, P[:, 1])
    # the probability peaks at a on the feeding center: species 2 is centered on η[2], and so
    # is species 4 on η[3]
    @test P[2, 2] == 0.8
    @test P[3, 4] == 0.8
    # one half range from the center the Gaussian kernel is a e^-1, and with p = 1 also a e^-1
    @test P[2, 4] ≈ 0.8 * exp(-1)
    @test probabilistic_niche_probabilities(η, c, r; a = 0.8, p = 1)[2, 4] ≈ 0.8 * exp(-1)
    # it falls off with distance from the center
    @test P[3, 5] > P[2, 5] > P[1, 5]

    # as p grows the kernel approaches the niche model's interval of width r centered on c
    Pbox = probabilistic_niche_probabilities(η, c, r; p = 1000)
    for i in 2:5, j in 1:5
        abs(abs(η[j] - c[i]) - r[i] / 2) > 0.01 || continue
        @test (Pbox[j, i] > 0.5) == (abs(η[j] - c[i]) < r[i] / 2)
    end

    @test_throws DimensionMismatch probabilistic_niche_probabilities(η, c, r[1:4])
    @test_throws ArgumentError probabilistic_niche_probabilities(η, c, -r)
    @test_throws ArgumentError probabilistic_niche_probabilities(η, c, r; a = 1.5)
    @test_throws ArgumentError probabilistic_niche_probabilities(η, c, r; p = 0)
end

@testset "probabilistic niche network" begin
    rng = Xoshiro(12)
    # webs from given parameters have sum(P) links on average
    S = 50
    η = sort(rand(rng, S))
    r = η .* rand(rng, S) .* 0.3
    c = [rand(rng) * (η[i] - r[i] / 2) + r[i] / 2 for i in 1:S]
    P = probabilistic_niche_probabilities(η, c, r)
    webs = [probabilistic_niche_network(rng, η, c, r) for _ in 1:2000]
    @test all(web isa PredationMatrix for web in webs)
    @test isapprox(mean(sum, webs), sum(P), rtol = 0.02)
    # and each link appears with its probability
    @test maximum(abs.(mean(webs) .- P)) < 0.05
    @test probabilistic_niche_network(Xoshiro(1), η, c, r) == probabilistic_niche_network(Xoshiro(1), η, c, r)

    # diets are close to intervals for sharp kernels and not for the Gaussian
    n_interval(web) = count(i -> is_interval(findall(web[:, i])), 1:60)
    gauss = mean(n_interval(probabilistic_niche_network(rng, 60, 0.1)) for _ in 1:20)
    sharp = mean(n_interval(probabilistic_niche_network(rng, 60, 0.1; p = 64)) for _ in 1:20)
    niche = mean(n_interval(niche_network(rng, 60, 0.1)) for _ in 1:20)
    @test gauss < sharp
    @test sharp >= 0.95 * niche

    # connectance stays on target for other a and p at small C, where little of the kernel
    # falls outside niche values 0 to 1
    for (a, p) in ((0.5, 2), (1, 1), (1, 8))
        mC = mean(connectance(probabilistic_niche_network(rng, 100, 0.05; a, p, tol = Inf)) for _ in 1:200)
        @test isapprox(mC, 0.05, rtol = 0.05)
    end

    @test_throws ArgumentError probabilistic_niche_network(rng, 20, 0.1; a = 0)
    @test_throws ArgumentError probabilistic_niche_network(rng, 20, 0.1; p = -1)
    # 0.2 / (0.2 Γ(3/2)) > 0.5
    @test_throws ArgumentError probabilistic_niche_network(rng, 20, 0.2; a = 0.2)
end
