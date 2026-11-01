using FoodWebs
using Test
using Random: Random, Xoshiro
using Statistics: mean, std, var
using Distributions: Normal, Uniform
using LinearAlgebra: Diagonal

@testset "modularity by hand" begin
    # two disconnected links, 1 → 2 and 3 → 4: each group holds one link, and the expected
    # number of links inside a group is 1/2, so Q = (1/2)(1/2 + 1/2) = 1/2 directed, and the same
    # undirected
    web = Bool[0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
    @test web_modularity(web, [1, 1, 2, 2]) ≈ 1 / 2
    @test web_modularity(web, [1, 1, 2, 2]; directed = false) ≈ 1 / 2
    # Q is linear in γ: Q = e − γ a with e = 1 and a = 1/2 here
    @test web_modularity(web, [1, 1, 2, 2]; γ = 0) ≈ 1
    @test web_modularity(web, [1, 1, 2, 2]; γ = 2) ≈ 0
    # splitting each link across groups
    @test web_modularity(web, [1, 2, 1, 2]) ≈ 0
    # a single group, and a web without links
    @test web_modularity(web, ones(Int, 4)) ≈ 0 atol = 1.0e-15
    @test web_modularity(falses(3, 3), [1, 2, 3]) == 0
end

@testset "modularity invariances" begin
    rng = Xoshiro(81)
    web = niche_network(rng, 30, 0.15)
    groups = rand(rng, 1:3, 30)
    W = rand(rng, 30, 30) .* web
    for directed in (true, false)
        # direction convention, scale of the weights, and the labels themselves do not matter
        @test web_modularity(web', groups; directed) ≈ web_modularity(web, groups; directed)
        @test web_modularity(W', groups; directed) ≈ web_modularity(W, groups; directed)
        @test web_modularity(7 .* W, groups; directed) ≈ web_modularity(W, groups; directed)
        @test web_modularity(W, [:a, :b, :c][groups]; directed) ≈ web_modularity(W, groups; directed)
        # unit weights on the links give the binary web
        @test web_modularity(Float64.(web), groups; directed) ≈ web_modularity(web, groups; directed)
        # cannibalistic links are left out
        Wc = copy(W)
        Wc[1, 1] = Wc[5, 5] = 3.0
        @test web_modularity(Wc, groups; directed) ≈ web_modularity(W, groups; directed)
        @test abs(web_modularity(W, groups; directed)) <= 1
    end

    @test_throws DimensionMismatch web_modularity(web, groups[1:29])
    @test_throws ArgumentError web_modularity(-W, groups)
    @test_throws ArgumentError web_modularity(web, groups; γ = -1)
end

@testset "modularity against NetworkX and igraph" begin
    # fixed partitions of every bundled web, from test/fixtures/make_modularity.py
    header, rows = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "modularity.csv"))
    @test header == ["web", "partition", "directed", "weighted", "gamma", "networkx", "igraph"]
    @test length(rows) == 26 * 16
    webs = Dict{String, Any}()
    for row in rows
        name, partition = row[1], row[2]
        directed, weighted, γ = row[3] == "True", row[4] == "True", parse(Float64, row[5])
        networkx, igraph = parse(Float64, row[6]), parse(Float64, row[7])
        web = get!(() -> foodweb(name), webs, name)
        S = size(web, 1)
        groups = partition == "thirds" ? [3 * (i - 1) ÷ S for i in 1:S] : [(i - 1) % 3 for i in 1:S]
        W = weighted ? [web[i, j] * (1 + (3 * (i - 1) + 7 * (j - 1)) % 5) for i in 1:S, j in 1:S] : web
        Q = web_modularity(W, groups; directed, γ)
        @test Q ≈ networkx atol = 1.0e-12
        @test Q ≈ igraph atol = 1.0e-12
    end
end

# three planted compartments of 20 species, linked inside with probability 0.3 and between with
# probability 0.01
function planted_web(rng)
    truth = repeat(1:3, inner = 20)
    web = [i != j && rand(rng) < (truth[i] == truth[j] ? 0.3 : 0.01) for i in 1:60, j in 1:60]
    return web, truth
end

@testset "partition search recovers planted compartments" begin
    web, truth = planted_web(Xoshiro(82))
    for method in (:louvain, :annealing), directed in (true, false)
        r = compartments(Xoshiro(1), web; method, directed, runs = 5)
        @test r.groups == truth
        @test r.modularity ≈ web_modularity(web, truth; directed)
        @test all(≈(r.modularity), r.run_modularity)
        @test all(x -> x == 0 || x == 1, r.coassignment)
    end
end

@testset "partition search properties" begin
    fw = foodweb("cheddar/tl84")
    S = size(fw, 1)
    linked(i, j) = fw[i, j] != 0 || fw[j, i] != 0
    for method in (:louvain, :annealing)
        r = compartments(Xoshiro(2), fw; method, runs = 4)
        @test r.modularity ≈ web_modularity(fw, r.groups)
        @test length(r.run_modularity) == 4 && all(q -> 0 < q <= 1, r.run_modularity)
        # labels 1, 2, … in order of first appearance
        @test r.groups[1] == 1 && all(maximum(r.groups[1:(i - 1)]) + 1 >= r.groups[i] for i in 2:S)
        # every group is connected
        for c in unique(r.groups)
            members = findall(==(c), r.groups)
            reached = Set([first(members)])
            frontier = [first(members)]
            while !isempty(frontier)
                v = pop!(frontier)
                for u in members
                    (u ∉ reached && linked(u, v)) && (push!(reached, u); push!(frontier, u))
                end
            end
            @test length(reached) == length(members)
        end
        D = r.coassignment
        @test D == D' && all(==(1), D[i, i] for i in 1:S) && all(x -> 0 <= x <= 1, D)
        # reproducible from a seeded rng
        @test compartments(Xoshiro(2), fw; method, runs = 4) == r
    end
    # a higher resolution gives at least as many groups
    @test maximum(compartments(Xoshiro(3), fw; γ = 2).groups) >= maximum(compartments(Xoshiro(3), fw; γ = 0.5).groups)
    # a web without links: every species alone
    r = compartments(Xoshiro(4), falses(4, 4))
    @test r.groups == 1:4 && r.modularity == 0

    @test_throws ArgumentError compartments(fw; method = :spectral)
    @test_throws ArgumentError compartments(fw; runs = 0)
    @test_throws ArgumentError compartments(fw; threshold = 0)
    @test_throws ArgumentError compartments(fw; γ = -1)
    @test_throws ArgumentError compartments(-ones(3, 3))
end

@testset "partition search against exact optima" begin
    # maximum modularity of the sparsest small webs, from test/fixtures/make_modularity_optimum.py
    header, rows = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "modularity_optimum.csv"))
    @test header == ["web", "directed", "optimum"]
    @test length(rows) == 10
    for (k, row) in enumerate(rows)
        fw, directed, optimum = foodweb(row[1]), row[2] == "True", parse(Float64, row[3])
        found = map((:louvain, :annealing)) do method
            r = compartments(Xoshiro(k), fw; directed, method, runs = 10)
            maximum(r.run_modularity)
        end
        @test all(<=(optimum + 1.0e-9), found)
        @test all(>=(optimum - 0.01), found)
        @test maximum(found) >= optimum - 0.001
    end
end

@testset "modularity z-scores" begin
    # webs drawn from the null model itself score z near 0 on average, for both kinds of null (the
    # mean of 15 z-scores has a standard error of about 0.27; the bound is about 3 of them)
    rng = Xoshiro(83)
    base = niche_network(rng, 30, 0.12)
    rewired_z = [modularity_zscore(rng, rewire_links(rng, base); samples = 30, runs = 3).z for _ in 1:15]
    @test abs(mean(rewired_z)) < 0.8
    niche_z = [modularity_zscore(rng, niche_network(rng, 30, 0.12); null = niche_network, samples = 30, runs = 3).z for _ in 1:15]
    @test abs(mean(niche_z)) < 0.8

    # planted compartments stand far above their rewired versions
    web, truth = planted_web(Xoshiro(84))
    r = modularity_zscore(Xoshiro(1), web; samples = 30, runs = 3)
    @test r.z > 5 && r.p == 0
    @test r.groups == truth
    @test r.observed ≈ web_modularity(web, truth)
    @test length(r.null_modularity) == 30
    @test r.mean ≈ mean(r.null_modularity) && r.sd ≈ std(r.null_modularity)
    # Tuesday Lake is more modular than its rewired versions (Guimerà et al. 2010)
    @test modularity_zscore(Xoshiro(2), foodweb("cheddar/tl84"); samples = 30, runs = 3).z > 3
    # reproducible from a seeded rng
    @test modularity_zscore(Xoshiro(3), base; samples = 10, runs = 2) == modularity_zscore(Xoshiro(3), base; samples = 10, runs = 2)

    @test_throws ArgumentError modularity_zscore(rand(Xoshiro(4), 5, 5))
    @test_throws ArgumentError modularity_zscore(base; null = :swap)
    @test_throws ArgumentError modularity_zscore(base; samples = 1)
end

@testset "species roles" begin
    # two three-species modules (1 → 2 → 3, 1 → 3, and the same for 4, 5, 6), and a top predator
    # 7, alone in its group, that eats 3 and 6
    web = falses(7, 7)
    for (prey, predator) in [(1, 2), (2, 3), (1, 3), (4, 5), (5, 6), (4, 6), (3, 7), (6, 7)]
        web[prey, predator] = true
    end
    groups = [1, 1, 1, 2, 2, 2, 3]
    r = species_roles(web, groups)
    # 7 splits its two links evenly over the two other groups; 3 has 2 of its 3 links at home
    @test r.P[7] ≈ 1 - 2 * (1 / 2)^2 && r.P[3] ≈ 1 - (2 / 3)^2 - (1 / 3)^2 && r.P[1] == 0
    # every member of a module has the same within-module degree, so z = 0
    @test all(iszero, r.z)
    @test r.role == [:peripheral for _ in 1:7]
    @test r.cartography[1] == :ultra_peripheral && r.cartography[7] == :peripheral
    # links to prey only: 7's prey are split over two groups, 1 has no prey
    prey = species_roles(web, groups; links = :prey)
    @test prey.P[7] ≈ 1 / 2 && prey.P[1] == 0 && prey.P[3] == 0
    # every species' predators are in a single group
    @test species_roles(web, groups; links = :predators).P == zeros(7)

    # a star module: species 1 linked to 2, 3, 4 and 5, so κ = (4, 1, 1, 1, 1) with mean 1.6 and
    # sample standard deviation √1.8
    star = falses(5, 5)
    star[1, 2:5] .= true
    @test species_roles(star, ones(Int, 5)).z ≈ [2.4, -0.6, -0.6, -0.6, -0.6] ./ sqrt(1.8)

    # weights: unit weights give the binary roles, and scale does not matter
    rng = Xoshiro(85)
    niche = niche_network(rng, 40, 0.12)
    g = compartments(rng, niche).groups
    W = rand(rng, 40, 40) .* niche
    @test species_roles(Float64.(niche), g) == species_roles(niche, g)
    R = species_roles(W, g)
    @test species_roles(5 .* W, g).P ≈ R.P && species_roles(5 .* W, g).z ≈ R.z
    @test all(p -> 0 <= p < 1, R.P)
    # the classes follow from z and P
    @test all(R.role[i] == (R.z[i] > 2.5 ? (R.P[i] > 0.62 ? :network_hub : :module_hub) : (R.P[i] > 0.62 ? :connector : :peripheral)) for i in 1:40)

    @test_throws ArgumentError species_roles(web, groups; links = :both)
    @test_throws DimensionMismatch species_roles(web, groups[1:6])
end

@testset "compartment homogeneity" begin
    # a chain 1 → 2 → 3 → 4 split into its lower and upper halves has one boundary, the fewest
    # possible, and more homogeneous than chance; alternating groups have three, the most
    chain = food_chain_network(4)
    @test FoodWebs.boundaries([1, 1, 2, 2], 1:4, [1, 2]) == 1
    @test compartment_homogeneity(chain, [1, 1, 2, 2]).trophic > 0
    @test compartment_homogeneity(chain, [1, 2, 1, 2]).trophic < 0

    # the exact null moments against random partitions with the same group sizes: shuffling the
    # labels samples that ensemble, and H_t is a z-score, so it has mean 0 and variance 1
    rng = Xoshiro(86)
    web = niche_network(rng, 30, 0.12; basal_paths = true)
    groups = compartments(rng, web).groups
    levels = trophic_levels(web) .+ 1.0e-6 .* (1:30)
    H = [compartment_homogeneity(web, Random.shuffle(rng, groups); levels) for _ in 1:20_000]
    trophic = [h.trophic for h in H]
    @test abs(mean(trophic)) < 0.03 && abs(var(trophic) - 1) < 0.05
    # H_s sums one z-score per species, so its mean is 0
    sink = [h.sink for h in H]
    @test abs(mean(sink)) < 4 * std(sink) / sqrt(length(sink))

    # with tied levels, the boundaries are the mean over random orderings of the ties
    tied = round.(trophic_levels(web); digits = 1)
    jittered = [FoodWebs.boundaries(groups, tied .+ 1.0e-6 .* rand(rng, 30), unique(groups)) for _ in 1:20_000]
    @test abs(mean(jittered) - FoodWebs.boundaries(groups, tied, unique(groups))) < 4 * std(jittered) / sqrt(length(jittered))

    @test_throws DimensionMismatch compartment_homogeneity(web, groups[1:29])
end

@testset "block structure" begin
    # constant blocks: group 1 = {1, 2} compete (−1 both ways), group 2 = {3, 4, 5} do not
    # interact, group 2 eats group 1 (−2 on group 1, +0.5 on group 2), diagonal −1
    M = [-1.0 -1 -2 -2 -2; -1 -1 -2 -2 -2; 0.5 0.5 -1 0 0; 0.5 0.5 0 -1 0; 0.5 0.5 0 0 -1]
    b = block_structure(M, [:a, :a, :b, :b, :b])
    @test b.labels == [:a, :b] && b.sizes == [2, 3]
    @test b.connectance == [1 1; 1 0]
    @test b.mean == [-1 -2; 0.5 0] && b.variance == zeros(2, 2)
    @test b.competition == [1 0; 0 0] && b.predator_prey == [0 6; 6 0]
    @test b.mutualism == zeros(Int, 2, 2) && b.one_way == zeros(Int, 2, 2)
    @test all(isnan, b.correlation)
    @test b.diagonal == [-1, -1]
    # species 1: 2 within its group and 3 × 2.5 with the other
    @test b.coupling[1] ≈ 7.5 / 9.5 && b.coupling[3] ≈ 1

    # one group gives the moments of elliptic_law
    rng = Xoshiro(86)
    J = random_predprey(rng, niche_network(rng, 40, 0.15), Normal(0, 1), 0.3)
    one = block_structure(J, ones(Int, 40))
    law = elliptic_law(J)
    @test one.mean[1, 1] ≈ law.μ && one.correlation[1, 1] ≈ law.ρ
    @test all(iszero, one.coupling)

    # block moments of a large matrix with known parameters: within blocks, pairs are present with
    # probability C and bivariate normal with mean μ, variance σ² and correlation ρ; between blocks,
    # entries are present with probability C_b and normal with mean μ_b
    S, C, μ, σ, ρ, C_b, μ_b = 400, 0.3, -0.5, 1.0, -0.6, 0.05, 1.0
    groups = repeat(1:2, inner = S ÷ 2)
    M = zeros(S, S)
    for i in 1:S, j in (i + 1):S
        if groups[i] == groups[j]
            rand(rng) < C || continue
            z1, z2 = randn(rng), randn(rng)
            M[i, j] = μ + σ * z1
            M[j, i] = μ + σ * (ρ * z1 + sqrt(1 - ρ^2) * z2)
        else
            rand(rng) < C_b && (M[i, j] = μ_b + σ * randn(rng))
            rand(rng) < C_b && (M[j, i] = μ_b + σ * randn(rng))
        end
    end
    b = block_structure(M, groups)
    variance(c, m) = c * (σ^2 + m^2) - c^2 * m^2
    @test b.connectance[1, 1] ≈ C atol = 0.01
    @test b.connectance[1, 2] ≈ C_b atol = 0.01
    @test b.mean[1, 1] ≈ C * μ atol = 0.01
    @test b.mean[2, 1] ≈ C_b * μ_b atol = 0.01
    @test b.variance[2, 2] ≈ variance(C, μ) atol = 0.02
    @test b.variance[1, 2] ≈ variance(C_b, μ_b) atol = 0.01
    @test b.correlation[1, 1] ≈ (C * (ρ * σ^2 + μ^2) - C^2 * μ^2) / variance(C, μ) atol = 0.03
    # independent entries between blocks: no pair correlation
    @test abs(b.correlation[1, 2]) < 0.03

    # every link of a cascade web is one predator–prey pair of random_predprey, and one one-way pair
    # of a flux matrix
    web = cascade_network(rng, 30, 0.15)
    groups = compartments(rng, web).groups
    L = count(!iszero, web)
    upper(X) = sum(X[r, s] for r in axes(X, 1) for s in axes(X, 2) if r <= s)
    @test upper(block_structure(random_predprey(rng, web, Uniform(0, 1), 0.3), groups).predator_prey) == L
    flux = rand(rng, 30, 30) .* web
    fb = block_structure(flux, groups)
    @test upper(fb.one_way) == L && upper(fb.predator_prey) == 0

    @test_throws DimensionMismatch block_structure(M, groups[1:29])
end

# the bipartite modularity of Guimerà, Sales-Pardo & Amaral (2007), from teams of members
function team_modularity(teams, S, groups; γ = 1)
    c = zeros(Int, S, S)
    for team in teams, i in team, j in team
        i != j && (c[i, j] += 1)
    end
    t = [count(team -> i in team, teams) for i in 1:S]
    K = sum(length(team) * (length(team) - 1) for team in teams)
    L = sum(length, teams)
    same = [(i, j) for i in 1:S for j in 1:S if i != j && groups[i] == groups[j]]
    return sum(c[i, j] for (i, j) in same; init = 0) / K - γ * sum(t[i] * t[j] for (i, j) in same; init = 0) / L^2
end

@testset "predator and prey views" begin
    # prey 1 and 2 both eaten by predators 3 and 4: predators 3 and 4 share two prey, so
    # M = (2 + 2)/4 − (2·2 + 2·2)/4² = 1/2 in the predator view, and the same in the prey view
    web = falses(4, 4)
    web[1:2, 3:4] .= true
    @test web_modularity(web, [1, 1, 2, 2]; view = :predator) ≈ 1 / 2
    @test web_modularity(web, [1, 1, 2, 2]; view = :prey) ≈ 1 / 2

    # against the definition from teams, on niche webs and random partitions
    rng = Xoshiro(84)
    for _ in 1:5
        web = niche_network(rng, 25, 0.15)
        groups = rand(rng, 1:4, 25)
        eats(i, j) = i != j && web[i, j] != 0
        prey_of = [[i for i in 1:25 if eats(i, j)] for j in 1:25]
        predators_of = [[j for j in 1:25 if eats(i, j)] for i in 1:25]
        # predator view: each prey's predators are a team; prey view: each predator's prey
        @test web_modularity(web, groups; view = :predator) ≈ team_modularity(filter(!isempty, predators_of), 25, groups)
        @test web_modularity(web, groups; view = :prey) ≈ team_modularity(filter(!isempty, prey_of), 25, groups)
        @test web_modularity(web, groups; view = :prey, γ = 0.5) ≈ team_modularity(filter(!isempty, prey_of), 25, groups; γ = 0.5)
        # the views agree for a symmetric web and differ for this one
        sym = web .| web'
        @test web_modularity(sym, groups; view = :predator) ≈ web_modularity(sym, groups; view = :prey)
        @test !(web_modularity(web, groups; view = :predator) ≈ web_modularity(web, groups; view = :prey))
    end

    # planted guilds: prey 1–10 eaten by predators 21–25 and prey 11–20 by predators 26–30, each
    # link present with probability 0.6; the prey view groups each set of prey, and the predator
    # view each set of predators
    web = falses(30, 30)
    for i in 1:20, j in 21:30
        ((i <= 10) == (j <= 25)) && rand(rng) < 0.6 && (web[i, j] = true)
    end
    # (a species without predators, or without prey, shares them with no one and stays alone)
    eaten = [i for i in 1:30 if any(web[i, :])]
    eating = [j for j in 1:30 if any(web[:, j])]
    prey = compartments(Xoshiro(1), web; view = :prey, runs = 5).groups
    A, B = intersect(eaten, 1:10), intersect(eaten, 11:20)
    @test allequal(prey[A]) && allequal(prey[B]) && prey[A[1]] != prey[B[1]]
    predators = compartments(Xoshiro(1), web; view = :predator, method = :annealing, runs = 3).groups
    A, B = intersect(eating, 21:25), intersect(eating, 26:30)
    @test allequal(predators[A]) && allequal(predators[B]) && predators[A[1]] != predators[B[1]]
    r = modularity_zscore(Xoshiro(2), web; view = :prey, samples = 20, runs = 3)
    @test r.z > 3

    @test_throws ArgumentError web_modularity(rand(Xoshiro(3), 4, 4), [1, 1, 2, 2]; view = :prey)
    @test_throws ArgumentError web_modularity(web, ones(Int, 30); view = :guild)
end

# a random topology with compartments in the strengths only: each pair of 60 species interacts
# with probability 0.25 as predator and prey, strongly (magnitude 1 to 2) within the planted halves
# and weakly (0 to 0.1) between them
function strength_planted(rng)
    truth = repeat(1:2, inner = 30)
    J, web = zeros(60, 60), falses(60, 60)
    for i in 1:60, j in (i + 1):60
        rand(rng) < 0.25 || continue
        m = truth[i] == truth[j] ? 1 + rand(rng) : 0.1 * rand(rng)
        prey, predator = rand(rng, Bool) ? (i, j) : (j, i)
        J[prey, predator], J[predator, prey] = -m, 0.3m
        web[prey, predator] = true
    end
    return J, web, truth
end

@testset "pair weights of community matrices" begin
    # a predator–prey pair (−2, 0.5), a one-way pair (−1, 0), and a competing pair (−1, −4)
    J = [-1.0 -2 -1 -1; 0.5 -1 0 0; 0 0 -1 0; -4 0 0 -1]
    G = pair_weights(J)
    @test G[1, 2] == G[2, 1] == 1 && G[1, 3] == 0 && G[1, 4] == 2 && all(iszero, G[i, i] for i in 1:4)
    T = pair_weights(J; method = :sum)
    @test T[1, 2] == 2.5 && T[1, 3] == 1 && T[1, 4] == 5 && T == T'

    rng = Xoshiro(89)
    web = niche_network(rng, 30, 0.15)
    J = random_predprey(rng, web, Uniform(0, 1), 0.3)
    # the geometric weight does not depend on the abundance units of each species
    D = Diagonal(exp.(randn(rng, 30)))
    @test pair_weights(D * J / D) ≈ pair_weights(J)
    # with equal magnitudes, the summed weights give the undirected modularity of a web without
    # two-way pairs (a pair that eat each other is two links of the binary web but one pair)
    web = cascade_network(rng, 30, 0.15)
    unit = [i == j ? -1.0 : web[i, j] ? -1.0 : web[j, i] ? 1.0 : 0.0 for i in 1:30, j in 1:30]
    groups = compartments(rng, web).groups
    @test web_modularity(pair_weights(unit; method = :sum), groups) ≈ web_modularity(web, groups; directed = false)

    # shuffled weights keep who interacts and the set of weights, and symmetry
    G = pair_weights(J)
    H = shuffle_weights(rng, G)
    @test H == H' && (H .!= 0) == (G .!= 0) && sort(vec(H)) == sort(vec(G))
    flux = rand(rng, 30, 30) .* web
    F = shuffle_weights(rng, flux)
    @test (F .!= 0) == (flux .!= 0) && sort(vec(F)) == sort(vec(flux))

    # compartments in the strengths only: the binary web is no more modular than its rewired
    # versions, while the weights are far more modular than their shuffled versions, and the
    # weighted groups follow the planted halves (a species with only weak links may join the other)
    J, web, truth = strength_planted(Xoshiro(90))
    @test modularity_zscore(Xoshiro(1), web; samples = 30, runs = 3).z < 2
    r = modularity_zscore(Xoshiro(1), pair_weights(J); null = :shuffle, samples = 30, runs = 3)
    @test r.z > 3
    purity = sum(maximum(count(==(h), truth[r.groups .== g]) for h in 1:2) for g in unique(r.groups)) / 60
    @test purity >= 0.85
    # weights drawn from the shuffle null itself score z near 0 on average
    W = pair_weights(J)
    z = [modularity_zscore(rng, shuffle_weights(rng, W); null = :shuffle, samples = 20, runs = 2).z for _ in 1:12]
    @test abs(mean(z)) < 0.7

    @test_throws ArgumentError pair_weights(J; method = :product)
    @test_throws ArgumentError modularity_zscore(pair_weights(J))
end

@testset "energy channels" begin
    # two chains 1 → 2 → 3 and 4 → 5 → 6, joined by 7, which takes three times as much from 3 as
    # from 6: 7 gets 3/4 of its energy from source 1 and joins its channel
    F = zeros(7, 7)
    F[1, 2] = F[2, 3] = F[4, 5] = F[5, 6] = 1
    F[3, 7], F[6, 7] = 3, 1
    r = energy_channels(F)
    @test r.reliance[7, :] ≈ [0.75, 0.25]
    @test r.reliance[1:6, :] == [1 0; 1 0; 1 0; 0 1; 0 1; 0 1]
    @test r.groups == [1, 1, 1, 2, 2, 2, 1] && r.dominant == [1, 1, 1, 2, 2, 2, 1]
    @test r.modularity ≈ web_modularity(F, r.groups)
    # a binary web splits each diet equally
    @test source_reliance(F .> 0)[7, :] ≈ [0.5, 0.5]

    # the reliance solves R = PᵀR + E on a flux web, and each row sums to 1
    rng = Xoshiro(90)
    web = niche_network(rng, 30, 0.15; basal_paths = true)
    W = rand(rng, 30, 30) .* web
    sources = [j for j in 1:30 if !any(i != j && web[i, j] for i in 1:30)]
    R = source_reliance(W)
    w = [i == j ? 0.0 : W[i, j] for i in 1:30, j in 1:30]
    for i in 1:30
        expected = i in sources ? float.(sources .== i) : sum(w[k, i] .* R[k, :] for k in 1:30) ./ sum(w[:, i])
        @test R[i, :] ≈ expected
    end
    # merging never lowers the modularity, and stops when no pair of channels would raise it
    merged, unmerged = energy_channels(W), energy_channels(W; merge = false)
    @test merged.modularity >= unmerged.modularity
    labels = unique(merged.groups)
    @test all(
        web_modularity(W, [g == b ? a : g for g in merged.groups]) <= merged.modularity + 1.0e-12
            for a in labels for b in labels if a != b
    )

    # detritus that eats (receives) from consumers is a source when named
    D = zeros(5, 5)
    D[1, 2] = D[2, 3] = 1   # producer 1 → 2 → 3
    D[3, 4] = 1             # 3 → detritus 4
    D[4, 5] = 1             # detritus 4 → detritivore 5
    @test source_reliance(D)[5, :] ≈ [1]
    @test source_reliance(D; sources = [1, 4])[5, :] ≈ [0, 1]

    # errors: a species that eats nothing but is not a source, and a loop no source enters
    @test_throws ArgumentError source_reliance(F; sources = [1])
    loop = zeros(3, 3)
    loop[2, 3] = loop[3, 2] = 1
    @test_throws ArgumentError source_reliance(loop; sources = [1])
end

@testset "energy channels of Chesapeake Bay (Zhao et al. 2018)" begin
    # the carbon flows of Baird & Ulanowicz (1989), from enaR, by test/fixtures/make_chesapeake_bay.py
    _, nodes = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "chesapeake_bay_nodes.csv"))
    _, flows = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "chesapeake_bay_flows.csv"))
    names = [row[2] for row in nodes]
    living = [row[3] == "True" for row in nodes]
    F = zeros(36, 36)
    for (from, to, flow) in flows
        F[parse(Int, from), parse(Int, to)] += parse(Float64, flow)
    end
    # the basal nodes: the non-living carbon pools and the species that eat nothing
    sources = [i for i in 1:36 if !living[i] || !any(F[k, i] > 0 for k in 1:36 if k != i)]
    @test names[sources] ==
        ["Phytoplankton", "Benthic Diatoms", "Dissolved Organic Carbon", "Suspended Particulate Carbon", "Sediment Partculate Carbon"]
    r = energy_channels(F; sources)
    # "the energy of 13 taxa was mainly supplied by sediment particulate nutrient, while suspended
    # particulate nutrient was the major supplier for another 12 taxa, phytoplankton was the major
    # supplier for four taxa, dissolved organic nutrient was the major supplier for two taxa, and
    # benthic diatoms were not a major supplier for any taxa"
    consumers = setdiff(1:36, sources)
    @test [count(==(k), r.dominant[consumers]) for k in 1:5] == [4, 0, 2, 12, 13]
    # "Three compartments were detected, one of which was mainly sustained by suspended particulate
    # nutrient and phytoplankton; another compartment was mainly sustained by sediment particulate
    # nutrient but also contained benthic diatoms, while the final compartment contained just three
    # taxa: free bacteria, heterotrophic microflagellates, and dissolved organic nutrient"
    @test maximum(r.groups) == 3
    member(name) = r.groups[findfirst(==(name), names)]
    @test member("Phytoplankton") == member("Suspended Particulate Carbon")
    @test member("Benthic Diatoms") == member("Sediment Partculate Carbon")
    @test sort(names[r.groups .== member("Dissolved Organic Carbon")]) ==
        ["Dissolved Organic Carbon", "Free Bacteria", "Heterotrophic Microflagelates"]
    # Q = 0.39 in their Table 1, the standard directed modularity (their 2W is the total flux)
    @test round(r.modularity; digits = 2) == 0.39
end

@testset "group model" begin
    # 1 eats nothing; 3 eats 1 and 2, 4 eats 1: three of the four links from group 1 = {1, 2} to
    # group 2 = {3, 4}, none elsewhere
    web = falses(4, 4)
    web[1, 3] = web[2, 3] = web[1, 4] = true
    m = group_model(web, [:a, :a, :b, :b])
    @test m.labels == [:a, :b] && m.p == [0 0.75; 0 0]
    ll = 3 * log(0.75) + log(0.25)
    @test m.loglik ≈ ll
    @test m.aic ≈ 2 * 4 + 2 * 4 - 2ll
    @test m.icl ≈ ll + 2 * 2 * log(1 / 2) - 2 * log(16) - log(4) / 2
    # one group: the connectance L / S²
    one = group_model(web, ones(Int, 4))
    @test only(one.p) == 3 / 16
    @test one.loglik ≈ 3 * log(3 / 16) + 13 * log(13 / 16)

    # three trophic layers of 15 species: each layer eaten by the next with probability 0.4,
    # other pairs with probability 0.02. The ICL selects three groups, the layers, whose
    # probabilities are dense above the diagonal and sparse on it; the AIC of the paper adds groups
    rng = Xoshiro(88)
    truth = repeat(1:3, inner = 15)
    layered = [i != j && rand(rng) < (truth[j] == truth[i] + 1 ? 0.4 : 0.02) for i in 1:45, j in 1:45]
    f = fit_group_model(Xoshiro(1), layered; k = 1:5)
    @test f.k == 3 && f.groups == truth
    @test f.p[1, 2] > 0.3 && f.p[2, 3] > 0.3 && all(f.p[r, r] < 0.1 for r in 1:3)
    @test f.loglik >= group_model(layered, truth).loglik - 1.0e-9
    @test f.criterion_by_k[3] == maximum(values(f.criterion_by_k))
    # a random web: one group
    random = [rand(rng) < 0.1 for _ in 1:45, _ in 1:45]
    @test fit_group_model(Xoshiro(2), random; k = 1:4).k == 1
    # reproducible from a seeded rng
    @test fit_group_model(Xoshiro(3), layered; k = 3) == fit_group_model(Xoshiro(3), layered; k = 3)

    @test_throws ArgumentError group_model(rand(Xoshiro(4), 4, 4), [1, 1, 2, 2])
    @test_throws ArgumentError fit_group_model(layered; k = 0:2)
    @test_throws ArgumentError fit_group_model(layered; criterion = :bic)
    @test_throws DimensionMismatch group_model(web, [1, 2])
end

@testset "signed modularity" begin
    rng = Xoshiro(91)
    # without negative weights, signed modularity is the modularity
    web = niche_network(rng, 30, 0.15)
    groups = compartments(rng, web).groups
    W = rand(rng, 30, 30) .* web
    @test web_modularity(W, groups; signed = true) ≈ web_modularity(W, groups)
    # with both signs, Q = (W⁺Q⁺ − W⁻Q⁻)/(W⁺ + W⁻) for the modularities of the two parts, which the
    # fixtures check against NetworkX and igraph
    M = randn(rng, 30, 30) .* (web .| web')
    P, N = max.(M, 0), max.(-M, 0)
    Wp, Wn = sum(P) - sum(P[i, i] for i in 1:30), sum(N) - sum(N[i, i] for i in 1:30)
    for directed in (true, false)
        expected = (Wp * web_modularity(P, groups; directed) - Wn * web_modularity(N, groups; directed)) / (Wp + Wn)
        @test web_modularity(M, groups; signed = true, directed) ≈ expected
    end

    # signed pair weights: predator–prey pairs negative, competing pairs positive
    J = [-1.0 -2 -1; 0.5 -1 0; -4 0 -1]
    @test pair_weights(J; method = :signed)[1, 2] ≈ -1 && pair_weights(J; method = :signed)[1, 3] ≈ 2

    # planted signed groups: positive weights within three groups of 15, negative between, on a
    # random topology; the signed search recovers them
    truth = repeat(1:3, inner = 15)
    M = zeros(45, 45)
    for i in 1:45, j in (i + 1):45
        rand(rng) < 0.3 || continue
        M[i, j] = M[j, i] = (truth[i] == truth[j] ? 1 : -1) * (0.5 + rand(rng))
    end
    # annealing finds them; Louvain, a faster heuristic, comes within 5%
    r = compartments(Xoshiro(1), M; signed = true, method = :annealing, runs = 3)
    @test r.groups == truth
    @test r.modularity ≈ web_modularity(M, truth; signed = true)
    r = compartments(Xoshiro(1), M; signed = true, runs = 3)
    @test maximum(r.run_modularity) >= 0.95 * web_modularity(M, truth; signed = true)

    @test_throws ArgumentError web_modularity(M, truth)
    @test_throws ArgumentError web_modularity(M .> 0, truth; view = :prey, signed = true)
    @test_throws ArgumentError pair_weights(J; method = :difference)
end

@testset "weighted block model" begin
    # two groups {1, 2} and {3, 4}: the entries from group 1 to group 2 are −1 and −3, those back
    # are 0.5 and 0.5, and the rest are 0 (diagonal ignored)
    M = [-9.0 0 -1 0; 0 -9 0 -3; 0.5 0 -9 0; 0 0.5 0 -9]
    b = block_model(M, [1, 1, 2, 2])
    @test b.p == [0 0.5; 0.5 0] && b.μ == [0 -2; 0.5 0]
    # shrunk variances: s₀² is the variance of all four nonzero entries
    x = [-1.0, -3, 0.5, 0.5]
    s02 = mean(abs2, x .- mean(x))
    @test b.σ[1, 2] ≈ sqrt((2 + s02) / 3) && b.σ[2, 1] ≈ sqrt((0 + s02) / 3) && b.σ[1, 1] ≈ sqrt(s02)
    normal(x, μ, σ) = -log(2π * σ^2) / 2 - (x - μ)^2 / (2σ^2)
    # each off-diagonal block has 2 nonzero entries of 4, and the diagonal blocks none of 2
    ll = 8 * log(0.5) + normal(-1, -2, b.σ[1, 2]) + normal(-3, -2, b.σ[1, 2]) + 2 * normal(0.5, 0.5, b.σ[2, 1])
    @test b.loglik ≈ ll
    # with the block means of block_structure: p μ, counting zero entries
    rng = Xoshiro(92)
    J = random_predprey(rng, niche_network(rng, 30, 0.15), Uniform(0, 1), 0.3)
    groups = compartments(rng, J .< 0).groups
    @test block_model(J, groups).p .* block_model(J, groups).μ ≈ block_structure(J, groups).mean

    # prey eaten by predators, with the same connectance within and between the groups: the
    # pattern of nonzero entries is random, and the signs and strengths set the two groups apart
    truth = [fill(1, 20); fill(2, 20)]
    M = block_random_matrix(rng, [20, 20]; connectance = 0.2, μ = [-0.5 -1.5; 1.0 -0.5], σ = 0.3, diagonal = -1)
    f = fit_block_model(Xoshiro(1), M; k = 1:4)
    @test f.k == 2 && f.groups == truth
    @test f.μ[1, 2] < 0 < f.μ[2, 1]
    @test f.loglik >= block_model(M, truth).loglik - 1.0e-9
    # a homogeneous matrix: one group
    @test fit_block_model(Xoshiro(2), block_random_matrix(rng, [40]; connectance = 0.2, μ = -0.5, σ = 0.5); k = 1:4).k == 1
    @test fit_block_model(Xoshiro(3), M; k = 2) == fit_block_model(Xoshiro(3), M; k = 2)

    @test_throws DimensionMismatch block_model(M, truth[1:3])
    @test_throws ArgumentError fit_block_model(M; k = 0:2)
end
