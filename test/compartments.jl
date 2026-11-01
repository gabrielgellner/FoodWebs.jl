using FoodWebs
using Test
using Random: Xoshiro

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
