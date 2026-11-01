using FoodWebs
using Test
using Random: Xoshiro

@testset "modularity by hand" begin
    # two disconnected links, 1 → 2 and 3 → 4: each group holds one link, and the expected
    # number of links inside a group is 1/2, so Q = (1/2)(1/2 + 1/2) = 1/2 directed, and the same
    # undirected
    web = Bool[0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
    @test modularity(web, [1, 1, 2, 2]) ≈ 1 / 2
    @test modularity(web, [1, 1, 2, 2]; directed = false) ≈ 1 / 2
    # Q is linear in γ: Q = e − γ a with e = 1 and a = 1/2 here
    @test modularity(web, [1, 1, 2, 2]; γ = 0) ≈ 1
    @test modularity(web, [1, 1, 2, 2]; γ = 2) ≈ 0
    # splitting each link across groups
    @test modularity(web, [1, 2, 1, 2]) ≈ 0
    # a single group, and a web without links
    @test modularity(web, ones(Int, 4)) ≈ 0 atol = 1.0e-15
    @test modularity(falses(3, 3), [1, 2, 3]) == 0
end

@testset "modularity invariances" begin
    rng = Xoshiro(81)
    web = niche_network(rng, 30, 0.15)
    groups = rand(rng, 1:3, 30)
    W = rand(rng, 30, 30) .* web
    for directed in (true, false)
        # direction convention, scale of the weights, and the labels themselves do not matter
        @test modularity(web', groups; directed) ≈ modularity(web, groups; directed)
        @test modularity(W', groups; directed) ≈ modularity(W, groups; directed)
        @test modularity(7 .* W, groups; directed) ≈ modularity(W, groups; directed)
        @test modularity(W, [:a, :b, :c][groups]; directed) ≈ modularity(W, groups; directed)
        # unit weights on the links give the binary web
        @test modularity(Float64.(web), groups; directed) ≈ modularity(web, groups; directed)
        # cannibalistic links are left out
        Wc = copy(W)
        Wc[1, 1] = Wc[5, 5] = 3.0
        @test modularity(Wc, groups; directed) ≈ modularity(W, groups; directed)
        @test abs(modularity(W, groups; directed)) <= 1
    end

    @test_throws DimensionMismatch modularity(web, groups[1:29])
    @test_throws ArgumentError modularity(-W, groups)
    @test_throws ArgumentError modularity(web, groups; γ = -1)
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
        Q = modularity(W, groups; directed, γ)
        @test Q ≈ networkx atol = 1.0e-12
        @test Q ≈ igraph atol = 1.0e-12
    end
end
