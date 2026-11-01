using FoodWebs
using Test
using Random: Xoshiro

# named counts of the motifs present
present(web; size = 3) = Dict(motif_name(i) => n for (i, n) in zip(motif_ids(size), motif_census(web; size)) if n > 0)

@testset "motif classes" begin
    # mfinder IDs of the 13 three-species motifs (Milo et al. 2002) and the count of
    # four-species classes, with the bi-fan and bi-parallel (diamond) IDs
    @test motif_ids(3) == [6, 12, 14, 36, 38, 46, 74, 78, 98, 102, 108, 110, 238]
    @test length(motif_ids(4)) == 199
    @test 204 in motif_ids(4) && 904 in motif_ids(4)
    @test motif_name(12) == "S1"
    @test motif_name(238) == "D6"
    @test motif_name(904) == "M904"
    @test_throws ArgumentError motif_ids(5)
end

@testset "motif census by hand" begin
    # 2 eats 1, 3 eats 1 and 2: omnivory
    @test present([0 1 1; 0 0 1; 0 0 0]) == Dict("S2" => 1)
    @test present(food_chain_network(4)) == Dict("S1" => 2)
    # a three-species loop
    @test present([0 0 1; 1 0 0; 0 1 0]) == Dict("S3" => 1)
    # 1 and 2 eat each other, and 3 eats 1 (D4) or 1 eats 3 (D3)
    @test present([0 1 1; 1 0 0; 0 0 0]) == Dict("D4" => 1)
    @test present([0 1 0; 1 0 0; 1 0 0]) == Dict("D3" => 1)
    # 2 and 3 eat 1 (exploitative competition), or 3 eats 1 and 2 (apparent competition)
    @test present([0 1 1; 0 0 0; 0 0 0]) == Dict("S4" => 1)
    @test present([0 0 1; 0 0 1; 0 0 0]) == Dict("S5" => 1)
    # the diamond: 2 and 3 eat 1, and 4 eats 2 and 3
    diamond = zeros(Bool, 4, 4)
    diamond[1, 2] = diamond[1, 3] = diamond[2, 4] = diamond[3, 4] = true
    @test present(diamond; size = 4) == Dict("M904" => 1)
    @test present(diamond) == Dict("S1" => 2, "S4" => 1, "S5" => 1)
    # cannibalism is left out, and disconnected sets are not counted
    @test present([1 1 1; 0 1 1; 0 0 1]) == Dict("S2" => 1)
    @test isempty(present(zeros(Bool, 4, 4)))
    @test_throws ArgumentError motif_census(diamond; size = 2)
end

@testset "motif census against igraph and NetworkX" begin
    # test/fixtures/make_motif_census.py counts the motifs of every bundled web with igraph
    # (three and four species) and NetworkX's triad census (three species)
    header, rows = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "motif_census.csv"))
    @test header == ["web", "size", "motif", "igraph", "networkx"]
    expected = Dict{Tuple{String, Int}, Dict{Int, Int}}()
    for row in rows
        web, size, motif, count = row[1], parse(Int, row[2]), parse(Int, row[3]), parse(Int, row[4])
        get!(expected, (web, size), Dict{Int, Int}())[motif] = count
        size == 3 && @test parse(Int, row[5]) == count
    end
    @test length(expected) == 2 * 26
    for ((web, size), counts) in expected
        fw = foodweb(web)
        @test Dict(zip(motif_ids(size), motif_census(fw; size))) == counts
    end
end

@testset "motif roles" begin
    # in the omnivory motif, 1 is eaten by two species, 2 is in the middle, 3 eats two
    roles = motif_roles([0 1 1; 0 0 1; 0 0 0])
    @test [FoodWebs.ROLE_KEYS[findfirst(>(0), roles[s, :])] for s in 1:3] ==
        [(38, 2, 0), (38, 1, 1), (38, 0, 2)]
    # top, middle and bottom of a food chain: positions 1, 2 and 3
    roles = motif_roles(food_chain_network(3))
    @test roles[3, 1] == roles[2, 2] == roles[1, 3] == 1
    @test sum(roles) == 3
    # each occurrence of a motif puts each of its three species in one of its positions
    for web in (foodweb("williams_martinez_2000/coachella_valley"), foodweb("cheddar/tl84"))
        roles = motif_roles(web)
        census = Dict(zip(motif_ids(3), motif_census(web)))
        for id in motif_ids(3)
            positions = [r for (r, key) in enumerate(FoodWebs.ROLE_KEYS) if key[1] == id]
            @test sum(roles[:, positions]) == 3 * census[id]
        end
    end
end

@testset "rewired webs" begin
    rng = Xoshiro(18)
    web = foodweb("williams_martinez_2000/coachella_valley")
    A = [i != j && web[j, i] for i in 1:29, j in 1:29]
    one_way(A, i) = (count(j -> A[i, j] && !A[j, i], 1:29), count(j -> A[j, i] && !A[i, j], 1:29))
    two_way(A, i) = count(j -> A[i, j] && A[j, i], 1:29)
    for _ in 1:5
        r = rewire_links(rng, web)
        B = [i != j && r[j, i] for i in 1:29, j in 1:29]
        # each species keeps its one-way prey and predators, its two-way links and cannibalism
        @test all(one_way(B, i) == one_way(A, i) for i in 1:29)
        @test all(two_way(B, i) == two_way(A, i) for i in 1:29)
        @test all(r[i, i] == web[i, i] for i in 1:29)
        @test r != web
    end
    @test rewire_links(Xoshiro(1), web) == rewire_links(Xoshiro(1), web)
    # a web with fewer than two links cannot be rewired
    @test rewire_links(rng, [0 1; 0 0]) == [0 1; 0 0]
end

@testset "motif z-scores" begin
    rng = Xoshiro(19)
    web = foodweb("williams_martinez_2000/st_martin_island")
    p = motif_zscores(rng, web; samples = 20)
    @test p.ids == motif_ids(3)
    @test p.observed == motif_census(web)
    @test length(p.z) == 13
    # St. Martin has no two-way links, so the rewired webs have none either and those motifs
    # have no variance
    @test all(isnan(p.z[k]) for (k, id) in enumerate(p.ids) if startswith(motif_name(id), "D"))
    @test all(isfinite, p.z[[findfirst(==(id), p.ids) for id in (12, 38, 36, 6)]])
    @test isequal(motif_zscores(Xoshiro(2), web; samples = 5).z, motif_zscores(Xoshiro(2), web; samples = 5).z)
    # the profile is the z-scores scaled to length 1, with the undefined ones as 0
    @test sum(abs2, p.profile) ≈ 1
    @test all(iszero(p.profile[k]) for k in eachindex(p.z) if isnan(p.z[k]))
    @test p.profile[1] ≈ p.z[1] / sqrt(sum(abs2, filter(!isnan, p.z)))
end

@testset "Stouffer et al. (2007) motif profiles" begin
    # Stouffer et al. (2007, Fig. 3) split empirical webs by the signs of their S2 (omnivory),
    # S4 and S5 (competition) z-scores against the same null model: (a) S2 over-represented and
    # S4, S5 under-represented, or (b) the reverse. Every bundled web they studied falls in its
    # published group. Across 10 seeds with 50 rewired webs the smallest |z| among these is
    # 0.66 (Benguela S2), so the signs are stable at 100.
    groups = Dict(
        "williams_martinez_2000/little_rock_lake" => :a, "williams_martinez_2000/coachella_valley" => :a,
        "williams_martinez_2000/chesapeake_bay" => :a, "williams_martinez_2000/st_martin_island" => :a,
        "cheddar/benguela" => :a, "williams_martinez_2000/bridge_brook_lake" => :b,
        "williams_martinez_2000/skipwith_pond" => :b, "williams_martinez_2000/ythan_estuary" => :b,
    )
    rng = Xoshiro(20)
    for (web, group) in groups
        p = motif_zscores(rng, foodweb(web); samples = 100)
        z(id) = p.z[findfirst(==(id), p.ids)]
        sign_a = z(38) > 0 && z(36) < 0 && z(6) < 0
        sign_b = z(38) < 0 && z(36) > 0 && z(6) > 0
        @test group === :a ? sign_a : sign_b
        # and three-species loops are under-represented in all seven Williams & Martinez webs
        startswith(web, "williams_martinez_2000/") && @test z(98) < 0
    end
end
