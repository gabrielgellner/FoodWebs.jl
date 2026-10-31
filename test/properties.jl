using FoodWebs
using Test
using Random: Xoshiro
using Statistics: mean

# `x` agrees with a published value to the precision it is printed with
function agrees(x, published::AbstractString)
    decimals = length(published) - something(findfirst('.', published), length(published))
    return abs(x - parse(Float64, published)) <= 0.5 * 10.0^(-decimals) + 1.0e-12
end

@testset "species types" begin
    # 1 basal, 2 eats 1, 3 eats 2 and itself, 4 isolated
    web = [0 1 0 0; 0 0 1 0; 0 0 1 0; 0 0 0 0]
    @test basal_species(web) == [1, 4]
    @test intermediate_species(web) == [2, 3]
    @test top_species(web) == [4]
    @test cannibals(web) == [3]
    # a cannibal with no other predators has itself as a predator, so it is not top
    @test top_species([0 1; 0 1]) == Int[]
    @test intermediate_species([0 1; 0 1]) == [2]
end

@testset "loops" begin
    # 1 basal, 2 and 3 eat each other and 2 eats 1, 4 eats 3, 5 eats 1
    web = zeros(Bool, 5, 5)
    web[1, 2] = web[2, 3] = web[3, 2] = web[3, 4] = web[1, 5] = true
    @test species_in_loops(web) == [2, 3, 4]
    # cannibalism is not a loop
    @test isempty(species_in_loops([1 1; 0 0]))
    @test isempty(species_in_loops(food_chain_network(5)))
end

@testset "omnivores" begin
    # 3 eats 1 and 2, so its chains have lengths 1 and 2
    @test omnivores([0 1 1; 0 0 1; 0 0 0]) == [3]
    # 3 eats two basal species: two prey, one chain length
    @test isempty(omnivores([0 0 1; 0 0 1; 0 0 0]))
    @test isempty(omnivores(food_chain_network(5)))
    # a cannibal eating one other species has one prey
    @test isempty(omnivores([0 1; 0 1]))
end

@testset "generality and vulnerability" begin
    web = [0 1 1; 0 0 1; 0 0 0]
    # S / L = 1
    @test generality(web) == [0, 1, 2]
    @test vulnerability(web) == [2, 1, 0]
    rng = Xoshiro(14)
    for _ in 1:10
        w = niche_network(rng, 40, 0.1)
        @test mean(generality(w)) ≈ 1
        @test mean(vulnerability(w)) ≈ 1
    end
end

@testset "trophic similarity" begin
    # 2 and 3 both eat 1, and 4 eats 3
    web = zeros(Bool, 4, 4)
    web[1, 2] = web[1, 3] = web[3, 4] = true
    s = trophic_similarity(web)
    @test s == s'
    # 2 and 3 share their one prey; 3 also has a predator: 1 shared of 2
    @test s[2, 3] == 1 / 2
    @test s[1, 4] == 0
    @test all(s[i, i] == 1 for i in 1:4)
    # maxima: 1 -> 0 (shares nothing), 2 -> 1/2, 3 -> 1/2, 4 -> 0
    @test mean_max_similarity(web) == 1 / 4
    @test trophic_similarity(zeros(Bool, 2, 2)) == zeros(2, 2)
end

@testset "trophic levels: shortest and short-weighted" begin
    web = [0 1 1; 0 0 1; 0 0 0]
    @test trophic_levels(web; method = :shortest) == [1, 2, 2]
    @test trophic_levels(web; method = :short_weighted) ≈ [1, 2, 2.25]
    @test trophic_levels(food_chain_network(4); method = :short_weighted) ≈ 1:4
    # cannibalism is ignored
    @test trophic_levels([1 1; 0 1]; method = :shortest) == [1, 2]
    @test_throws ArgumentError trophic_levels([0 1; 1 0]; method = :shortest)
    @test_throws ArgumentError trophic_levels(web; method = :longest)
    # the shortest level is never above the prey-averaged one
    rng = Xoshiro(15)
    for _ in 1:20
        w = generalized_cascade_network(rng, 40, 0.15)
        @test all(trophic_levels(w; method = :shortest) .<= trophic_levels(w) .+ 1.0e-9)
    end
end

@testset "characteristic path length" begin
    # undirected distances 1, 1 and 2, each counted twice, over 3^2 pairs
    @test characteristic_path_length(food_chain_network(3)) ≈ 8 / 9
    # link direction is ignored
    @test characteristic_path_length([0 0; 1 0]) == characteristic_path_length([0 1; 0 0])
    @test_throws ArgumentError characteristic_path_length(zeros(Bool, 3, 3))
end

@testset "structural properties" begin
    p = structural_properties([0 1 1; 0 0 1; 0 0 0])
    @test p.T + p.I + p.B ≈ 1
    @test p.Omniv == 1 / 3
    @test p.TL ≈ (1 + 2 + 2.25) / 3
    @test p.MaxTL ≈ 2.25
    # trophic levels are undefined for species cut off from basal species: 1 and 2 eat each
    # other, and 3 eats 1
    @test isnan(structural_properties([0 1 1; 1 0 0; 0 0 0]).TL)
    # and path lengths for disconnected webs
    @test isnan(structural_properties([0 1 0; 0 0 0; 0 0 0]).PathLength)
end

@testset "Williams & Martinez (2000) Table 2" begin
    # Table 2 of Williams & Martinez (2000): T, I, B, GenSD, VulSD, MxSim, Loop, Cannib, Omniv.
    # Little Rock Lake has no Omniv. Their Ythan Estuary has 78 species and ours 83 (Dunne et
    # al. 2002), so it is left out.
    table = Dict(
        "skipwith_pond" => ("0.040", "0.92", "0.040", "0.92", "0.54", "0.76", "0.40", "0.32", "0.60"),
        "little_rock_lake" => ("0.011", "0.86", "0.13", "1.42", "0.61", "0.75", "0.33", "0.14", ""),
        "bridge_brook_lake" => ("0", "0.68", "0.32", "1.09", "0.61", "0.74", "0.32", "0.12", "0.40"),
        "chesapeake_bay" => ("0.32", "0.52", "0.16", "0.78", "1.12", "0.50", "0", "0.032", "0.52"),
        "coachella_valley" => ("0", "0.90", "0.10", "0.73", "0.60", "0.72", "0.62", "0.66", "0.76"),
        "st_martin_island" => ("0.17", "0.69", "0.14", "1.02", "0.78", "0.54", "0", "0", "0.60"),
    )
    names = (:T, :I, :B, :GenSD, :VulSD, :MxSim, :Loop, :Cannib, :Omniv)
    # Values that differ from Table 2 beyond its rounding, with the largest difference found:
    # Chesapeake Bay has one link fewer than their version (C = 0.071 against 0.072); our MxSim
    # is lower for three webs, and Coachella's loops leave Omniv at 0.79.
    tolerance = Dict(
        ("chesapeake_bay", :GenSD) => 0.03, ("chesapeake_bay", :VulSD) => 0.02,
        ("little_rock_lake", :GenSD) => 0.01, ("little_rock_lake", :MxSim) => 0.01,
        ("bridge_brook_lake", :MxSim) => 0.04, ("skipwith_pond", :MxSim) => 0.04,
        ("coachella_valley", :Omniv) => 0.04,
    )
    for (web, published) in table
        p = structural_properties(foodweb("williams_martinez_2000/$web"))
        for (name, value) in zip(names, published)
            isempty(value) && continue
            if haskey(tolerance, (web, name))
                @test abs(p[name] - parse(Float64, value)) <= tolerance[(web, name)]
            else
                @test agrees(p[name], value)
            end
        end
    end
end

@testset "Dunne, Williams & Martinez (2002) path lengths" begin
    # characteristic path length D from their Table 1, including the 83-species Ythan Estuary
    published = Dict(
        "skipwith_pond" => "1.33", "little_rock_lake" => "1.89", "bridge_brook_lake" => "1.85",
        "chesapeake_bay" => "2.65", "ythan_estuary" => "2.20", "coachella_valley" => "1.42",
        "st_martin_island" => "1.88",
    )
    for (web, D) in published
        @test agrees(characteristic_path_length(foodweb("williams_martinez_2000/$web")), D)
    end
end
