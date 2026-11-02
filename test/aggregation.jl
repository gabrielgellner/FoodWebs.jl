using FoodWebs
using Test

@testset "trophic species by hand" begin
    # 1 and 2 are eaten by both 3 and 4, which eat nothing else and have no predators; 5 eats 3
    web = zeros(Int, 5, 5)
    web[1, 3] = web[1, 4] = web[2, 3] = web[2, 4] = 1
    web[3, 5] = 1
    @test trophic_species(web) == [1, 1, 2, 3, 4]
    # two identical consumers lump to one
    web[4, 5] = 1
    @test trophic_species(web) == [1, 1, 2, 2, 3]
    lumped = lump_species(web, trophic_species(web))
    @test lumped == PredationMatrix([0 1 0; 0 0 1; 0 0 0])
    # lumping is idempotent: the groups of a lumped web are its species
    @test trophic_species(lumped) == 1:3
    # cannibals eating only themselves are different trophic species
    @test trophic_species([1 0; 0 1]) == [1, 2]
    # a group eats itself when one member eats another
    @test lump_species([0 1; 0 0], [:a, :a]) == PredationMatrix(fill(1, 1, 1))
    # any labels, grouped in order of first appearance
    @test lump_species(web, ["b", "b", "a", "a", "c"]) == lumped
    @test_throws DimensionMismatch lump_species(web, [1, 2])
end

@testset "trophic species of the Williams & Martinez webs" begin
    # these webs are trophic species already (Williams & Martinez 2000, Table 1)
    for id in filter(startswith("williams_martinez_2000"), foodweb_ids())
        fw = foodweb(id)
        @test trophic_species(fw) == 1:size(fw, 1)
        @test lump_species(fw, trophic_species(fw)) == fw
    end
end

@testset "lumping a FoodWeb" begin
    fw = FoodWeb(
        [0 0 1; 0 0 1; 0 0 0]; species = ["a", "b", "c"], name = "web",
        traits = Dict("M" => [1.0, 4.0, 10.0], "N" => [3.0, 1.0, 2.0], "category" => ["producer", "producer", "invertebrate"])
    )
    lumped = lump_species(fw, trophic_species(fw))
    @test lumped.species == ["a, b", "c"]
    @test lumped.name == "web"
    @test lumped.traits["N"] == [4.0, 2.0]
    # the mass weighted by abundance keeps the biomass: (1 × 3 + 4 × 1) / 4
    @test lumped.traits["M"] == [7 / 4, 10.0]
    @test lumped.traits["M"] .* lumped.traits["N"] == [7.0, 20.0]
    @test lumped.traits["category"] == ["producer", "invertebrate"]
    custom = lump_species(fw, [1, 1, 2]; names = ["plants", "grazer"], combine = Dict("M" => maximum))
    @test custom.species == ["plants", "grazer"] && custom.traits["M"] == [4.0, 10.0]
    # without abundances the mass is the plain mean, and missing values give missing
    plain = FoodWeb([0 0 1; 0 0 1; 0 0 0]; traits = Dict("M" => [1.0, 4.0, missing]))
    @test isequal(lump_species(plain, [1, 1, 2]).traits["M"], [2.5, missing])
    @test_throws DimensionMismatch lump_species(fw, [1, 1, 2]; names = ["one"])
end

@testset "trophic species against cheddar" begin
    # the trophic species and lumped webs of every cheddar web, from
    # test/fixtures/make_trophic_species.R, which sorts cheddar's links so that it compares prey
    # and predators as sets
    header, rows = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "trophic_species.csv"))
    col = Dict(h => k for (k, h) in enumerate(header))
    expected = Dict{String, Vector{Int}}()
    for row in rows
        push!(get!(expected, row[col["web"]], Int[]), parse(Int, row[col["trophic_species"]]))
    end
    header, rows = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "trophic_species_lumped.csv"))
    lcol = Dict(h => k for (k, h) in enumerate(header))
    value(x) = x == "missing" ? missing : parse(Float64, x)
    same(a, b) = ismissing(a) || ismissing(b) ? isequal(a, b) : isapprox(a, b; rtol = 1.0e-12)
    @test length(expected) == 19
    for (id, labels) in expected
        fw = foodweb(id)
        @test trophic_species(fw) == labels
        lumped = lump_species(fw, labels)
        groups = filter(r -> r[lcol["web"]] == id, rows)
        @test length(groups) == size(lumped, 1)
        for r in groups
            g = parse(Int, r[lcol["group"]])
            prey = isempty(r[lcol["resources"]]) ? Int[] : parse.(Int, split(r[lcol["resources"]]))
            @test findall(lumped[:, g]) == prey
            if haskey(fw.traits, "M")
                @test same(lumped.traits["M"][g], value(r[lcol["M"]]))
            end
            if haskey(fw.traits, "N")
                # cheddar keeps the mean abundance, the sum over the group's members here
                n = count(==(g), labels)
                @test same(lumped.traits["N"][g], n * value(r[lcol["N"]]))
            end
        end
        @test trophic_species(lumped) == 1:size(lumped, 1)
    end
end
