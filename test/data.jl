using FoodWebs
using Test

@testset "FoodWeb" begin
    fw = FoodWeb([0 1; 0 0]; species = ["grass", "rabbit"], name = "meadow")
    @test fw isa AbstractMatrix{Bool}
    @test size(fw) == (2, 2)
    @test fw[1, 2] && !fw[2, 1]
    @test fw.species == ["grass", "rabbit"]
    @test connectance(fw) == 1 / 4
    @test basal_species(fw) == [1]
    @test sprint(show, MIME("text/plain"), fw) == "FoodWeb \"meadow\" with 2 species and 1 links"
    # species are named by index by default
    @test FoodWeb(food_chain_network(3)).species == ["1", "2", "3"]
    @test_throws DimensionMismatch FoodWeb([0 1; 0 0]; species = ["a"])
    @test_throws ArgumentError FoodWeb([0 1; 0 0]; species = ["a", "a"])
    @test_throws ArgumentError FoodWeb([0 2; 0 0])
end

@testset "CSV lines" begin
    split = FoodWebs.split_csv_line
    @test split("a,b,c") == ["a", "b", "c"]
    @test split("a,,c") == ["a", "", "c"]
    @test split("\"Daphnia, large\",b") == ["Daphnia, large", "b"]
    @test split("\"say \"\"hi\"\"\",x") == ["say \"hi\"", "x"]
    @test split("\"\"") == [""]
end

@testset "read_foodweb" begin
    dir = mktempdir()
    write(joinpath(dir, "nodes.csv"), "node,category\nalgae,producer\n\"Daphnia, large\",consumer\nfish,consumer\n")
    write(joinpath(dir, "trophic.links.csv"), "resource,consumer\nalgae,\"Daphnia, large\"\n\"Daphnia, large\",fish\nfish,fish\n")
    fw = read_foodweb(dir; name = "pond")
    @test fw.species == ["algae", "Daphnia, large", "fish"]
    @test fw.name == "pond"
    @test fw == Bool[0 1 0; 0 0 1; 0 0 1]
    # columns are found by name
    write(joinpath(dir, "trophic.links.csv"), "consumer,resource\nfish,algae\n")
    @test read_foodweb(dir)[1, 3]
    write(joinpath(dir, "trophic.links.csv"), "resource,consumer\nalgae,heron\n")
    @test_throws ArgumentError read_foodweb(dir)
    write(joinpath(dir, "trophic.links.csv"), "from,to\nalgae,fish\n")
    @test_throws ArgumentError read_foodweb(dir)
end

@testset "bundled food webs" begin
    ids = foodweb_ids()
    @test count(startswith("cheddar/"), ids) == 19
    @test count(startswith("williams_martinez_2000/"), ids) == 7
    for id in filter(id -> startswith(id, "cheddar/") || startswith(id, "williams_martinez_2000/"), ids)
        fw = foodweb(id)
        meta = fw.metadata
        @test size(fw, 1) == meta["species"]
        @test count(fw.web.links) == meta["links"]
        @test !isempty(meta["citation"]) && !isempty(meta["license"]) && !isempty(meta["source"])
    end
    @test_throws ArgumentError foodweb("cheddar/no_such_web")

    # each bundled collection's datapackage.json lists both files of every web in its catalog
    for collection in ("cheddar", "williams_martinez_2000")
        dir = joinpath(pkgdir(FoodWebs), "data", collection)
        package = read(joinpath(dir, "datapackage.json"), String)
        for id in keys(FoodWebs.TOML.parsefile(joinpath(dir, "catalog.toml")))
            @test occursin("\"$id/nodes.csv\"", package)
            @test occursin("\"$id/trophic.links.csv\"", package)
        end
    end

    # cheddar tags primary producers, which must be basal: this checks the orientation
    for id in ("cheddar/tl84", "cheddar/benguela", "cheddar/broadstonestream")
        fw = foodweb(id)
        header, rows = FoodWebs.read_csv(joinpath(pkgdir(FoodWebs), "data", id, "nodes.csv"))
        category = findfirst(==("category"), header)
        producers = [k for (k, row) in enumerate(rows) if row[category] == "producer"]
        @test !isempty(producers)
        @test all(in(basal_species(fw)), producers)
    end

    # species and connectance of the Williams & Martinez (2000) webs, as listed in Dunne,
    # Williams & Martinez (2002) Table 1
    published = Dict(
        "skipwith_pond" => (25, 0.315), "coachella_valley" => (29, 0.312),
        "chesapeake_bay" => (31, 0.071), "st_martin_island" => (42, 0.116),
        "bridge_brook_lake" => (25, 0.171), "little_rock_lake" => (92, 0.118),
        "ythan_estuary" => (83, 0.057),
    )
    for (name, (S, C)) in published
        fw = foodweb("williams_martinez_2000/$name")
        @test size(fw, 1) == S
        @test round(connectance(fw), digits = 3) == C
    end
end

@testset "downloaded food webs" begin
    cache = mktempdir()
    dir = mkpath(joinpath(cache, "mine", "pond"))
    write(joinpath(dir, "nodes.csv"), "node\nalgae\nfish\n")
    write(joinpath(dir, "trophic.links.csv"), "resource,consumer\nalgae,fish\n")
    write(
        joinpath(cache, "mine", "catalog.toml"),
        "[pond]\ntitle = \"A pond\"\nspecies = 2\nlinks = 1\ncitation = \"none\"\n"
    )
    withenv("FOODWEBS_CACHE" => cache) do
        @test FoodWebs.foodweb_cache() == cache
        @test "mine/pond" in foodweb_ids()
        fw = foodweb("mine/pond")
        @test fw.name == "A pond"
        @test fw == Bool[0 1; 0 0]
    end
end
