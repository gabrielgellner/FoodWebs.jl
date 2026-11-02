using FoodWebs
using Test
using Random: Xoshiro
using Statistics: mean

@testset "node traits" begin
    fw = foodweb("cheddar/tl84")
    @test fw.traits["M"] isa Vector{Union{Missing, Float64}}
    @test fw.traits["category"][1] == "producer"
    # detritus has no measured mass
    @test count(ismissing, foodweb("cheddar/broadstonestream").traits["M"]) == 5
    # webs whose nodes.csv has only names have no traits
    @test isempty(foodweb("williams_martinez_2000/skipwith_pond").traits)
    # traits survive simplification
    @test remove_cannibalism(fw)[1].traits["M"] === fw.traits["M"]
    @test_throws DimensionMismatch FoodWeb([0 1; 0 0]; traits = Dict("M" => [1.0]))
end

@testset "observed biomass" begin
    fw = foodweb("cheddar/tl84")
    @test isequal(biomass(fw), fw.traits["M"] .* fw.traits["N"])
    @test all(>(0), biomass(fw))
    @test count(ismissing, biomass(foodweb("cheddar/ythanestuary"))) == 1
    @test_throws ArgumentError biomass(foodweb("cheddar/benguela"))
    @test_throws ArgumentError biomass(foodweb("williams_martinez_2000/skipwith_pond"))
end

@testset "steady-state biomass" begin
    # a food chain with efficiency TE and turnover falling by k per level changes biomass by
    # TE k per level
    chain = food_chain_network(5)
    x = steady_state_biomass(chain; efficiency = 0.1, turnover = 2.0 .^ -(0:4))
    @test all(x[k + 1] / x[k] ≈ 0.2 for k in 1:4)
    @test pyramid_slope(chain, x) ≈ log(0.2)
    @test is_pyramid(chain, x) && !is_inverted_pyramid(chain, x)
    # allometric turnover with a constant mass ratio R gives TE R^(1/4) per level, so the
    # pyramid inverts when TE R^(1/4) > 1
    for (R, inverted) in ((1.0e2, false), (1.0e6, true))
        fw = FoodWeb(chain; traits = Dict("M" => R .^ (0:4)))
        x = steady_state_biomass(fw; efficiency = 0.1, turnover = :allometric)
        @test all(x[k + 1] / x[k] ≈ 0.1 * R^(1 / 4) for k in 1:4)
        @test is_inverted_pyramid(fw, x) == inverted
    end

    # flow balance: P = P_basal + diag(TE) Dᵀ P, with x = P / τ
    rng = Xoshiro(24)
    for gen in (cascade_network, niche_network, nested_hierarchy_network)
        web = gen(rng, 30, 0.12; basal_paths = true)
        τ = rand(rng, 30) .+ 0.5
        x = steady_state_biomass(web; efficiency = 0.2, turnover = τ)
        @test all(>(0), x)
        P = x .* τ
        D = FoodWebs.production_shares(web, :equal)
        basal = [!any(i != j && web[i, j] for i in 1:30) for j in 1:30]
        @test P ≈ basal .+ 0.2 .* (D' * P)
        # consumers produce less in total than the basal species supply
        @test sum(P[.!basal]) < sum(P[basal])
    end

    # shares as a matrix: a predator taking all of one prey's production
    web = Bool[0 1 1; 0 0 0; 0 0 0]
    x = steady_state_biomass(web; efficiency = 0.5, shares = [0 1 0; 0 0 0; 0 0 0])
    @test x ≈ [1, 0.5, 0]
    @test_throws ArgumentError steady_state_biomass(web; shares = [0 0.7 0.7; 0 0 0; 0 0 0])
    @test_throws ArgumentError steady_state_biomass(web; shares = [0 1 0; 1 0 0; 0 0 0])

    # noise is reproducible and multiplicative
    web = cascade_network(Xoshiro(25), 20, 0.15)
    x0 = steady_state_biomass(web)
    x1 = steady_state_biomass(Xoshiro(1), web; noise = 0.3)
    @test x1 == steady_state_biomass(Xoshiro(1), web; noise = 0.3)
    @test x1 != x0 && all(>(0), x1)

    @test_throws ArgumentError steady_state_biomass([0 1; 1 0])
    @test_throws ArgumentError steady_state_biomass(chain; efficiency = 1.0)
    @test_throws ArgumentError steady_state_biomass(chain; turnover = :allometric)
    @test_throws ArgumentError steady_state_biomass(foodweb("cheddar/broadstonestream"); turnover = :allometric)
end

@testset "pyramid measures" begin
    # 2 and 3 eat 1, 3 eats 2
    web = Bool[0 1 1; 0 0 1; 0 0 0]
    r = biomass_ratios(web, [4.0, 2.0, 1.0])
    @test Set((t.prey, t.predator, t.ratio) for t in r) == Set([(1, 2, 0.5), (1, 3, 0.25), (2, 3, 0.5)])
    @test is_pyramid(web, [4.0, 2.0, 1.0])
    @test !is_pyramid(web, [4.0, 1.0, 2.0])
    @test is_inverted_pyramid(web, [1.0, 2.0, 4.0])
    @test pyramid_slope(web, [1.0, 1.0, 1.0]) == 0
    @test_throws ArgumentError pyramid_slope(web, [1.0, 0.0, 1.0])
end

@testset "body mass ratios" begin
    # 1 → 2 → 3 with 3 also eating 1, and 2 a cannibal
    web = [0 1 1; 0 1 1; 0 0 0]
    M = [1.0, 10.0, 1000.0]
    @test mass_ratios(web, M) == [
        (prey = 1, predator = 2, ratio = 10.0), (prey = 1, predator = 3, ratio = 1000.0),
        (prey = 2, predator = 3, ratio = 100.0),
    ]
    # Cohen, Schittler, Raffaelli & Reuman (2009), Table S3: links with the prey lighter than,
    # as heavy as and heavier than the predator, without cannibalism and Ythan's detritus
    for (id, counts) in (("cheddar/tl84", (262, 0, 2)), ("cheddar/tl86", (232, 0, 4)), ("cheddar/ythanestuary", (368, 2, 9)))
        ratios = collect(skipmissing(r.ratio for r in mass_ratios(foodweb(id))))
        @test (count(>(1), ratios), count(==(1), ratios), count(<(1), ratios)) == counts
    end
    @test_throws ArgumentError mass_ratios(foodweb("williams_martinez_2000/skipwith_pond"))
end

@testset "mass-abundance relationship" begin
    # N = 100 M^(-3/4) exactly
    M = [1.0, 10.0, 100.0, 1000.0]
    fit = mass_abundance_slope(M, 100 .* M .^ -0.75)
    @test fit.slope ≈ -0.75 && fit.intercept ≈ 2
    # Cohen, Jonsson & Carpenter (2003) report -0.83 for Tuesday Lake in 1984; the coefficients
    # are those of cheddar's tests
    tl84 = mass_abundance_slope(foodweb("cheddar/tl84"))
    @test tl84.slope ≈ -0.82711453844476423569
    @test tl84.intercept ≈ -2.6862753 atol = 1.0e-7
    broadstone = mass_abundance_slope(foodweb("cheddar/broadstonestream"))
    @test round.((broadstone.intercept, broadstone.slope), digits = 2) == (1.08, -0.74)
    @test_throws DimensionMismatch mass_abundance_slope(M, M[1:3])
    @test_throws ArgumentError mass_abundance_slope(M, -M)
    @test_throws ArgumentError mass_abundance_slope([1.0], [1.0])
    @test_throws ArgumentError mass_abundance_slope(foodweb("cheddar/benguela"))
end

@testset "links in the mass-abundance plane" begin
    # a resource of mass 1 and abundance 100 eaten by a consumer of mass 10 and abundance 1
    link = only(mass_abundance_links([0 1; 0 0], [1.0, 10.0], [100.0, 1.0]))
    @test link.prey == 1 && link.predator == 2
    @test link.length ≈ 3
    @test link.slope ≈ -2
    @test link.angle ≈ -63.43494882292201
    # a lighter, rarer predator points straight back: -180, never 180
    @test only(mass_abundance_links([0 1; 0 0], [10.0, 1.0], [1.0, 1.0])).angle == -180
    # Cohen, Schittler, Raffaelli & Reuman (2009), Table 1: the number of links and the mean link
    # length, without cannibalism and species missing a mass or abundance
    for (id, L, mean_length) in (("cheddar/tl84", 264, 6.33), ("cheddar/tl86", 236, 5.9), ("cheddar/ythanestuary", 379, 7.29))
        links = mass_abundance_links(foodweb(id))
        @test length(links) == L
        @test round(mean(l.length for l in links), digits = 2) == mean_length
        @test all(l -> -180 <= l.angle < 180, links)
    end
    @test_throws DimensionMismatch mass_abundance_links([0 1; 0 0], [1.0], [1.0, 1.0])
    @test_throws ArgumentError mass_abundance_links([0 1; 0 0], [1.0, 0.0], [1.0, 1.0])
end

@testset "trophic level against body mass" begin
    # a chain whose masses grow tenfold per level gains one level per decade of mass
    fit = mass_level_slope(food_chain_network(4), [1.0, 10.0, 100.0, 1000.0])
    @test fit.slope ≈ 1 && fit.intercept ≈ 1
    @test mass_level_slope(food_chain_network(3), [1.0, 1.0e2, 1.0e4]; levels = [1, 2, 3]).slope ≈ 0.5
    # Ythan's detritus has no mass and is left out
    fw = foodweb("cheddar/ythanestuary")
    levels = trophic_levels(fw)
    keep = .!ismissing.(fw.traits["M"])
    @test mass_level_slope(fw).slope ≈ mass_level_slope(fw.web[keep, keep], Float64.(fw.traits["M"][keep]); levels = levels[keep]).slope
    @test mass_level_slope(foodweb("cheddar/tl84")).slope > 0
    @test_throws DimensionMismatch mass_level_slope(food_chain_network(3), [1.0, 2.0])
    @test_throws ArgumentError mass_level_slope(food_chain_network(2), [1.0, -1.0])
end
