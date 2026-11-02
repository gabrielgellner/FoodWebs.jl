using FoodWebs
using Test
using Random: Xoshiro

@testset "diet breadth diets by hand" begin
    # four prey and three predators, named below by their masses: prey 1, 2, 4, 8 and predators
    # 10, 20, 50. With ai = aj = ni = 0 every encounter rate is a, E_i = M_i, and
    # H_ij = 1 / (b - M_i / M_j) for M_i / M_j < b.
    M = [1.0, 2, 4, 8, 10, 20, 50]
    consumers = [false, false, false, false, true, true, true]
    web = diet_breadth_network(M; handling = RatioHandling(; b = 0.45), a = 0.3, ni = 0, consumers)
    # predator 10 handles 1, 2 and 4, with profitabilities 0.35, 0.5 and 0.2. The rate on {2} is
    # 0.6 / 2.2 = 0.27 < 0.35, so it adds 1; the rate on {1, 2} is 0.9 / 3.06 = 0.29 > 0.2.
    @test M[web[:, 5]] == [1, 2]
    # predator 20 handles 1, 2, 4 and 8, with profitabilities 0.4, 0.7, 1.0 and 0.4. The rate on
    # {4} is 1.2 / 2.2 = 0.55 < 0.7; the rate on {2, 4} is 1.8 / 3.06 = 0.59 > 0.4.
    @test M[web[:, 6]] == [2, 4]
    # predator 50 handles every species up to 20, most profitable 10 (2.5) then 8 (2.32) then
    # 4 (1.48). The rate on {10} is 3 / 2.2 = 1.36; the rate on {8, 10} is 5.4 / 3.23 = 1.67.
    @test M[web[:, 7]] == [8, 10]
    @test !any(web[:, 1:4])
    # smaller encounter rates widen every diet: at a = 0.1 predator 10 eats all it can handle
    wide = diet_breadth_network(M; handling = RatioHandling(; b = 0.45), a = 0.1, ni = 0, consumers)
    @test M[wide[:, 5]] == [1, 2, 4]
    @test all(wide .>= web)
end

@testset "diet breadth model invariants" begin
    for seed in 1:5
        rng = Xoshiro(seed)
        M = exp10.(6 .* rand(rng, 40) .- 3)
        ratio = RatioHandling(; b = 0.3)
        web = diet_breadth_network(M; handling = ratio, a = 0.01, ai = 0.5, aj = 0.2)
        # diets are intervals of the resources ordered by mass, below b M_j
        order = sortperm(M)
        for j in 1:40
            ranks = findall(web[order, j])
            @test isempty(ranks) || ranks == first(ranks):last(ranks)
            @test all(M[web[:, j]] .< 0.3 * M[j])
        end
        # the diets depend on a n h, ai + ni and aj only, and not on e
        @test diet_breadth_network(M; handling = ratio, a = 0.01, ai = 0.5, aj = 0.2, e = 7) == web
        @test diet_breadth_network(M; handling = RatioHandling(; b = 0.3, h = 4), a = 0.0025, ai = 0.5, aj = 0.2) == web
        @test diet_breadth_network(M; handling = ratio, a = 0.005, n = 2, ai = 0.5, aj = 0.2) == web
        @test diet_breadth_network(M; handling = ratio, a = 0.01, ai = 0.25, ni = -0.5, aj = 0.2) == web
        @test diet_breadth_network(M; handling = ratio, a = 0.01, ai = 0.5, aj = 0.2, N = M .^ -0.75) == web
        # larger a narrows every diet
        narrow = diet_breadth_network(M; handling = ratio, a = 1, ai = 0.5, aj = 0.2)
        @test all(narrow .<= web)
        @test count(narrow) < count(web)
    end
end

@testset "diet breadth grows with consumer mass" begin
    # with power handling every consumer ranks resources the same way, and the diet of consumer
    # j depends on its mass only through a h M_j^(aj + hj), which narrows the diet as it grows
    for seed in 1:5
        M = sort(exp10.(4 .* rand(Xoshiro(seed), 30)))
        for (aj, hj) in ((0.3, -0.75), (0.5, 0.25))
            web = diet_breadth_network(M; handling = PowerHandling(; hi = 1.5, hj), a = 1.0e-3, ai = 0.2, aj)
            breadth = vec(sum(web; dims = 1))
            @test aj + hj < 0 ? issorted(breadth) : issorted(breadth; rev = true)
            @test breadth[1] != breadth[end]
            # the diets are nested, smallest resources first
            @test all(j -> findall(web[:, j]) == 1:breadth[j], 1:30)
        end
    end
end

@testset "diet breadth for empirical webs" begin
    fw = foodweb("cheddar/benguela")
    web = diet_breadth_network(fw; handling = RatioHandling(; b = 0.1), a = 1.0e-3)
    # only species with prey in the observed web eat
    @test !any(web[:, 1]) && !any(web[:, 2])
    @test any(web[:, 4])
    @test web == diet_breadth_network(
        Float64.(fw.traits["M"]); handling = RatioHandling(; b = 0.1), a = 1.0e-3,
        consumers = vec(any(fw; dims = 1))
    )
    @test_throws ArgumentError diet_breadth_network(foodweb("cheddar/broadstonestream"); handling = RatioHandling(; b = 0.1))
    @test_throws ArgumentError diet_breadth_network(foodweb("cheddar/skipwithpond"); handling = RatioHandling(; b = 0.1))
end

@testset "fitting the diet breadth model" begin
    # Petchey et al. (2008) predict 0.57 of the Benguela links with ratio handling times and 0.42
    # with power handling times (their Fig. 3)
    fw = foodweb("cheddar/benguela")
    for (handling, published) in ((:ratio, 0.57), (:power, 0.42))
        for seed in 1:3
            fit = fit_diet_breadth(Xoshiro(seed), fw; handling, starts = 5, steps = 200)
            @test count(fit.web) == count(fw)
            @test fit.correct == count(fit.web .& fw)
            @test fit.fraction == fit.correct / count(fw)
            @test fit.fraction > published
            @test diet_breadth_network(Float64.(fw.traits["M"]); fit.parameters...) == fit.web
        end
    end
    fit = fit_diet_breadth(Xoshiro(1), fw; starts = 2, steps = 20)
    @test fit_diet_breadth(Xoshiro(1), fw; starts = 2, steps = 20).web == fit.web
    @test fit.parameters.handling isa RatioHandling
end

@testset "diet breadth argument errors" begin
    M = [1.0, 2.0, 4.0]
    @test_throws ArgumentError RatioHandling(; b = 0)
    @test_throws ArgumentError RatioHandling(; b = 1, h = -1)
    @test_throws ArgumentError PowerHandling(; hi = 1, hj = 0, h = 0)
    @test_throws ArgumentError diet_breadth_network([1.0, 0.0]; handling = RatioHandling(; b = 1))
    @test_throws ArgumentError diet_breadth_network(M; handling = RatioHandling(; b = 1), a = 0)
    @test_throws ArgumentError diet_breadth_network(M; handling = RatioHandling(; b = 1), N = [1, 1, 0])
    @test_throws DimensionMismatch diet_breadth_network(M; handling = RatioHandling(; b = 1), N = [1, 1])
    @test_throws DimensionMismatch diet_breadth_network(M; handling = RatioHandling(; b = 1), consumers = [true])
    web = [0 1 0; 0 0 1; 0 0 0]
    @test_throws ArgumentError fit_diet_breadth(web, M; handling = :linear)
    @test_throws ArgumentError fit_diet_breadth(web, M; starts = 0)
    @test_throws ArgumentError fit_diet_breadth(zeros(Int, 3, 3), M)
    @test_throws DimensionMismatch fit_diet_breadth(web, [1.0, 2.0])
end
