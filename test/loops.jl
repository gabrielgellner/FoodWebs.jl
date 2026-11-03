using FoodWebs
using Test
using CairoMakie: CairoMakie
using Distributions: Uniform
using Graphs: DiGraph, add_edge!, simplecycles_limited_length
using LinearAlgebra: Diagonal, diag
using Random: Xoshiro

# Species 1 is eaten by 2 and 3, and 3 also eats 2: an omnivory triangle
omnivory_matrix() = [
    -1.0 -2.0 -3.0
    0.5 0.0 -4.0
    0.6 0.7 0.0
]

# A cycle as a sorted rotation starting from its lowest species, for comparison with Graphs.jl
function canonical(cycle)
    k = argmin(cycle)
    return [cycle[k:end]; cycle[1:(k - 1)]]
end

@testset "feeding loops by hand" begin
    J = omnivory_matrix()
    loops = feeding_loops(J)
    @test length(loops) == 5
    @test issorted([l.weight for l in loops]; rev = true)
    pairs = filter(l -> l.length == 2, loops)
    @test Set(l.species for l in pairs) == Set([[1, 2], [1, 3], [2, 3]])
    @test all(l -> l.sign == -1, pairs)
    @test only(l.weight for l in pairs if l.species == [2, 3]) ≈ sqrt(4 * 0.7)
    # 1 → 2 → 3 → 1 runs through two top-down effects and one bottom-up effect: two negative
    # entries and one positive, so the loop is positive
    up = only(l for l in loops if l.species == [1, 2, 3])
    @test up.links == [(1, 2), (2, 3), (3, 1)]
    @test up.sign == 1
    @test up.weight ≈ cbrt(2 * 4 * 0.6)
    down = only(l for l in loops if l.species == [1, 3, 2])
    @test down.sign == -1
    @test down.weight ≈ cbrt(3 * 0.7 * 0.5)
    # the diagonal takes no part
    @test feeding_loops(J - Diagonal(diag(J))) == loops
    # min_length, max_length and one direction per loop
    @test length(feeding_loops(J; min_length = 3)) == 2
    @test length(feeding_loops(J; max_length = 2)) == 3
    @test feeding_loops(J; min_length = 3, both_directions = false) == [up]
    @test heaviest_loop(J) == up
    @test heaviest_loop(J; min_length = 2) == loops[1]
    # weights scale with J
    @test [l.weight for l in feeding_loops(-3J)] ≈ 3 .* [l.weight for l in loops]

    # a food chain has only predator–prey loops
    chain = [
        -1.0 -1.0 0.0 0.0
        1.0 0.0 -1.0 0.0
        0.0 1.0 0.0 -1.0
        0.0 0.0 1.0 0.0
    ]
    @test length(feeding_loops(chain)) == 3
    @test isempty(feeding_loops(chain; min_length = 3))
    @test heaviest_loop(chain) === nothing

    @test_throws ArgumentError feeding_loops(J; min_length = 1)
    @test_throws ArgumentError feeding_loops(J; min_length = 4, max_length = 3)
    @test_throws ArgumentError feeding_loops(J; limit = 4)
    @test_throws ArgumentError heaviest_loop(J; max_length = 2)
    @test_throws DimensionMismatch feeding_loops(ones(2, 3))
end

@testset "feeding loops of random webs" begin
    for seed in 1:4
        web = niche_network(Xoshiro(seed + 900), 10, 0.2)
        J = random_predprey(Xoshiro(seed), web, Uniform(0.1, 2), 0.5)
        S = size(J, 1)
        g = DiGraph(S)
        for i in 1:S, j in 1:S
            i != j && J[i, j] != 0 && add_edge!(g, i, j)
        end
        L = 5
        loops = feeding_loops(J; max_length = L)
        @test sort([l.species for l in loops]) == sort(canonical.(simplecycles_limited_length(g, L)))
        @test all(l -> l.weight ≈ prod(abs(J[i, j]) for (i, j) in l.links)^(1 / l.length), loops)
        @test all(l -> l.sign == sign(prod(J[i, j] for (i, j) in l.links)), loops)
        # each loop of three or more species is listed in both directions, since every link
        # gives J two nonzero entries
        one = feeding_loops(J; max_length = L, both_directions = false)
        @test count(l -> l.length > 2, loops) == 2count(l -> l.length > 2, one)
        @test count(l -> l.length == 2, loops) == count(l -> l.length == 2, one)
        # the pruned search agrees with the list
        for m in 2:3
            @test heaviest_loop(J; min_length = m, max_length = L) ==
                first(feeding_loops(J; min_length = m, max_length = L))
        end
        # and finds loops at least as heavy when it may go longer
        @test heaviest_loop(J; max_length = 8).weight >= heaviest_loop(J; max_length = L).weight
    end
end

@testset "maximum mean cycle" begin
    # in the omnivory triangle the heaviest loop over all lengths is 1 → 2 → 3 → 1, so Karp's
    # algorithm finds it without a search
    J = omnivory_matrix()
    logw = [J[i, j] == 0 || i == j ? -Inf : log(abs(J[i, j])) for i in 1:3, j in 1:3]
    out = [[j for j in 1:3 if j != i && J[i, j] != 0] for i in 1:3]
    bound, cycle = FoodWebs.max_mean_cycle(logw, out)
    @test exp(bound) ≈ cbrt(2 * 4 * 0.6)
    @test cycle in ([1, 2, 3], [2, 3, 1], [3, 1, 2])
    @test FoodWebs.max_mean_cycle(fill(-Inf, 2, 2), [Int[], Int[]]) == (-Inf, nothing)
    for seed in 1:20
        web = niche_network(Xoshiro(seed + 990), 7, 0.25)
        J = random_predprey(Xoshiro(seed), web, Uniform(0.1, 2), 0.5)
        S = size(J, 1)
        logw = [J[i, j] == 0 || i == j ? -Inf : log(abs(J[i, j])) for i in 1:S, j in 1:S]
        out = [[j for j in 1:S if j != i && J[i, j] != 0] for i in 1:S]
        bound, cycle = FoodWebs.max_mean_cycle(logw, out)
        # the bound is the heaviest of all loops, of any length
        all_loops = feeding_loops(J; max_length = S)
        @test exp(bound) ≈ all_loops[1].weight
        @test cycle !== nothing
        @test exp(sum(logw[cycle[m], cycle[mod1(m + 1, end)]] for m in eachindex(cycle)) / length(cycle)) ≈
            all_loops[1].weight
        @test heaviest_loop(J; min_length = 2, max_length = S) == all_loops[1]
        # when the heaviest loop is a predator–prey pair the search finds the heaviest longer one
        for (m, L) in ((3, S), (3, 4), (4, 5))
            loops = feeding_loops(J; min_length = m, max_length = L)
            @test heaviest_loop(J; min_length = m, max_length = L) == (isempty(loops) ? nothing : loops[1])
        end
    end
end

@testset "loop membership" begin
    J = omnivory_matrix()
    loops = feeding_loops(J; min_length = 3)
    @test loop_membership(loops, 4) == Bool[1 1; 1 1; 1 1; 0 0]
    # each species is in two of the three predator–prey loops
    @test vec(sum(loop_membership(feeding_loops(J; max_length = 2), 3); dims = 2)) == [2, 2, 2]
    mask = loop_mask(loops[1:1], 3)
    @test findall(mask) == sort(CartesianIndex.(loops[1].links))
    # every off-diagonal entry lies on some predator–prey loop
    @test loop_mask(feeding_loops(J; max_length = 2), 3) == (J .!= 0) .& .!Bool[1 0 0; 0 1 0; 0 0 1]
end

@testset "minimum self-regulation" begin
    for seed in 1:4
        web = niche_network(Xoshiro(seed + 950), 12, 0.15)
        J = random_predprey(Xoshiro(seed), web, Uniform(0.1, 2), 0.5)
        offdiag = J - Diagonal(diag(J))
        # with equal scale the diagonal shifts every eigenvalue
        @test min_self_regulation(J) ≈ spectral_abscissa(offdiag) rtol = 1.0e-8
        @test min_self_regulation(J; scale = fill(2.0, 12)) ≈ min_self_regulation(J) / 2 rtol = 1.0e-8
        d = rand(Xoshiro(seed + 50), 12) .+ 0.1
        s = min_self_regulation(J; scale = d)
        @test spectral_abscissa(offdiag - (1 + 1.0e-6) * s * Diagonal(d)) < 0
        @test spectral_abscissa(offdiag - (1 - 1.0e-6) * s * Diagonal(d)) >= 0
    end
    # the off-diagonal part has trace 0, so its spectral abscissa is at least 0; at 0, any
    # self-regulation stabilizes a predator–prey pair and a chain
    @test min_self_regulation([0.0 -1.0; 1.0 0.0]) == 0
    @test min_self_regulation([0.0 -1.0; 0.0 0.0]) == 0
    @test min_self_regulation(omnivory_matrix()) > 0
    @test_throws DimensionMismatch min_self_regulation(ones(2, 2); scale = ones(3))
    @test_throws ArgumentError min_self_regulation(ones(2, 2); scale = [1.0, 0.0])
end

@testset "loop plots" begin
    web = niche_network(Xoshiro(3), 12, 0.2; basal_paths = true)
    J = random_predprey(Xoshiro(3), web, Uniform(0.1, 2), 0.5)
    loops = feeding_loops(J; min_length = 3, both_directions = false)[1:4]
    fig, ax = loopplot(web, loops; title = "loops")
    @test fig isa CairoMakie.Figure
    path = tempname() * ".png"
    CairoMakie.save(path, fig)
    @test filesize(path) > 0
    fig = CairoMakie.Figure()
    @test loopplot!(CairoMakie.Axis(fig[1, 1]), web, loops[1:0]) isa CairoMakie.Axis
    @test_throws ArgumentError loopplot(web, [(species = [1, 13], links = [(1, 13), (13, 1)], length = 2, sign = -1, weight = 1.0)])
    # arrows move to the right of their direction
    ext = Base.get_extension(FoodWebs, :FoodWebsMakieExt)
    o, d = ext.offset_arrow_geometry(
        [CairoMakie.Point2f(0, 0)], [CairoMakie.Point2f(1, 0)], [0], [0], CairoMakie.Vec2f(100, 100), 5
    )
    @test o[1] ≈ CairoMakie.Point2f(0, -0.05)
    @test d[1] ≈ CairoMakie.Vec2f(1, 0)
end
