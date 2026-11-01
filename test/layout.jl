using FoodWebs
using FoodWebs: edge_crossings
using Test
using Random: Xoshiro
using CairoMakie: CairoMakie
using Distributions: Uniform

chain(S) = (
    web = zeros(Bool, S, S); for i in 1:(S - 1)
        web[i, i + 1] = true
    end; web
)

@testset "trophic layout" begin
    @testset "heights and bands" begin
        web = niche_network(Xoshiro(3), 15, 0.15; basal_paths = true)
        levels = trophic_levels(web)
        positions = trophic_layout(web)
        @test last.(positions) == levels
        @test all(p -> all(isfinite, p), positions)
        # any vector of levels can be given
        @test last.(trophic_layout(web; levels = 2 .* levels)) == 2 .* levels
        @test_throws DimensionMismatch trophic_layout(web; levels = levels[1:3])
        # a food chain is a vertical line
        @test all(p -> p[1] == 0, trophic_layout(chain(5)))
        # two species on the same level sit side by side, one unit apart
        diamond = zeros(Bool, 4, 4)
        diamond[1, 2] = diamond[1, 3] = diamond[2, 4] = diamond[3, 4] = true
        p = trophic_layout(diamond)
        @test sort([p[2][1], p[3][1]]) == [-0.5, 0.5]
        @test p[1][1] == p[4][1] == 0
        # evenly spaced levels split into bands no taller than `band`, rather than chaining
        # into one band
        @test FoodWebs.trophic_bands([1.0, 1.2, 1.4, 1.6, 1.8], 0.3) == [[1, 2], [3, 4], [5]]
        # two single species in close bands interleave instead of sitting on top of each other
        x = FoodWebs.place!(zeros(2), [[1], [2]], [2.0, 2.2], 0.3)
        @test x == [0.0, 0.5]
        # bands far apart stay centred
        @test FoodWebs.place!(zeros(2), [[1], [2]], [1.0, 2.0], 0.3) == [0.0, 0.0]
    end

    @testset "deterministic and fewer crossings" begin
        for seed in 1:20
            web = niche_network(Xoshiro(seed), 20, 0.12; basal_paths = true)
            @test trophic_layout(web) == trophic_layout(web)
            @test edge_crossings(web, trophic_layout(web)) <= edge_crossings(web, trophic_layout(web; sweeps = 0))
        end
    end

    @testset "crossings" begin
        # two links that form an X cross once, two parallel links not at all
        web = zeros(Bool, 4, 4)
        web[1, 4] = web[2, 3] = true
        @test edge_crossings(web, [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]) == 1
        @test edge_crossings(web, [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]) == 0
    end
end

@testset "layout by group" begin
    # two food chains 1 → 3 → 5 and 2 → 4 → 6, numbered so that the plain layout interleaves them
    web = zeros(Bool, 6, 6)
    web[1, 3] = web[3, 5] = web[2, 4] = web[4, 6] = true
    groups = [1, 2, 1, 2, 1, 2]
    p = trophic_layout(web; groups)
    @test last.(p) == trophic_levels(web)
    # each chain keeps one side in every band, as a vertical channel
    side(i) = sign(p[i][1])
    @test side(1) == side(3) == side(5) && side(2) == side(4) == side(6) && side(1) != side(2)
    # in each band of a larger web, every group's species are contiguous along x
    web = niche_network(Xoshiro(5), 30, 0.12; basal_paths = true)
    groups = compartments(Xoshiro(6), web).groups
    p = trophic_layout(web; groups)
    for band in FoodWebs.trophic_bands(trophic_levels(web), 0.3)
        order = band[sortperm([p[i][1] for i in band])]
        g = groups[order]
        @test all(count(k -> g[k] != g[k + 1] && g[k + 1] == c, 1:(length(g) - 1)) + (g[1] == c) <= 1 for c in unique(g))
    end
    @test_throws DimensionMismatch trophic_layout(web; groups = groups[1:3])
end

@testset "Makie extension" begin
    for web in (
            niche_network(Xoshiro(1), 12, 0.15; basal_paths = true),
            cascade_network(Xoshiro(1), 12, 0.15),
            chain(4),
        )
        fig, ax = foodwebplot(web; title = "web")
        @test fig isa CairoMakie.Figure
        path = tempname() * ".png"
        CairoMakie.save(path, fig)
        @test filesize(path) > 0
    end
    # into an existing axis, with custom colours and no labels
    fig = CairoMakie.Figure()
    ax = CairoMakie.Axis(fig[1, 1])
    @test foodwebplot!(ax, chain(3); color = :tomato, labels = nothing) === ax
    # an empirical web with named species, and a web too large for labels by default
    fig, ax = foodwebplot(foodweb("williams_martinez_2000/chesapeake_bay"))
    @test fig isa CairoMakie.Figure
    fig, ax = foodwebplot(foodweb("williams_martinez_2000/little_rock_lake"))
    @test !any(p -> p isa CairoMakie.Text, ax.scene.plots)

    ext = Base.get_extension(FoodWebs, :FoodWebsMakieExt)
    # arrows end at the node edge in pixels, whatever the aspect ratio: an arrow from (0, 0) to
    # (1, 1) stops 10 pixels short of (1, 1) on a wide and on a tall axis
    for scale in (CairoMakie.Vec2f(100, 10), CairoMakie.Vec2f(10, 100))
        o, d = ext.arrow_geometry([CairoMakie.Point2f(0, 0)], [CairoMakie.Point2f(1, 1)], [5], [10], scale)
        tip = (o[1] .+ d[1]) .* scale
        @test sqrt(sum(abs2, tip .- scale)) ≈ 10 atol = 1.0e-3
        @test sqrt(sum(abs2, o[1] .* scale)) ≈ 5 atol = 1.0e-3
    end
    # nodes closer than their radii get no arrow
    o, _ = ext.arrow_geometry([CairoMakie.Point2f(0, 0)], [CairoMakie.Point2f(0.1, 0)], [5], [10], CairoMakie.Vec2f(100, 100))
    @test isempty(o)
    # labels are white on dark nodes and black on light ones
    @test ext.label_color(:black) == CairoMakie.RGBf(1, 1, 1)
    @test ext.label_color(:yellow) == CairoMakie.RGBf(0, 0, 0)
    # thinner drawing for larger webs
    @test ext.default_markersize(100) < ext.default_markersize(20)
    @test ext.default_edgealpha(1000) < ext.default_edgealpha(50)

    # compartments: groups, links between groups only, roles and a legend
    web = niche_network(Xoshiro(7), 25, 0.15; basal_paths = true)
    groups = compartments(Xoshiro(8), web).groups
    roles = species_roles(web, groups).role
    for kwargs in ((; groups), (; groups, show_within = false), (; groups, roles), (; roles))
        fig, ax = foodwebplot(web; kwargs...)
        path = tempname() * ".png"
        CairoMakie.save(path, fig)
        @test filesize(path) > 0
    end
    @test any(b -> b isa CairoMakie.Legend, foodwebplot(web; groups)[1].content)
    @test !any(b -> b isa CairoMakie.Legend, foodwebplot(web; groups, legend = false)[1].content)
    # single-species groups are grey, the others take the palette in order
    colors = ext.group_colors([1, 1, 2, 3, 3], [:red, :blue])
    @test colors[1] == :red && colors[3] == :blue && colors[2] == ext.SINGLE_COLOR
    @test_throws DimensionMismatch foodwebplot(web; groups = groups[1:3])

    # block-ordered heatmaps of a signed community matrix and of a co-assignment matrix
    J = random_predprey(Xoshiro(9), PredationMatrix(web), Uniform(0, 1), 0.3)
    fig, ax = blockplot(J, groups)
    @test fig isa CairoMakie.Figure
    r = compartments(Xoshiro(10), web; runs = 4)
    fig, ax = blockplot(r.coassignment, r.groups)
    path = tempname() * ".png"
    CairoMakie.save(path, fig)
    @test filesize(path) > 0
    @test_throws DimensionMismatch blockplot(J, groups[1:3])

    # eigenvalues with the predictions of block_spectrum
    M = modular_random_matrix(Xoshiro(11), 200; α = 0.25, C = 0.2, Q = -0.3, μ = -1.0)
    fig, ax = eigenplot(M; groups = [fill(1, 50); fill(2, 150)])
    @test fig isa CairoMakie.Figure
end
