using FoodWebs
using FoodWebs: edge_crossings
using Test
using Random: Xoshiro
using CairoMakie: CairoMakie

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
end
