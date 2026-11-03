using FoodWebs
using Test
using LinearAlgebra: Diagonal
using Random: Xoshiro

# The example webs of fluxweb with its inputs and fluxes, from test/fixtures/make_fluxweb.R
function fluxweb_example(id)
    _, species = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "fluxweb_species.csv"))
    _, links = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "fluxweb_links.csv"))
    rows = filter(r -> r[1] == id, species)
    S = length(rows)
    x, M, loss, efficiency = (parse.(Float64, getindex.(rows, k)) for k in 3:6)
    web, preference, flux = zeros(Bool, S, S), zeros(S, S), zeros(S, S)
    for r in links
        r[1] == id || continue
        i, j = parse(Int, r[2]), parse(Int, r[3])
        web[i, j] = true
        preference[i, j], flux[i, j] = parse(Float64, r[4]), parse(Float64, r[5])
    end
    return (; web, x, M, loss, efficiency, preference, flux)
end

# assimilated intake less predation and metabolic losses, for each consumer
function flux_imbalance(web, F, x, X, e)
    S = size(web, 1)
    E = e isa Real ? fill(e, S, S) : repeat(e, 1, S)
    consumers = [j for j in 1:S if any(web[:, j])]
    return [sum(E[:, j] .* F[:, j]) - X[j] * x[j] - sum(F[j, :]) for j in consumers]
end

@testset "fluxes by hand" begin
    # a chain 1 → 2 → 3 with losses 0.5 per unit biomass and efficiency 0.5: the top predator
    # assimilates its losses, F23 = 0.5 x3 / 0.5, and the middle species its losses and F23
    chain = food_chain_network(3)
    x = [10.0, 4.0, 1.0]
    F = fluxes_from_metabolism(chain, x, ones(3); efficiency = 0.5, intercept = 0.5, exponent = 0.0)
    @test F ≈ [0 6 0; 0 0 1; 0 0 0]
    # allometric losses 0.71 M^-0.25 with M = 16 for the top predator
    F = fluxes_from_metabolism(chain, x, [missing, 1.0, 16.0]; efficiency = 0.5)
    @test F[2, 3] ≈ 0.71 / 2 * x[3] / 0.5
    @test F[1, 2] ≈ (0.71 * x[2] + F[2, 3]) / 0.5

    # an omnivore eating 1 and 2 splits its intake by prey biomass, 10 : 4, or equally
    omnivory = Bool[0 1 1; 0 0 1; 0 0 0]
    F = fluxes_from_metabolism(omnivory, x, ones(3); efficiency = 0.5, intercept = 0.5, exponent = 0.0)
    @test F[1, 3] ≈ 10 / 14 && F[2, 3] ≈ 4 / 14
    @test F[1, 2] ≈ (2 + 4 / 14) / 0.5
    F = fluxes_from_metabolism(omnivory, x, ones(3); efficiency = 0.5, intercept = 0.5, exponent = 0.0, preferences = :equal)
    @test F[1, 3] ≈ 0.5 && F[2, 3] ≈ 0.5
    # efficiencies by prey: the omnivore assimilates 0.5 of species 1 and 0.8 of species 2
    F = fluxes_from_metabolism(omnivory, x, ones(3); efficiency = [0.5, 0.8, 0.9], intercept = 0.5, exponent = 0.0, preferences = :equal)
    @test F[1, 3] ≈ 0.5 / 0.65 / 2 && F[2, 3] ≈ F[1, 3]
    # a preference matrix gives the shares, here 3 : 1
    F = fluxes_from_metabolism(omnivory, x, ones(3); efficiency = 0.5, intercept = 0.5, exponent = 0.0, preferences = [0 1 3; 0 0 1; 0 0 0])
    @test F[1, 3] / F[2, 3] ≈ 3
    # link-specific efficiencies as a matrix match the same efficiencies by prey
    e = [0.5, 0.8, 0.9]
    @test fluxes_from_metabolism(omnivory, x, ones(3); efficiency = repeat(e, 1, 3)) ≈
        fluxes_from_metabolism(omnivory, x, ones(3); efficiency = e)
end

@testset "flux balance" begin
    rng = Xoshiro(73)
    for gen in (cascade_network, generalized_cascade_network, niche_network), _ in 1:5
        # with random biomasses most niche webs have a cannibal that is most of its own diet and
        # cannot assimilate enough to balance; the fluxweb webs test cannibalism
        web = first(remove_cannibalism(gen(rng, 25, 0.12; basal_paths = true)))
        x = exp10.(2 .* rand(rng, 25))
        M = exp10.(4 .* rand(rng, 25))
        e = 0.3 .+ 0.6 .* rand(rng, 25)
        F = try
            fluxes_from_metabolism(web, x, M; efficiency = e)
        catch err
            # about 1 niche web in 30 has a feeding loop that no non-negative fluxes balance
            err isa ArgumentError && gen === niche_network && continue
            rethrow()
        end
        @test all(>=(0), F) && all(F[i, j] == 0 for i in 1:25, j in 1:25 if !web[i, j])
        @test maximum(abs, flux_imbalance(web, F, x, 0.71 .* M .^ -0.25, e)) < 1.0e-8 * maximum(F)
        # each consumer's intake splits by prey biomass
        for j in 1:25
            prey = findall(web[:, j])
            isempty(prey) || @test F[prey, j] ./ sum(F[prey, j]) ≈ x[prey] ./ sum(x[prey])
        end
    end
    # acyclic webs always balance: the system is triangular with a positive diagonal
    for _ in 1:20
        web = cascade_network(rng, 30, 0.15; basal_paths = true)
        F = fluxes_from_metabolism(web, rand(rng, 30) .+ 0.1, rand(rng, 30) .+ 0.1; efficiency = 0.2)
        @test all(>=(0), F)
    end
end

@testset "fluxes against fluxweb" begin
    for id in ("simple_case", "species_level", "groups_level")
        ex = fluxweb_example(id)
        S = length(ex.x)
        # fluxweb weights its preference matrix by prey biomass; simple_case gives losses per unit
        # biomass directly
        F = fluxes_from_metabolism(
            ex.web, ex.x, ones(S); efficiency = ex.efficiency, intercept = ex.loss, exponent = 0.0,
            preferences = ex.preference .* ex.x
        )
        @test F ≈ ex.flux rtol = 1.0e-12
    end
    # the allometric losses of species_level from its body masses, with fluxweb's binary matrix
    ex = fluxweb_example("species_level")
    @test fluxes_from_metabolism(ex.web, ex.x, ex.M; efficiency = ex.efficiency) ≈ ex.flux rtol = 1.0e-12
end

@testset "fixed fluxes in feasible communities" begin
    # species_level has 26 cannibals: with self-regulation only from cannibalism, the community
    # matrix at the fixed fluxes is fluxweb's create.jacob
    ex = fluxweb_example("species_level")
    S = length(ex.x)
    e = ex.efficiency
    fc0 = FeasibleCommunities(ex.web, ex.x; fixed_fluxes = ex.flux, efficiency = e, self_regulation = :none, rates = zeros(S))
    θ = map(coordinates(fc0)) do c
        i, j = c.i, c.j
        c.kind === :loss ? ex.flux[i, j] / (ex.x[i] * ex.x[j]) :
            c.kind === :gain ? e[j] * ex.flux[j, i] / (ex.x[i] * ex.x[j]) :
            (1 - e[i]) * ex.flux[i, i] / ex.x[i]^2
    end
    r = -interaction_matrix(fc0, θ) * ex.x
    fc = FeasibleCommunities(ex.web, ex.x; fixed_fluxes = ex.flux, efficiency = e, self_regulation = :none, rates = r)
    P = halfspaces(fc)
    @test P.C * θ ≈ P.d
    @test fluxes(fc, θ) ≈ ex.flux .* .!Diagonal(trues(S))
    # with self-regulation only from cannibalism, each consumer's growth rate is minus its metabolic
    # loss rate
    consumers = [j for j in 1:S if any(ex.web[:, j])]
    @test r[consumers] ≈ -0.71 .* ex.M[consumers] .^ -0.25
    # against create.jacob, with the abscissa as row 0
    _, rows = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "fluxweb_jacobian.csv"))
    Jr = zeros(S, S)
    abscissa = 0.0
    for row in rows
        i, j, v = parse(Int, row[1]), parse(Int, row[2]), parse(Float64, row[3])
        i == 0 ? (abscissa = v) : (Jr[i, j] = v)
    end
    J = community_matrix(fc, θ)
    @test J ≈ Jr rtol = 1.0e-10
    @test spectral_abscissa(J) ≈ abscissa rtol = 1.0e-6

    # with rate intervals each species keeps its growth rate and self-regulation, tied by its row
    fc = FeasibleCommunities(ex.web, ex.x; fixed_fluxes = ex.flux, efficiency = e, rates = fill((-10.0, 10.0), S))
    @test dimension(fc) == S
    @test isbounded(fc)
    fc = FeasibleCommunities(ex.web, ex.x; fixed_fluxes = ex.flux, efficiency = e, rates = fill((-10.0, Inf), S))
    @test !isbounded(fc)

    # without efficiencies the gains are free: a self-regulated consumer is unbounded unless its
    # self-regulation share is below 1
    chain = food_chain_network(3)
    x = [10.0, 4.0, 1.0]
    F = fluxes_from_metabolism(chain, x, ones(3); efficiency = 0.5)
    @test !isbounded(FeasibleCommunities(chain, x; fixed_fluxes = F))
    @test isbounded(FeasibleCommunities(chain, x; fixed_fluxes = F, self_share = (0.0, 0.5)))
    @test isbounded(FeasibleCommunities(chain, x; fixed_fluxes = F, self_regulation = :basal))

    @test_throws DimensionMismatch FeasibleCommunities(chain, x; fixed_fluxes = zeros(2, 2))
    @test_throws ArgumentError FeasibleCommunities(chain, x; fixed_fluxes = -F)
    @test_throws ArgumentError FeasibleCommunities(chain, x; fixed_fluxes = F .+ 1)
    @test_throws DimensionMismatch FeasibleCommunities(chain, x; efficiency = [0.5, 0.5])
end

@testset "assimilation efficiencies" begin
    @test assimilation_efficiencies([:plant, :animal, :detritus]) == [0.545, 0.906, 0.158]
    @test assimilation_efficiencies([:plant, :animal]; plant = 0.4) == [0.4, 0.906]
    @test_throws ArgumentError assimilation_efficiencies([:fungus])
    @test_throws ArgumentError assimilation_efficiencies([:plant]; plant = 0.0)

    fw = foodweb("cheddar/tl84")
    e = assimilation_efficiencies(fw)
    @test e == [c == "producer" ? 0.545 : 0.906 for c in fw.traits["category"]]
    # CPOM and FPOM have no category and no mass
    e = assimilation_efficiencies(foodweb("cheddar/broadstonestream"))
    @test count(==(0.158), e) == 2
    # the Millstream fungi have no category but a mass
    @test_throws ArgumentError assimilation_efficiencies(foodweb("cheddar/millstream_c4"))
    @test_throws ArgumentError assimilation_efficiencies(foodweb("williams_martinez_2000/skipwith_pond"))

    # fluxes of Tuesday Lake from its biomasses, body masses and categories
    F = fluxes_from_metabolism(fw)
    x, M = biomass(fw), fw.traits["M"]
    @test all(>=(0), F)
    @test maximum(abs, flux_imbalance(fw.web, F, x, 0.71 .* M .^ -0.25, assimilation_efficiencies(fw))) < 1.0e-8 * maximum(F)
end

@testset "flux argument errors" begin
    chain = food_chain_network(3)
    x = [10.0, 4.0, 1.0]
    @test_throws UndefKeywordError fluxes_from_metabolism(chain, x, ones(3))
    @test_throws DimensionMismatch fluxes_from_metabolism(chain, x[1:2], ones(3); efficiency = 0.5)
    @test_throws DimensionMismatch fluxes_from_metabolism(chain, x, ones(2); efficiency = 0.5)
    @test_throws ArgumentError fluxes_from_metabolism(chain, [10.0, 0.0, 1.0], ones(3); efficiency = 0.5)
    @test_throws ArgumentError fluxes_from_metabolism(chain, x, [1.0, missing, 1.0]; efficiency = 0.5)
    @test_throws ArgumentError fluxes_from_metabolism(chain, x, ones(3); efficiency = 1.5)
    @test_throws DimensionMismatch fluxes_from_metabolism(chain, x, ones(3); efficiency = [0.5, 0.5])
    @test_throws DimensionMismatch fluxes_from_metabolism(chain, x, ones(3); efficiency = 0.5, intercept = [1.0])
    @test_throws ArgumentError fluxes_from_metabolism(chain, x, ones(3); efficiency = 0.5, preferences = :size)
    @test_throws ArgumentError fluxes_from_metabolism(chain, x, ones(3); efficiency = 0.5, preferences = ones(3, 3))
    @test_throws ArgumentError fluxes_from_metabolism(chain, x, ones(3); efficiency = 0.5, preferences = zeros(3, 3))
    # a cannibal taking most of its own intake cannot assimilate enough to pay for its losses
    cannibal = Bool[0 1; 0 1]
    @test_throws ArgumentError fluxes_from_metabolism(cannibal, [1.0, 100.0], ones(2); efficiency = 0.5)
    @test_throws ArgumentError fluxes_from_metabolism(foodweb("cheddar/benguela"))
end
