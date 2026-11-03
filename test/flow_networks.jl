using FoodWebs
using Test
using LinearAlgebra: I

# A cycle 1 → 2 → 3 → 1 through the detritus compartment 3: compartment 1 takes an input of 100
# and the 20 that 3 returns, 2 respires half its throughflow, and 3 exports 40
cycle_network() = FlowNetwork(
    [0 120 0; 0 0 60; 20 0 0]; inputs = [100, 0, 0], exports = [0, 0, 40], respiration = [0, 60, 0],
    living = [true, true, false]
)

@testset "flow indices by hand" begin
    net = cycle_network()
    @test throughflow(net) == [120, 120, 60]
    @test total_system_throughflow(net) == 300
    @test total_system_throughput(net) == 400
    @test average_path_length(net) == 3
    # a unit leaving any compartment of the cycle returns with probability 1 · 1/2 · 1/3 = 1/6, so
    # each returns 1 - 1/n_ii = 1/6 of its throughflow, and so does the whole system
    @test finn_cycling_index(net) ≈ 1 / 6

    # the extended matrix has row sums 120, 120, 60 and input 100, and column sums 120, 120, 60,
    # export 40 and respiration 60, of a total of 400
    a = ascendency(net)
    @test a.ascendency ≈ 340 * log2(10 / 3) + 20 * log2(10 / 9) + 40 * log2(20 / 3)
    @test a.capacity ≈ -sum(t * log2(t / 400) for t in (120, 60, 20, 40, 60, 100))
    @test a.overhead ≈ a.capacity - a.ascendency
    @test a.diversity ≈ a.capacity / 400 && a.mutual_information ≈ a.ascendency / 400
    @test 0 < a.ascendency < a.capacity

    # detritus has level 1, so 1 feeding on it is at 2 and 2 at 3
    @test trophic_levels(net) ≈ [2, 3, 1]
    @test omnivory_index(net) ≈ zeros(3) atol = 1.0e-12

    # direct impacts: 1 → 2 benefits 2 by its whole diet and costs 1 all of its production; 3 is
    # dead, so 2's flow to it costs 2 nothing; 3 → 1 is a sixth of 1's diet and a third of 3's
    # production
    Q = [0 1 -1 / 3; -1 0 1; 1 / 6 0 0]
    @test mixed_trophic_impact(net) ≈ inv(I - Q) - I

    # an omnivore 3 takes 10 from the basal 1 and 30 from 2, so its level is 1 + (1 + 3 · 2) / 4
    # and its omnivory index the variance of 1 and 2 with weights 1/4 and 3/4
    net = FlowNetwork(
        [0 40 10; 0 0 30; 0 0 0]; inputs = [100, 0, 0], respiration = [50, 10, 40]
    )
    @test trophic_levels(net) ≈ [1, 2, 2.75]
    @test omnivory_index(net) ≈ [0, 0, 0.1875]
    @test finn_cycling_index(net) ≈ 0 atol = 1.0e-12
    @test average_path_length(net) ≈ 180 / 100
end

@testset "flow indices against enaR and NetIndices" begin
    # the Chesapeake Bay carbon flows (Baird & Ulanowicz 1989), from test/fixtures/make_flow_indices.R
    _, nodes = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "flow_indices_nodes.csv"))
    _, links = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "flow_indices_flows.csv"))
    _, values = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "flow_indices.csv"))
    _, mti = FoodWebs.read_csv(joinpath(@__DIR__, "fixtures", "flow_indices_mti.csv"))
    n = length(nodes)
    column(k) = parse.(Float64, getindex.(nodes, k))
    flows = zeros(n, n)
    for (i, j, f) in links
        flows[parse(Int, i), parse(Int, j)] = parse(Float64, f)
    end
    net = FlowNetwork(
        flows; inputs = column(4), exports = column(5), respiration = column(6),
        living = getindex.(nodes, 3) .== "TRUE", names = getindex.(nodes, 2)
    )
    expected = Dict(k => parse(Float64, v) for (k, v) in values)
    @test total_system_throughflow(net) ≈ expected["TST"] rtol = 1.0e-12
    @test total_system_throughput(net) ≈ expected["TSTp"] rtol = 1.0e-12
    @test average_path_length(net) ≈ expected["APL"] rtol = 1.0e-12
    @test finn_cycling_index(net) ≈ expected["FCI"] rtol = 1.0e-9
    a = ascendency(net)
    @test a.diversity ≈ expected["H"] rtol = 1.0e-12
    @test a.mutual_information ≈ expected["AMI"] rtol = 1.0e-12
    @test a.capacity ≈ expected["CAP"] rtol = 1.0e-12
    @test a.ascendency ≈ expected["ASC"] rtol = 1.0e-12
    @test a.overhead ≈ expected["OH"] rtol = 1.0e-12
    # NetIndices keeps the blue crab's flow to itself in its diet but not in its level, which
    # changes the blue crab (19) and the striped bass (33) that eats it; the levels here leave
    # cannibalism out, so the blue crab is 1 plus the mean level of its other prey
    TL, OI = trophic_levels(net), omnivory_index(net)
    others = setdiff(1:n, [19, 33])
    @test TL[others] ≈ column(7)[others] rtol = 1.0e-10
    @test OI[others] ≈ column(8)[others] atol = 1.0e-10
    prey = [j for j in 1:n if j != 19 && flows[j, 19] > 0]
    @test TL[19] ≈ 1 + sum(flows[prey, 19] .* TL[prey]) / sum(flows[prey, 19])
    @test round(TL[19]; digits = 2) == 3.5 && round(column(7)[19]; digits = 2) == 3.41
    M = zeros(n, n)
    for (i, j, v) in mti
        M[parse(Int, i), parse(Int, j)] = parse(Float64, v)
    end
    @test mixed_trophic_impact(net) ≈ M atol = 1.0e-9
end

@testset "flow network from fluxes" begin
    # the budget of fluxes_from_metabolism: basal species take in what their predators eat, and each
    # consumer respires its metabolic losses and exports what it does not assimilate
    web = Bool[0 1 1 0; 0 0 1 1; 0 0 0 1; 0 0 0 0]
    x, M, e = [100.0, 20.0, 5.0, 1.0], [1.0, 10.0, 100.0, 1000.0], [0.5, 0.8, 0.9, 0.9]
    F = fluxes_from_metabolism(web, x, M; efficiency = e)
    intake = vec(sum(F; dims = 1))
    assimilated = vec(sum(e .* F; dims = 1))
    basal = intake .== 0
    net = FlowNetwork(
        F; inputs = basal .* vec(sum(F; dims = 2)), exports = intake .- assimilated,
        respiration = 0.71 .* M .^ -0.25 .* x .* .!basal, tol = 1.0e-10
    )
    @test throughflow(net) ≈ vec(sum(F; dims = 2)) .+ net.exports .+ net.respiration
    @test finn_cycling_index(net) ≈ 0 atol = 1.0e-12
    # the levels weight each predator's diet by its fluxes, here by prey biomass
    @test trophic_levels(net)[3] ≈ 1 + (100 * 1 + 20 * 2) / 120
end

@testset "flow network argument errors" begin
    @test_throws ArgumentError FlowNetwork([0 1; 0 0]; inputs = [1, 0])
    net = FlowNetwork([0 1; 0 0]; inputs = [1, 0], tol = Inf)
    @test throughflow(net) == [1, 1]
    @test_throws ArgumentError FlowNetwork([0 -1; 0 0]; inputs = [1, 0], tol = Inf)
    @test_throws DimensionMismatch FlowNetwork([0 1; 0 0]; inputs = [1])
    @test_throws DimensionMismatch FlowNetwork([0 1; 0 0]; inputs = 1, living = [true])
    @test_throws DimensionMismatch FlowNetwork([0 1; 0 0]; inputs = 1, names = ["a"])
    @test_throws DimensionMismatch FlowNetwork([0 1 0; 0 0 1]; inputs = 1)
    # within tol of its outflows a compartment is at steady state
    @test FlowNetwork([0 1; 0 0]; inputs = [1, 0], respiration = [0, 1.04]) isa FlowNetwork
    @test_throws ArgumentError FlowNetwork([0 1; 0 0]; inputs = [1, 0], respiration = [0, 1.1])
end
