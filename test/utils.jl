using FoodWebs
using Test
using Random: Xoshiro

@testset "connectance" begin
    @test connectance([0 1; 0 0]) == 1 / 4
    @test connectance(PredationMatrix([0 1; 0 0])) == 1 / 4
    # cannibalistic links count once
    @test connectance([1 1; 0 0]) == 2 / 4
    @test connectance(zeros(Int, 3, 3)) == 0
    @test_throws DimensionMismatch connectance([0 1 0; 0 0 1])
end

@testset "basal and top species" begin
    # 1 and 2 basal, 3 eats both, 4 eats 3
    web = [0 0 1 0; 0 0 1 0; 0 0 0 1; 0 0 0 0]
    @test basal_species(web) == [1, 2]
    @test top_species(web) == [4]
    # a cannibal is not basal
    @test basal_species([1 1; 0 0]) == Int[]
end

@testset "disconnected from basal" begin
    @test disconnected_from_basal([0 1 0; 0 0 1; 0 0 0]) == Int[]
    # a loop with access to a basal species is fine
    @test disconnected_from_basal([0 1 0; 0 0 1; 0 1 0]) == Int[]
    # 2 and 3 only eat each other, and 1 eats nothing
    @test disconnected_from_basal([0 0 0; 0 0 1; 0 1 0]) == [2, 3]
    # 4 eats the closed loop, so it is cut off too
    @test disconnected_from_basal([0 0 0 0; 0 0 1 1; 0 1 0 0; 0 0 0 0]) == [2, 3, 4]
    @test disconnected_from_basal([0 1; 1 0]) == [1, 2]
    # a species that only eats itself counts as basal
    @test disconnected_from_basal([1 1; 0 0]) == Int[]
    @test disconnected_from_basal(food_chain_network(5)) == Int[]
end

@testset "trophic levels" begin
    linchain = [0 1 0; 0 0 1; 0 0 0]
    @test trophic_levels(linchain) ≈ [1.0, 2.0, 3.0]

    omnchain = [0 1 1; 0 0 1; 0 0 0]
    @test trophic_levels(omnchain) ≈ [1.0, 2.0, 2.5]

    omn_relabel = [0 0 1; 1 0 1; 0 0 0]
    @test trophic_levels(omn_relabel) ≈ [2.0, 1.0, 2.5]

    # cannibalism is ignored
    cannibal = [1 1 0; 0 1 1; 0 0 1]
    @test trophic_levels(cannibal) ≈ [1.0, 2.0, 3.0]
    @test trophic_levels([1 1; 0 0]) ≈ [1.0, 2.0]

    coupled_chain = [0 0 0 0 0; 0 0 1 0 0; 1 0 0 0 0; 1 0 0 0 0; 0 0 0 1 0]
    @test trophic_levels(coupled_chain) ≈ [3, 1, 2, 2, 1]

    coupled_chain_omn = [0 0 0 0 0; 1 0 1 0 0; 1 0 0 0 0; 1 0 0 0 0; 0 0 0 1 0]
    @test trophic_levels(coupled_chain_omn) ≈ [8 / 3, 1, 2, 2, 1]

    # a feeding loop with access to a basal species: 2 eats 1 and 3, 3 eats 2, so
    # t2 = 1 + (1 + t3) / 2 and t3 = 1 + t2
    loop = [0 1 0; 0 0 1; 0 1 0]
    @test trophic_levels(loop) ≈ [1.0, 4.0, 5.0]

    # closed loops with no basal species have undefined trophic levels
    @test_throws ArgumentError trophic_levels([0 1; 1 0])
    @test_throws ArgumentError trophic_levels([0 0 0; 0 0 1; 0 1 0])

    # every species in a generated web is at least level 1, and a predator is always above
    # its lowest prey. Cascade webs are acyclic, so trophic levels are always defined.
    rng = Xoshiro(7)
    for _ in 1:20
        web = generalized_cascade_network(rng, 40, 0.15)
        t = trophic_levels(web)
        @test all(>=(1 - 1.0e-12), t)
        @test all(t[b] ≈ 1 for b in basal_species(web))
        for j in 1:40
            prey = [i for i in 1:40 if web[i, j] == 1 && i != j]
            isempty(prey) || @test t[j] >= 1 + minimum(t[prey]) - 1.0e-9
        end
    end
end
