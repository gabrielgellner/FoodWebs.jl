using FoodWebs
using Test
using Random
using Aqua
using Distributions: Normal
using LinearAlgebra: diag
using Statistics: mean, std

Random.seed!(20261030)

@testset "FoodWebs.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(FoodWebs)
    end

    @testset "trophic levels" begin
        linchain = [0 1 0; 0 0 1; 0 0 0]
        @test trophic_levels(linchain) ≈ [1.0, 2.0, 3.0]

        omnchain = [0 1 1; 0 0 1; 0 0 0]
        @test trophic_levels(omnchain) ≈ [1.0, 2.0, 2.5]

        omn_relabel = [0 0 1; 1 0 1; 0 0 0]
        @test trophic_levels(omn_relabel) ≈ [2.0, 1.0, 2.5]

        cannibal = [1 1 0; 0 1 1; 0 0 1]
        @test trophic_levels(cannibal) ≈ [1.0, 2.0, 3.0]

        coupled_chain = [0 0 0 0 0; 0 0 1 0 0; 1 0 0 0 0; 1 0 0 0 0; 0 0 0 1 0]
        @test trophic_levels(coupled_chain) ≈ [3, 1, 2, 2, 1]

        coupled_chain_omn = [0 0 0 0 0; 1 0 1 0 0; 1 0 0 0 0; 1 0 0 0 0; 0 0 0 1 0]
        @test trophic_levels(coupled_chain_omn) ≈ [8 / 3, 1, 2, 2, 1]
    end

    @testset "connectance" begin
        for C in (0.25, 0.53, 0.72)
            @test isapprox(connectance(may_network(500, C)), C, rtol = 2e-2)
        end
        for C in (0.25, 0.57, 0.91)
            @test isapprox(connectance(cascade_network(500, C)), C, rtol = 2e-2)
        end
        # at high C niche ranges are truncated at 1, so the realized connectance falls below
        # the target; only check where the model is well behaved
        # single niche webs vary a lot in connectance, so check the mean over many webs
        for C in (0.11, 0.2, 0.3)
            @test isapprox(mean(connectance(niche_network(200, C)) for _ in 1:50), C, rtol = 3e-2)
        end
    end

    @testset "niche diets are intervals" begin
        # species are sorted by niche value, so each diet is a contiguous run of row indices
        for _ in 1:20
            web = niche_network(50, 0.2)
            for i in 1:50
                prey = findall(==(1), web[:, i])
                isempty(prey) || @test prey == first(prey):last(prey)
            end
        end
    end

    @testset "generalized cascade" begin
        # one Beta(1, β) draw per predator gives diet fractions with far more spread than the
        # plain cascade's binomial ones
        S, C = 200, 0.2
        frac(web) = [sum(web[1:(j - 1), j]) / (j - 1) for j in (S ÷ 2):S]
        @test std(frac(generalized_cascade_network(S, C))) > 3 * std(frac(cascade_network(S, C)))
        @test isapprox(connectance(generalized_cascade_network(500, C)), C, rtol = 0.1)
    end

    @testset "generalized niche" begin
        for diet in (0.0, 0.5, 1.0)
            web = generalized_niche_network(60, 0.2, diet)
            @test all(iszero, diag(web))
            @test !isempty(basal_species(web))
        end
    end

    @testset "random_predprey" begin
        # species 2 (higher index) is eaten by species 1: a lower triangle link
        web = FoodWebs.PredationMatrix([0 0; 1 0])
        cmat = random_predprey(web, Normal(), 0.5)
        @test cmat[1, 2] != 0 && cmat[2, 1] != 0
        @test sign(cmat[1, 2]) == -sign(cmat[2, 1])
    end

    @testset "resilience" begin
        A1 = [-1.0 0; 0 -10]
        A2 = [-1.0 15; 0 -10]
        @test λ_stability(A1) == λ_stability(A2)
        @test ν_stability(A1) ≈ -1.0
        @test isapprox(ν_stability(A2), 3.25, rtol = 1e-2)
    end
end
