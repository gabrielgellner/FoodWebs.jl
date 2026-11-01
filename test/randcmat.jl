using FoodWebs
using Test
using Random: Xoshiro
using Distributions: Normal, Dirac, Exponential
using LinearAlgebra: I, diag

@testset "random_predprey" begin
    rng = Xoshiro(8)

    # 2 eats 1: the predator's effect on the prey is negative (J[prey, pred] < 0) and the prey's
    # effect on the predator is positive (J[pred, prey] > 0)
    J = random_predprey(rng, PredationMatrix([0 1; 0 0]), Normal(), 0.5)
    @test J[1, 2] < 0 && J[2, 1] > 0
    @test diag(J) == [-1, -1]

    # links below the diagonal are kept: 1 eats 2
    J = random_predprey(rng, PredationMatrix([0 0; 1 0]), Normal(), 0.5)
    @test J[2, 1] < 0 && J[1, 2] > 0

    for _ in 1:10
        web = niche_network(rng, 40, 0.15)
        J = random_predprey(rng, web, Normal(), 0.3)
        for i in 1:40, j in 1:40
            i == j && continue
            if web[i, j] == 1 && web[j, i] == 0
                @test J[i, j] < 0 && J[j, i] > 0
            elseif web[i, j] == 0 && web[j, i] == 0
                @test J[i, j] == 0
            end
        end
        @test all(==(-1), diag(J))

        # shuffling keeps a +/- pair on every link, but not who is the predator
        Js = random_predprey(rng, web, Normal(), 0.3; shuffle_prey = true)
        links = [(i, j) for i in 1:40 for j in 1:40 if i != j && web[i, j] == 1 && web[j, i] == 0]
        @test all(sign(Js[i, j]) == -sign(Js[j, i]) for (i, j) in links)
        @test count(Js[i, j] > 0 for (i, j) in links) > 0
    end

    # f scales the prey -> predator effects: with |dist| == 1 they are exactly f
    J = random_predprey(rng, food_chain_network(4), Dirac(1.0), 0.25)
    @test all(J[i + 1, i] == 0.25 for i in 1:3)
    @test all(J[i, i + 1] == -1 for i in 1:3)

    # two species that eat each other: each entry is a gain minus a loss, f - 1 with unit
    # magnitudes, whichever order the pair is visited in
    J = random_predprey(rng, PredationMatrix([0 1; 1 0]), Dirac(1.0), 0.25)
    @test J[1, 2] == J[2, 1] == 0.25 - 1
    # in niche webs both orientations of a two-way pair now occur, where the old code always
    # kept the lower-indexed species as the predator
    signs = Int[]
    for _ in 1:20
        w = niche_network(rng, 40, 0.2)
        J = random_predprey(rng, w, Normal(), 1.0)
        for i in 1:40, j in (i + 1):40
            w[i, j] && w[j, i] && push!(signs, sign(J[j, i]))
        end
    end
    @test !isempty(signs) && any(==(1), signs) && any(==(-1), signs)
    # webs without two-way pairs are drawn exactly as before
    web = cascade_network(Xoshiro(10), 20, 0.15)
    J = random_predprey(Xoshiro(11), web, Normal(), 0.5)
    expected = -Matrix(1.0I, 20, 20)
    draws = Xoshiro(11)
    for i in 1:20, j in 1:20
        if i != j && web[i, j]
            expected[i, j] = -abs(randn(draws))
            expected[j, i] = 0.5 * abs(randn(draws))
        end
    end
    @test J == expected

    web = niche_network(Xoshiro(9), 30, 0.15)
    @test random_predprey(Xoshiro(1), web, Normal(), 0.5) ==
        random_predprey(Xoshiro(1), web, Normal(), 0.5)
end

@testset "expected_community_matrix" begin
    # 2 eats 1: the predator's effect is -E|a|, the prey's f E|a|; with Dirac magnitudes the
    # expected matrix is the matrix itself
    web = PredationMatrix([0 1 1; 0 0 1; 0 0 0])
    @test expected_community_matrix(web, Dirac(1.0), 0.5) == random_predprey(Xoshiro(1), web, Dirac(1.0), 0.5)
    m = sqrt(2 / π)
    E = expected_community_matrix(PredationMatrix([0 1; 0 0]), Normal(), 0.5)
    @test E ≈ [-1 -m; 0.5m -1] rtol = 1.0e-4
    # species that eat each other get both links' effects; shuffled pairs average them
    @test expected_community_matrix(PredationMatrix([0 1; 1 0]), Normal(), 0.5) ≈ [-1 -0.5m; -0.5m -1] rtol = 1.0e-4
    @test expected_community_matrix(PredationMatrix([0 1; 0 0]), Normal(), 0.5; shuffle_prey = true) ≈
        [-1 -0.25m; -0.25m -1] rtol = 1.0e-4
    # the mean of many draws
    rng = Xoshiro(23)
    web = niche_network(rng, 20, 0.15)
    draws = [random_predprey(rng, web, Exponential(2.0), 0.3) for _ in 1:4000]
    @test maximum(abs.(sum(draws) / length(draws) .- expected_community_matrix(web, Exponential(2.0), 0.3))) < 0.15
    # with f = 1 the expected matrix plus the identity is skew-symmetric
    E = expected_community_matrix(niche_network(rng, 40, 0.15), Normal(), 1.0)
    @test E + I ≈ -(E + I)'
end
