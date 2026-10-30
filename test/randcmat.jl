using FoodWebs
using Test
using Random: Xoshiro
using Distributions: Normal, Dirac
using LinearAlgebra: diag

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

    web = niche_network(Xoshiro(9), 30, 0.15)
    @test random_predprey(Xoshiro(1), web, Normal(), 0.5) ==
        random_predprey(Xoshiro(1), web, Normal(), 0.5)
end
