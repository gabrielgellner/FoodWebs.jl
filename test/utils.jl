using FoodWebs
using Test
using Random: Xoshiro, randperm

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

# a web where species at level k > 1 eat only species at level k - 1, so every link spans
# exactly one level: no loops and no omnivory
function layered_web(rng, nlevels, width)
    level = vcat(1, 2:nlevels, rand(rng, 1:nlevels, width))
    S = length(level)
    web = zeros(Bool, S, S)
    for j in 1:S, i in 1:S
        level[i] == level[j] - 1 && rand(rng) < 0.4 && (web[i, j] = true)
    end
    # every consumer needs at least one prey a level below
    for j in 1:S
        level[j] > 1 && !any(web[:, j]) && (web[findfirst(==(level[j] - 1), level), j] = true)
    end
    return web, level
end

@testset "MacKay trophic levels" begin
    # feed-forward motif, MacKay, Johnson & Sansom (2020) eq. 6.4: h = [-2/3, 0, 2/3] up to
    # a shift, and F₀ = 1/9
    ffl = [0 1 1; 0 0 1; 0 0 0]
    @test trophic_levels(ffl; method = :mackay) ≈ [1, 5 / 3, 7 / 3]
    @test trophic_incoherence(ffl) ≈ 1 / 9

    @test trophic_levels(food_chain_network(6); method = :mackay) ≈ 1:6
    @test trophic_incoherence(food_chain_network(6)) ≈ 0 atol = 1.0e-12

    # with no loops and no omnivory the two methods agree, and the web is coherent
    rng = Xoshiro(8)
    for _ in 1:20
        web, level = layered_web(rng, 4, 20)
        @test trophic_levels(web; method = :mackay) ≈ level
        @test trophic_levels(web) ≈ level
        @test trophic_incoherence(web) ≈ 0 atol = 1.0e-10
    end

    # each weakly connected component is shifted so its lowest level is 1, and an isolated
    # species has level 1
    two_chains = [0 1 0 0 0; 0 0 0 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 0 0 0]
    @test trophic_levels(two_chains; method = :mackay) ≈ [1, 2, 1, 2, 3]
    @test trophic_levels(zeros(Bool, 3, 3); method = :mackay) == ones(3)

    # a feeding loop has no net direction: equal levels and F₀ = 1
    @test trophic_levels([0 1; 1 0]; method = :mackay) ≈ [1, 1]
    @test trophic_incoherence([0 1; 1 0]) ≈ 1
    @test trophic_incoherence([0 1 0; 0 0 1; 1 0 0]) ≈ 1

    # cannibalism leaves the levels unchanged and adds 1 to the numerator of F₀
    @test trophic_levels([1 1; 0 0]; method = :mackay) ≈ [1, 2]
    @test trophic_incoherence([1 1; 0 0]) ≈ 1 / 2

    # scaling every weight changes nothing
    @test trophic_levels(3.5 * ffl; method = :mackay) ≈ trophic_levels(ffl; method = :mackay)
    @test trophic_incoherence(3.5 * ffl) ≈ trophic_incoherence(ffl)
    # weight a on the 1 -> 3 link: by symmetry both short links span x, and minimizing
    # 2(x - 1)^2 + a(2x - 1)^2 gives x = (1 + a) / (1 + 2a)
    a = 10
    x = (1 + a) / (1 + 2a)
    @test trophic_levels([0 1 a; 0 0 1; 0 0 0]; method = :mackay) ≈ [1, 1 + x, 1 + 2x]
    @test trophic_incoherence([0 1 a; 0 0 1; 0 0 0]) ≈ (2(x - 1)^2 + a * (2x - 1)^2) / (2 + a)

    # defined on niche webs with loops cut off from basal species, and on May networks, which
    # often have no basal species at all
    rng = Xoshiro(9)
    webs = [niche_network(rng, 40, 0.25) for _ in 1:500]
    unfed = filter(w -> !isempty(disconnected_from_basal(w)), webs)
    @test !isempty(unfed)
    for web in vcat(unfed, [may_network(rng, 40, 0.2) for _ in 1:20])
        @test all(isfinite, trophic_levels(web; method = :mackay))
        @test 0 <= trophic_incoherence(web) <= 1
    end

    # invariant to relabeling species
    for web in (niche_network(rng, 30, 0.15), may_network(rng, 30, 0.15))
        p = randperm(rng, 30)
        @test trophic_levels(web[p, p]; method = :mackay) ≈ trophic_levels(web; method = :mackay)[p]
        @test trophic_incoherence(web[p, p]) ≈ trophic_incoherence(web)
    end

    @test_throws ArgumentError trophic_levels([0 -1; 0 0]; method = :mackay)
    @test_throws ArgumentError trophic_levels(ffl; method = :levine)
    @test_throws ArgumentError trophic_incoherence(zeros(Bool, 3, 3))
    @test_throws DimensionMismatch trophic_levels([0 1 0; 0 0 1]; method = :mackay)
end

@testset "remove cannibalism and orient two-way links" begin
    # 1 and 2 eat each other, 3 eats 1 and itself
    web = Bool[0 1 1; 1 0 0; 0 0 1]
    w, removed = remove_cannibalism(web)
    @test removed == [(3, 3)]
    @test w == Bool[0 1 1; 1 0 0; 0 0 0]
    @test remove_cannibalism(w)[2] == Tuple{Int, Int}[]

    # by index, 2 eats 1 and the link where 1 eats 2 is removed; the cannibal link stays
    w, removed = orient_two_way_links(web; keep = :index)
    @test w == Bool[0 1 1; 0 0 0; 0 0 1]
    @test removed == [(2, 1)]
    w, removed = orient_two_way_links(web; keep = :drop)
    @test w == Bool[0 0 1; 0 0 0; 0 0 1]
    @test sort(removed) == [(1, 2), (2, 1)]
    @test_throws ArgumentError orient_two_way_links(web; keep = :mass)
    @test_throws UndefKeywordError orient_two_way_links(web)

    # up the hierarchy: 1 and 2 eat each other, 3 eats 2 and 4 eats 3, so 2 sits above 1 and
    # keeps eating it
    up = zeros(Bool, 4, 4)
    up[1, 2] = up[2, 1] = up[2, 3] = up[3, 4] = true
    w, removed = orient_two_way_links(up; keep = :up)
    @test removed == [(2, 1)]
    @test w[1, 2] && !w[2, 1]

    rng = Xoshiro(21)
    for _ in 1:20
        nw = niche_network(rng, 40, 0.2)
        pairs = [(i, j) for i in 1:40 for j in (i + 1):40 if nw[i, j] && nw[j, i]]
        oneway = [(i, j) for i in 1:40, j in 1:40 if i != j && nw[i, j] && !nw[j, i]]
        for keep in (:up, :index, :random, :drop)
            w, removed = orient_two_way_links(rng, nw; keep)
            @test w isa PredationMatrix
            # no two-way pairs remain, every one-way link and cannibal link is kept
            @test !any(w[i, j] && w[j, i] for i in 1:40, j in 1:40 if i < j)
            @test all(w[i, j] for (i, j) in oneway)
            @test all(w[i, i] == nw[i, i] for i in 1:40)
            @test length(removed) == (keep === :drop ? 2 : 1) * length(pairs)
            # applying the rule again changes nothing
            @test orient_two_way_links(rng, w; keep)[1] == w
        end
        # in niche order the higher-niche species is kept as the predator
        w, _ = orient_two_way_links(nw; keep = :index)
        @test all(w[i, j] && !w[j, i] for (i, j) in pairs)
    end
    nw = niche_network(Xoshiro(22), 40, 0.2)
    @test orient_two_way_links(Xoshiro(1), nw; keep = :random) == orient_two_way_links(Xoshiro(1), nw; keep = :random)

    # a FoodWeb keeps its species, name and metadata
    fw = foodweb("williams_martinez_2000/coachella_valley")
    w, removed = remove_cannibalism(fw)
    @test w isa FoodWeb && w.species == fw.species && w.name == fw.name
    @test length(removed) == length(cannibals(fw))
    @test isempty(cannibals(w))
end
