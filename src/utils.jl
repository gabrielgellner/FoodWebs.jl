"""
    connectance(web)

The connectance `L / S^2` of the `S` species food web `web`, where `L` is the number of links.
"""
function connectance(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    return count(!iszero, web) / S^2
end

"""
    basal_species(web)

The species that eat nothing, where `web[i, j] != 0` means species `j` eats species `i`.
"""
basal_species(web::AbstractMatrix) = [j for j in axes(web, 2) if all(iszero, view(web, :, j))]

"""
    top_species(web)

The species that nothing eats, where `web[i, j] != 0` means species `j` eats species `i`.
"""
top_species(web::AbstractMatrix) = [i for i in axes(web, 1) if all(iszero, view(web, i, :))]

"""
    disconnected_from_basal(web)

The species in `web` with no path of feeding links up from a basal species, where
`web[i, j] != 0` means species `j` eats species `i`. These are species in closed loops that
only eat each other, so no energy can flow to them. Cannibalistic links are ignored, so a
species that only eats itself counts as basal.
"""
function disconnected_from_basal(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    eats(pred, prey) = pred != prey && web[prey, pred] != 0
    basal = [j for j in 1:S if !any(eats(j, i) for i in 1:S)]
    # search up the web from the basal species along prey -> predator links
    reached = falses(S)
    reached[basal] .= true
    frontier = copy(basal)
    while !isempty(frontier)
        prey = pop!(frontier)
        for pred in 1:S
            if eats(pred, prey) && !reached[pred]
                reached[pred] = true
                push!(frontier, pred)
            end
        end
    end
    return findall(!, reached)
end

"""
    trophic_levels(web; method = :prey_averaged)

The trophic levels of the species in `web`, where `web[i, j] != 0` means species `j` eats
species `i`.

With `method = :prey_averaged`, basal species have trophic level 1 and every other species
has 1 plus the mean trophic level of its prey (Levine 1980), that is, the solution of
`(I - P) t = 1` where `P` is the diet matrix with predators as rows normalized to sum to 1.
Links are unweighted and cannibalistic links are ignored. Throws an `ArgumentError` if some
species have no path of feeding links up from a basal species (see
[`disconnected_from_basal`](@ref)), as their trophic levels are then undefined.

With `method = :mackay`, the levels `h` solve `Λ h = v` (MacKay, Johnson & Sansom 2020), where
`W = web` holds the link weights in the direction of energy flow (prey to predator),
`Λ = Diagonal(u) - W - W'`, `u = w_in + w_out` and `v = w_in - w_out`, with `w_in` and `w_out`
the column and row sums of `W`. This is a least-squares fit of `h[j] - h[i] = 1` over the
links, so it is defined for every web, with or without basal species. The levels are unique
up to a constant within each weakly connected component; each component is shifted so its
lowest level is 1. With this shift the two methods agree on webs with no loops and no
omnivory, where every link spans exactly one level. `web` may be any non-negative weighted matrix, and
cannibalistic links do not change the levels.

With `method = :shortest`, a species' level is 1 plus the number of links in its shortest food
chain down to a basal species, and with `method = :short_weighted` it is the mean of the
shortest and prey-averaged levels (Williams & Martinez 2004). Like `:prey_averaged`, these
ignore cannibalistic links and throw an `ArgumentError` for species with no feeding path up
from a basal species.

MacKay, R. S., Johnson, S. & Sansom, B. How directed is a directed network? R. Soc. Open
Sci. 7, 201138 (2020).

Williams, R. J. & Martinez, N. D. Limits to trophic levels and omnivory in complex food webs:
theory and data. Am. Nat. 163, 458–468 (2004).
"""
function trophic_levels(web::AbstractMatrix; method::Symbol = :prey_averaged)
    method === :prey_averaged && return prey_averaged_trophic_levels(web)
    method === :mackay && return mackay_trophic_levels(web)
    method === :shortest && return shortest_trophic_levels(web)
    method === :short_weighted && return (shortest_trophic_levels(web) + prey_averaged_trophic_levels(web)) / 2
    throw(
        ArgumentError(
            "method must be :prey_averaged, :short_weighted, :shortest or :mackay, given :$method"
        )
    )
end

function check_fed(web)
    unfed = disconnected_from_basal(web)
    isempty(unfed) || throw(
        ArgumentError(
            "species $unfed have no feeding path up from a basal species, so their trophic " *
                "levels are undefined"
        )
    )
    return nothing
end

function shortest_trophic_levels(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    check_fed(web)
    eats(pred, prey) = pred != prey && web[prey, pred] != 0
    level = zeros(Int, S)
    frontier = [j for j in 1:S if !any(eats(j, i) for i in 1:S)]
    level[frontier] .= 1
    # breadth-first search up the web from the basal species
    while !isempty(frontier)
        next = Int[]
        for prey in frontier, pred in 1:S
            if level[pred] == 0 && eats(pred, prey)
                level[pred] = level[prey] + 1
                push!(next, pred)
            end
        end
        frontier = next
    end
    return Float64.(level)
end

function prey_averaged_trophic_levels(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    check_fed(web)
    # diet matrix with predators as rows and without cannibalism
    P = Matrix{Float64}(transpose(web) .!= 0)
    P[diagind(P)] .= 0
    P ./= max.(sum(P, dims = 2), 1)
    return (I - P) \ ones(S)
end

function mackay_trophic_levels(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    all(>=(0), web) || throw(ArgumentError("link weights must be non-negative"))
    W = Matrix{Float64}(web)
    w_in = vec(sum(W, dims = 1))
    w_out = vec(sum(W, dims = 2))
    Λ = Diagonal(w_in + w_out) - W - W'
    v = w_in - w_out
    h = ones(S)
    # Λ is singular, with the constant vector on each weakly connected component in its null
    # space. Fix the first species of each component at 0, solve for the rest, then shift.
    for comp in weakly_connected_components(DiGraph(W .!= 0))
        rest = comp[2:end]
        h[comp[1]] = 0
        h[rest] = Λ[rest, rest] \ v[rest]
        h[comp] .+= 1 - minimum(h[comp])
    end
    return h
end

"""
    trophic_incoherence(web)

The trophic incoherence `F₀` of `web` (MacKay, Johnson & Sansom 2020), the weighted mean of
`(h[j] - h[i] - 1)^2` over the links `web[i, j] != 0`, with `h` the levels from
`trophic_levels(web; method = :mackay)`. `F₀` lies in `[0, 1]`: it is 0 when every link spans
exactly one level, as in a food chain, and 1 when every link joins species at the same level,
as in a web with no net direction such as a feeding loop. Each cannibalistic link adds its
weight to the numerator. `web` may be any non-negative weighted matrix with at least one link.

MacKay, R. S., Johnson, S. & Sansom, B. How directed is a directed network? R. Soc. Open
Sci. 7, 201138 (2020).
"""
function trophic_incoherence(web::AbstractMatrix)
    h = mackay_trophic_levels(web)
    total = sum(web)
    total > 0 || throw(ArgumentError("the web has no links, so its trophic incoherence is undefined"))
    return sum(web[i, j] * (h[j] - h[i] - 1)^2 for i in axes(web, 1), j in axes(web, 2)) / total
end

# A web of the same type as `web` with the links of the Bool matrix `A`
rebuild_web(web::PredationMatrix, A) = PredationMatrix(A)
rebuild_web(web::AbstractMatrix, A) = Matrix{Bool}(A)

"""
    remove_cannibalism(web)

`web` without its cannibalistic links, and the removed links as `(i, i)` pairs. The result has
the type of `web`: a [`PredationMatrix`](@ref), a [`FoodWeb`](@ref) keeping its species and
metadata, or a `Matrix{Bool}`. Removing links lowers connectance.
"""
function remove_cannibalism(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    A = Matrix{Bool}(web .!= 0)
    removed = [(i, i) for i in 1:S if A[i, i]]
    foreach(((i, _),) -> A[i, i] = false, removed)
    return rebuild_web(web, A), removed
end

"""
    orient_two_way_links([rng], web; keep)

`web` with one link kept from each pair of species that eat each other, and the removed links
as `(prey, predator)` pairs, where `web[i, j] != 0` means species `j` eats species `i`. `keep`
chooses the link:

- `:up`: the species at the higher trophic level eats the other, by the MacKay levels of `web`,
  which are defined for every web; ties go to the higher index;
- `:index`: the higher-indexed species eats the other, the order of cascade and niche webs;
- `:random`: a fair coin for each pair, from `rng`;
- `:drop`: both links are removed.

No rule for this is established in the literature. Two-way links arise from life stages that
eat each other (intraguild predation), from aggregating taxa into trophic species, and by design
in the niche model, so each rule changes the web; comparing results under several rules, and
against keeping both links, shows how much the choice matters. Removing links lowers
connectance, and can disconnect species or cut them off from basal species (see
[`disconnected_from_basal`](@ref)). Cannibalistic links are left as they are (see
[`remove_cannibalism`](@ref)). The result has the type of `web`.
"""
function orient_two_way_links(rng::AbstractRNG, web::AbstractMatrix; keep::Symbol)
    keep in (:up, :index, :random, :drop) ||
        throw(ArgumentError("keep must be :up, :index, :random or :drop, given :$keep"))
    S = LinearAlgebra.checksquare(web)
    A = Matrix{Bool}(web .!= 0)
    levels = keep === :up ? trophic_levels(web; method = :mackay) : nothing
    removed = Tuple{Int, Int}[]
    for i in 1:S, j in (i + 1):S
        A[i, j] && A[j, i] || continue
        # whether j, the higher-indexed species, is kept as the predator of i
        j_eats_i = if keep === :up
            levels[j] >= levels[i]
        elseif keep === :index
            true
        elseif keep === :random
            rand(rng, Bool)
        end
        if keep === :drop
            A[i, j] = A[j, i] = false
            push!(removed, (i, j), (j, i))
        elseif j_eats_i
            A[j, i] = false
            push!(removed, (j, i))
        else
            A[i, j] = false
            push!(removed, (i, j))
        end
    end
    return rebuild_web(web, A), removed
end
function orient_two_way_links(web::AbstractMatrix; keep::Symbol)
    return orient_two_way_links(Random.default_rng(), web; keep)
end
