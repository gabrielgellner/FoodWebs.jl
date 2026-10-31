# # Structural properties
# The food-web properties of Williams & Martinez (2000), with the conventions that reproduce
# their Table 2 on the same webs, and the characteristic path length of Dunne, Williams &
# Martinez (2002). As elsewhere, `web[i, j] != 0` means species `j` eats species `i`.

"""
    intermediate_species(web)

The species that have both prey and predators, where `web[i, j] != 0` means species `j` eats
species `i`. A cannibalistic link counts as both, so every species is exactly one of
[`top_species`](@ref), intermediate or [`basal_species`](@ref), except a species that eats
nothing and is eaten by nothing.
"""
function intermediate_species(web::AbstractMatrix)
    return [i for i in axes(web, 1) if any(!iszero, view(web, i, :)) && any(!iszero, view(web, :, i))]
end

"""
    cannibals(web)

The species that eat themselves, where `web[i, j] != 0` means species `j` eats species `i`.
"""
cannibals(web::AbstractMatrix) = [i for i in axes(web, 1) if web[i, i] != 0]

"""
    species_in_loops(web)

The species with a food chain that passes through a feeding loop longer than a cannibalistic
link: the species in such a loop and the species that eat down into one (Williams & Martinez
2000). A food chain runs from a species down its prey links to a basal species.
"""
function species_in_loops(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    # edges from consumer to prey, without cannibalism
    g = DiGraph(S)
    for i in 1:S, j in 1:S
        i != j && web[i, j] != 0 && add_edge!(g, j, i)
    end
    looped = falses(S)
    for c in strongly_connected_components(g)
        length(c) > 1 && (looped[c] .= true)
    end
    # a species whose chains reach a loop: anything that reaches a looped species
    reaches = [any(looped[v] for v in vertices(g) if has_path(g, s, v)) for s in 1:S]
    return findall(reaches)
end

# Whether the loop-free food chains from `s` down to basal species have lengths other than
# `seen`, stopping at the first second length
function second_chain_length(web, s, onpath, depth, basal, seen)
    if basal[s]
        seen[] == -1 && (seen[] = depth)
        return seen[] != depth
    end
    onpath[s] = true
    for p in axes(web, 1)
        (p == s || onpath[p] || web[p, s] == 0) && continue
        second_chain_length(web, p, onpath, depth + 1, basal, seen) && (onpath[s] = false; return true)
    end
    onpath[s] = false
    return false
end

"""
    omnivores(web)

The species that eat two or more species and have food chains of different lengths
(Williams & Martinez 2000), where a food chain is a loop-free path down the prey links to a
basal species. Cannibalistic links are ignored.
"""
function omnivores(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    eats(pred, prey) = pred != prey && web[prey, pred] != 0
    basal = [!any(eats(j, i) for i in 1:S) for j in 1:S]
    off = [eats(j, i) for i in 1:S, j in 1:S]
    return [
        s for s in 1:S if count(i -> eats(s, i), 1:S) >= 2 &&
            second_chain_length(off, s, falses(S), 0, basal, Ref(-1))
    ]
end

"""
    generality(web)

The normalized number of prey of each species, `G[j] = (S / L) Σ_i web[i, j]`, which has mean 1
(Williams & Martinez 2000).
"""
function generality(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    return vec(count(!iszero, web; dims = 1)) .* (S / count(!iszero, web))
end

"""
    vulnerability(web)

The normalized number of predators of each species, `V[i] = (S / L) Σ_j web[i, j]`, which has
mean 1 (Williams & Martinez 2000).
"""
function vulnerability(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    return vec(count(!iszero, web; dims = 2)) .* (S / count(!iszero, web))
end

"""
    trophic_similarity(web)

The matrix of trophic similarities `s[i, j]`: the number of predators and prey species `i` and
`j` share, divided by their total number of distinct predators plus distinct prey (Williams &
Martinez 2000). Two species with no predators or prey have similarity 0.
"""
function trophic_similarity(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    A = web .!= 0
    s = zeros(S, S)
    for i in 1:S, j in 1:S
        shared = count(view(A, :, i) .& view(A, :, j)) + count(view(A, i, :) .& view(A, j, :))
        total = count(view(A, :, i) .| view(A, :, j)) + count(view(A, i, :) .| view(A, j, :))
        s[i, j] = total == 0 ? 0.0 : shared / total
    end
    return s
end

"""
    mean_max_similarity(web)

The mean over species of each species' largest [`trophic_similarity`](@ref) to another species
(MxSim, Williams & Martinez 2000).
"""
function mean_max_similarity(web::AbstractMatrix)
    s = trophic_similarity(web)
    S = size(s, 1)
    return mean(maximum(s[i, j] for j in 1:S if j != i) for i in 1:S)
end

"""
    characteristic_path_length(web)

The mean number of links in the shortest path between two species, with links taken as
undirected and the mean taken over all `S^2` ordered pairs, including each species with itself
at distance 0 (Dunne, Williams & Martinez 2002). This reproduces their values for the
Williams & Martinez (2000) webs. Throws an `ArgumentError` if the web is disconnected.
"""
function characteristic_path_length(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    g = Graph(S)
    for i in 1:S, j in 1:S
        i != j && web[i, j] != 0 && add_edge!(g, i, j)
    end
    is_connected(g) || throw(ArgumentError("the web is disconnected, so path lengths are undefined"))
    return sum(sum(gdistances(g, s)) for s in 1:S) / S^2
end

"""
    structural_properties(web)

The structural properties of Williams & Martinez (2000) and the characteristic path length of
Dunne, Williams & Martinez (2002), as a named tuple:

- `T`, `I`, `B`: the fractions of top, intermediate and basal species
- `GenSD`, `VulSD`: the sample standard deviations of [`generality`](@ref) and
  [`vulnerability`](@ref)
- `MxSim`: [`mean_max_similarity`](@ref)
- `Loop`, `Cannib`, `Omniv`: the fractions of [`species_in_loops`](@ref), [`cannibals`](@ref)
  and [`omnivores`](@ref)
- `TL`, `MaxTL`: the mean and maximum short-weighted [`trophic_levels`](@ref), or `NaN` when
  some species have no feeding path from a basal species
- `PathLength`: [`characteristic_path_length`](@ref), or `NaN` when the web is disconnected

Species with no links at all, which some empirical webs include, count as basal and top.

On the seven Williams & Martinez webs these reproduce their Table 2 to its two significant
figures for `T`, `I`, `B`, `Loop` and `Cannib`, and for most webs `GenSD`, `VulSD`, `MxSim`
and `Omniv` (see the tests for the exceptions).

Williams, R. J. & Martinez, N. D. Simple rules yield complex food webs. Nature 404, 180–183
(2000).

Dunne, J. A., Williams, R. J. & Martinez, N. D. Food-web structure and network theory: the
role of connectance and size. PNAS 99, 12917–12922 (2002).
"""
function structural_properties(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    tl = isempty(disconnected_from_basal(web)) ? trophic_levels(web; method = :short_weighted) : [NaN]
    return (
        T = length(top_species(web)) / S,
        I = length(intermediate_species(web)) / S,
        B = length(basal_species(web)) / S,
        GenSD = std(generality(web)),
        VulSD = std(vulnerability(web)),
        MxSim = mean_max_similarity(web),
        Loop = length(species_in_loops(web)) / S,
        Cannib = length(cannibals(web)) / S,
        Omniv = length(omnivores(web)) / S,
        TL = mean(tl),
        MaxTL = maximum(tl),
        PathLength = is_connected(DiGraph(web .!= 0)) ? characteristic_path_length(web) : NaN,
    )
end
