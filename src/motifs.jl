# # Network motifs
# Motifs are the connected subgraphs of `k` species, counted up to relabeling. A motif is
# identified by the mfinder ID of Milo et al. (2002), used by mfinder and pymfinder: write the
# motif's `k`×`k` adjacency `A`, with `A[i, j] = 1` when species `i` eats species `j` (the
# transpose of `web`), as a binary number whose bit `(i - 1)k + (j - 1)` is `A[i, j]`, and take
# the smallest such number over all relabelings. Cannibalistic links are left out.

# The motif ID of each `k`×`k` adjacency code, or 0 for a disconnected subgraph, built once
const MOTIF_TABLES = Dict{Int, Vector{Int}}()

permutations_of(n) = n == 1 ? [[1]] : [vcat(p[1:(k - 1)], n, p[k:end]) for p in permutations_of(n - 1) for k in 1:n]

adjacency_code(A, k) = sum(Int(A[i, j]) << ((i - 1) * k + (j - 1)) for i in 1:k, j in 1:k)

function weakly_connected_pattern(A, k)
    seen = falses(k)
    seen[1] = true
    stack = [1]
    while !isempty(stack)
        v = pop!(stack)
        for w in 1:k
            if !seen[w] && (A[v, w] != 0 || A[w, v] != 0)
                seen[w] = true
                push!(stack, w)
            end
        end
    end
    return all(seen)
end

function motif_table(k::Int)
    k in (3, 4) || throw(ArgumentError("motif size must be 3 or 4, given $k"))
    return get!(MOTIF_TABLES, k) do
        P = permutations_of(k)
        slots = [(i, j) for i in 1:k for j in 1:k if i != j]
        table = zeros(Int, 2^(k^2))
        for pattern in 0:(2^length(slots) - 1)
            A = zeros(Int, k, k)
            for (b, (i, j)) in enumerate(slots)
                A[i, j] = (pattern >> (b - 1)) & 1
            end
            weakly_connected_pattern(A, k) || continue
            table[adjacency_code(A, k) + 1] = minimum(adjacency_code(A[p, p], k) for p in P)
        end
        table
    end
end

"""
    motif_ids(size)

The mfinder IDs of the connected motifs of `size` species, in increasing order: 13 for
`size = 3` and 199 for `size = 4`. [`motif_census`](@ref) counts motifs in this order.
"""
motif_ids(size::Int) = sort!(unique!(filter(!iszero, motif_table(size))))

# Stouffer et al. (2007) names of the three-species motifs: S1 food chain, S2 omnivory, S3
# loop, S4 exploitative competition (two predators sharing a prey), S5 apparent competition (a
# predator with two prey), and D1-D8 with two-way links
const STOUFFER_NAMES = Dict(
    12 => "S1", 38 => "S2", 98 => "S3", 36 => "S4", 6 => "S5", 46 => "D1", 108 => "D2",
    14 => "D3", 74 => "D4", 102 => "D5", 238 => "D6", 110 => "D7", 78 => "D8",
)

"""
    motif_name(id)

The name of the motif with mfinder ID `id`: the Stouffer et al. (2007) names for the 13
three-species motifs (`"S1"` food chain, `"S2"` omnivory, `"S3"` loop, `"S4"` exploitative
competition, `"S5"` apparent competition, `"D1"` to `"D8"` with two-way links), and `"M\$id"`
for the others. The four-species diamond, two consumers sharing a prey and a predator, is
`"M904"`.
"""
motif_name(id::Integer) = get(STOUFFER_NAMES, Int(id), "M$id")

# predator -> prey adjacency without cannibalism: A[i, j] when species i eats species j
function eats_matrix(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    return [i != j && web[j, i] != 0 for i in 1:S, j in 1:S]
end

# Grow `members` by species from `extension` until it has `k`, calling `visit` on each full set.
# `covered[u]` counts how many members `u` is or is next to; exclusive neighbours of a new species
# are those above `root` with no count (Wernicke 2006).
function esu_extend!(visit, members, extension, root, k, neighbours, covered)
    if length(members) == k
        visit(members)
        return nothing
    end
    extension = copy(extension)
    while !isempty(extension)
        w = pop!(extension)
        added = [u for u in neighbours[w] if u > root && covered[u] == 0]
        push!(members, w)
        covered[w] += 1
        foreach(u -> covered[u] += 1, neighbours[w])
        esu_extend!(visit, members, vcat(extension, added), root, k, neighbours, covered)
        foreach(u -> covered[u] -= 1, neighbours[w])
        covered[w] -= 1
        pop!(members)
    end
    return nothing
end

# Call `visit(members)` for every weakly connected set of `k` species, each once, with the ESU
# algorithm of Wernicke (2006). `members` is reused between calls.
function each_connected_subset(visit, A::AbstractMatrix{Bool}, k::Int)
    S = size(A, 1)
    neighbours = [[u for u in 1:S if u != v && (A[v, u] || A[u, v])] for v in 1:S]
    members = Int[]
    covered = zeros(Int, S)
    for v in 1:S
        push!(members, v)
        covered[v] += 1
        foreach(u -> covered[u] += 1, neighbours[v])
        esu_extend!(visit, members, [u for u in neighbours[v] if u > v], v, k, neighbours, covered)
        foreach(u -> covered[u] -= 1, neighbours[v])
        covered[v] -= 1
        pop!(members)
    end
    return nothing
end

function subset_id(A, members, k, table)
    code = 0
    for (a, i) in enumerate(members), (b, j) in enumerate(members)
        A[i, j] && (code |= 1 << ((a - 1) * k + (b - 1)))
    end
    return table[code + 1]
end

"""
    motif_census(web; size = 3)

The number of times each connected motif of `size` species (3 or 4) occurs in `web` as an
induced subgraph, in the order of [`motif_ids`](@ref)`(size)`, where `web[i, j] != 0` means
species `j` eats species `i`. Cannibalistic links are left out. The counts sum to the number of
weakly connected sets of `size` species.
"""
function motif_census(web::AbstractMatrix; size::Int = 3)
    table = motif_table(size)
    ids = motif_ids(size)
    index = Dict(id => n for (n, id) in enumerate(ids))
    A = eats_matrix(web)
    counts = zeros(Int, length(ids))
    each_connected_subset(A, size) do members
        counts[index[subset_id(A, members, size, table)]] += 1
    end
    return counts
end

# The 30 positions of species in the three-species motifs (Stouffer et al. 2012), as
# (motif ID, predators, prey) within the motif, in pymfinder's order
const ROLE_KEYS = [
    (12, 0, 1), (12, 1, 1), (12, 1, 0), (38, 0, 2), (38, 1, 1), (38, 2, 0), (98, 1, 1),
    (36, 0, 1), (36, 2, 0), (6, 0, 2), (6, 1, 0), (46, 1, 2), (46, 2, 0), (108, 0, 2),
    (108, 2, 1), (14, 1, 2), (14, 1, 1), (14, 1, 0), (74, 0, 1), (74, 1, 1), (74, 2, 1),
    (102, 1, 1), (102, 2, 1), (102, 1, 2), (238, 2, 2), (110, 1, 2), (110, 2, 1), (110, 2, 2),
    (78, 1, 1), (78, 2, 2),
]

"""
    motif_roles(web)

The number of times each species occupies each of the 30 positions in the three-species motifs
(Stouffer et al. 2012), as an `S`×30 matrix. Position `r` is `FoodWebs.ROLE_KEYS[r]`: the
motif's mfinder ID and the species' numbers of predators and prey within the motif, in the order
used by pymfinder. For example position 1 is the top of a food chain (`S1`, no predators, one
prey) and position 3 its bottom. Each motif occurrence adds one count to each of its three
species, so the matrix sums to three times the number of connected sets of three species.
"""
function motif_roles(web::AbstractMatrix)
    table = motif_table(3)
    A = eats_matrix(web)
    position = Dict(key => r for (r, key) in enumerate(ROLE_KEYS))
    roles = zeros(Int, size(A, 1), length(ROLE_KEYS))
    each_connected_subset(A, 3) do members
        id = subset_id(A, members, 3, table)
        for m in members
            predators = count(o -> o != m && A[o, m], members)
            prey = count(o -> o != m && A[m, o], members)
            roles[m, position[(id, predators, prey)]] += 1
        end
    end
    return roles
end

"""
    rewire_links([rng], web; swaps = 10 * L)

A randomized copy of `web` with the same number of prey and predators for every species and the
same numbers of one-way and two-way links per species (Milo et al. 2002; Stouffer et al. 2007),
for the null model of [`motif_zscores`](@ref). Each of `swaps` attempts picks two one-way links
`a → b` and `c → d` (`a` eats `b`) and rewires them to `a → d` and `c → b`, or two two-way links
and rewires them likewise, rejecting swaps that would duplicate a link or join two linked
species. Cannibalistic links are kept as they are. `L` is the number of links.
"""
function rewire_links(rng::AbstractRNG, web::AbstractMatrix; swaps::Integer = 10 * count(!iszero, web))
    S = LinearAlgebra.checksquare(web)
    A = eats_matrix(web)
    single = [(a, b) for a in 1:S for b in 1:S if A[a, b] && !A[b, a]]
    double = [(a, b) for a in 1:S for b in (a + 1):S if A[a, b] && A[b, a]]
    linked(x, y) = A[x, y] || A[y, x]
    for _ in 1:swaps
        n = length(single) + length(double)
        n < 2 && break
        two_way = rand(rng) < length(double) / n
        links = two_way ? double : single
        length(links) < 2 && continue
        p, q = rand(rng, 1:length(links)), rand(rng, 1:length(links))
        p == q && continue
        (a, b), (c, d) = links[p], links[q]
        # a two-way link has no direction, so either pairing is a valid swap
        two_way && rand(rng, Bool) && ((c, d) = (d, c))
        length(Set((a, b, c, d))) == 4 || continue
        (linked(a, d) || linked(c, b)) && continue
        A[a, b] = A[c, d] = false
        A[a, d] = A[c, b] = true
        if two_way
            A[b, a] = A[d, c] = false
            A[d, a] = A[b, c] = true
            links[p], links[q] = minmax(a, d), minmax(c, b)
        else
            links[p], links[q] = (a, d), (c, b)
        end
    end
    # back to web[prey, predator], with the original cannibalistic links
    rewired = [A[j, i] for i in 1:S, j in 1:S]
    for i in 1:S
        rewired[i, i] = web[i, i] != 0
    end
    return PredationMatrix(rewired)
end
rewire_links(web::AbstractMatrix; kwargs...) = rewire_links(Random.default_rng(), web; kwargs...)

"""
    motif_zscores([rng], web; size = 3, samples = 100, swaps = 10 * L)

How over- or under-represented each motif of `size` species is in `web` compared with
`samples` randomized webs from [`rewire_links`](@ref). Returns a named tuple with the motif
`ids`, the `observed` counts, the `mean` and `sd` of the randomized counts, the z-scores
`z = (observed - mean) / sd`, which are `NaN` for motifs whose randomized count never varies,
and the normalized `profile`: the z-scores scaled to length 1, with `NaN` taken as 0 (Milo et al.
2004). Profiles compare webs of different sizes, as in Stouffer et al. (2007), and the dot
product of two profiles is their uncentred correlation.

Milo, R. et al. Superfamilies of evolved and designed networks. Science 303, 1538–1542 (2004).

Stouffer, D. B., Camacho, J., Jiang, W. & Amaral, L. A. N. Evidence for the existence of a
robust pattern of prey selection in food webs. Proc. R. Soc. B 274, 1931–1940 (2007).
"""
function motif_zscores(
        rng::AbstractRNG, web::AbstractMatrix; size::Int = 3, samples::Integer = 100,
        swaps::Integer = 10 * count(!iszero, web)
    )
    observed = motif_census(web; size)
    random = reduce(hcat, (motif_census(rewire_links(rng, web; swaps); size) for _ in 1:samples))
    μ = vec(mean(random; dims = 2))
    σ = vec(std(random; dims = 2))
    z = [s > 0 ? (o - m) / s : NaN for (o, m, s) in zip(observed, μ, σ)]
    finite = replace(z, NaN => 0.0)
    length_z = sqrt(sum(abs2, finite))
    profile = length_z > 0 ? finite ./ length_z : finite
    return (ids = motif_ids(size), observed, mean = μ, sd = σ, z, profile)
end
motif_zscores(web::AbstractMatrix; kwargs...) = motif_zscores(Random.default_rng(), web; kwargs...)
