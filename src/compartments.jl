# # Compartments
# Compartments are groups of species that interact more among themselves than with others. As
# elsewhere, `W[i, j]` is the link from prey `i` to predator `j`: 1 or 0 in a web, or a
# non-negative weight such as an energy flux. Cannibalistic links are left out.

"""
    web_modularity(W, groups; directed = true, γ = 1)

The modularity of the partition `groups` (one label per species) of the web or non-negative
weight matrix `W`, where `W[i, j]` is the link from prey `i` to predator `j`, such as an energy
flux. Cannibalistic links (the diagonal) are left out.

With `directed = true`, the directed modularity of Leicht & Newman (2008), weighted as in
Arenas et al. (2007),

    Q = (1/W) Σ_ij [W_ij − γ s_i^out s_j^in / W] δ(g_i, g_j),

where `s_i^out = Σ_j W_ij`, `s_j^in = Σ_i W_ij` and `W = Σ_ij W_ij`; `Q` is the same for `W` and its
transpose, so the direction convention does not matter. With `directed = false`, the modularity
of Newman & Girvan (2004) on the undirected web `W + Wᵀ`,

    Q = (1/2m) Σ_ij [A_ij − γ k_i k_j / 2m] δ(g_i, g_j),

in which each link counts once, so that two species that eat each other are joined by two links.
`γ` is the resolution of Reichardt & Bornholdt (2006): larger values favour smaller groups. `Q`
is 0 for a single group and for a web without links.

Leicht, E. A. & Newman, M. E. J. Community structure in directed networks. Phys. Rev. Lett. 100,
118703 (2008). Arenas, A., Duch, J., Fernández, A. & Gómez, S. Size reduction of complex networks
preserving modularity. New J. Phys. 9, 176 (2007). Newman, M. E. J. & Girvan, M. Finding and
evaluating community structure in networks. Phys. Rev. E 69, 026113 (2004). Reichardt, J. &
Bornholdt, S. Statistical mechanics of community detection. Phys. Rev. E 74, 016110 (2006).
"""
function web_modularity(W::AbstractMatrix, groups::AbstractVector; directed::Bool = true, γ::Real = 1)
    S = LinearAlgebra.checksquare(W)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S, got $(length(groups))"))
    γ >= 0 || throw(ArgumentError("the resolution γ must be non-negative, got $γ"))
    return partition_modularity(link_weights(W, directed), groups, γ)
end

# The links of W without the diagonal, as `w + wᵀ` when undirected
function link_weights(W::AbstractMatrix, directed::Bool)
    S = LinearAlgebra.checksquare(W)
    w = [i == j ? 0.0 : float(W[i, j]) for i in 1:S, j in 1:S]
    all(>=(0), w) || throw(ArgumentError("weights must be non-negative"))
    return directed ? w : w + w'
end

# Q of a partition of the link weights `w`, in either form
function partition_modularity(w::AbstractMatrix, groups::AbstractVector, γ::Real)
    S = size(w, 1)
    total = sum(w)
    total == 0 && return 0.0
    s_out = vec(sum(w; dims = 2))
    s_in = vec(sum(w; dims = 1))
    Q = 0.0
    for j in 1:S, i in 1:S
        groups[i] == groups[j] && (Q += w[i, j] - γ * s_out[i] * s_in[j] / total)
    end
    return Q / total
end

"""
    compartments([rng], W; directed = true, γ = 1, method = :louvain, runs = 10, threshold = 0.5)

A partition of the species of the web or non-negative weight matrix `W` (`W[i, j]` the link from
prey `i` to predator `j`) of high modularity, as scored by [`web_modularity`](@ref) with the same
`directed` and `γ`. Returns a named tuple:

- `groups`: the consensus partition, labelled `1, 2, …` in order of first appearance;
- `modularity`: its modularity;
- `best`: the partition of highest modularity among the searches, labelled the same way;
- `run_modularity`: the modularity of each of the `runs` independent searches;
- `coassignment`: the fraction of those searches that put each pair of species in one group.

Each search uses `method`:

- `:louvain`, the Louvain method (Blondel et al. 2008), in its directed form (Dugué & Perez 2015)
  when `directed = true`, fast enough for any web;
- `:annealing`, simulated annealing (Guimerà & Amaral 2005), with moves of single species and
  merges of groups, followed by a greedy pass; slower, and a check on the Louvain result.

After each search, groups that fall apart into unlinked pieces are split, which never lowers the
modularity, so that every group is connected. Many partitions have nearly the highest modularity
(Good, de Montjoye & Clauset 2010), so `compartments` combines the searches by consensus
clustering (Lancichinetti & Fortunato 2012): it keeps the pairs that at least `threshold` of the
searches put together, searches that co-assignment matrix again, and repeats until every search
agrees. If they still disagree after 20 rounds, the partition of highest modularity among the last
searches is returned. The spread of `run_modularity` and the entries of `coassignment` between 0
and 1 show how well defined the compartments are. On weakly modular webs the searches disagree,
and the consensus, which keeps what they share, can score well below `best`. On the bundled webs
small enough to solve exactly (10 cases), the best of 10 searches came within 0.005 of the maximum
modularity by either method, and within 0.001 by the better of the two.

Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. Fast unfolding of communities in
large networks. J. Stat. Mech. P10008 (2008). Dugué, N. & Perez, A. Directed Louvain: maximizing
modularity in directed networks. hal-01231784 (2015). Guimerà, R. & Amaral, L. A. N. Functional
cartography of complex metabolic networks. Nature 433, 895–900 (2005). Good, B. H., de
Montjoye, Y.-A. & Clauset, A. Performance of modularity maximization in practical contexts. Phys.
Rev. E 81, 046106 (2010). Lancichinetti, A. & Fortunato, S. Consensus clustering in complex
networks. Sci. Rep. 2, 336 (2012).
"""
function compartments(
        rng::AbstractRNG, W::AbstractMatrix; directed::Bool = true, γ::Real = 1, method::Symbol = :louvain,
        runs::Integer = 10, threshold::Real = 0.5
    )
    γ >= 0 || throw(ArgumentError("the resolution γ must be non-negative, got $γ"))
    method in (:louvain, :annealing) || throw(ArgumentError("method must be :louvain or :annealing, given :$method"))
    runs >= 1 || throw(ArgumentError("runs must be at least 1"))
    0 < threshold <= 1 || throw(ArgumentError("threshold must be in (0, 1], got $threshold"))
    w = link_weights(W, directed)
    search(x, γx) = connected_groups(x, method === :louvain ? louvain_groups(rng, x, γx) : annealed_groups(rng, x, γx))
    first_round = [search(w, γ) for _ in 1:runs]
    run_modularity = [partition_modularity(w, g, γ) for g in first_round]
    coassignment = coassignment_matrix(first_round)
    partitions = first_round
    groups = nothing
    for _ in 1:20
        D = coassignment_matrix(partitions)
        S = size(D, 1)
        kept = [i != j && D[i, j] >= threshold ? D[i, j] : 0.0 for i in 1:S, j in 1:S]
        if all(x -> x == 0 || x == 1, kept)
            groups = connected_groups(kept, collect(1:S); merge = true)
            break
        end
        partitions = [search(kept, 1) for _ in 1:runs]
    end
    groups === nothing && (groups = partitions[argmax([partition_modularity(w, g, γ) for g in partitions])])
    # the searches of the co-assignment matrix can join species that no link of `w` connects
    groups = relabel(connected_groups(w, groups))
    best = relabel(first_round[argmax(run_modularity)])
    return (; groups, modularity = partition_modularity(w, groups, γ), best, run_modularity, coassignment)
end
compartments(W::AbstractMatrix; kwargs...) = compartments(Random.default_rng(), W; kwargs...)

# the fraction of partitions that put each pair of species in one group
function coassignment_matrix(partitions)
    S = length(first(partitions))
    return [count(g -> g[i] == g[j], partitions) / length(partitions) for i in 1:S, j in 1:S]
end

# labels 1, 2, … in order of first appearance
function relabel(groups)
    labels = Dict{eltype(groups), Int}()
    return [get!(labels, g, length(labels) + 1) for g in groups]
end

# Split each group into the pieces joined by links; with `merge = true`, start from singletons
# and join every linked pair instead, giving the connected components of `w`
function connected_groups(w::AbstractMatrix, groups::AbstractVector; merge::Bool = false)
    S = size(w, 1)
    label = zeros(Int, S)
    next = 0
    for start in 1:S
        label[start] == 0 || continue
        next += 1
        label[start] = next
        stack = [start]
        while !isempty(stack)
            v = pop!(stack)
            for u in 1:S
                (label[u] == 0 && (w[v, u] > 0 || w[u, v] > 0)) || continue
                (merge || groups[u] == groups[v]) || continue
                label[u] = next
                push!(stack, u)
            end
        end
    end
    return label
end

# The change in Q from moving species `v` to group `b`. Moving `v` from group `a` to `b` changes
# the weight inside groups by `L_v(b) − L_v(a)`, with `L_v(c)` the weight between `v` and the
# other members of `c` in both directions (a weight on the diagonal moves with `v`), and the
# expected weight by `γ/W [so_v (G_in(b) − G_in(a) + si_v) + si_v (G_out(b) − G_out(a) + so_v)]`,
# where `G_out` and `G_in` are the groups' total out- and in-weights.
mutable struct ModularityState{M <: AbstractMatrix{Float64}}
    const w::M
    const γ::Float64
    const total::Float64
    const so::Vector{Float64}
    const si::Vector{Float64}
    const groups::Vector{Int}
    const members::Vector{Set{Int}}
    const Gout::Vector{Float64}
    const Gin::Vector{Float64}
end
function ModularityState(w::AbstractMatrix{Float64}, γ::Real, groups::AbstractVector{<:Integer})
    S = size(w, 1)
    so, si = vec(sum(w; dims = 2)), vec(sum(w; dims = 1))
    members = [Set{Int}() for _ in 1:S]
    Gout, Gin = zeros(S), zeros(S)
    for v in 1:S
        push!(members[groups[v]], v)
        Gout[groups[v]] += so[v]
        Gin[groups[v]] += si[v]
    end
    return ModularityState(w, float(γ), sum(w), so, si, collect(groups), members, Gout, Gin)
end

toward(st::ModularityState, v, c) = sum((st.w[v, u] + st.w[u, v] for u in st.members[c] if u != v); init = 0.0)
function move_gain(st::ModularityState, v, b)
    a = st.groups[v]
    a == b && return 0.0
    links = toward(st, v, b) - toward(st, v, a)
    expected = st.γ / st.total *
        (st.so[v] * (st.Gin[b] - st.Gin[a] + st.si[v]) + st.si[v] * (st.Gout[b] - st.Gout[a] + st.so[v]))
    return (links - expected) / st.total
end
function move!(st::ModularityState, v, b)
    a = st.groups[v]
    delete!(st.members[a], v)
    push!(st.members[b], v)
    st.Gout[a] -= st.so[v]
    st.Gin[a] -= st.si[v]
    st.Gout[b] += st.so[v]
    st.Gin[b] += st.si[v]
    st.groups[v] = b
    return nothing
end
empty_group(st::ModularityState) = findfirst(isempty, st.members)

# Greedy moves of single species, in random order, to the group of a linked species (or, with
# `new_groups`, to a new group) that raises Q most, until no move raises it
function greedy_moves!(rng::AbstractRNG, st::ModularityState; new_groups::Bool = false)
    S = length(st.groups)
    linked = [[u for u in 1:S if u != v && (st.w[v, u] > 0 || st.w[u, v] > 0)] for v in 1:S]
    improved = true
    while improved
        improved = false
        for v in Random.shuffle(rng, 1:S)
            targets = unique(st.groups[u] for u in linked[v])
            if new_groups
                c = empty_group(st)
                c === nothing || push!(targets, c)
            end
            isempty(targets) && continue
            gains = [move_gain(st, v, b) for b in targets]
            k = argmax(gains)
            if gains[k] > 1.0e-12
                move!(st, v, targets[k])
                improved = true
            end
        end
    end
    return st
end

# The Louvain method (Blondel et al. 2008), in the directed form of Dugué & Perez (2015) when `w`
# is asymmetric: greedy moves of single species, then the groups become the species of a smaller
# weight matrix that keeps the weight inside each group on its diagonal, until no move is made
function louvain_groups(rng::AbstractRNG, w::AbstractMatrix, γ::Real)
    S = size(w, 1)
    all(iszero, w) && return collect(1:S)
    assignment = collect(1:S)
    current = Matrix{Float64}(w)
    while true
        n = size(current, 1)
        st = greedy_moves!(rng, ModularityState(current, γ, collect(1:n)))
        groups = relabel(st.groups)
        k = maximum(groups)
        k == n && break
        assignment = groups[assignment]
        aggregated = zeros(k, k)
        for i in 1:n, j in 1:n
            aggregated[groups[i], groups[j]] += current[i, j]
        end
        current = aggregated
    end
    return assignment
end

# Simulated annealing of modularity (Guimerà & Amaral 2005): Metropolis moves of single species
# to the group of a linked species or to a new group, and merges of two groups, under geometric
# cooling from a temperature of about the change in Q of one link, then greedy moves from the
# best partition found
function annealed_groups(rng::AbstractRNG, w::AbstractMatrix, γ::Real; cooling::Real = 0.98, sweeps::Integer = 20)
    S = size(w, 1)
    total = sum(w)
    total == 0 && return collect(1:S)
    st = ModularityState(Matrix{Float64}(w), γ, collect(1:S))
    linked = [[u for u in 1:S if u != v && (w[v, u] > 0 || w[u, v] > 0)] for v in 1:S]
    Q = partition_modularity(st.w, st.groups, γ)
    best, bestQ = copy(st.groups), Q
    T = 1 / total
    while T > 1.0e-4 / total
        for _ in 1:(sweeps * S)
            v = rand(rng, 1:S)
            b = if !isempty(linked[v]) && rand(rng) < 0.8
                st.groups[rand(rng, linked[v])]
            else
                something(empty_group(st), st.groups[v])
            end
            ΔQ = move_gain(st, v, b)
            if ΔQ > 0 || rand(rng) < exp(ΔQ / T)
                move!(st, v, b)
                Q += ΔQ
            end
        end
        for _ in 1:S
            nonempty = [c for c in 1:S if !isempty(st.members[c])]
            length(nonempty) < 2 && break
            a, b = rand(rng, nonempty), rand(rng, nonempty)
            a == b && continue
            links = sum((st.w[u, v] + st.w[v, u] for u in st.members[a] for v in st.members[b]); init = 0.0)
            ΔQ = (links - st.γ / total * (st.Gout[a] * st.Gin[b] + st.Gout[b] * st.Gin[a])) / total
            if ΔQ > 0 || rand(rng) < exp(ΔQ / T)
                foreach(v -> move!(st, v, a), collect(st.members[b]))
                Q += ΔQ
            end
        end
        Q > bestQ && ((best, bestQ) = (copy(st.groups), Q))
        T *= cooling
    end
    return greedy_moves!(rng, ModularityState(st.w, γ, best); new_groups = true).groups
end
