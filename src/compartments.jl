# # Compartments
# Compartments are groups of species that interact more among themselves than with others. As
# elsewhere, `W[i, j]` is the link from prey `i` to predator `j`: 1 or 0 in a web, or a
# non-negative weight such as an energy flux. Cannibalistic links are left out.

"""
    web_modularity(W, groups; view = :density, directed = true, γ = 1, signed = false)

The modularity of the partition `groups` (one label per species) of the web or non-negative
weight matrix `W`, where `W[i, j]` is the link from prey `i` to predator `j`, such as an energy
flux. Cannibalistic links (the diagonal) are left out. Guimerà et al. (2010) define compartments
in three ways, chosen by `view`.

`view = :density`, species densely linked to each other. With `directed = true`, the directed
modularity of Leicht & Newman (2008), weighted as in Arenas et al. (2007),

    Q = (1/W) Σ_ij [W_ij − γ s_i^out s_j^in / W] δ(g_i, g_j),

where `s_i^out = Σ_j W_ij`, `s_j^in = Σ_i W_ij` and `W = Σ_ij W_ij`; `Q` is the same for `W` and its
transpose, so the direction convention does not matter. With `directed = false`, the modularity
of Newman & Girvan (2004) on the undirected web `W + Wᵀ`,

    Q = (1/2m) Σ_ij [A_ij − γ k_i k_j / 2m] δ(g_i, g_j),

in which each link counts once, so that two species that eat each other are joined by two links.

`view = :predator`, species that share many prey, and `view = :prey`, species that share many
predators: the bipartite modularity of Guimerà, Sales-Pardo & Amaral (2007), with each prey's set
of predators (or each predator's set of prey) as a team,

    M = Σ_s [ Σ_{i ≠ j ∈ s} c_ij / Σ_a m_a(m_a − 1) − γ Σ_{i ≠ j ∈ s} t_i t_j / L² ],

where `c_ij` counts the prey (predators) that `i` and `j` share, `t_i` the prey (predators) of
`i`, `m_a` the predators of prey `a` (the prey of predator `a`) and `L` the links. The two views
differ for a web and agree for a symmetric one. These views need a binary web and ignore
`directed`.

`signed = true` (density view) allows negative weights, such as the signed pair weights of
[`pair_weights`](@ref), with the signed modularity of Gómez, Jensen & Arenas (2009),

    Q = (1/(W⁺ + W⁻)) Σ_ij [W_ij − γ (s⁺_i^out s⁺_j^in / W⁺ − s⁻_i^out s⁻_j^in / W⁻)] δ(g_i, g_j),

with `s⁺`, `W⁺` and `s⁻`, `W⁻` the strengths and totals of the positive and negative parts. It
rewards positive weights within groups and negative weights between them, and equals
`(W⁺ Q⁺ − W⁻ Q⁻) / (W⁺ + W⁻)` for the modularities `Q⁺` and `Q⁻` of the two parts; without negative
weights it is the modularity above.

`γ` is the resolution of Reichardt & Bornholdt (2006): larger values favour smaller groups. `Q`
is 0 for a single group in the density view and for a web without links.

Leicht, E. A. & Newman, M. E. J. Community structure in directed networks. Phys. Rev. Lett. 100,
118703 (2008). Arenas, A., Duch, J., Fernández, A. & Gómez, S. Size reduction of complex networks
preserving modularity. New J. Phys. 9, 176 (2007). Newman, M. E. J. & Girvan, M. Finding and
evaluating community structure in networks. Phys. Rev. E 69, 026113 (2004). Gómez, S., Jensen, P.
& Arenas, A. Analysis of community structure in networks of correlated data. Phys. Rev. E 80,
016114 (2009). Guimerà, R.,
Sales-Pardo, M. & Amaral, L. A. N. Module identification in bipartite and directed networks.
Phys. Rev. E 76, 036102 (2007). Guimerà, R. et al. Origin of compartmentalization in food webs.
Ecology 91, 2941–2951 (2010). Reichardt, J. & Bornholdt, S. Statistical mechanics of community
detection. Phys. Rev. E 74, 016110 (2006).
"""
function web_modularity(
        W::AbstractMatrix, groups::AbstractVector; view::Symbol = :density, directed::Bool = true, γ::Real = 1,
        signed::Bool = false
    )
    S = LinearAlgebra.checksquare(W)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S, got $(length(groups))"))
    γ >= 0 || throw(ArgumentError("the resolution γ must be non-negative, got $γ"))
    return partition_modularity(modularity_objective(W, view, directed, γ; signed), groups)
end

# The links of W without the diagonal, as `w + wᵀ` when undirected
function link_weights(W::AbstractMatrix, directed::Bool)
    S = LinearAlgebra.checksquare(W)
    w = [i == j ? 0.0 : float(W[i, j]) for i in 1:S, j in 1:S]
    all(>=(0), w) || throw(ArgumentError("weights must be non-negative"))
    return directed ? w : w + w'
end

# Every modularity here has the form Q = (1/T) Σ_ij [w_ij − Σ_t κ_t so_t,i si_t,j] δ(g_i, g_j) + offset,
# a sum of null terms t with coefficients κ_t, and T the total weight (Σ w, or W⁺ + W⁻ when signed). The density view uses the link weights and their row and column sums. The role
# views use the shared prey or predators c_ij and the prey or predator counts t_i: with
# K = Σ_a m_a(m_a − 1) = Σ c_ij, the bipartite modularity is this form with γ K² / L² in place of
# γ, plus the terms i = j of the null that it leaves out, γ Σ t_i² / L², which do not depend on the
# partition. Searches use (w, so, si, γ); the offset only corrects the value reported.
function modularity_objective(W::AbstractMatrix, view::Symbol, directed::Bool, γ::Real; signed::Bool = false)
    if signed
        view === :density || throw(ArgumentError("signed modularity is defined for the density view only"))
        return signed_objective(W, directed, γ)
    end
    view === :density && return density_objective(link_weights(W, directed), γ)
    view in (:predator, :prey) || throw(ArgumentError("view must be :density, :predator or :prey, given :$view"))
    a = link_weights(W, true)
    all(x -> x == 0 || x == 1, a) || throw(ArgumentError("the $view view needs a binary web"))
    # teams are the rows of `members`, members the columns: for the predator view each prey's
    # predators, for the prey view each predator's prey
    members = view === :predator ? a : Matrix(a')
    S = size(a, 1)
    c = [i == j ? 0.0 : sum(members[:, i] .* members[:, j]) for i in 1:S, j in 1:S]
    t = vec(sum(members; dims = 1))
    L, K = sum(t), sum(c)
    (L == 0 || K == 0) && return (; w = c, total = K, terms = Tuple{Float64, Vector{Float64}, Vector{Float64}}[], offset = 0.0, refine = true)
    return (; w = c, total = K, terms = [(γ * K / L^2, t, t)], offset = γ * sum(abs2, t) / L^2, refine = true)
end

# The density view: one null term γ s_i^out s_j^in / W
function density_objective(w::AbstractMatrix, γ::Real)
    W = sum(w)
    terms = W > 0 ? [(γ / W, vec(sum(w; dims = 2)), vec(sum(w; dims = 1)))] : Tuple{Float64, Vector{Float64}, Vector{Float64}}[]
    return (; w, total = W, terms, offset = 0.0, refine = true)
end

# Signed modularity (Gómez, Jensen & Arenas 2009): with w⁺ and w⁻ the positive and negative parts,
# Q = (1/(W⁺ + W⁻)) Σ [w_ij − γ (s⁺_i^out s⁺_j^in / W⁺ − s⁻_i^out s⁻_j^in / W⁻)] δ(g_i, g_j). Splitting a
# group into unlinked pieces can lower it, through the negative term, so groups are not refined.
function signed_objective(W::AbstractMatrix, directed::Bool, γ::Real)
    S = LinearAlgebra.checksquare(W)
    w = [i == j ? 0.0 : float(W[i, j]) for i in 1:S, j in 1:S]
    # the positive and negative parts, each made undirected on its own so that a pair with
    # entries of opposite sign keeps both
    plus, minus = max.(w, 0.0), max.(-w, 0.0)
    directed || (plus = plus + plus'; minus = minus + minus')
    w = plus - minus
    terms = Tuple{Float64, Vector{Float64}, Vector{Float64}}[]
    Wp, Wm = sum(plus), sum(minus)
    Wp > 0 && push!(terms, (γ / Wp, vec(sum(plus; dims = 2)), vec(sum(plus; dims = 1))))
    Wm > 0 && push!(terms, (-γ / Wm, vec(sum(minus; dims = 2)), vec(sum(minus; dims = 1))))
    return (; w, total = Wp + Wm, terms, offset = 0.0, refine = false)
end

# Q of a partition for an objective
function partition_modularity(obj, groups::AbstractVector)
    (; w, total, terms, offset) = obj
    S = size(w, 1)
    total == 0 && return 0.0
    Q = 0.0
    for j in 1:S, i in 1:S
        groups[i] == groups[j] || continue
        Q += w[i, j] - sum((κ * so[i] * si[j] for (κ, so, si) in terms); init = 0.0)
    end
    return Q / total + offset
end

"""
    compartments([rng], W; view = :density, directed = true, γ = 1, signed = false, method = :louvain,
                 runs = 10, threshold = 0.5)

A partition of the species of the web or non-negative weight matrix `W` (`W[i, j]` the link from
prey `i` to predator `j`) of high modularity, as scored by [`web_modularity`](@ref) with the same
`view`, `directed`, `γ` and `signed`. Returns a named tuple:

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
modularity, so that every group is connected (in the predator and prey views, through shared prey
or predators); signed modularity can fall when a group is split, so its groups are not refined.
For signed modularity, Louvain can stop a few percent short of the best partition; use
`method = :annealing`. Many partitions have nearly the highest modularity
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
        rng::AbstractRNG, W::AbstractMatrix; view::Symbol = :density, directed::Bool = true, γ::Real = 1,
        signed::Bool = false, method::Symbol = :louvain, runs::Integer = 10, threshold::Real = 0.5
    )
    γ >= 0 || throw(ArgumentError("the resolution γ must be non-negative, got $γ"))
    method in (:louvain, :annealing) || throw(ArgumentError("method must be :louvain or :annealing, given :$method"))
    runs >= 1 || throw(ArgumentError("runs must be at least 1"))
    0 < threshold <= 1 || throw(ArgumentError("threshold must be in (0, 1], got $threshold"))
    obj = modularity_objective(W, view, directed, γ; signed)
    w = obj.w
    refine(o, g) = o.refine ? connected_groups(o.w, g) : g
    search(o) = refine(o, method === :louvain ? louvain_groups(rng, o) : annealed_groups(rng, o))
    first_round = [search(obj) for _ in 1:runs]
    run_modularity = [partition_modularity(obj, g) for g in first_round]
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
        partitions = [search(density_objective(kept, 1)) for _ in 1:runs]
    end
    groups === nothing && (groups = partitions[argmax([partition_modularity(obj, g) for g in partitions])])
    # the searches of the co-assignment matrix can join species that no link of `w` connects
    groups = relabel(refine(obj, groups))
    best = relabel(first_round[argmax(run_modularity)])
    return (; groups, modularity = partition_modularity(obj, groups), best, run_modularity, coassignment)
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
    const total::Float64
    # each null term t: its coefficient κ[t] and species vectors so[:, t], si[:, t]
    const κ::Vector{Float64}
    const so::Matrix{Float64}
    const si::Matrix{Float64}
    const groups::Vector{Int}
    const members::Vector{Set{Int}}
    # group sums of so and si for each term
    const Gout::Matrix{Float64}
    const Gin::Matrix{Float64}
end
function ModularityState(w::AbstractMatrix{Float64}, terms, groups::AbstractVector{<:Integer}; total::Real = sum(abs, w))
    S = size(w, 1)
    κ = Float64[t[1] for t in terms]
    so = isempty(terms) ? zeros(S, 0) : reduce(hcat, (Vector{Float64}(t[2]) for t in terms))
    si = isempty(terms) ? zeros(S, 0) : reduce(hcat, (Vector{Float64}(t[3]) for t in terms))
    members = [Set{Int}() for _ in 1:S]
    Gout, Gin = zeros(S, length(κ)), zeros(S, length(κ))
    for v in 1:S
        push!(members[groups[v]], v)
        Gout[groups[v], :] .+= so[v, :]
        Gin[groups[v], :] .+= si[v, :]
    end
    return ModularityState(w, float(total), κ, so, si, collect(groups), members, Gout, Gin)
end
state_terms(st::ModularityState) = [(st.κ[t], st.so[:, t], st.si[:, t]) for t in eachindex(st.κ)]

toward(st::ModularityState, v, c) = sum((st.w[v, u] + st.w[u, v] for u in st.members[c] if u != v); init = 0.0)
function move_gain(st::ModularityState, v, b)
    a = st.groups[v]
    a == b && return 0.0
    links = toward(st, v, b) - toward(st, v, a)
    expected = 0.0
    for t in eachindex(st.κ)
        so, si = st.so[v, t], st.si[v, t]
        expected += st.κ[t] * (so * (st.Gin[b, t] - st.Gin[a, t] + si) + si * (st.Gout[b, t] - st.Gout[a, t] + so))
    end
    return (links - expected) / st.total
end
function move!(st::ModularityState, v, b)
    a = st.groups[v]
    delete!(st.members[a], v)
    push!(st.members[b], v)
    for t in eachindex(st.κ)
        st.Gout[a, t] -= st.so[v, t]
        st.Gin[a, t] -= st.si[v, t]
        st.Gout[b, t] += st.so[v, t]
        st.Gin[b, t] += st.si[v, t]
    end
    st.groups[v] = b
    return nothing
end
empty_group(st::ModularityState) = findfirst(isempty, st.members)

# Greedy moves of single species, in random order, to the group of a linked species (or, with
# `new_groups`, to a new group) that raises Q most, until no move raises it
function greedy_moves!(rng::AbstractRNG, st::ModularityState; new_groups::Bool = false)
    S = length(st.groups)
    linked = [[u for u in 1:S if u != v && (st.w[v, u] != 0 || st.w[u, v] != 0)] for v in 1:S]
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
# weight matrix that keeps the weight inside each group on its diagonal, until no move is made;
# then greedy moves of single species from the final partition
function louvain_groups(rng::AbstractRNG, obj)
    S = size(obj.w, 1)
    all(iszero, obj.w) && return collect(1:S)
    assignment = collect(1:S)
    current, terms = Matrix{Float64}(obj.w), obj.terms
    while true
        n = size(current, 1)
        st = greedy_moves!(rng, ModularityState(current, terms, collect(1:n); obj.total))
        groups = relabel(st.groups)
        k = maximum(groups)
        k == n && break
        assignment = groups[assignment]
        aggregated = zeros(k, k)
        for i in 1:n, j in 1:n
            aggregated[groups[i], groups[j]] += current[i, j]
        end
        terms = [(κ, [sum(so[groups .== r]) for r in 1:k], [sum(si[groups .== r]) for r in 1:k]) for (κ, so, si) in terms]
        current = aggregated
    end
    # moves of single species from the final partition, which the moves of aggregated groups
    # cannot make
    return greedy_moves!(rng, ModularityState(Matrix{Float64}(obj.w), obj.terms, assignment; obj.total); new_groups = true).groups
end

# Simulated annealing of modularity (Guimerà & Amaral 2005): Metropolis moves of single species
# to the group of a linked species or to a new group, and merges of two groups, under geometric
# cooling from a temperature of about the change in Q of one link, then greedy moves from the
# best partition found
function annealed_groups(rng::AbstractRNG, obj; cooling::Real = 0.98, sweeps::Integer = 20)
    w, total = obj.w, obj.total
    S = size(w, 1)
    total == 0 && return collect(1:S)
    st = ModularityState(Matrix{Float64}(w), obj.terms, collect(1:S); total)
    linked = [[u for u in 1:S if u != v && (w[v, u] != 0 || w[u, v] != 0)] for v in 1:S]
    Q = partition_modularity((; w = st.w, total, terms = obj.terms, offset = 0.0), st.groups)
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
            expected = sum(st.κ[t] * (st.Gout[a, t] * st.Gin[b, t] + st.Gout[b, t] * st.Gin[a, t]) for t in eachindex(st.κ); init = 0.0)
            ΔQ = (links - expected) / total
            if ΔQ > 0 || rand(rng) < exp(ΔQ / T)
                foreach(v -> move!(st, v, a), collect(st.members[b]))
                Q += ΔQ
            end
        end
        Q > bestQ && ((best, bestQ) = (copy(st.groups), Q))
        T *= cooling
    end
    return greedy_moves!(rng, ModularityState(st.w, obj.terms, best; total); new_groups = true).groups
end

"""
    modularity_zscore([rng], web; null = :rewire, samples = 100, swaps = 10 * L, view = :density,
                      directed = true, γ = 1, method = :louvain, runs = 5)

How much more modular `web` is than webs from a null model, with the modularity of
each web, in the `view` of [`web_modularity`](@ref), the highest found by [`compartments`](@ref)
in `runs` searches by `method`. Random webs
have high modularity from fluctuations alone (Guimerà, Sales-Pardo & Amaral 2004), so the same
search runs on every null web: comparing null webs on the partition found for `web` would bias
the test towards significance.

- `null = :rewire`: [`rewire_links`](@ref) with `swaps` swaps, keeping each species' numbers of
  prey and predators and of one-way and two-way links;
- `null = :shuffle`: the weights of a weight matrix, such as [`pair_weights`](@ref) of a
  community matrix, permuted at random among its links (among pairs when it is symmetric),
  keeping the topology; this tests whether the strengths are arranged in compartments beyond
  what the topology gives;
- `null = f`, a function `f(rng, S, C)` returning a web with `S` species and connectance about `C`,
  such as `niche_network`, `cascade_network` or `(rng, S, C) -> generalized_niche_network(rng, S,
  C, 0.6)`: a structural model at the species richness and connectance of `web`. Guimerà et al.
  (2010) found that empirical food webs are no more modular than the generalized niche model.

Returns a named tuple: the `observed` modularity and its partition `groups`, the modularity of
each null web `null_modularity`, their `mean` and `sd`, `z = (observed − mean) / sd` (`NaN` when
`sd` is 0), and `p`, the fraction of null webs at least as modular as `web`.

Guimerà, R., Sales-Pardo, M. & Amaral, L. A. N. Modularity from fluctuations in random graphs and
complex networks. Phys. Rev. E 70, 025101 (2004). Guimerà, R. et al. Origin of
compartmentalization in food webs. Ecology 91, 2941–2951 (2010).
"""
function modularity_zscore(
        rng::AbstractRNG, web::AbstractMatrix; null = :rewire, samples::Integer = 100,
        swaps::Integer = 10 * count(!iszero, web), view::Symbol = :density, directed::Bool = true, γ::Real = 1,
        method::Symbol = :louvain, runs::Integer = 5
    )
    samples >= 2 || throw(ArgumentError("samples must be at least 2"))
    S = LinearAlgebra.checksquare(web)
    binary = all(x -> x == 0 || x == 1, web)
    null === :shuffle || binary ||
        throw(ArgumentError("the $null null generates binary webs, so it needs a binary web; use null = :shuffle for weights"))
    draw = if null === :rewire
        () -> rewire_links(rng, web; swaps)
    elseif null === :shuffle
        () -> shuffle_weights(rng, web)
    elseif null isa Function
        C = connectance(web)
        () -> null(rng, S, C)
    else
        throw(ArgumentError("null must be :rewire, :shuffle or a function (rng, S, C) -> web, given $null"))
    end
    best(W) = compartments(rng, W; view, directed, γ, method, runs)
    found = best(web)
    observed = maximum(found.run_modularity)
    null_modularity = [maximum(best(draw()).run_modularity) for _ in 1:samples]
    μ, σ = mean(null_modularity), std(null_modularity)
    z = σ > 0 ? (observed - μ) / σ : NaN
    p = count(>=(observed), null_modularity) / samples
    return (; observed, groups = found.best, null_modularity, mean = μ, sd = σ, z, p)
end
modularity_zscore(web::AbstractMatrix; kwargs...) = modularity_zscore(Random.default_rng(), web; kwargs...)

"""
    species_roles(W, groups; links = :all)

The role of each species in the partition `groups` of the web or non-negative weight matrix `W`
(`W[i, j]` the link from prey `i` to predator `j`), from its links: all of them (`links = :all`,
counting each link once), those to its prey (`:prey`) or those to its predators (`:predators`).
With weights, the counts below are sums of weights (strengths). Cannibalistic links are left out.
Returns a named tuple:

- `z`: the within-module degree, `z_i = (κ_i − κ̄) / σ_κ`, where `κ_i` counts the links of `i` to
  other members of its group and `κ̄` and `σ_κ` are the mean and standard deviation of `κ` in the
  group, with the sample standard deviation as in `czvalues` of the R package bipartite; 0 when
  `σ_κ` is 0 or undefined, as in groups of one species;
- `P`: the participation coefficient, `P_i = 1 − Σ_s (κ_is / k_i)²`, where `κ_is` counts the links
  of `i` to group `s` and `k_i` all its links; 0 near the group's own links only and close to 1 when
  spread evenly over all groups; 0 for a species without links;
- `role`: the classes of Olesen et al. (2007): `:peripheral` (`z ≤ 2.5`, `P ≤ 0.62`),
  `:connector` (`z ≤ 2.5`, `P > 0.62`), `:module_hub` (`z > 2.5`, `P ≤ 0.62`) and `:network_hub`
  (`z > 2.5`, `P > 0.62`);
- `cartography`: the seven roles of Guimerà & Amaral (2005), for non-hubs (`z < 2.5`)
  `:ultra_peripheral` (`P ≤ 0.05`), `:peripheral` (`P ≤ 0.62`), `:connector` (`P ≤ 0.80`) and
  `:kinless` (`P > 0.80`), and for hubs (`z ≥ 2.5`) `:provincial_hub` (`P ≤ 0.30`),
  `:connector_hub` (`P ≤ 0.75`) and `:kinless_hub` (`P > 0.75`).

Connectors, the species with high `P`, are the species that link compartments.

Guimerà, R. & Amaral, L. A. N. Functional cartography of complex metabolic networks. Nature 433,
895–900 (2005). Olesen, J. M., Bascompte, J., Dupont, Y. L. & Jordano, P. The modularity of
pollination networks. Proc. Natl Acad. Sci. USA 104, 19891–19896 (2007).
"""
function species_roles(W::AbstractMatrix, groups::AbstractVector; links::Symbol = :all)
    S = LinearAlgebra.checksquare(W)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S, got $(length(groups))"))
    w = link_weights(W, true)
    # a[i, j]: the link of species i to species j counted for i's role
    a = if links === :all
        w + w'
    elseif links === :prey
        Matrix(w')
    elseif links === :predators
        w
    else
        throw(ArgumentError("links must be :all, :prey or :predators, given :$links"))
    end
    labels = unique(groups)
    toward = [sum((a[i, j] for j in 1:S if groups[j] == s); init = 0.0) for i in 1:S, s in labels]
    own = [findfirst(==(groups[i]), labels) for i in 1:S]
    κ = [toward[i, own[i]] for i in 1:S]
    k = vec(sum(toward; dims = 2))
    z = zeros(S)
    for s in eachindex(labels)
        members = findall(==(s), own)
        μ, σ = mean(κ[members]), length(members) > 1 ? std(κ[members]) : 0.0
        σ > 0 && (z[members] .= (κ[members] .- μ) ./ σ)
    end
    P = [k[i] > 0 ? 1 - sum(abs2, toward[i, :] ./ k[i]) : 0.0 for i in 1:S]
    role = [
        z[i] > 2.5 ? (P[i] > 0.62 ? :network_hub : :module_hub) : (P[i] > 0.62 ? :connector : :peripheral)
            for i in 1:S
    ]
    cartography = map(1:S) do i
        if z[i] >= 2.5
            P[i] <= 0.3 ? :provincial_hub : P[i] <= 0.75 ? :connector_hub : :kinless_hub
        else
            P[i] <= 0.05 ? :ultra_peripheral : P[i] <= 0.62 ? :peripheral : P[i] <= 0.8 ? :connector : :kinless
        end
    end
    return (; z, P, role, cartography)
end

"""
    compartment_homogeneity(web, groups; levels = trophic_levels(web))

How far the compartments `groups` of `web` follow trophic structure, as the trophic-level and sink
homogeneity of Guimerà et al. (2010), each a z-score against random partitions with the same
group sizes. Positive values mean more homogeneous than chance. Returns a named tuple:

- `trophic`: `H_t = (⟨B⟩ − B) / σ_B`, where `B` counts the boundaries between groups along the
  species sorted by `levels`; compartments of contiguous trophic levels have few boundaries.
  Species with equal levels are ordered at random, and `B` is its mean over those orderings.
- `sink`: `H_s = Σ_i (o_i − ⟨o_i⟩) / σ_{o_i}`, where `o_i` counts the other members of species
  `i`'s group that are in its sink web (its prey, their prey and so on); species whose `o_i`
  cannot vary are left out.

The means and standard deviations over random partitions are exact.

Guimerà, R. et al. Origin of compartmentalization in food webs. Ecology 91, 2941–2951 (2010).
"""
function compartment_homogeneity(web::AbstractMatrix, groups::AbstractVector; levels::AbstractVector = trophic_levels(web))
    S = LinearAlgebra.checksquare(web)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S, got $(length(groups))"))
    length(levels) == S || throw(DimensionMismatch("levels must have length $S, got $(length(levels))"))
    labels = unique(groups)
    sizes = [count(==(s), groups) for s in labels]
    return (; trophic = trophic_homogeneity(groups, levels, labels, sizes), sink = sink_homogeneity(web, groups, sizes))
end

# The expected number of boundaries along the species sorted by level, over random orderings
# of species with equal levels: within a block of `m` tied species with group counts `m_c`,
# neighbours differ with probability 1 − Σ m_c (m_c − 1) / (m (m − 1)); across consecutive blocks
# with counts `a_c` and `b_c`, with probability 1 − Σ a_c b_c / (m_a m_b).
function boundaries(groups, levels, labels)
    order = sortperm(levels)
    blocks = Vector{Vector{Int}}()
    for i in order
        if isempty(blocks) || !isapprox(levels[i], levels[last(blocks)[1]]; atol = 1.0e-9)
            push!(blocks, [i])
        else
            push!(last(blocks), i)
        end
    end
    counts = [[count(i -> groups[i] == s, block) for s in labels] for block in blocks]
    B = 0.0
    for (block, c) in zip(blocks, counts)
        m = length(block)
        m > 1 && (B += (m - 1) * (1 - sum(c .* (c .- 1)) / (m * (m - 1))))
    end
    for k in 2:length(blocks)
        B += 1 - sum(counts[k - 1] .* counts[k]) / (length(blocks[k - 1]) * length(blocks[k]))
    end
    return B
end

# H_t. Under a random arrangement of the labels, with `n_c` species in group `c`, the indicator
# `X_t` that neighbours `t` and `t + 1` share a group has mean `p2 = Σ n_c (n_c − 1) / (S)_2`, two
# overlapping indicators are both 1 with probability `p3 = Σ (n_c)_3 / (S)_3`, and two disjoint
# ones with probability `p4 = [Σ (n_c)_4 + Σ_{c ≠ d} (n_c)_2 (n_d)_2] / (S)_4`, where `(n)_k` is the
# falling factorial. `B = (S − 1) − Σ X_t`.
function trophic_homogeneity(groups, levels, labels, sizes)
    S = length(groups)
    S < 4 && return NaN
    falling(n, k) = prod(n - j for j in 0:(k - 1); init = 1.0)
    p2 = sum(falling(n, 2) for n in sizes) / falling(S, 2)
    p3 = sum(falling(n, 3) for n in sizes) / falling(S, 3)
    pairs2 = [falling(n, 2) for n in sizes]
    p4 = (sum(falling(n, 4) for n in sizes) + sum(pairs2)^2 - sum(abs2, pairs2)) / falling(S, 4)
    adjacent = S - 2
    disjoint = (S - 1) * (S - 2) ÷ 2 - adjacent
    mean_B = (S - 1) * (1 - p2)
    var_B = (S - 1) * p2 * (1 - p2) + 2 * adjacent * (p3 - p2^2) + 2 * disjoint * (p4 - p2^2)
    var_B > 0 || return NaN
    return (mean_B - boundaries(groups, levels, labels)) / sqrt(var_B)
end

# H_s. In a random partition with the observed group sizes, species `i` falls in a group of size
# `n` with probability `n · (number of groups of size n) / S`, and the other `n − 1` members are a
# random draw from the other `S − 1` species, so `o_i` is hypergeometric given `n`: with `t_i`
# species in the sink web of `i`, mean `(n − 1) t_i / (S − 1)` and variance
# `(n − 1) t_i (S − 1 − t_i) (S − n) / ((S − 1)² (S − 2))`.
function sink_homogeneity(web, groups, sizes)
    S = length(groups)
    S < 3 && return NaN
    eats = [i != j && web[i, j] != 0 for i in 1:S, j in 1:S]
    weights = Dict(n => n * count(==(n), sizes) / S for n in unique(sizes))
    H = 0.0
    for i in 1:S
        # the sink web of i, by search down its prey
        sink = falses(S)
        stack = [i]
        while !isempty(stack)
            v = pop!(stack)
            for u in 1:S
                (eats[u, v] && !sink[u] && u != i) || continue
                sink[u] = true
                push!(stack, u)
            end
        end
        t = count(sink)
        o = count(u -> u != i && sink[u] && groups[u] == groups[i], 1:S)
        μ = sum(p * (n - 1) * t / (S - 1) for (n, p) in weights)
        second = sum(
            p * ((n - 1) * t * (S - 1 - t) * (S - n) / ((S - 1)^2 * (S - 2)) + ((n - 1) * t / (S - 1))^2)
                for (n, p) in weights
        )
        σ = sqrt(max(second - μ^2, 0.0))
        σ > 1.0e-12 && (H += (o - μ) / σ)
    end
    return H
end

"""
    block_structure(M, groups)

The within- and between-group structure of a community (Jacobian) matrix or flux matrix `M` for
the partition `groups`, the quantities that set the eigenvalues of block-structured random
matrices: the mean of the entries between two groups sets the outlier eigenvalues, and their
variance and pair correlation the bulk (Grilli, Rogers & Allesina 2016). Groups are ordered by
first appearance in `groups`. Returns a named tuple of `k × k` matrices over pairs of groups
`(r, s)`, for the off-diagonal entries `M[i, j]` with `i` in `r` and `j` in `s`:

- `connectance`: the fraction of nonzero entries;
- `mean`, `variance`: their mean and population variance, counting zero entries, as in
  [`elliptic_law`](@ref), so that a block with nonzero entries of mean `μ` and variance `σ²` at
  connectance `C` has mean `Cμ` and variance `C(σ² + μ²) − C²μ²`;
- `correlation`: the correlation of the pairs `(M[i, j], M[j, i])`, with the geometric mean of
  the two variances in the denominator as in `elliptic_law` (`NaN` when either is 0);
- `predator_prey`, `competition`, `mutualism`, `one_way`: the numbers of pairs with entries of
  opposite sign, both negative, both positive, and one zero, symmetric in `r` and `s`;

and, for the groups and the species:

- `labels`, `sizes`: the groups and their numbers of species;
- `diagonal`: the mean diagonal entry of each group;
- `coupling`: for each species, the share of its interaction magnitude, `Σ_j |M[i, j]| + |M[j, i]|`
  over `j ≠ i`, that is with other groups; the coupling species are those with `coupling > 0`.

With one group, `mean` and `correlation` are those of `elliptic_law`.

Grilli, J., Rogers, T. & Allesina, S. Modularity and stability in ecological communities. Nat.
Commun. 7, 12031 (2016).
"""
function block_structure(M::AbstractMatrix{<:Real}, groups::AbstractVector)
    S = LinearAlgebra.checksquare(M)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S, got $(length(groups))"))
    labels = unique(groups)
    k = length(labels)
    g = [findfirst(==(groups[i]), labels) for i in 1:S]
    members = [findall(==(r), g) for r in 1:k]
    entries(r, s) = [float(M[i, j]) for i in members[r] for j in members[s] if i != j]
    connectance, μ, variance, correlation = (fill(NaN, k, k) for _ in 1:4)
    predator_prey, competition, mutualism, one_way = (zeros(Int, k, k) for _ in 1:4)
    for r in 1:k, s in 1:k
        x = entries(r, s)
        isempty(x) && continue
        connectance[r, s] = count(!iszero, x) / length(x)
        μ[r, s] = mean(x)
        variance[r, s] = mean(abs2, x .- μ[r, s])
        r <= s || continue
        pairs = [(float(M[i, j]), float(M[j, i])) for i in members[r] for j in members[s] if (r == s ? i < j : true)]
        a, b = first.(pairs), last.(pairs)
        va, vb = mean(abs2, a .- mean(a)), mean(abs2, b .- mean(b))
        ρ = va > 0 && vb > 0 ? mean((a .- mean(a)) .* (b .- mean(b))) / sqrt(va * vb) : NaN
        correlation[r, s] = correlation[s, r] = ρ
        counts = (
            count(p -> p[1] * p[2] < 0, pairs), count(p -> p[1] < 0 && p[2] < 0, pairs),
            count(p -> p[1] > 0 && p[2] > 0, pairs), count(p -> xor(p[1] == 0, p[2] == 0), pairs),
        )
        for (matrix, n) in zip((predator_prey, competition, mutualism, one_way), counts)
            matrix[r, s] = matrix[s, r] = n
        end
    end
    diagonal = [mean(float(M[i, i]) for i in members[r]) for r in 1:k]
    coupling = map(1:S) do i
        total = sum((abs(M[i, j]) + abs(M[j, i]) for j in 1:S if j != i); init = 0.0)
        outside = sum((abs(M[i, j]) + abs(M[j, i]) for j in 1:S if g[j] != g[i]); init = 0.0)
        total > 0 ? outside / total : 0.0
    end
    return (;
        labels, sizes = length.(members), connectance, mean = μ, variance, correlation,
        predator_prey, competition, mutualism, one_way, diagonal, coupling,
    )
end


"""
    shuffle_weights([rng], W)

A copy of the weight matrix `W` with its nonzero off-diagonal weights permuted at random among
the same positions, keeping which species interact. When `W` is symmetric, as for
[`pair_weights`](@ref), the weights are permuted among pairs so that the copy stays symmetric.
The diagonal is kept. The null model of `modularity_zscore(...; null = :shuffle)`.
"""
function shuffle_weights(rng::AbstractRNG, W::AbstractMatrix)
    S = LinearAlgebra.checksquare(W)
    shuffled = Matrix{Float64}(W)
    if W == W'
        positions = [(i, j) for i in 1:S for j in (i + 1):S if W[i, j] != 0]
        values = Random.shuffle(rng, [float(W[i, j]) for (i, j) in positions])
        for ((i, j), v) in zip(positions, values)
            shuffled[i, j] = shuffled[j, i] = v
        end
    else
        positions = [(i, j) for i in 1:S for j in 1:S if i != j && W[i, j] != 0]
        values = Random.shuffle(rng, [float(W[i, j]) for (i, j) in positions])
        for ((i, j), v) in zip(positions, values)
            shuffled[i, j] = v
        end
    end
    return shuffled
end
shuffle_weights(W::AbstractMatrix) = shuffle_weights(Random.default_rng(), W)

"""
    pair_weights(J; method = :geometric)

Symmetric weights for the pairs of species of a community (Jacobian) matrix `J`,
for [`web_modularity`](@ref), [`compartments`](@ref) and [`modularity_zscore`](@ref). Modularity
needs non-negative weights, and the signed entries of `J` do not give them: a predator–prey pair
has one negative and one positive entry, and their sum `J_ij + J_ji` cancels strong trophic
links. The weight of the pair `(i, j)` is

- `method = :geometric`: `√|J_ij J_ji|`, unchanged by rescaling the abundance units of any
  species (`J → D J D⁻¹` for diagonal `D`), as for feasible communities `J = diag(x) A`; a pair
  with only one nonzero entry, such as a donor-controlled link, gets weight 0;
- `method = :sum`: `|J_ij| + |J_ji|`, which keeps one-way pairs and reduces to the binary web
  when every magnitude is equal;
- `method = :signed`: `sign(J_ij J_ji) √|J_ij J_ji|`, negative for predator–prey pairs and positive
  for competition and mutualism, for signed modularity (`signed = true`), which then groups
  competitors and mutualists and separates predators from their prey.

The sign of `J_ij J_ji` tells predator–prey pairs (negative) from competition and mutualism
(positive); [`block_structure`](@ref) counts them by block. Neither weight has been evaluated
in the literature as a basis for compartments.
"""
function pair_weights(J::AbstractMatrix{<:Real}; method::Symbol = :geometric)
    S = LinearAlgebra.checksquare(J)
    weight = if method === :geometric
        (a, b) -> sqrt(abs(a * b))
    elseif method === :sum
        (a, b) -> abs(a) + abs(b)
    elseif method === :signed
        (a, b) -> sign(a * b) * sqrt(abs(a * b))
    else
        throw(ArgumentError("method must be :geometric, :sum or :signed, given :$method"))
    end
    return [i == j ? 0.0 : float(weight(J[i, j], J[j, i])) for i in 1:S, j in 1:S]
end

"""
    source_reliance(W; sources = basal species)

The fraction of the energy of each species that comes from each basal source, the reliance of
Rooney et al. (2008) used by Zhao et al. (2018), for the web or flux matrix `W` (`W[i, j]` the
flux from prey `i` to predator `j`). Returns an `S × B` matrix, one row per species and one
column per source in `sources` (by default the species that eat nothing; give them to include
detritus that receives inputs, for example). A source relies only on itself, and every other
species on its prey in proportion to its diet,

    R_i = Σ_k p_ki R_k,   p_ki = W_ki / Σ_l W_li,

so `R = (I − Pᵀ)⁻¹ E` with `E` the indicator of the sources. In a binary web the diet is split
equally among the prey. Rows sum to 1. Cannibalistic links are left out. Throws an `ArgumentError`
if some species get no energy from any source, or if the system is singular, as for a feeding
loop that no source enters.

Rooney, N., McCann, K. S. & Moore, J. C. A landscape theory for food web architecture. Ecol.
Lett. 11, 867–881 (2008). Zhao, L., Zhang, H., Tian, W. & Xu, X. Identifying compartments in
ecological networks based on energy channels. Ecol. Evol. 8, 309–318 (2018).
"""
function source_reliance(W::AbstractMatrix; sources::AbstractVector{<:Integer} = default_sources(W))
    S = LinearAlgebra.checksquare(W)
    w = link_weights(W, true)
    isempty(sources) && throw(ArgumentError("no sources"))
    all(s -> 1 <= s <= S, sources) || throw(ArgumentError("sources must be species indices"))
    intake = vec(sum(w; dims = 1))
    is_source = falses(S)
    is_source[sources] .= true
    for i in 1:S
        is_source[i] || intake[i] > 0 || throw(ArgumentError("species $i eats nothing but is not a source"))
    end
    # P[k, i]: the share of species i's diet from species k, for species that are not sources
    P = [is_source[i] ? 0.0 : w[k, i] / intake[i] for k in 1:S, i in 1:S]
    E = [float(sources[b] == i) for i in 1:S, b in eachindex(sources)]
    A = I - P'
    LinearAlgebra.rank(A) == S ||
        throw(ArgumentError("the reliance is undefined: some species are on a feeding loop that no source enters"))
    R = A \ E
    all(r -> isapprox(sum(r), 1; atol = 1.0e-8), eachrow(R)) ||
        throw(ArgumentError("some species get no energy from any source"))
    return R
end
default_sources(W) = [j for j in axes(W, 2) if !any(i != j && W[i, j] != 0 for i in axes(W, 1))]

"""
    energy_channels(W; sources = basal species, merge = true, γ = 1)

Compartments as energy channels (Zhao et al. 2018): each species joins the channel of the source
it relies on most ([`source_reliance`](@ref); ties go to the first source), and then, with
`merge = true`, the two channels whose merging raises the directed modularity of `W`
([`web_modularity`](@ref), resolution `γ`) the most are merged, until no merging raises it.
Returns a named tuple: the channels `groups`, labelled `1, 2, …` in order of first appearance,
their `modularity`, the `reliance` matrix and the `dominant` source of each species (an index
into `sources`).

The channels follow the flow of energy from the sources, an ecological definition to compare
with the partitions of highest modularity from [`compartments`](@ref). With its non-living
carbon pools and the species that eat nothing as sources, the Chesapeake Bay network reproduces
the worked example of Zhao et al. (2018): the same main source for every consumer, the same three
channels and `Q = 0.39`.
"""
function energy_channels(
        W::AbstractMatrix; sources::AbstractVector{<:Integer} = default_sources(W), merge::Bool = true, γ::Real = 1
    )
    R = source_reliance(W; sources)
    dominant = [argmax(r) for r in eachrow(R)]
    groups = copy(dominant)
    Q = web_modularity(W, groups; γ)
    while merge
        labels = unique(groups)
        length(labels) > 1 || break
        best, pair = 0.0, nothing
        for (k, a) in enumerate(labels), b in labels[(k + 1):end]
            merged = [g == b ? a : g for g in groups]
            ΔQ = web_modularity(W, merged; γ) - Q
            ΔQ > best && ((best, pair) = (ΔQ, (a, b)))
        end
        pair === nothing && break
        a, b = pair
        groups = [g == b ? a : g for g in groups]
        Q += best
    end
    groups = relabel(groups)
    return (; groups, modularity = web_modularity(W, groups; γ), reliance = R, dominant)
end

"""
    group_model(web, groups)

The group model of Allesina & Pascual (2009) for the partition `groups` of the binary `web`
(`web[i, j]` the link from prey `i` to predator `j`): a directed stochastic block model in which a
species of group `r` is eaten by a species of group `s` with probability `p[r, s]`, independently
for each of the `S²` ordered pairs (cannibalism included, as in the connectance `L / S²`). Returns
a named tuple with the maximum likelihood probabilities `p[r, s] = L_rs / (S_r S_s)`, groups in
order of first appearance in `groups` (`labels`), the log-likelihood

    log L = Σ_rs [L_rs log p_rs + (S_r S_s − L_rs) log(1 − p_rs)],

with `0 log 0 = 0`, the Akaike information criterion of the paper, `AIC = 2k² + 2S − 2 log L`
for `k` groups, and the integrated classification likelihood of Daudin, Picard & Robin (2008),

    ICL = log L + Σ_r S_r log(S_r / S) − (k² / 2) log S² − ((k − 1) / 2) log S,

which also counts the assignment of species to groups. Dense blocks on the diagonal of `p` are compartments; dense blocks off it are
groups of species that eat the same prey or are eaten by the same predators, such as trophic
guilds, which modularity does not detect. Degree correction (Karrer & Newman 2011) is not
included.

Allesina, S. & Pascual, M. Food web models: a plea for groups. Ecol. Lett. 12, 652–662 (2009).
Daudin, J.-J., Picard, F. & Robin, S. A mixture model for random graphs. Stat. Comput. 18,
173–183 (2008).
"""
function group_model(web::AbstractMatrix, groups::AbstractVector)
    S = LinearAlgebra.checksquare(web)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S, got $(length(groups))"))
    all(x -> x == 0 || x == 1, web) || throw(ArgumentError("the group model needs a binary web"))
    labels = unique(groups)
    g = [findfirst(==(groups[i]), labels) for i in 1:S]
    k = length(labels)
    n = [count(==(r), g) for r in 1:k]
    L = zeros(Int, k, k)
    for i in 1:S, j in 1:S
        web[i, j] != 0 && (L[g[i], g[j]] += 1)
    end
    p = [L[r, s] / (n[r] * n[s]) for r in 1:k, s in 1:k]
    loglik = sum(block_loglik(L[r, s], n[r] * n[s]) for r in 1:k, s in 1:k)
    # the complete-data log-likelihood adds the group assignments, with proportions n_r / S
    icl = loglik + sum(xlogx(n[r]) for r in 1:k) - xlogx(S) - k^2 / 2 * log(S^2) - (k - 1) / 2 * log(S)
    return (; labels, p, loglik, aic = 2k^2 + 2S - 2loglik, icl)
end

# The maximized log-likelihood of `l` links among `m` ordered pairs, with 0 log 0 = 0
xlogx(x) = x > 0 ? x * log(x) : 0.0
block_loglik(l, m) = m == 0 ? 0.0 : xlogx(l) + xlogx(m - l) - xlogx(m)

"""
    fit_group_model([rng], web; k = 1:min(S, 10), criterion = :icl, runs = 5, cooling = 0.98, sweeps = 20)

The partition of the binary `web` into `k` groups of highest likelihood under the
[`group_model`](@ref), or, for a range of `k`, the one chosen by `criterion`: the highest ICL
(`:icl`, the default) or the lowest AIC (`:aic`, the criterion of Allesina & Pascual 2009). The
AIC of the paper charges `2k²` for the probabilities but nothing for choosing the partition, whose
`2S` term is the same for every `k`, so a search over partitions fits noise: on random webs it
keeps adding groups. The ICL counts the assignment of species to groups and selects one group for
a random web, as the R package sbm does. For each `k`, `runs`
searches by simulated annealing, with moves of single species between groups, are followed by
moves that raise the likelihood until none does; the best is kept. Allesina & Pascual (2009) used
a genetic algorithm with hill climbing for the same objective. Returns the named tuple of
`group_model` for the selected partition, with its `groups` (labelled `1, 2, …` in order of first
appearance), its `k`, and `criterion_by_k`, the value of the criterion for the best partition
found for each `k` tried.
"""
function fit_group_model(
        rng::AbstractRNG, web::AbstractMatrix; k = 1:min(size(web, 1), 10), criterion::Symbol = :icl,
        runs::Integer = 5, cooling::Real = 0.98, sweeps::Integer = 20
    )
    criterion in (:icl, :aic) || throw(ArgumentError("criterion must be :icl or :aic, given :$criterion"))
    S = LinearAlgebra.checksquare(web)
    all(x -> x == 0 || x == 1, web) || throw(ArgumentError("the group model needs a binary web"))
    ks = collect(k)
    all(x -> 1 <= x <= S, ks) || throw(ArgumentError("k must be between 1 and $S"))
    runs >= 1 || throw(ArgumentError("runs must be at least 1"))
    A = Matrix{Int}(web .!= 0)
    best_by_k = map(ks) do kk
        candidates = [anneal_groups_likelihood(rng, A, kk; cooling, sweeps) for _ in 1:runs]
        candidates[argmax([group_model(A, c).loglik for c in candidates])]
    end
    fits = [group_model(A, g) for g in best_by_k]
    values = criterion === :icl ? [f.icl for f in fits] : [f.aic for f in fits]
    b = criterion === :icl ? argmax(values) : argmin(values)
    groups = relabel(best_by_k[b])
    fit = group_model(A, groups)
    return (; groups, k = maximum(groups), fit..., criterion_by_k = Dict(zip(ks, values)))
end
fit_group_model(web::AbstractMatrix; kwargs...) = fit_group_model(Random.default_rng(), web; kwargs...)

# Simulated annealing of the group model log-likelihood for k groups, then moves that raise it
# until none does. The block link counts L and group sizes n are updated by each move; the change
# in log-likelihood comes from the blocks in the rows and columns of the two groups involved.
function anneal_groups_likelihood(rng::AbstractRNG, A::AbstractMatrix{<:Integer}, k::Integer; cooling, sweeps)
    S = size(A, 1)
    k == 1 && return ones(Int, S)
    g = rand(rng, 1:k, S)
    n = [count(==(r), g) for r in 1:k]
    L = zeros(Int, k, k)
    for i in 1:S, j in 1:S
        A[i, j] != 0 && (L[g[i], g[j]] += 1)
    end
    total(L, n) = sum(block_loglik(L[r, s], n[r] * n[s]) for r in 1:k, s in 1:k)
    # links of v to (out) and from (in) each group, not counting its self-loop
    function counts(v)
        out, inn = zeros(Int, k), zeros(Int, k)
        for u in 1:S
            u == v && continue
            A[v, u] != 0 && (out[g[u]] += 1)
            A[u, v] != 0 && (inn[g[u]] += 1)
        end
        return out, inn
    end
    function move_counts(v, b)
        a = g[v]
        out, inn = counts(v)
        self = A[v, v]
        L2, n2 = copy(L), copy(n)
        for s in 1:k
            L2[a, s] -= out[s]; L2[b, s] += out[s]
            L2[s, a] -= inn[s]; L2[s, b] += inn[s]
        end
        L2[a, a] -= self; L2[b, b] += self
        n2[a] -= 1; n2[b] += 1
        return L2, n2
    end
    current = total(L, n)
    best, bestll = copy(g), current
    T = 1.0
    while T > 1.0e-3
        for _ in 1:(sweeps * S)
            v = rand(rng, 1:S)
            b = rand(rng, 1:k)
            b == g[v] && continue
            L2, n2 = move_counts(v, b)
            new = total(L2, n2)
            if new >= current || rand(rng) < exp((new - current) / T)
                g[v], L, n, current = b, L2, n2, new
                current > bestll && ((best, bestll) = (copy(g), current))
            end
        end
        T *= cooling
    end
    # moves that raise the likelihood, from the best partition found
    g .= best
    n = [count(==(r), g) for r in 1:k]
    L = zeros(Int, k, k)
    for i in 1:S, j in 1:S
        A[i, j] != 0 && (L[g[i], g[j]] += 1)
    end
    current = total(L, n)
    improved = true
    while improved
        improved = false
        for v in Random.shuffle(rng, 1:S), b in 1:k
            b == g[v] && continue
            L2, n2 = move_counts(v, b)
            new = total(L2, n2)
            if new > current + 1.0e-10
                g[v], L, n, current = b, L2, n2, new
                improved = true
            end
        end
    end
    return g
end
