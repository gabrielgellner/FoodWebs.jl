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

"""
    modularity_zscore([rng], web; null = :rewire, samples = 100, swaps = 10 * L, directed = true,
                      γ = 1, method = :louvain, runs = 5)

How much more modular the binary `web` is than webs from a null model, with the modularity of
each web the highest found by [`compartments`](@ref) in `runs` searches by `method`. Random webs
have high modularity from fluctuations alone (Guimerà, Sales-Pardo & Amaral 2004), so the same
search runs on every null web: comparing null webs on the partition found for `web` would bias
the test towards significance.

- `null = :rewire`: [`rewire_links`](@ref) with `swaps` swaps, keeping each species' numbers of
  prey and predators and of one-way and two-way links;
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
        swaps::Integer = 10 * count(!iszero, web), directed::Bool = true, γ::Real = 1, method::Symbol = :louvain,
        runs::Integer = 5
    )
    all(x -> x == 0 || x == 1, web) ||
        throw(ArgumentError("modularity_zscore needs a binary web; its null models generate binary webs"))
    samples >= 2 || throw(ArgumentError("samples must be at least 2"))
    S = LinearAlgebra.checksquare(web)
    draw = if null === :rewire
        () -> rewire_links(rng, web; swaps)
    elseif null isa Function
        C = connectance(web)
        () -> null(rng, S, C)
    else
        throw(ArgumentError("null must be :rewire or a function (rng, S, C) -> web, given $null"))
    end
    best(W) = compartments(rng, W; directed, γ, method, runs)
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
