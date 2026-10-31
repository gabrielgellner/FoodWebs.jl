# # Conventions
# A food web on `S` species is an `S`×`S` binary matrix where `web[i, j] == 1` means species `j`
# eats species `i` (prey are rows, predators are columns). Connectance is `C = L / S^2` where `L`
# is the number of links, following Williams & Martinez (2000).

"""
    PredationMatrix(links::AbstractMatrix{<:Integer})

A square binary food web matrix where `web[i, j] == 1` means species `j` eats species `i`.
"""
struct PredationMatrix <: AbstractMatrix{Bool}
    links::Matrix{Bool}
    function PredationMatrix(links::AbstractMatrix{<:Integer})
        S1, S2 = size(links)
        S1 == S2 || throw(DimensionMismatch("PredationMatrix must be square, given ($S1, $S2)"))
        all(x -> x == 0 || x == 1, links) || throw(ArgumentError("PredationMatrix must be binary"))
        return new(Matrix{Bool}(links))
    end
end

Base.size(web::PredationMatrix) = size(web.links)
Base.getindex(web::PredationMatrix, i::Int, j::Int) = web.links[i, j]

"""
    food_chain_network(S)

Generate a linear food chain of length `S`, where species `i + 1` eats species `i`.
"""
function food_chain_network(S::Int)
    S > 0 || throw(ArgumentError("S must be positive, given $S"))
    return PredationMatrix(diagm(1 => trues(S - 1)))
end

# Williams & Martinez (2000) style acceptance: keep generating webs until one is (weakly)
# connected, so there are no isolated species or separate components, and has a number of links
# `L` within a relative tolerance `tol` of the target `C S^2`. The window is at least ±1/2 so that
# it always holds the nearest whole number of links; for small webs `tol C S^2` can fall below
# 1/2, leaving no achievable `L`. With `basal_paths = true` webs are also rejected if some species
# have no feeding path up from a basal species.
function accept_web(generate, S, C, tol, maxiter, basal_paths)
    target = C * S^2
    halfwidth = max(tol * target, 0.5)
    # a connected web on S species needs at least S - 1 links
    target + halfwidth >= S - 1 || throw(
        ArgumentError(
            "A connected web on $S species needs at least $(S - 1) links, but connectance $C " *
                "with tolerance $tol allows at most $(floor(Int, target + halfwidth))"
        )
    )
    for _ in 1:maxiter
        adj = generate()
        if is_connected(DiGraph(adj)) && abs(count(adj) - target) <= halfwidth &&
                (!basal_paths || isempty(disconnected_from_basal(adj)))
            return PredationMatrix(adj)
        end
    end
    links = isfinite(halfwidth) ?
        "$(max(ceil(Int, target - halfwidth), 0)) to $(floor(Int, target + halfwidth)) links " *
        "(connectance $C, tolerance $tol) " : "any number of links "
    error(
        "Could not generate a connected web on $S species with $links" *
            (basal_paths ? "and feeding paths from basal species to all species " : "") *
            "in $maxiter attempts"
    )
end

# Upper bound on C for models where only one direction of each species pair can be a link
function check_pairwise_connectance(S, C)
    S > 1 || throw(ArgumentError("S must be greater than 1, given $S"))
    return 0 < C <= (S - 1) / (2S) ||
        throw(ArgumentError("C must be in (0, (S - 1) / 2S] = (0, $((S - 1) / (2S))], given $C"))
end

function check_niche_connectance(S, C)
    S > 1 || throw(ArgumentError("S must be greater than 1, given $S"))
    return 0 < C < 0.5 || throw(ArgumentError("C must be in (0, 0.5), given $C"))
end

"""
    may_network([rng], S, C; tol = 0.03, maxiter = 10_000, basal_paths = false)

Generate a random `S` species web with connectance `C`, following:
May, R. M. Will a large complex system be stable? Nature 238, 413–414 (1972).

Each pair of species is linked with probability `2CS / (S - 1)`, and the direction of the
link (who eats whom) is chosen at random. There are no self links or mutual predation. As
these are random graphs rather than food web models, all but sparse webs usually have no
basal species, so `basal_paths = true` is rarely satisfiable.

Webs are resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` of `C S^2`, widened
for small webs to include the whole number of links nearest `C S^2`. With `basal_paths = true`
they are also resampled until every species has a feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)).
"""
function may_network(
        rng::AbstractRNG, S::Int, C::Real; tol = 0.03, maxiter = 10_000,
        basal_paths = false
    )
    check_pairwise_connectance(S, C)
    p = 2C * S / (S - 1)
    return accept_web(S, C, tol, maxiter, basal_paths) do
        adj = zeros(Bool, S, S)
        for i in 1:(S - 1), j in (i + 1):S
            if rand(rng) < p
                rand(rng, Bool) ? (adj[i, j] = true) : (adj[j, i] = true)
            end
        end
        adj
    end
end
may_network(S::Int, C::Real; kwargs...) = may_network(Random.default_rng(), S, C; kwargs...)

"""
    cascade_network([rng], S, C; tol = 0.03, maxiter = 10_000, basal_paths = false)

Generate a cascade model web with `S` species and connectance `C`, following:
Cohen, J. E. & Newman, C. M. A stochastic theory of community food webs I. Models and
aggregated data. Proc. R. Soc. Lond. B 224, 421–448 (1985).

Each species eats each lower ranked species with probability `2CS / (S - 1)`.

Webs are resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` of `C S^2`, widened
for small webs to include the whole number of links nearest `C S^2`. With `basal_paths = true`
they are also resampled until every species has a feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)).
"""
function cascade_network(
        rng::AbstractRNG, S::Int, C::Real; tol = 0.03, maxiter = 10_000,
        basal_paths = false
    )
    check_pairwise_connectance(S, C)
    p = 2C * S / (S - 1)
    return accept_web(S, C, tol, maxiter, basal_paths) do
        adj = zeros(Bool, S, S)
        for i in 1:(S - 1), j in (i + 1):S
            if rand(rng) < p
                adj[i, j] = true
            end
        end
        adj
    end
end
cascade_network(S::Int, C::Real; kwargs...) = cascade_network(Random.default_rng(), S, C; kwargs...)

"""
    generalized_cascade_network([rng], S, C; tol = 0.03, maxiter = 10_000, basal_paths = false)

Generate a generalized cascade model web with `S` species and connectance `C`, following:
Stouffer, D. B., Camacho, J., Guimerà, R., Ng, C. A. & Amaral, L. A. N. Quantitative patterns
in the structure of model and empirical food webs. Ecology 86, 1301–1311 (2005).

Each species `j` draws a single probability `x_j ~ Beta(1, β)` of eating each lower ranked
species, with `β` chosen so that `E[x] = 2CS / (S - 1)`.

Webs are resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` of `C S^2`, widened
for small webs to include the whole number of links nearest `C S^2`. With `basal_paths = true`
they are also resampled until every species has a feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)).
"""
function generalized_cascade_network(
        rng::AbstractRNG, S::Int, C::Real; tol = 0.03,
        maxiter = 10_000, basal_paths = false
    )
    check_pairwise_connectance(S, C)
    β = (S - 1) / (2C * S) - 1
    return accept_web(S, C, tol, maxiter, basal_paths) do
        # drawing x per pair instead would marginalize to the plain cascade model
        x = rand(rng, Beta(1, β), S)
        adj = zeros(Bool, S, S)
        for i in 1:(S - 1), j in (i + 1):S
            if rand(rng) < x[j]
                adj[i, j] = true
            end
        end
        adj
    end
end
function generalized_cascade_network(S::Int, C::Real; kwargs...)
    return generalized_cascade_network(Random.default_rng(), S, C; kwargs...)
end

# Niche values, ranges and centers for the niche model. Species are sorted by niche value and
# the species with the smallest niche value has range 0, so every web has a basal species.
function niche_parameters(rng::AbstractRNG, S::Int, C::Real)
    β = 1 / (2C) - 1
    η = sort(rand(rng, S))
    r = η .* rand(rng, Beta(1, β), S)
    r[1] = 0.0
    centers = [rand(rng, Uniform(r[i] / 2, η[i])) for i in 1:S]
    return η, r, centers
end

"""
    niche_network([rng], S, C; tol = 0.03, maxiter = 10_000, basal_paths = false)

Generate a niche model web with `S` species and connectance `C`, following:
Williams, R. J. & Martinez, N. D. Simple rules yield complex food webs. Nature 404, 180–183
(2000).

Each species `i` has a niche value `η_i ~ U(0, 1)` and eats every species whose niche value
lies in the interval of width `r_i = x η_i`, `x ~ Beta(1, 1/2C - 1)`, centered at
`c_i ~ U(r_i / 2, η_i)`. Species are ordered by niche value.

Webs are resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` of `C S^2`, widened
for small webs to include the whole number of links nearest `C S^2`. With `basal_paths = true`
they are also resampled until every species has a feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)).
"""
function niche_network(
        rng::AbstractRNG, S::Int, C::Real; tol = 0.03, maxiter = 10_000,
        basal_paths = false
    )
    check_niche_connectance(S, C)
    return accept_web(S, C, tol, maxiter, basal_paths) do
        η, r, centers = niche_parameters(rng, S, C)
        adj = zeros(Bool, S, S)
        for i in 1:S, j in 1:S
            if centers[i] - r[i] / 2 <= η[j] <= centers[i] + r[i] / 2
                adj[j, i] = true
            end
        end
        adj
    end
end
niche_network(S::Int, C::Real; kwargs...) = niche_network(Random.default_rng(), S, C; kwargs...)

"""
    generalized_niche_network([rng], S, C, diet; tol = 0.03, maxiter = 10_000, basal_paths = false)

Generate a generalized niche model web with `S` species, connectance `C`, and diet
contiguity `diet ∈ [0, 1]`, following:
Stouffer, D. B., Camacho, J. & Amaral, L. A. N. A robust measure of food web intervality.
PNAS 103, 19015–19020 (2006).

Each species keeps a fraction `diet` of its niche model range as a contiguous interval, and
makes up the missing links with cascade style links to random lower ranked species. With
`diet = 1` this is the niche model without cannibalism, and with `diet = 0` it is close to
the generalized cascade model.

Webs are resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` of `C S^2`, widened
for small webs to include the whole number of links nearest `C S^2`. With `basal_paths = true`
they are also resampled until every species has a feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)).
"""
function generalized_niche_network(
        rng::AbstractRNG, S::Int, C::Real, diet::Real; tol = 0.03,
        maxiter = 10_000, basal_paths = false
    )
    check_niche_connectance(S, C)
    0 <= diet <= 1 || throw(ArgumentError("diet must be in [0, 1], given $diet"))
    return accept_web(S, C, tol, maxiter, basal_paths) do
        η, r, _ = niche_parameters(rng, S, C)
        reduced_r = diet * r
        centers = [rand(rng, Uniform(reduced_r[i] / 2, η[i])) for i in 1:S]
        n_missing_prey = round.(Int, (r - reduced_r) * S)
        adj = zeros(Bool, S, S)
        for i in 1:S
            for j in 1:S
                if i != j && centers[i] - reduced_r[i] / 2 <= η[j] <= centers[i] + reduced_r[i] / 2
                    adj[j, i] = true
                end
            end
            candidates = [j for j in 1:(i - 1) if adj[j, i] == 0]
            k = min(length(candidates), n_missing_prey[i])
            for j in sample(rng, candidates, k; replace = false)
                adj[j, i] = true
            end
        end
        adj
    end
end
function generalized_niche_network(S::Int, C::Real, diet::Real; kwargs...)
    return generalized_niche_network(Random.default_rng(), S, C, diet; kwargs...)
end

"""
    nested_hierarchy_network([rng], S, C; tol = 0.03, maxiter = 10_000, basal_paths = false)

Generate a nested-hierarchy model web with `S` species and connectance `C`, following:
Cattin, M.-F., Bersier, L.-F., Banašek-Richter, C., Baltensperger, R. & Gabriel, J.-P.
Phylogenetic constraints and adaptation explain food-web structure. Nature 427, 835–839
(2004), as described by Stouffer, D. B., Camacho, J., Guimerà, R., Ng, C. A. & Amaral, L. A. N.
Quantitative patterns in the structure of model and empirical food webs. Ecology 86,
1301–1311 (2005), Appendix B.

Species are ordered by niche value, and each species `i` eats `round(S r_i)` species, with
`r_i` drawn as the niche model's range. Species pick their prey in order of niche value:

1. The first prey `j` is a random species with a lower niche value.
2. If another species also eats `j`, further prey come from the prey of `j`'s consumer group:
   the consumers of `j` and the consumers sharing a prey with them. Only species with lower
   niche values than `i` are drawn.
3. Then species with lower niche values that have no predators yet.
4. Then species with niche values greater than or equal to that of `i`, including `i` itself.
5. Then any remaining species with lower niche values, so every species gets its full diet.

Step 5 is not in the published description and is rarely needed. With steps 1–4, about 8% of
links go to species with niche values at least that of the predator, as Stouffer et al.
report (their Fig. 4). Diets are not intervals, unlike the niche model.

Webs are resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` of `C S^2`, widened
for small webs to include the whole number of links nearest `C S^2`. With `basal_paths = true`
they are also resampled until every species has a feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)).
"""
function nested_hierarchy_network(
        rng::AbstractRNG, S::Int, C::Real; tol = 0.03, maxiter = 10_000,
        basal_paths = false
    )
    check_niche_connectance(S, C)
    return accept_web(S, C, tol, maxiter, basal_paths) do
        _, r, _ = niche_parameters(rng, S, C)
        # r ≤ η ≤ 1, so no species needs more than S prey, and the lowest species has r = 0
        n_prey = round.(Int, r * S)
        adj = zeros(Bool, S, S)
        for i in 2:S
            n_prey[i] == 0 && continue
            j = rand(rng, 1:(i - 1))
            adj[j, i] = true
            consumers = [c for c in 1:S if c != i && adj[j, c]]
            if !isempty(consumers)
                group = [c for c in 1:S if c != i && any(adj[p, c] && any(view(adj, p, consumers)) for p in 1:S)]
                pool = [p for p in 1:(i - 1) if any(view(adj, p, group))]
                add_prey!(rng, adj, i, pool, n_prey[i])
            end
            add_prey!(rng, adj, i, [p for p in 1:(i - 1) if !any(view(adj, p, :))], n_prey[i])
            add_prey!(rng, adj, i, i:S, n_prey[i])
            add_prey!(rng, adj, i, 1:(i - 1), n_prey[i])
        end
        adj
    end
end
function nested_hierarchy_network(S::Int, C::Real; kwargs...)
    return nested_hierarchy_network(Random.default_rng(), S, C; kwargs...)
end

# Add random prey of `i` from `candidates` until it has `n_prey` prey or no candidates are left
function add_prey!(rng, adj, i, candidates, n_prey)
    new = [p for p in candidates if !adj[p, i]]
    k = min(length(new), n_prey - count(view(adj, :, i)))
    k > 0 && (adj[sample(rng, new, k; replace = false), i] .= true)
    return adj
end
