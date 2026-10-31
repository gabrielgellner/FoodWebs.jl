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

# Expected links of the minimum potential niche model whose niche ranges are drawn at
# connectance `Cp`: the `Cp S^2` potential links of the niche model, of which the `E` extremes
# are kept and the rest forbidden with probability `f`. A species with range `r = x η` has about
# `k ~ Binomial(S, r)` species in range and keeps `min(k, 2)` extremes. Averaging over
# `η ~ U(0, 1)` in closed form gives `E[min(k, 2) | x] = 2 - 3A + (1 - x)^S` with
# `A = (1 - (1 - x)^(S + 1)) / ((S + 1) x)`; the average over `x ~ Beta(1, β)` is a midpoint
# rule on its quantiles.
function expected_mpn_links(S, Cp, f; n = 200)
    β = 1 / (2Cp) - 1
    function kept_extremes(x)
        A = -expm1((S + 1) * log1p(-x)) / ((S + 1) * x)
        return 2 - 3A + (1 - x)^S
    end
    E = (S - 1) * sum(kept_extremes(1 - (1 - (k - 0.5) / n)^(1 / β)) for k in 1:n) / n
    return (1 - f) * Cp * S^2 + f * E
end

# The connectance at which to draw the niche ranges so that the minimum potential niche model
# has `C S^2` links on average after forbidding links, found by bisection
function potential_connectance(S, C, f)
    f == 0 && return C
    lo, hi = float(C), 0.5 - 1.0e-9
    expected_mpn_links(S, hi, f) >= C * S^2 || throw(
        ArgumentError(
            "connectance $C is out of reach with f = $f: the niche ranges would need a " *
                "connectance of 0.5 or more"
        )
    )
    for _ in 1:60
        mid = (lo + hi) / 2
        expected_mpn_links(S, mid, f) < C * S^2 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    minimum_potential_niche_network([rng], S, C, f; tol = 0.03, maxiter = 10_000, basal_paths = false)

Generate a minimum potential niche model web with `S` species, connectance `C` and fraction
`f ∈ [0, 1)` of forbidden links, following:
Allesina, S., Alonso, D. & Pascual, M. A general model for food web structure. Science 320,
658–661 (2008).

Each species has a niche interval as in [`niche_network`](@ref). The species at the two ends
of the interval are always prey, and each species inside it is a prey with probability
`1 - f`, so diets need not be intervals. With `f = 0` this is the niche model, and the same rng
gives the same web.

The niche intervals are drawn at a potential connectance above `C`, chosen so that the
expected connectance after forbidding links is `C`. It accounts for the kept ends of each
interval (Johnson et al. 2014, SI 1.5), using a binomial approximation for the number of
species in an interval. Mean connectance is within a few percent of `C` while the potential
connectance stays below about 0.3, beyond which the niche model's intervals lose range past
niche value 1. Throws an `ArgumentError` if `C` needs a potential connectance of 0.5 or more.

Webs are resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` of `C S^2`, widened
for small webs to include the whole number of links nearest `C S^2`. With `basal_paths = true`
they are also resampled until every species has a feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)).

Johnson, S., Domínguez-García, V., Donetti, L. & Muñoz, M. A. Trophic coherence determines
food-web stability. PNAS 111, 17923–17928 (2014).
"""
function minimum_potential_niche_network(
        rng::AbstractRNG, S::Int, C::Real, f::Real; tol = 0.03, maxiter = 10_000,
        basal_paths = false
    )
    check_niche_connectance(S, C)
    0 <= f < 1 || throw(ArgumentError("f must be in [0, 1), given $f"))
    Cp = potential_connectance(S, C, f)
    return accept_web(S, C, tol, maxiter, basal_paths) do
        η, r, centers = niche_parameters(rng, S, Cp)
        adj = zeros(Bool, S, S)
        for i in 1:S
            # species are sorted by niche value, so the interval is a run of indices
            prey = [j for j in 1:S if centers[i] - r[i] / 2 <= η[j] <= centers[i] + r[i] / 2]
            isempty(prey) && continue
            adj[first(prey), i] = true
            adj[last(prey), i] = true
            for j in prey[2:(end - 1)]
                adj[j, i] = f == 0 || rand(rng) >= f
            end
        end
        adj
    end
end
function minimum_potential_niche_network(S::Int, C::Real, f::Real; kwargs...)
    return minimum_potential_niche_network(Random.default_rng(), S, C, f; kwargs...)
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

"""
    probabilistic_niche_probabilities(η, c, r; a = 1, p = 2)

The link probabilities of the probabilistic niche model, following:
Williams, R. J., Anandanadesan, A. & Purves, D. The probabilistic niche model reveals the
niche structure and role of body size in a complex food web. PLoS ONE 5, e12092 (2010).

Returns the `S`×`S` matrix `P` where `P[j, i]` is the probability that species `i` eats species
`j` (prey are rows, predators are columns),

    P[j, i] = a exp(-|(η[j] - c[i]) / (r[i] / 2)|^p),

for niche values `η`, feeding centers `c` and feeding ranges `r`. The paper uses the Gaussian
kernel `p = 2`. The probability peaks at `a` when `η[j] == c[i]`, and as `p` grows the kernel
approaches the niche model's interval of width `r[i]` centered on `c[i]`. A species with
`r[i] == 0` eats nothing. Williams et al. fix `a = 0.9999` when fitting, as `a = 1` gives
links outside the data a log-likelihood of `-Inf`.
"""
function probabilistic_niche_probabilities(
        η::AbstractVector{<:Real}, c::AbstractVector{<:Real}, r::AbstractVector{<:Real};
        a::Real = 1, p::Real = 2
    )
    S = length(η)
    length(c) == length(r) == S ||
        throw(DimensionMismatch("η, c and r must have the same length, given $S, $(length(c)), $(length(r))"))
    all(>=(0), r) || throw(ArgumentError("feeding ranges r must be non-negative"))
    0 <= a <= 1 || throw(ArgumentError("a must be in [0, 1], given $a"))
    p > 0 || throw(ArgumentError("p must be positive, given $p"))
    return [r[i] == 0 ? 0.0 : a * exp(-abs((η[j] - c[i]) / (r[i] / 2))^p) for j in 1:S, i in 1:S]
end

"""
    probabilistic_niche_network([rng], η, c, r; a = 1, p = 2)
    probabilistic_niche_network([rng], S, C; a = 1, p = 2, tol = 0.03, maxiter = 10_000, basal_paths = false)

Generate a probabilistic niche model web (Williams, Anandanadesan & Purves 2010), where each
species `i` eats each species `j` independently with probability `P[j, i]` from
[`probabilistic_niche_probabilities`](@ref).

Given niche values `η`, feeding centers `c` and feeding ranges `r`, it draws one web from
these probabilities.

Given `S` and `C`, it draws the parameters as in [`niche_network`](@ref). The kernel
integrates to `a Γ(1 + 1/p) r[i]` (`a √π r[i] / 2` for the Gaussian), against `r[i]` for the
niche model's interval, so the ranges are drawn at the niche model connectance
`C / (a Γ(1 + 1/p))` to give `C S^2` links on average. This must be below 0.5. Webs are
resampled until they are connected and their number of links `L` satisfies
`abs(L - C S^2) <= max(tol C S^2, 1/2)`: within a relative tolerance `tol` of `C S^2`, widened
for small webs to include the whole number of links nearest `C S^2`. With `basal_paths = true`
they are also resampled until every species has a feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)).
"""
function probabilistic_niche_network(
        rng::AbstractRNG, η::AbstractVector{<:Real}, c::AbstractVector{<:Real},
        r::AbstractVector{<:Real}; a::Real = 1, p::Real = 2
    )
    P = probabilistic_niche_probabilities(η, c, r; a, p)
    return PredationMatrix(rand(rng, size(P)...) .< P)
end
function probabilistic_niche_network(
        rng::AbstractRNG, S::Int, C::Real; a::Real = 1, p::Real = 2, tol = 0.03,
        maxiter = 10_000, basal_paths = false
    )
    check_niche_connectance(S, C)
    0 < a <= 1 || throw(ArgumentError("a must be in (0, 1], given $a"))
    p > 0 || throw(ArgumentError("p must be positive, given $p"))
    Cn = C / (a * gamma(1 + 1 / p))
    Cn < 0.5 || throw(
        ArgumentError(
            "connectance $C is out of reach with a = $a, p = $p: the niche ranges would " *
                "need a connectance of $Cn, at least 0.5"
        )
    )
    return accept_web(S, C, tol, maxiter, basal_paths) do
        η, r, centers = niche_parameters(rng, S, Cn)
        rand(rng, S, S) .< probabilistic_niche_probabilities(η, centers, r; a, p)
    end
end
function probabilistic_niche_network(η::AbstractVector{<:Real}, c::AbstractVector{<:Real}, r::AbstractVector{<:Real}; kwargs...)
    return probabilistic_niche_network(Random.default_rng(), η, c, r; kwargs...)
end
function probabilistic_niche_network(S::Int, C::Real; kwargs...)
    return probabilistic_niche_network(Random.default_rng(), S, C; kwargs...)
end
