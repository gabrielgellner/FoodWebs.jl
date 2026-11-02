# # The allometric diet breadth model
# Petchey, Beckerman, Riede & Warren (2008) predict each consumer's diet from body masses by
# optimal foraging. Consumer `j` ranks its resources by profitability `E_i / H_ij` and eats the
# `k` most profitable, with `k` maximizing the rate of energy intake of the contingency model,
# `Σ λ_ij E_i / (1 + Σ λ_ij H_ij)` over the diet. Encounter rates `λ_ij = A_ij N_i`, energy
# contents `E_i`, abundances `N_i` and handling times `H_ij` scale with the masses.
#
# The diets depend on fewer parameters than the model has. Scaling every encounter rate by `c`
# and every handling time by `1 / c` scales every intake rate by `c` and leaves the diets
# unchanged, so only the product `a n h` matters; `e` scales every intake rate and does not
# matter; and `λ_ij` depends on the resource mass through `M_i^(ai + ni)`, so only `ai + ni`
# matters.

"""
    RatioHandling(; b, h = 1)

The ratio handling time function of Petchey et al. (2008): consumer `j` handles an individual
of resource `i` in time `H_ij = h / (b - M_i / M_j)` when `M_i / M_j < b` and cannot eat it
otherwise. Handling time rises from `h / b` for tiny resources to infinity at the mass ratio
`b`, so the most profitable resource of consumer `j` has mass `b M_j / 2`.
"""
struct RatioHandling
    h::Float64
    b::Float64
    function RatioHandling(; b::Real, h::Real = 1)
        b > 0 || throw(ArgumentError("b must be positive, given $b"))
        h > 0 || throw(ArgumentError("h must be positive, given $h"))
        return new(h, b)
    end
end

"""
    PowerHandling(; hi, hj, h = 1)

The power handling time function of Petchey et al. (2008): consumer `j` handles an individual
of resource `i` in time `H_ij = h M_i^hi M_j^hj`. Every consumer finds the smallest resource
most profitable when `hi > 1` and the largest when `hi < 1`.
"""
struct PowerHandling
    h::Float64
    hi::Float64
    hj::Float64
    function PowerHandling(; hi::Real, hj::Real, h::Real = 1)
        h > 0 || throw(ArgumentError("h must be positive, given $h"))
        return new(h, hi, hj)
    end
end

handling_time(H::RatioHandling, Mi, Mj) = Mi / Mj < H.b ? H.h / (H.b - Mi / Mj) : Inf
handling_time(H::PowerHandling, Mi, Mj) = H.h * Mi^H.hi * Mj^H.hj

"""
    diet_breadth_network(M; handling, a = 1, ai = 0, aj = 0, e = 1, n = 1, ni = -0.75,
                         N = nothing, consumers = trues(length(M)))
    diet_breadth_network(fw::FoodWeb; kwargs...)

The food web of the allometric diet breadth model, following:
Petchey, O. L., Beckerman, A. P., Riede, J. O. & Warren, P. H. Size, foraging, and food web
structure. PNAS 105, 4191–4196 (2008).

Species `i` has body mass `M[i]`; for a `FoodWeb`, the masses are its `"M"` trait. The model
scales the foraging variables of consumer `j` and resource `i` with their masses:

- abundance `N_i = n M_i^ni`, unless the abundances `N` are given;
- energy content `E_i = e M_i`;
- attack rate `A_ij = a M_i^ai M_j^aj`, and encounter rate `λ_ij = A_ij N_i`;
- handling time `H_ij` from `handling`, a [`RatioHandling`](@ref) or [`PowerHandling`](@ref).

Each consumer ranks the resources it can handle (finite `H_ij`) by profitability `E_i / H_ij`
and eats the `k` most profitable, choosing `k` to maximize its rate of energy intake
`Σ λ_ij E_i / (1 + Σ λ_ij H_ij)` (the contingency model of optimal foraging). Equivalently, it
adds resources in order of profitability while the next one is more profitable than the intake
rate of the diet so far. Ties in profitability go to the lower index.

Every consumer eats at least one resource when it can handle any. Only the species marked in
`consumers` eat, so producers can be left without a diet; for a `FoodWeb` the default is the
species with prey in the web.

The diets depend on `a n h`, `ai + ni` and `aj` and the handling exponents, and not on `e`:
scaling encounter rates up by a factor and handling times down by the same factor scales every
intake rate and leaves the diets unchanged. Larger `a` narrows every diet.
"""
function diet_breadth_network(
        M::AbstractVector{<:Real}; handling::Union{RatioHandling, PowerHandling}, a::Real = 1,
        ai::Real = 0, aj::Real = 0, e::Real = 1, n::Real = 1, ni::Real = -0.75, N = nothing,
        consumers::AbstractVector{Bool} = trues(length(M))
    )
    S = length(M)
    all(m -> m > 0, M) || throw(ArgumentError("body masses must be positive"))
    a > 0 || throw(ArgumentError("a must be positive, given $a"))
    e > 0 || throw(ArgumentError("e must be positive, given $e"))
    length(consumers) == S ||
        throw(DimensionMismatch("$(length(consumers)) consumer flags for $S species"))
    if N === nothing
        n > 0 || throw(ArgumentError("n must be positive, given $n"))
        N = n .* M .^ ni
    else
        length(N) == S || throw(DimensionMismatch("$(length(N)) abundances for $S species"))
        all(x -> x > 0, N) || throw(ArgumentError("abundances must be positive"))
    end
    return diet_web(rank_resources(M, handling, ai, aj, e .* M, N, consumers), a)
end
function diet_breadth_network(fw::FoodWeb; kwargs...)
    return diet_breadth_network(
        body_masses(fw); consumers = [any(view(fw.web.links, :, j)) for j in axes(fw, 2)], kwargs...
    )
end

function body_masses(fw::FoodWeb)
    haskey(fw.traits, "M") || throw(ArgumentError("the web has no body mass trait \"M\""))
    M = fw.traits["M"]
    missing_mass = [fw.species[i] for i in eachindex(M) if ismissing(M[i])]
    isempty(missing_mass) ||
        throw(ArgumentError("species without a body mass: $(join(missing_mass, ", "))"))
    return Float64.(M)
end

# For each consumer, the resources it can handle in order of profitability (ties to the lower
# index), with their encounter rates over `a`, energy contents and handling times. Empty for
# species that do not eat.
function rank_resources(M, handling, ai, aj, E, N, consumers)
    S = length(M)
    return map(1:S) do j
        consumers[j] || return (; index = Int[], λ = Float64[], E = Float64[], H = Float64[])
        H = [handling_time(handling, M[i], M[j]) for i in 1:S]
        index = [i for i in 1:S if isfinite(H[i])]
        sort!(index; by = i -> E[i] / H[i], rev = true, alg = Base.Sort.DEFAULT_STABLE)
        (; index, λ = [M[i]^ai * M[j]^aj * N[i] for i in index], E = E[index], H = H[index])
    end
end

# The length of the diet maximizing Σ λ E / (1 + Σ λ H) over the most profitable resources. The
# intake rate rises while the next resource is more profitable than the current rate and falls
# after, so the diet stops at the first resource no more profitable than the rate.
function diet_length(r, a)
    gain, time = 0.0, 1.0
    k = 0
    for (λ, E, H) in zip(r.λ, r.E, r.H)
        k > 0 && !(E / H > gain / time) && break
        gain += a * λ * E
        time += a * λ * H
        k += 1
    end
    return k
end

diet_links(ranked, a) = sum(r -> diet_length(r, a), ranked)

function diet_web(ranked, a)
    S = length(ranked)
    adj = zeros(Bool, S, S)
    for (j, r) in enumerate(ranked)
        adj[r.index[1:diet_length(r, a)], j] .= true
    end
    return PredationMatrix(adj)
end

"""
    fit_diet_breadth([rng], fw::FoodWeb; handling = :ratio, kwargs...)
    fit_diet_breadth([rng], web, M; handling = :ratio, consumers, starts = 20,
                     steps = 300, ni = -0.75)

Fit the allometric diet breadth model ([`diet_breadth_network`](@ref)) to the observed `web`
with body masses `M`, maximizing the fraction of links predicted correctly, as Petchey et al.
(2008) did. Only the species in `consumers` eat (default: the species with prey in `web`).

The model's number of links falls as the attack rate constant `a` grows, so for every value of
the other parameters the fit sets `a` by bisection to give the number of links `L` nearest the
observed number; the fraction correct is then the number of predicted links present in `web`
over `L`, or over the observed number of links when the prediction has more. The fit searches
the remaining parameters, `ai` and `aj` and then `b` for `handling = :ratio` or `hi` and `hj`
for `handling = :power`, by hill climbing from `starts` random starting points, each for
`steps` random moves. It holds `h = n = e = 1` and the abundance exponent `ni`, since the diets
depend on `a n h` and `ai + ni` only.

Returns a named tuple with the predicted `web`, the fitted `parameters` (the keyword arguments
of `diet_breadth_network`, including `handling`), the number of predicted links `correct` in
`web`, and the `fraction` correct.
"""
function fit_diet_breadth(
        rng::AbstractRNG, web::AbstractMatrix, M::AbstractVector{<:Real};
        handling::Symbol = :ratio, consumers::Union{Nothing, AbstractVector{Bool}} = nothing,
        starts::Int = 20, steps::Int = 300, ni::Real = -0.75
    )
    observed = PredationMatrix(web).links
    S = size(observed, 1)
    length(M) == S || throw(DimensionMismatch("$(length(M)) body masses for $S species"))
    consumers = something(consumers, vec(any(observed; dims = 1)))
    length(consumers) == S ||
        throw(DimensionMismatch("$(length(consumers)) consumer flags for $S species"))
    handling in (:ratio, :power) ||
        throw(ArgumentError("handling must be :ratio or :power, given :$handling"))
    starts > 0 || throw(ArgumentError("starts must be positive, given $starts"))
    L = count(observed)
    L > 0 || throw(ArgumentError("the web has no links"))
    # ai, aj, then log10 b or hi, hj; the boxes for random starts
    box = handling === :ratio ? [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)] :
        [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 2.0), (-1.0, 1.0)]
    model(θ) = handling === :ratio ? RatioHandling(; b = exp10(θ[3])) : PowerHandling(; hi = θ[3], hj = θ[4])
    function evaluate(θ)
        kw = (; handling = model(θ), ai = θ[1], aj = θ[2], ni, consumers)
        ranked = rank_resources(M, kw.handling, kw.ai, kw.aj, M, M .^ ni, consumers)
        loga = match_links(log10a -> diet_links(ranked, exp10(log10a)), L)
        predicted = diet_web(ranked, exp10(loga))
        Lp = count(predicted.links)
        correct = count(predicted.links .& observed)
        return (;
            fraction = correct / max(L, Lp), gap = abs(Lp - L), correct, web = predicted,
            parameters = (; a = exp10(loga), kw...),
        )
    end
    better(x, y) = x.fraction > y.fraction || (x.fraction == y.fraction && x.gap < y.gap)
    # hill climbing that accepts equal moves, so it drifts across plateaus of the fraction
    function climb()
        θ = [lo + (hi - lo) * rand(rng) for (lo, hi) in box]
        current = evaluate(θ)
        for step in 1:steps
            width = 0.5 * (1 - (step - 1) / steps) + 0.01
            θ′ = θ .+ width .* randn(rng, length(θ))
            candidate = evaluate(θ′)
            if !better(current, candidate)
                θ, current = θ′, candidate
            end
        end
        return current
    end
    best = reduce((x, y) -> better(y, x) ? y : x, [climb() for _ in 1:starts])
    return (; best.web, best.parameters, best.correct, best.fraction)
end
function fit_diet_breadth(rng::AbstractRNG, fw::FoodWeb; kwargs...)
    return fit_diet_breadth(rng, fw.web, body_masses(fw); kwargs...)
end
fit_diet_breadth(fw::FoodWeb; kwargs...) = fit_diet_breadth(Random.default_rng(), fw; kwargs...)
function fit_diet_breadth(web::AbstractMatrix, M::AbstractVector{<:Real}; kwargs...)
    return fit_diet_breadth(Random.default_rng(), web, M; kwargs...)
end

# The log10 a whose number of links is nearest L. Links fall as a grows, so bisect for the
# smallest log10 a with at most L links and compare it with the value just below it.
function match_links(links, L; lo = -60.0, hi = 60.0, iterations = 40)
    links(lo) <= L && return lo
    links(hi) > L && return hi
    for _ in 1:iterations
        mid = (lo + hi) / 2
        links(mid) <= L ? (hi = mid) : (lo = mid)
    end
    return L - links(hi) <= links(lo) - L ? hi : lo
end
