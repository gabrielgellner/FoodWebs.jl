module FoodWebsPolytopeSamplersExt

using Distributions: Exponential, gradlogpdf, rate
using FoodWebs: FoodWebs, CommunityPrior, ExponentialTilt, FeasibleCommunities, community_matrix, coordinates,
    draw_growth_rates, growth_rates, halfspaces, interaction_matrix, isbounded
using PolytopeSamplers: PolytopeSamplers, BilliardWalk, CoordinateHitAndRun, ExponentialTarget, Polytope,
    Target, Warmup, sample_until_converged, warmup
using Statistics: mean
using Random: AbstractRNG, Random

function PolytopeSamplers.Polytope(fc::FeasibleCommunities)
    P = halfspaces(fc)
    return Polytope(P.lower, P.upper; P.A, P.b, P.C, P.d)
end

"""
    isempty(fc::FeasibleCommunities; kwargs...)

Whether the polytope of `fc` has no point, by a linear program through PolytopeSamplers.jl.
"""
Base.isempty(fc::FeasibleCommunities; kwargs...) = isempty(Polytope(fc); kwargs...)

# The log density of a prior and its gradient on the sampled coordinates φ, with θ = scale .* φ,
# without the feasibility check of `logdensity`, as the sampler only proposes points of the
# polytope. The map is diagonal and constant, so the density in φ is that of θ up to a constant.
function prior_target(w::Warmup, prior::CommunityPrior, fc::FeasibleCommunities, scale)
    dists = [FoodWebs.coordinate_distribution(prior, fc, c) for c in coordinates(fc)]
    active = [k for k in eachindex(dists) if dists[k] !== nothing]
    positions = FoodWebs.kind_positions(fc)
    ℓ(φ) = sum((FoodWebs.coordinate_logpdf(dists[k], scale[k] * φ[k], positions[k]) for k in active); init = 0.0)
    gradient(d::ExponentialTilt, x, j) = -FoodWebs.tilt_rate(d, j)
    gradient(d, x, j) = gradlogpdf(d, x)
    function ∇ℓ(g, φ)
        fill!(g, 0)
        for k in active
            g[k] = scale[k] * gradient(dists[k], scale[k] * φ[k], positions[k])
        end
        return g
    end
    return Target(w, ℓ; gradient = ∇ℓ)
end

# The rates c of a prior whose density is exp(-cᵀθ) on the polytope, or `nothing` if it is not:
# uniform coordinates have rate 0, tilts their rate, and exponentials, on coordinates that cannot
# be negative, the rate of the distribution
function exponential_rates(prior::CommunityPrior, fc::FeasibleCommunities)
    lower = halfspaces(fc).lower
    positions = FoodWebs.kind_positions(fc)
    c = zeros(length(coordinates(fc)))
    for (k, coordinate) in enumerate(coordinates(fc))
        d = FoodWebs.coordinate_distribution(prior, fc, coordinate)
        if d isa ExponentialTilt
            c[k] = FoodWebs.tilt_rate(d, positions[k])
        elseif d isa Exponential && lower[k] >= 0
            c[k] = rate(d)
        elseif d !== nothing
            return nothing
        end
    end
    return c
end

# Tilts that do not decrease must have a bounded polytope to normalize
function check_tilts(prior::CommunityPrior, fc::FeasibleCommunities)
    tilts = [prior.loss, prior.gain, prior.self_regulation, prior.growth]
    prior.growth isa NamedTuple && append!(tilts, values(prior.growth))
    any(d -> d isa ExponentialTilt && any(<=(0), d.rate), tilts) && !isbounded(fc) &&
        throw(ArgumentError("an ExponentialTilt with rate ≤ 0 does not normalize on an unbounded polytope"))
    return nothing
end

# A well rounded polytope holds a ball of radius near 1; below this the walk mixes slowly
const MIN_ROUNDED_RADIUS = 0.5

# The warmup of the polytope of `fc` in per capita coordinates (`true`) or in those of `fc`
# (`false`), the scale from its coordinates to those of `fc`, and the coordinates used. The
# ellipsoid rounding fails on some of these polytopes, whose coordinate widths span many orders
# of magnitude: the warmup throws or leaves an inscribed ball far smaller than 1. Each system is
# then tried again with `extremes = true`, whose 2n linear programs fix directions narrower than
# its `atol` and so remove the thinnest ones, and `:auto` falls back to the coordinates of `fc`
# only when both per capita attempts fail, as chains in them can leave the polytope. When every
# attempt is poorly rounded it keeps the roundest; when every attempt throws it rethrows the
# first error, which is also the error of an empty or unbounded polytope.
function feasible_warmup(rng, fc::FeasibleCommunities, per_capita)
    systems = per_capita === :auto ? (true, false) : (per_capita,)
    first_error, best = nothing, nothing
    for pc in systems, extremes in (false, true)
        H = halfspaces(fc; per_capita = pc)
        P = Polytope(H.lower, H.upper; H.A, H.b, H.C, H.d)
        w = try
            warmup(rng, P; extremes)
        catch err
            err isa Union{ArgumentError, ErrorException} || rethrow()
            first_error = something(first_error, err)
            continue
        end
        found = (w, pc ? H.scale : ones(length(H.lower)), pc)
        w.radius >= MIN_ROUNDED_RADIUS && return found
        (best === nothing || w.radius > best[1].radius) && (best = found)
    end
    best === nothing || return best
    throw(first_error)
end

function FoodWebs.sample_communities(
        rng::AbstractRNG, fc::FeasibleCommunities, prior::Union{Nothing, CommunityPrior} = nothing;
        target_ess::Real = 1000, sampler = nothing, exact::Bool = true,
        per_capita::Union{Bool, Symbol} = :auto, kwargs...
    )
    per_capita in (true, false, :auto) ||
        throw(ArgumentError("per_capita must be true, false or :auto, given $per_capita"))
    if prior !== nothing
        FoodWebs.check_tilt_lengths(prior, fc)
        check_tilts(prior, fc)
    end
    w, scale, used = feasible_warmup(rng, fc, per_capita)
    model = if prior === nothing
        w
    else
        c = exact ? exponential_rates(prior, fc) : nothing
        c === nothing ? prior_target(w, prior, fc, scale) : ExponentialTarget(w, c .* scale)
    end
    method = something(sampler, prior === nothing ? BilliardWalk() : CoordinateHitAndRun())
    result = sample_until_converged(rng, model, method; target_ess, kwargs...)
    points = scale .* reduce(hcat, result.chains)
    communities = [
        (A = interaction_matrix(fc, θ), r = growth_rates(fc, θ), J = community_matrix(fc, θ))
            for θ in eachcol(points)
    ]
    return (; points, communities, result.diagnostics, result.converged, per_capita = used)
end
function FoodWebs.sample_communities(fc::FeasibleCommunities, args...; kwargs...)
    return FoodWebs.sample_communities(Random.default_rng(), fc, args...; kwargs...)
end

function FoodWebs.mean_coordinates(rng::AbstractRNG, fc::FeasibleCommunities; target_ess::Real = 1000, kwargs...)
    points = FoodWebs.sample_communities(rng, fc; target_ess, kwargs...).points
    coordinate = vec(mean(points; dims = 2))
    of_kind(kind) = coordinate[findall(c -> c.kind === kind, coordinates(fc))]
    return (; coordinate, loss = of_kind(:loss), gain = of_kind(:gain), self = of_kind(:self), growth = of_kind(:growth))
end
function FoodWebs.mean_coordinates(fc::FeasibleCommunities; kwargs...)
    return FoodWebs.mean_coordinates(Random.default_rng(), fc; kwargs...)
end

function FoodWebs.feasible_growth_rates(
        rng::AbstractRNG, web::AbstractMatrix, x::AbstractVector{<:Real}; maxiter::Integer = 1000,
        basal = nothing, consumer = nothing, kwargs...
    )
    draw_kwargs = (; (basal === nothing ? () : (; basal))..., (consumer === nothing ? () : (; consumer))...)
    for draws in 1:maxiter
        r = draw_growth_rates(rng, web; draw_kwargs...)
        fc = FeasibleCommunities(web, x; rates = r, kwargs...)
        isempty(fc) || return (; r, communities = fc, draws)
    end
    return error("no feasible growth rates in $maxiter draws")
end
function FoodWebs.feasible_growth_rates(web::AbstractMatrix, x::AbstractVector{<:Real}; kwargs...)
    return FoodWebs.feasible_growth_rates(Random.default_rng(), web, x; kwargs...)
end

end
