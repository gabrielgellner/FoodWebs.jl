module FoodWebsPolytopeSamplersExt

using Distributions: gradlogpdf, logpdf
using FoodWebs: FoodWebs, CommunityPrior, FeasibleCommunities, community_matrix, coordinates,
    draw_growth_rates, growth_rates, halfspaces, interaction_matrix
using PolytopeSamplers: PolytopeSamplers, BilliardWalk, CoordinateHitAndRun, Polytope, Target,
    sample_until_converged, warmup
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

# The log density of a prior and its gradient on the coordinates, without the feasibility check of
# `logdensity`, as the sampler only proposes points of the polytope
function prior_target(P::Polytope, prior::CommunityPrior, fc::FeasibleCommunities)
    dists = [FoodWebs.coordinate_distribution(prior, fc, c) for c in coordinates(fc)]
    active = [k for k in eachindex(dists) if dists[k] !== nothing]
    ℓ(θ) = sum((logpdf(dists[k], θ[k]) for k in active); init = 0.0)
    function ∇ℓ(g, θ)
        fill!(g, 0)
        for k in active
            g[k] = gradlogpdf(dists[k], θ[k])
        end
        return g
    end
    return Target(P, ℓ; gradient = ∇ℓ)
end

function FoodWebs.sample_communities(
        rng::AbstractRNG, fc::FeasibleCommunities, prior::Union{Nothing, CommunityPrior} = nothing;
        target_ess::Real = 1000, sampler = nothing, kwargs...
    )
    P = Polytope(fc)
    model = prior === nothing ? warmup(rng, P) : prior_target(P, prior, fc)
    method = something(sampler, prior === nothing ? BilliardWalk() : CoordinateHitAndRun())
    result = sample_until_converged(rng, model, method; target_ess, kwargs...)
    points = reduce(hcat, result.chains)
    communities = [
        (A = interaction_matrix(fc, θ), r = growth_rates(fc, θ), J = community_matrix(fc, θ))
            for θ in eachcol(points)
    ]
    return (; points, communities, result.diagnostics, result.converged)
end
function FoodWebs.sample_communities(fc::FeasibleCommunities, args...; kwargs...)
    return FoodWebs.sample_communities(Random.default_rng(), fc, args...; kwargs...)
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
