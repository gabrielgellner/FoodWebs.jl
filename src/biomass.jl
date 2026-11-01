# # Equilibrium biomass
# Biomasses for equilibrium-based community models: observed from body mass and abundance, or
# computed from steady-state energy flow through the web, and measures of the resulting
# biomass pyramids. As elsewhere, `web[i, j] != 0` means species `j` eats species `i`.

"""
    biomass(fw::FoodWeb)

The observed biomass of each species, body mass times numerical abundance, from the traits
`"M"` and `"N"` of `fw`, in the units of the web's catalog entry (`fw.metadata["properties"]`).
Species without a measured mass or abundance, such as detritus, have `missing`. Throws an
`ArgumentError` if the web has no `"M"` or `"N"` trait.
"""
function biomass(fw::FoodWeb)
    (haskey(fw.traits, "M") && haskey(fw.traits, "N")) ||
        throw(ArgumentError("$(fw.name) has no body mass M and abundance N traits"))
    return fw.traits["M"] .* fw.traits["N"]
end

# The share D[i, j] of prey i's production taken by predator j, without cannibalism
function production_shares(web, shares)
    S = size(web, 1)
    eats(pred, prey) = pred != prey && web[prey, pred] != 0
    if shares === :equal
        D = zeros(S, S)
        for i in 1:S
            predators = [j for j in 1:S if eats(j, i)]
            isempty(predators) || (D[i, predators] .= 1 / length(predators))
        end
        return D
    end
    shares isa AbstractMatrix || throw(ArgumentError("shares must be :equal or a matrix, given $shares"))
    size(shares) == (S, S) || throw(DimensionMismatch("shares must be $S×$S, got $(size(shares))"))
    all(>=(0), shares) || throw(ArgumentError("shares must be non-negative"))
    all(shares[i, j] == 0 || eats(j, i) for i in 1:S, j in 1:S) ||
        throw(ArgumentError("shares must be zero except on feeding links"))
    all(<=(1 + 1.0e-12), sum(shares; dims = 2)) ||
        throw(ArgumentError("each prey's shares must sum to at most 1"))
    return Matrix{Float64}(shares)
end

"""
    steady_state_biomass([rng], web; efficiency = 0.1, turnover = 1.0, basal_production = 1.0,
                         shares = :equal, turnover_exponent = -1/4, noise = 0.0)

Equilibrium biomasses from steady-state energy flow through `web`. Each basal species produces
`basal_production`; each consumer `j` produces `efficiency[j]` times the production it takes from
its prey, a share `D[i, j]` of each prey `i`'s production, so that

    P = P_basal + diag(efficiency) Dᵀ P,

and each species' biomass is its production divided by its turnover (production per unit
biomass), `x = P ./ turnover`. No trophic levels are needed, and webs with loops are allowed.

- `efficiency`: the transfer efficiency, production of a consumer per unit of production it
  eats, a number or one per species, in `(0, 1)`.
- `turnover`: a number, one per species, or `:allometric` for `M^turnover_exponent` from the
  body mass trait `"M"` of a [`FoodWeb`](@ref) (metabolic scaling, exponent `-1/4` by default).
- `basal_production`: a number, or one per species (used for basal species).
- `shares`: `:equal` to split each prey's production equally among its predators, or a matrix
  `D[prey, predator]` that is zero off the feeding links and sums to at most 1 over each prey's
  predators.
- `noise`: each biomass is multiplied by `exp(noise * z)` with `z` standard normal from `rng`.

Cannibalistic links are left out. With a constant efficiency `TE` along a food chain and turnover
falling by a factor `k` per level, biomass changes by `TE * k` per level: with allometric
turnover and a constant predator–prey mass ratio `R`, by `TE * R^(1/4)`, so the pyramid inverts
when `TE * R^(1/4) > 1`, the conditions McCauley et al. (2018) give for top-heavy communities.
Throws an `ArgumentError` for species with no feeding path up from a basal species (see
[`disconnected_from_basal`](@ref)), which receive no production.

McCauley, D. J. et al. On the prevalence and dynamics of inverted trophic pyramids and otherwise
top-heavy communities. Ecol. Lett. 21, 439–454 (2018).
"""
function steady_state_biomass(
        rng::AbstractRNG, web::AbstractMatrix; efficiency = 0.1, turnover = 1.0,
        basal_production = 1.0, shares = :equal, turnover_exponent::Real = -1 / 4, noise::Real = 0.0
    )
    S = LinearAlgebra.checksquare(web)
    unfed = disconnected_from_basal(web)
    isempty(unfed) || throw(
        ArgumentError("species $unfed have no feeding path up from a basal species, so they receive no production")
    )
    TE = efficiency isa Real ? fill(float(efficiency), S) : float.(efficiency)
    length(TE) == S || throw(DimensionMismatch("efficiency must have length $S"))
    all(e -> 0 < e < 1, TE) || throw(ArgumentError("efficiencies must be in (0, 1)"))
    τ = if turnover === :allometric
        (web isa FoodWeb && haskey(web.traits, "M")) ||
            throw(ArgumentError("allometric turnover needs a FoodWeb with a body mass trait M"))
        M = web.traits["M"]
        any(ismissing, M) && throw(ArgumentError("allometric turnover needs a body mass for every species"))
        float.(M) .^ turnover_exponent
    else
        turnover isa Real ? fill(float(turnover), S) : float.(turnover)
    end
    length(τ) == S || throw(DimensionMismatch("turnover must have length $S"))
    all(>(0), τ) || throw(ArgumentError("turnover must be positive"))
    basal = Set(j for j in 1:S if !any(i != j && web[i, j] != 0 for i in 1:S))
    production = basal_production isa Real ? fill(float(basal_production), S) : float.(basal_production)
    length(production) == S || throw(DimensionMismatch("basal_production must have length $S"))
    Pb = [j in basal ? production[j] : 0.0 for j in 1:S]
    D = production_shares(web, shares)
    P = (I - Diagonal(TE) * D') \ Pb
    x = P ./ τ
    noise > 0 && (x .*= exp.(noise .* randn(rng, S)))
    return x
end
steady_state_biomass(web::AbstractMatrix; kwargs...) = steady_state_biomass(Random.default_rng(), web; kwargs...)

"""
    biomass_ratios(web, x)

The consumer-to-resource biomass ratio `x[predator] / x[prey]` on each feeding link of `web`, as
named tuples `(prey, predator, ratio)`, leaving out cannibalism. Ratios near or above 1 mark
top-heavy links (McCauley et al. 2018).
"""
function biomass_ratios(web::AbstractMatrix, x::AbstractVector)
    S = LinearAlgebra.checksquare(web)
    length(x) == S || throw(DimensionMismatch("x must have length $S"))
    return [(prey = i, predator = j, ratio = x[j] / x[i]) for j in 1:S for i in 1:S if i != j && web[i, j] != 0]
end

"""
    pyramid_slope(web, x; levels = trophic_levels(web))

The least-squares slope of `log.(x)` against the trophic `levels`, so that `exp` of the slope is
the mean biomass ratio per trophic level: below 1 for an Eltonian pyramid, above 1 for an
inverted one. Biomasses must be positive.
"""
function pyramid_slope(web::AbstractMatrix, x::AbstractVector; levels::AbstractVector = trophic_levels(web))
    length(x) == length(levels) || throw(DimensionMismatch("x and levels must have the same length"))
    all(>(0), x) || throw(ArgumentError("biomasses must be positive"))
    y = log.(x)
    l̄, ȳ = mean(levels), mean(y)
    return sum((levels .- l̄) .* (y .- ȳ)) / sum(abs2, levels .- l̄)
end

"""
    is_pyramid(web, x)

Whether every consumer's biomass is at most each of its prey's, along the feeding links of `web`
(leaving out cannibalism): biomass never increases up a link. Webs with feeding loops cannot
satisfy this strictly unless the species on a loop have equal biomass.
"""
function is_pyramid(web::AbstractMatrix, x::AbstractVector)
    return all(r.ratio <= 1 for r in biomass_ratios(web, x))
end

"""
    is_inverted_pyramid(web, x)

Whether every consumer's biomass is at least each of its prey's, along the feeding links of
`web` (leaving out cannibalism): biomass never decreases up a link.
"""
function is_inverted_pyramid(web::AbstractMatrix, x::AbstractVector)
    return all(r.ratio >= 1 for r in biomass_ratios(web, x))
end
