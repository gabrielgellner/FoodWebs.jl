# # Energetics
# Energy fluxes on the links of a web, estimated top-down from metabolism at equilibrium
# (fluxweb, Gauzens et al. 2019): each consumer's assimilated intake pays for its metabolic
# losses and the biomass its predators take, and splits among its prey by preference. As
# elsewhere, `web[i, j] != 0` means species `j` eats species `i`, and fluxes are over
# `[prey, predator]`.

# Assimilation efficiencies of plant, animal and detritus resources (Lang et al. 2017), the values
# of the fluxweb vignette
const PLANT_EFFICIENCY = 0.545
const ANIMAL_EFFICIENCY = 0.906
const DETRITUS_EFFICIENCY = 0.158

"""
    assimilation_efficiencies(types; plant = 0.545, animal = 0.906, detritus = 0.158)
    assimilation_efficiencies(fw::FoodWeb; plant, animal, detritus)

The efficiency with which each species is assimilated when eaten, from its resource type,
`:plant`, `:animal` or `:detritus`: the fraction of the eaten biomass that pays for the
consumer's metabolism and production. The defaults are the means Lang et al. (2017) estimate
across consumers.

For a `FoodWeb` the types come from the cheddar `"category"` trait: `producer` is a plant,
`invertebrate`, `vert.ecto` and `vert.endo` are animals, and a species with an empty category
and no body mass, such as CPOM and FPOM, is detritus. Throws an `ArgumentError` for any other
species, such as the fungi of the Millstream webs, whose types must then be given.

Lang, B., Ehnes, R. B., Brose, U. & Rall, B. C. Temperature and consumer type dependencies of
energy flows in natural communities. Oikos 126, 1717–1725 (2017).
"""
function assimilation_efficiencies(
        types::AbstractVector{Symbol}; plant::Real = PLANT_EFFICIENCY, animal::Real = ANIMAL_EFFICIENCY,
        detritus::Real = DETRITUS_EFFICIENCY
    )
    efficiency = Dict(:plant => float(plant), :animal => float(animal), :detritus => float(detritus))
    all(e -> 0 < e <= 1, values(efficiency)) || throw(ArgumentError("efficiencies must be in (0, 1]"))
    return map(types) do t
        haskey(efficiency, t) || throw(ArgumentError("resource types are :plant, :animal or :detritus, given :$t"))
        efficiency[t]
    end
end
function assimilation_efficiencies(fw::FoodWeb; kwargs...)
    category = species_trait(fw, "category")
    M = get(fw.traits, "M", fill(missing, size(fw, 1)))
    types = map(eachindex(category)) do i
        c = category[i]
        c == "producer" && return :plant
        c in ("invertebrate", "vert.ecto", "vert.endo") && return :animal
        (ismissing(c) || c == "") && ismissing(M[i]) && return :detritus
        return :unknown
    end
    unknown = fw.species[types .=== :unknown]
    isempty(unknown) || throw(
        ArgumentError("species $(join(unknown, ", ")) have no plant, animal or detritus category; pass their types to assimilation_efficiencies")
    )
    return assimilation_efficiencies(types; kwargs...)
end

# The share W[i, j] of predator j's intake taken from prey i
function diet_shares(web, x, preferences)
    S = size(web, 1)
    w = if preferences === :biomass
        [web[i, j] != 0 ? float(x[i]) : 0.0 for i in 1:S, j in 1:S]
    elseif preferences === :equal
        [web[i, j] != 0 ? 1.0 : 0.0 for i in 1:S, j in 1:S]
    elseif preferences isa AbstractMatrix
        size(preferences) == (S, S) || throw(DimensionMismatch("preferences must be $S×$S, got $(size(preferences))"))
        all(>=(0), preferences) || throw(ArgumentError("preferences must be non-negative"))
        all(preferences[i, j] == 0 || web[i, j] != 0 for i in 1:S, j in 1:S) ||
            throw(ArgumentError("preferences must be zero except on feeding links"))
        Matrix{Float64}(preferences)
    else
        throw(ArgumentError("preferences must be :biomass, :equal or a matrix, given $preferences"))
    end
    total = sum(w; dims = 1)
    for j in 1:S
        any(web[i, j] != 0 for i in 1:S) || continue
        total[j] > 0 || throw(ArgumentError("species $j has prey but no positive preference for any of them"))
        w[:, j] ./= total[j]
    end
    return w
end

"""
    fluxes_from_metabolism(web, x, M; efficiency, intercept = 0.71, exponent = -0.25,
                           preferences = :biomass)
    fluxes_from_metabolism(fw::FoodWeb; efficiency = assimilation_efficiencies(fw), kwargs...)

The biomass flux on each feeding link of `web` at equilibrium, `F[prey, predator]`, the rate at
which the prey loses biomass to the predator, from the biomasses `x` and body masses `M`
(fluxweb, Gauzens et al. 2019). Each consumer `j` loses `X_j = intercept * M_j^exponent` per unit
biomass to metabolism, and its assimilated intake pays for those losses and for what its
predators take:

    Σ_i e[i, j] F[i, j] = X_j x_j + Σ_k F[j, k],

with its intake split among its prey by `F[i, j] = W[i, j] Σ_k F[k, j]`. This is one linear
equation per consumer, solved top-down; basal species need no body mass.

- `efficiency`: the assimilation efficiency, a number, a vector over species (the efficiency with
  which each species is assimilated when eaten, as from [`assimilation_efficiencies`](@ref)), or
  a matrix over `[prey, predator]`.
- `intercept`, `exponent`: the allometric metabolic rate per unit biomass, numbers or one per
  species. The defaults are the average over metabolic types in the fluxweb vignette; Brown et
  al. (2004) give constants by type. The intercept depends on the units of `M` and of time, so
  rescale it for other units.
- `preferences`: the diet shares `W`, `:biomass` (proportional to prey biomass, fluxweb's
  default), `:equal`, or a matrix of weights over `[prey, predator]`, zero off the feeding links,
  that each predator's shares are proportional to; multiply the weights by `x` to also weight by
  biomass.

For a `FoodWeb` the biomasses are [`biomass`](@ref) and the body masses its `"M"` trait. Throws an
`ArgumentError` when the balance has no non-negative solution, as when predators take more than
their prey's assimilated intake can supply.

Passing `F` with the same efficiency to [`FeasibleCommunities`](@ref) as `fixed_fluxes` fixes the
losses and gains of the Lotka–Volterra model. A cannibal's self-regulation `σ_i` then holds its
net cannibalism `(1 - e[i, i]) F[i, i] / x_i²`, and each consumer's growth rate is `-X_j` plus its
self-regulation beyond that. With self-regulation only from cannibalism, the community matrix is
fluxweb's `create.jacob`, `J[i, j] = (e[j, i] F[j, i] - F[i, j]) / x_j`, for webs without
detritus, to which `create.jacob` gives a positive diagonal.

Gauzens, B. et al. fluxweb: An R package to easily estimate energy fluxes in food webs. Methods
Ecol. Evol. 10, 270–279 (2019).

Brown, J. H., Gillooly, J. F., Allen, A. P., Savage, V. M. & West, G. B. Toward a metabolic theory
of ecology. Ecology 85, 1771–1789 (2004).
"""
function fluxes_from_metabolism(
        web::AbstractMatrix, x::AbstractVector, M::AbstractVector; efficiency,
        intercept = 0.71, exponent = -0.25, preferences = :biomass
    )
    S = LinearAlgebra.checksquare(web)
    length(x) == S || throw(DimensionMismatch("x must have length $S"))
    length(M) == S || throw(DimensionMismatch("M must have length $S"))
    all(xi -> !ismissing(xi) && xi > 0, x) || throw(ArgumentError("biomasses must be positive"))
    consumers = [j for j in 1:S if any(web[i, j] != 0 for i in 1:S)]
    all(j -> !ismissing(M[j]) && M[j] > 0, consumers) || throw(ArgumentError("consumers need a positive body mass"))
    E = if efficiency isa Real
        fill(float(efficiency), S, S)
    elseif efficiency isa AbstractVector
        length(efficiency) == S || throw(DimensionMismatch("efficiency must have length $S"))
        repeat(float.(efficiency), 1, S)
    else
        size(efficiency) == (S, S) || throw(DimensionMismatch("efficiency must be $S×$S"))
        Matrix{Float64}(efficiency)
    end
    all(0 < E[i, j] <= 1 for i in 1:S, j in 1:S if web[i, j] != 0) ||
        throw(ArgumentError("efficiencies on feeding links must be in (0, 1]"))
    a = intercept isa Real ? fill(float(intercept), S) : float.(intercept)
    b = exponent isa Real ? fill(float(exponent), S) : float.(exponent)
    length(a) == S && length(b) == S || throw(DimensionMismatch("intercept and exponent must have length $S"))
    W = diet_shares(web, x, preferences)
    # the intake of each consumer, from its balance: assimilated intake less predation equals losses
    assimilated = [sum(W[i, j] * E[i, j] for i in 1:S) for j in consumers]
    losses = [a[j] * float(M[j])^b[j] * x[j] for j in consumers]
    intake = try
        (Diagonal(assimilated) - W[consumers, consumers]) \ losses
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        throw(ArgumentError("the flux balance is singular, so no fluxes meet these losses"))
    end
    F = zeros(S, S)
    F[:, consumers] = W[:, consumers] .* intake'
    all(>=(0), F) || throw(ArgumentError("no non-negative fluxes balance these losses and efficiencies"))
    return F
end
function fluxes_from_metabolism(fw::FoodWeb; efficiency = assimilation_efficiencies(fw), kwargs...)
    return fluxes_from_metabolism(fw.web, biomass(fw), species_trait(fw, "M"); efficiency, kwargs...)
end
