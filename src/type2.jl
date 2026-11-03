# # Holling type II models
# A multi-prey Holling type II model with the equilibrium and link fluxes of a point of a
# `FeasibleCommunities`. Predator `j` takes prey `i` at the rate
#
#     F_ij(x) = a_ij x_i x_j / D_j(x),  D_j(x) = 1 + Σ_k a_kj h_kj x_k,
#
# and each species' per capita growth rate is
#
#     g_i(x) = r_i - σ_i x_i + Σ_j e_ji a_ji x_j / D_i(x) - Σ_j a_ij x_j / D_j(x),
#
# gains from its prey `j` and losses to its predators `j`. At the equilibrium `x` of the
# Lotka–Volterra point, with loss rates `λ_ij` and handling fractions `τ_j = Σ_k λ_kj h_kj x_k`,
# the attack rates `a_ij = λ_ij / (1 - τ_j)` give `D_j = 1 / (1 - τ_j)` and `F_ij = λ_ij x_i x_j`,
# the Lotka–Volterra fluxes, so the growth rates, self-regulation and efficiencies carry over.

"""
    TypeIIModel

A multi-prey Holling type II model with an equilibrium, built by [`type2_model`](@ref). Its fields
are the equilibrium biomasses `x`, growth rates `r`, self-regulation `self_regulation` (`σ`), and
matrices over `[prey, predator]` of attack rates `attack` (`a`), handling times `handling` (`h`)
and efficiencies `efficiency` (`e`). Predator `j` takes prey `i` at the rate
`F_ij(x) = a_ij x_i x_j / (1 + Σ_k a_kj h_kj x_k)`, and

    dx_i/dt = x_i (r_i - σ_i x_i) + Σ_j e_ji F_ji(x) - Σ_j F_ij(x).

[`community_matrix`](@ref community_matrix(::TypeIIModel)) gives its Jacobian at the equilibrium,
[`handling_fractions`](@ref) how saturated each predator is, and `type2_dynamics` and
`pulse_response` integrate it (load SciMLBase.jl).
"""
struct TypeIIModel
    x::Vector{Float64}
    r::Vector{Float64}
    self_regulation::Vector{Float64}
    attack::Matrix{Float64}
    handling::Matrix{Float64}
    efficiency::Matrix{Float64}
end

Base.show(io::IO, m::TypeIIModel) = print(
    io, "TypeIIModel with ", length(m.x), " species and ", count(>(0), m.attack), " links"
)

# D_j(x) = 1 + Σ_k a_kj h_kj x_k over the present species
denominators(m::TypeIIModel, x, present) = [
    1 + sum(m.attack[k, j] * m.handling[k, j] * x[k] for k in eachindex(x) if present[k]; init = 0.0)
        for j in eachindex(x)
]

# The per capita growth rates g, with absent species neither feeding nor eaten
function type2_per_capita!(g, x, m::TypeIIModel, present)
    D = denominators(m, x, present)
    S = length(x)
    for i in 1:S
        gi = m.r[i] - m.self_regulation[i] * x[i]
        for j in 1:S
            present[j] || continue
            m.attack[j, i] > 0 && (gi += m.efficiency[j, i] * m.attack[j, i] * x[j] / D[i])
            m.attack[i, j] > 0 && (gi -= m.attack[i, j] * x[j] / D[j])
        end
        g[i] = gi
    end
    return g
end

# The derivatives G[i, l] = ∂g_i/∂x_l over the present species
function type2_per_capita_jacobian!(G, x, m::TypeIIModel, present)
    D = denominators(m, x, present)
    S = length(x)
    fill!(G, 0)
    for i in 1:S
        G[i, i] -= m.self_regulation[i]
        for j in 1:S
            present[j] || continue
            # gains of predator i from prey j: e_ji a_ji x_j / D_i
            if m.attack[j, i] > 0
                c = m.efficiency[j, i] * m.attack[j, i]
                G[i, j] += c / D[i]
                for l in 1:S
                    present[l] && (G[i, l] -= c * x[j] * m.attack[l, i] * m.handling[l, i] / D[i]^2)
                end
            end
            # losses of prey i to predator j: a_ij x_j / D_j
            if m.attack[i, j] > 0
                G[i, j] -= m.attack[i, j] / D[j]
                for l in 1:S
                    present[l] && (G[i, l] += m.attack[i, j] * x[j] * m.attack[l, j] * m.handling[l, j] / D[j]^2)
                end
            end
        end
    end
    return G
end

"""
    community_matrix(m::TypeIIModel)

The Jacobian of the type II model `m` at its equilibrium, `diag(x) ∂g/∂x` with `g` the per capita
growth rates. Against the Lotka–Volterra point it was built from, saturation adds `F_ij λ_kj h_kj`
to the effect of each prey `k` on each prey `i` of a predator `j` they share (`k = i` included),
and lowers the effect of each prey on its predators; the effect of a predator on its prey is
unchanged.
"""
function community_matrix(m::TypeIIModel)
    S = length(m.x)
    present = trues(S)
    G = type2_per_capita_jacobian!(zeros(S, S), m.x, m, present)
    g = type2_per_capita!(zeros(S), m.x, m, present)
    return Diagonal(g) + Diagonal(m.x) * G
end

"""
    handling_fractions(m::TypeIIModel)

The fraction of time each predator spends handling prey at the equilibrium, `τ_j = 1 - 1 / D_j`
with `D_j = 1 + Σ_k a_kj h_kj x_k`; 0 for species that eat nothing. For a predator with one prey,
`1 - τ_j` is the elasticity of its consumption to its prey's density, the parameter of
generalized models (Gross et al. 2009) that is 1 in Lotka–Volterra models.

Gross, T., Rudolf, L., Levin, S. A. & Dieckmann, U. Generalized models reveal stabilizing factors
in food webs. Science 325, 747–750 (2009).
"""
handling_fractions(m::TypeIIModel) = 1 .- 1 ./ denominators(m, m.x, trues(length(m.x)))

# One value per species from a number or a vector
per_species(v::Real, S, name) = fill(float(v), S)
function per_species(v::AbstractVector, S, name)
    length(v) == S || throw(DimensionMismatch("$name must have length $S"))
    return float.(v)
end

"""
    type2_model(fc::FeasibleCommunities, θ; handling, handling_fraction, max_consumption,
                match_jacobian = false)

The multi-prey Holling type II model ([`TypeIIModel`](@ref)) with the equilibrium biomasses of
`fc` and the link fluxes of its point `θ`. Each link's loss rate `λ_ij` at `θ` gives the attack
rate `a_ij = λ_ij / (1 - τ_j)`, where `τ_j = Σ_k λ_kj h_kj x_k` is the fraction of time predator
`j` spends handling prey at the equilibrium. This reproduces every flux `λ_ij x_i x_j`, so the
growth rates, self-regulation and efficiencies `γ / λ` of `θ` carry over, and handling times
going to 0 recover the Lotka–Volterra model. Give the handling times by exactly one of:

- `handling`: a matrix of handling times over `[prey, predator]`;
- `handling_fraction`: `τ_j`, a number or one per species, with the same handling time for each
  of a predator's prey;
- `max_consumption`: each predator's maximum consumption rate per unit biomass, a number or one per
  species, so `h_kj = 1 / max_consumption[j]`; an allometric rate such as `y M^-0.25` (Yodzis &
  Innes 1992) gives handling from body mass.

Saturation changes the Jacobian (see [`community_matrix`](@ref community_matrix(::TypeIIModel))).
With `match_jacobian = true` the self-regulation and efficiencies are raised to cancel those
changes, with the growth rates that keep the equilibrium: each prey's self-regulation rises by
`Σ_j λ_ij² h_ij x_j`, its positive density dependence from saturation, and predator `j`'s
efficiencies rise to `e_kj + h_kj G_j / (1 - τ_j)`, with `G_j = Σ_i e_ij λ_ij x_i`, which
restores its response to each prey. The Jacobian then differs from the Lotka–Volterra one at `θ`
only by the effects `F_ij λ_kj h_kj` between distinct prey `i` and `k` of a shared predator `j`,
which no self-regulation or efficiency cancels. Where a predator also eats a prey it shares a
predator with, least squares would lower that link's efficiency to absorb the prey–prey effect,
often below 0; the efficiencies here only rise. Throws an `ArgumentError` when some `τ_j ≥ 1`,
when `θ` has a gain on a link with no loss, which no flux can carry, or when a raised efficiency
exceeds 1.

Cannibalism stays in the self-regulation, without saturation.

Yodzis, P. & Innes, S. Body size and consumer-resource dynamics. Am. Nat. 139, 1151–1175 (1992).
"""
function type2_model(
        fc::FeasibleCommunities, θ::AbstractVector{<:Real}; handling = nothing, handling_fraction = nothing,
        max_consumption = nothing, match_jacobian::Bool = false
    )
    length(θ) == length(fc.coordinates) ||
        throw(DimensionMismatch("expected $(length(fc.coordinates)) coordinates, got $(length(θ))"))
    count(!isnothing, (handling, handling_fraction, max_consumption)) == 1 ||
        throw(ArgumentError("give exactly one of handling, handling_fraction and max_consumption"))
    S = length(fc.x)
    x = fc.x
    λ, γ, σ = zeros(S, S), zeros(S, S), zeros(S)
    for (k, c) in enumerate(fc.coordinates)
        c.kind === :loss && (λ[c.i, c.j] = θ[k])
        c.kind === :gain && (γ[c.j, c.i] = θ[k])
        c.kind === :self && (σ[c.i] = θ[k])
    end
    for i in 1:S, j in 1:S
        γ[i, j] > 0 && λ[i, j] == 0 &&
            throw(ArgumentError("the link from $i to $j has a gain and no loss, which no flux can carry"))
    end
    e = [λ[i, j] > 0 ? γ[i, j] / λ[i, j] : 0.0 for i in 1:S, j in 1:S]
    # each predator's consumption per unit biomass at the equilibrium
    intake = [sum(λ[k, j] * x[k] for k in 1:S) for j in 1:S]
    h = if handling !== nothing
        size(handling) == (S, S) || throw(DimensionMismatch("handling must be $S×$S"))
        all(>=(0), handling) || throw(ArgumentError("handling times must be non-negative"))
        [fc.eats[i, j] ? float(handling[i, j]) : 0.0 for i in 1:S, j in 1:S]
    elseif handling_fraction !== nothing
        τ = per_species(handling_fraction, S, "handling_fraction")
        all(t -> 0 <= t < 1, τ) || throw(ArgumentError("handling fractions must be in [0, 1)"))
        [fc.eats[i, j] && intake[j] > 0 ? τ[j] / intake[j] : 0.0 for i in 1:S, j in 1:S]
    else
        y = per_species(max_consumption, S, "max_consumption")
        all(>(0), y) || throw(ArgumentError("maximum consumption rates must be positive"))
        [fc.eats[i, j] ? 1 / y[j] : 0.0 for i in 1:S, j in 1:S]
    end
    τ = [sum(λ[k, j] * h[k, j] * x[k] for k in 1:S) for j in 1:S]
    saturated = findall(>=(1), τ)
    isempty(saturated) || throw(
        ArgumentError("predators $saturated would spend all their time handling prey at the equilibrium (τ ≥ 1)")
    )
    a = λ ./ (1 .- τ')
    r = growth_rates(fc, θ)
    m = TypeIIModel(Vector{Float64}(x), r, σ, a, h, e)
    return match_jacobian ? compensated_type2(m) : m
end

# The model with raised self-regulation and efficiencies that cancel the effects of saturation on
# the diagonal and on each predator's response to its prey, and the growth rates that keep the
# equilibrium. At the equilibrium each prey i gains F_ij λ_ij h_ij / x_i of positive density
# dependence from each predator j, which σ_i absorbs. Predator j's response to prey k,
# x_j λ_kj (e_kj - h_kj G_j) with G_j = Σ_i e_ij λ_ij x_i, returns to x_j λ_kj e_kj with
# e'_kj = e_kj + h_kj G_j / (1 - τ_j), as then G'_j = G_j / (1 - τ_j).
function compensated_type2(m::TypeIIModel)
    S = length(m.x)
    x = m.x
    τ = handling_fractions(m)
    λ = m.attack .* (1 .- τ')
    σ = copy(m.self_regulation)
    for i in 1:S, j in 1:S
        σ[i] += λ[i, j]^2 * m.handling[i, j] * x[j]
    end
    G = [sum(m.efficiency[i, j] * λ[i, j] * x[i] for i in 1:S) for j in 1:S]
    e = [m.attack[k, j] > 0 ? m.efficiency[k, j] + m.handling[k, j] * G[j] / (1 - τ[j]) : 0.0 for k in 1:S, j in 1:S]
    for k in 1:S, j in 1:S
        e[k, j] <= 1 || throw(
            ArgumentError("matching the Jacobian needs an efficiency of $(e[k, j]) on the link from $k to $j, above 1; lower the handling times")
        )
    end
    g = type2_per_capita!(zeros(S), x, TypeIIModel(x, zeros(S), σ, m.attack, m.handling, e), trues(S))
    return TypeIIModel(x, -g, σ, m.attack, m.handling, e)
end

"""
    type2_problem(m::TypeIIModel, x0, tspan)

The `ODEProblem` of the type II model `m` from `x0` over `tspan`, with its analytic Jacobian
(load SciMLBase.jl to enable this method).
"""
function type2_problem end

"""
    type2_dynamics(alg, m::TypeIIModel, x0, tspan; extinction_threshold = 1e-9, kwargs...)

Integrate the type II model `m` from `x0` over `tspan` with the SciML solver `alg`, with
extinctions as in [`lotka_volterra`](@ref), and the same named tuple as the result (load
SciMLBase.jl to enable this method).
"""
function type2_dynamics end
