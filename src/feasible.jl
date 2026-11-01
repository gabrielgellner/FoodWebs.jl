# # Feasible communities
# Lotka–Volterra dynamics `dx_i/dt = x_i (r_i + Σ_j A_ij x_j)` have the equilibrium `x > 0` when
# `r + A x = 0`. At a fixed `x` and web this is linear in the interaction strengths and in `r`, so
# the parameters consistent with the equilibrium, together with linear ecological constraints,
# form a polytope. This follows the inverse formulation of ecological-reactivity (`inverse.jl`),
# generalized to typed coordinates, growth rates as variables, cannibalism and species that eat
# each other.
#
# Coordinates. Each feeding link, prey `i` eaten by predator `j`, has two non-negative
# coordinates: the prey's loss `λ`, at `A[i, j] -= λ`, and the predator's gain `γ`, at
# `A[j, i] += γ`. Species that eat each other put a loss and a gain at the same entry, which add.
# A self-regulated species has a self-regulation coordinate `σ`, at `A[i, i] -= σ`; cannibalism
# merges into it, as a cannibal's loss `a x²` and gain `e a x²` net to a self-limitation term.
# With rates given as intervals, each species also has its growth rate as a coordinate.
#
# Equilibrium. Row `i` of `r + A x = 0` holds the coordinates in row `i` of `A`, each times the
# density of the species in its column, plus `r_i`. Every interaction coordinate appears in one
# row only.
#
# Efficiency. With densities in biomass units, prey `i` loses biomass to predator `j` at rate
# `λ x_i x_j` and the predator gains `e λ x_i x_j`, so `γ = e λ`; a band `e_min ≤ e ≤ e_max` is
# `e_min λ ≤ γ ≤ e_max λ`.
#
# Self-regulation share. A species' turnover `T_i` is its per capita gains at equilibrium: its
# growth rate if basal, `Σ_prey x_k γ` if a consumer. A share band `s_min T_i ≤ x_i σ_i ≤ s_max T_i`
# is linear.
#
# Boundedness. Without efficiency bands every coordinate sits in one row, so the polytope is
# bounded exactly when no row has unbounded coordinates of both signs that can grow together;
# a self-regulated top predator with a share below 1 is bounded, as its self-regulation cannot
# outgrow its gains. With efficiency bands, the budget `x ⋅ (A x) = -x ⋅ r` gives
# `Σ_links (λ - γ) x_i x_j + Σ_i σ_i x_i² = x ⋅ r`, every term non-negative, so flux is bounded
# on every link whose band stays below 1; the polytope is unbounded exactly when links whose band
# reaches 1 form a feeding loop not drained by a positive minimum self-regulation share.

"""
    FeasibleCommunities(web, x; rates, self_regulation = :all, efficiency = nothing,
                        self_share = nothing, bounds = (0.0, Inf), per_capita_bounds = nothing)

The Lotka–Volterra parameters on the food web `web` (`web[i, j] != 0` when species `j` eats
species `i`) that have the equilibrium biomasses `x > 0`: the interaction matrix `A` and growth
rates `r` with `r + A x = 0`, under linear ecological constraints. They form a polytope in the
coordinates listed by [`coordinates`](@ref): a loss and a gain magnitude for each feeding link, a
self-regulation magnitude for each self-regulated species, and the growth rates when they are
variables.

- `rates`: a vector of fixed growth rates, or intervals that make the rates coordinates: a named
  tuple `(basal = (lo, hi), consumer = (lo, hi))` by role (basal species eat nothing), or a vector
  of `(lo, hi)` pairs. The default, `(basal = (0.0, 5.0), consumer = (-1.0, 0.0))`, matches the
  draws of the published constrained food web studies.
- `self_regulation`: the species with a self-regulation term: `:all` (the default; a magnitude of
  0 is allowed, so this includes no self-regulation, but each term adds a dimension and changes a
  uniform prior), `:all_but_top` (species no other eats, as in Gellner, McCann & Hastings),
  `:basal`, `:none`, or a vector of booleans. Cannibals always have one, which holds their net
  cannibalism.
- `efficiency`: `nothing`, a number `e`, a band `(e_min, e_max)`, or a matrix or pair of
  matrices over `[prey, predator]`, constraining each link's gain to `e_min λ ≤ γ ≤ e_max λ` with
  `0 ≤ e_min ≤ e_max ≤ 1`.
- `self_share`: `nothing`, a share `s` or a band `(s_min, s_max)` in `[0, 1]` of each
  self-regulated species' turnover taken by its self-regulation, `s_min T_i ≤ x_i σ_i ≤ s_max T_i`.
- `bounds`: `(lo, hi)` for every loss, gain and self-regulation magnitude.
- `per_capita_bounds`: `nothing`, or `(lo, hi)` for the per capita rate of every loss, gain and
  self-regulation term: a coordinate `θ` at `A[i, j]` contributes `θ x[j]` to species `i`'s per
  capita growth rate, so the bound is `lo / x[j] ≤ θ ≤ hi / x[j]`. It applies together with
  `bounds`.

The community matrix `J = diag(x) A` has the eigenvalues of `A diag(x)`, the matrix of per capita
rates, and the equilibrium, the rate intervals, `self_share` and `per_capita_bounds` constrain per
capita rates only. Rescaling `x` by a constant therefore leaves the spectra unchanged, and
without efficiency bands so does any change of `x`. A finite `bounds` caps the per capita rate of
a coordinate at `A[i, j]` at `hi x[j]`, which depends on the scale and the shape of `x`. Use
`per_capita_bounds` to bound the interactions on the scale of the growth rates.

Species that eat each other get a loss and a gain at each of their two entries of `A`; only the
difference enters the equilibrium, so each such pair adds directions that the equilibrium leaves
to the prior. See [`halfspaces`](@ref) for the constraints, [`dimension`](@ref),
[`isbounded`](@ref), [`structurally_singular`](@ref) and [`diagnose`](@ref) for their structure,
and [`interaction_matrix`](@ref), [`growth_rates`](@ref) and [`community_matrix`](@ref) for the
model at a point.
"""
struct FeasibleCommunities
    eats::Matrix{Bool}
    x::Vector{Float64}
    fixed_rates::Union{Nothing, Vector{Float64}}
    rate_bounds::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}
    self_regulated::Vector{Bool}
    efficiency::Union{Nothing, Tuple{Matrix{Float64}, Matrix{Float64}}}
    self_share::Union{Nothing, Tuple{Float64, Float64}}
    bounds::Tuple{Float64, Float64}
    per_capita_bounds::Union{Nothing, Tuple{Float64, Float64}}
    coordinates::Vector{@NamedTuple{kind::Symbol, i::Int, j::Int}}
end

function FeasibleCommunities(
        web::AbstractMatrix, x::AbstractVector{<:Real};
        rates = (basal = (0.0, 5.0), consumer = (-1.0, 0.0)), self_regulation = :all,
        efficiency = nothing, self_share = nothing, bounds = (0.0, Inf), per_capita_bounds = nothing
    )
    S = LinearAlgebra.checksquare(web)
    length(x) == S || throw(DimensionMismatch("x must have length $S, got $(length(x))"))
    all(>(0), x) || throw(ArgumentError("biomasses must be positive"))
    eats = [i != j && web[i, j] != 0 for i in 1:S, j in 1:S]
    cannibal = [web[i, i] != 0 for i in 1:S]
    basal = [!any(eats[:, j]) for j in 1:S]
    top = [!any(eats[i, :]) for i in 1:S]
    fixed_rates, rate_bounds = if rates isa AbstractVector{<:Real}
        length(rates) == S || throw(DimensionMismatch("rates must have length $S"))
        Vector{Float64}(rates), nothing
    elseif rates isa NamedTuple
        (haskey(rates, :basal) && haskey(rates, :consumer)) ||
            throw(ArgumentError("rate intervals by role need basal and consumer"))
        lo = [basal[i] ? first(rates.basal) : first(rates.consumer) for i in 1:S]
        hi = [basal[i] ? last(rates.basal) : last(rates.consumer) for i in 1:S]
        nothing, (float.(lo), float.(hi))
    elseif rates isa AbstractVector
        length(rates) == S || throw(DimensionMismatch("rates must have length $S"))
        nothing, (float.(first.(rates)), float.(last.(rates)))
    else
        throw(ArgumentError("rates must be a vector of rates, a vector of intervals or a named tuple of intervals by role"))
    end
    rate_bounds === nothing || all(rate_bounds[1] .<= rate_bounds[2]) ||
        throw(ArgumentError("rate intervals must have lo ≤ hi"))
    self_regulated = if self_regulation isa AbstractVector
        length(self_regulation) == S || throw(DimensionMismatch("self_regulation must have length $S"))
        Vector{Bool}(self_regulation)
    elseif self_regulation === :all
        fill(true, S)
    elseif self_regulation === :all_but_top
        .!top
    elseif self_regulation === :basal
        basal
    elseif self_regulation === :none
        fill(false, S)
    else
        throw(ArgumentError("self_regulation must be :all, :all_but_top, :basal, :none or a vector of booleans"))
    end
    self_regulated .|= cannibal
    band = efficiency === nothing ? nothing : efficiency_bands(eats, efficiency)
    share = if self_share === nothing
        nothing
    else
        lo, hi = self_share isa Tuple ? self_share : (self_share, self_share)
        0 <= lo <= hi <= 1 || throw(ArgumentError("self_share must satisfy 0 ≤ s_min ≤ s_max ≤ 1, got $self_share"))
        (float(lo), float(hi))
    end
    first(bounds) <= last(bounds) || throw(ArgumentError("bounds must have lo ≤ hi"))
    first(bounds) >= 0 || throw(ArgumentError("magnitudes are non-negative, so the lower bound must be at least 0"))
    per_capita = if per_capita_bounds === nothing
        nothing
    else
        lo, hi = per_capita_bounds
        0 <= lo <= hi || throw(ArgumentError("per_capita_bounds must satisfy 0 ≤ lo ≤ hi, got $per_capita_bounds"))
        (float(lo), float(hi))
    end
    coords = @NamedTuple{kind::Symbol, i::Int, j::Int}[]
    for i in 1:S
        self_regulated[i] && push!(coords, (kind = :self, i = i, j = i))
        for j in 1:S
            # i is the prey of j: its loss; i is the predator of j: its gain
            eats[i, j] && push!(coords, (kind = :loss, i = i, j = j))
            eats[j, i] && push!(coords, (kind = :gain, i = i, j = j))
        end
    end
    rate_bounds === nothing || append!(coords, [(kind = :growth, i = i, j = i) for i in 1:S])
    return FeasibleCommunities(
        eats, Vector{Float64}(x), fixed_rates, rate_bounds, self_regulated, band, share,
        (float(first(bounds)), float(last(bounds))), per_capita, coords
    )
end

# an efficiency on every link, from a number or a matrix over [prey, predator]
link_efficiencies(eats, e::Real) = [eats[i, j] ? float(e) : 0.0 for i in axes(eats, 1), j in axes(eats, 2)]
function link_efficiencies(eats, e::AbstractMatrix)
    size(e) == size(eats) || throw(DimensionMismatch("efficiency matrices must be $(size(eats, 1))×$(size(eats, 1))"))
    return [eats[i, j] ? float(e[i, j]) : 0.0 for i in axes(eats, 1), j in axes(eats, 2)]
end

function efficiency_bands(eats, efficiency)
    S = size(eats, 1)
    lo, hi = efficiency isa Tuple ? (link_efficiencies(eats, efficiency[1]), link_efficiencies(eats, efficiency[2])) :
        (link_efficiencies(eats, efficiency), link_efficiencies(eats, efficiency))
    for i in 1:S, j in 1:S
        eats[i, j] || continue
        0 <= lo[i, j] <= hi[i, j] <= 1 ||
            throw(ArgumentError("efficiency on the link from $i to $j must satisfy 0 ≤ e_min ≤ e_max ≤ 1"))
    end
    return (lo, hi)
end

Base.show(io::IO, fc::FeasibleCommunities) = print(
    io, "FeasibleCommunities on ", length(fc.x), " species with ", length(fc.coordinates),
    " coordinates", fc.rate_bounds === nothing ? " (fixed rates)" : ""
)

"""
    coordinates(fc::FeasibleCommunities)

One named tuple `(kind, i, j)` per coordinate of the polytope, in order: `kind` is `:loss` (prey
`i`'s loss to predator `j`, at `A[i, j]`), `:gain` (predator `i`'s gain from prey `j`, at
`A[i, j]`), `:self` (species `i`'s self-regulation, at `A[i, i]`) or `:growth` (species `i`'s growth
rate, with `j = i`).
"""
coordinates(fc::FeasibleCommunities) = fc.coordinates

# The sign of each interaction coordinate in A
coordinate_sign(c) = c.kind === :gain ? 1 : -1

"""
    halfspaces(fc::FeasibleCommunities)

The polytope of `fc` as `(lower, upper, A, b, C, d)`: `lower .≤ θ .≤ upper`, `A θ .≤ b` and
`C θ .== d`, the form of PolytopeSamplers.jl, with `θ` in the order of [`coordinates`](@ref). The
equilibrium rows `r + A x = 0` come first in `C`, then fixed efficiencies and shares; efficiency
and share bands are in `A`.
"""
function halfspaces(fc::FeasibleCommunities)
    coords = fc.coordinates
    n = length(coords)
    S = length(fc.x)
    index = Dict((c.kind, c.i, c.j) => k for (k, c) in enumerate(coords))
    r_lo, r_hi = something(fc.rate_bounds, (Float64[], Float64[]))
    fixed = fc.fixed_rates
    lower = fill(first(fc.bounds), n)
    upper = fill(last(fc.bounds), n)
    C = zeros(S, n)
    d = zeros(S)
    for (k, c) in enumerate(coords)
        if c.kind === :growth
            C[c.i, k] = 1
            lower[k], upper[k] = r_lo[c.i], r_hi[c.i]
        else
            C[c.i, k] = coordinate_sign(c) * fc.x[c.j]
            if fc.per_capita_bounds !== nothing
                lo, hi = fc.per_capita_bounds
                lower[k] = max(lower[k], lo / fc.x[c.j])
                upper[k] = min(upper[k], hi / fc.x[c.j])
            end
        end
    end
    fixed === nothing || (d .= .-fixed)
    inequalities, equalities = Vector{Float64}[], Vector{Float64}[]
    rhs_in, rhs_eq = Float64[], Float64[]
    if fc.efficiency !== nothing
        lo, hi = fc.efficiency
        for i in 1:S, j in 1:S
            fc.eats[i, j] || continue
            loss, gain = index[(:loss, i, j)], index[(:gain, j, i)]
            row(e) = (r = zeros(n); r[gain] = 1; r[loss] = -e; r)
            if lo[i, j] == hi[i, j]
                push!(equalities, row(hi[i, j]))
                push!(rhs_eq, 0)
            else
                push!(inequalities, row(hi[i, j]))
                push!(rhs_in, 0)
                lo[i, j] > 0 && (push!(inequalities, -row(lo[i, j])); push!(rhs_in, 0))
            end
        end
    end
    if fc.self_share !== nothing
        s_min, s_max = fc.self_share
        for i in 1:S
            fc.self_regulated[i] || continue
            if s_min == s_max
                row, constant = share_row(fc, index, fixed, i, s_max)
                push!(equalities, row)
                push!(rhs_eq, constant)
            else
                row, constant = share_row(fc, index, fixed, i, s_max)
                push!(inequalities, row)
                push!(rhs_in, constant)
                if s_min > 0
                    row, constant = share_row(fc, index, fixed, i, s_min)
                    push!(inequalities, -row)
                    push!(rhs_in, -constant)
                end
            end
        end
    end
    A = isempty(inequalities) ? zeros(0, n) : permutedims(reduce(hcat, inequalities))
    if !isempty(equalities)
        C = vcat(C, permutedims(reduce(hcat, equalities)))
        d = vcat(d, rhs_eq)
    end
    return (; lower, upper, A, b = rhs_in, C, d)
end

# x_i σ_i - s T_i for species i, as coefficients and a constant: T_i is the growth rate of a basal
# species and the gains Σ_prey x_k γ of a consumer
function share_row(fc::FeasibleCommunities, index, fixed, i, s)
    S = length(fc.x)
    row = zeros(length(fc.coordinates))
    row[index[(:self, i, i)]] = fc.x[i]
    constant = 0.0
    if !any(fc.eats[:, i])
        if fixed === nothing
            row[index[(:growth, i, i)]] -= s
        else
            constant = s * fixed[i]
        end
    else
        for k in 1:S
            fc.eats[k, i] && (row[index[(:gain, i, k)]] -= s * fc.x[k])
        end
    end
    return row, constant
end

"""
    interaction_matrix(fc::FeasibleCommunities, θ)

The interaction matrix `A` at the point `θ` of `fc`: each entry the sum of its gain coordinates
minus its loss coordinates, and minus the self-regulation on the diagonal.
"""
function interaction_matrix(fc::FeasibleCommunities, θ::AbstractVector{<:Real})
    length(θ) == length(fc.coordinates) ||
        throw(DimensionMismatch("expected $(length(fc.coordinates)) coordinates, got $(length(θ))"))
    S = length(fc.x)
    A = zeros(S, S)
    for (k, c) in enumerate(fc.coordinates)
        c.kind === :growth || (A[c.i, c.j] += coordinate_sign(c) * θ[k])
    end
    return A
end

"""
    growth_rates(fc::FeasibleCommunities, θ)

The growth rates at the point `θ` of `fc`: the fixed rates, or the growth coordinates of `θ`.
"""
function growth_rates(fc::FeasibleCommunities, θ::AbstractVector{<:Real})
    fc.fixed_rates === nothing || return copy(fc.fixed_rates)
    return [θ[k] for (k, c) in enumerate(fc.coordinates) if c.kind === :growth]
end

"""
    community_matrix(fc::FeasibleCommunities, θ)

The Jacobian of the Lotka–Volterra dynamics at the equilibrium, `diag(x) A`, at the point `θ`.
"""
community_matrix(fc::FeasibleCommunities, θ::AbstractVector{<:Real}) = Diagonal(fc.x) * interaction_matrix(fc, θ)

"""
    fluxes(fc::FeasibleCommunities, θ)

The biomass flux on each feeding link at the point `θ`, `F[prey, predator] = λ x_prey x_predator`,
the rate at which the prey loses biomass to the predator.
"""
function fluxes(fc::FeasibleCommunities, θ::AbstractVector{<:Real})
    S = length(fc.x)
    F = zeros(S, S)
    for (k, c) in enumerate(fc.coordinates)
        c.kind === :loss && (F[c.i, c.j] = θ[k] * fc.x[c.i] * fc.x[c.j])
    end
    return F
end

"""
    dimension(fc::FeasibleCommunities)

The dimension of the affine hull of the equality constraints of `fc`: the number of coordinates
less the rank of the equalities and of coordinates fixed by equal bounds. This is the dimension
of the polytope when it has a point strictly inside its inequalities.
"""
function dimension(fc::FeasibleCommunities)
    P = halfspaces(fc)
    fixed = [P.lower[k] == P.upper[k] for k in eachindex(P.lower)]
    E = vcat(P.C, Matrix{Float64}(I, length(fixed), length(fixed))[fixed, :])
    return length(P.lower) - (isempty(E) ? 0 : LinearAlgebra.rank(E))
end

"""
    isbounded(fc::FeasibleCommunities)

Whether the polytope of `fc` is bounded, decided from its structure without a solver. With finite
`bounds` or `per_capita_bounds` it always is. Without efficiency bands, each coordinate sits in one species' equilibrium
row, and the polytope is bounded unless some row has unbounded coordinates of both signs that can
grow together; a self-regulated top predator whose self-regulation share is below 1 cannot. With
efficiency bands, flux is bounded on every link whose band stays below 1, and the polytope is
unbounded exactly when links whose band reaches 1 form a feeding loop with no species on it
draining a positive minimum self-regulation share.
"""
function isbounded(fc::FeasibleCommunities)
    isfinite(last(fc.bounds)) && return true
    fc.per_capita_bounds !== nothing && isfinite(last(fc.per_capita_bounds)) && return true
    S = length(fc.x)
    s_min, s_max = fc.self_share === nothing ? (0.0, Inf) : fc.self_share
    if fc.efficiency === nothing
        return all(1:S) do i
            kinds = [c.kind for c in fc.coordinates if c.i == i && c.kind !== :growth]
            positive = any(==(:gain), kinds)
            negative = any(k -> k === :loss || k === :self, kinds)
            positive && negative || return true
            # only gains and self-regulation, with a share below 1: a self-regulated top predator
            return !any(==(:loss), kinds) && s_max < 1
        end
    end
    drains = [fc.self_regulated[i] && s_min > 0 for i in 1:S]
    open_links = [fc.eats[i, j] && last(fc.efficiency)[i, j] == 1 && !drains[i] && !drains[j] for i in 1:S, j in 1:S]
    return !has_feeding_loop(open_links)
end

# Kahn's algorithm: a directed graph has a cycle when its species cannot all be ordered with
# every prey before its predators
function has_feeding_loop(eats::AbstractMatrix{Bool})
    S = size(eats, 1)
    prey_left = vec(sum(eats; dims = 1))
    ready = [j for j in 1:S if prey_left[j] == 0]
    ordered = 0
    while !isempty(ready)
        i = pop!(ready)
        ordered += 1
        for j in 1:S
            eats[i, j] || continue
            prey_left[j] -= 1
            prey_left[j] == 0 && push!(ready, j)
        end
    end
    return ordered < S
end

"""
    structurally_singular(fc::FeasibleCommunities)

Whether the interaction matrix of `fc` is singular for every point: whether its pattern of
entries with a coordinate has no perfect matching of rows to columns. Two top predators without
self-regulation that share a single prey are an example: their rows are multiples of the same
unit vector.
"""
function structurally_singular(fc::FeasibleCommunities)
    S = length(fc.x)
    pattern = falses(S, S)
    for c in fc.coordinates
        c.kind === :growth || (pattern[c.i, c.j] = true)
    end
    return !has_perfect_matching(pattern)
end

# Try to match row i to a column through nonzero entries, rematching earlier rows along an
# augmenting path
function augment_matching!(match, nonzero, i, visited)
    for j in axes(nonzero, 2)
        (nonzero[i, j] && !visited[j]) || continue
        visited[j] = true
        if match[j] == 0 || augment_matching!(match, nonzero, match[j], visited)
            match[j] = i
            return true
        end
    end
    return false
end

# a perfect matching of rows to columns through nonzero entries, by augmenting paths
function has_perfect_matching(nonzero::AbstractMatrix{Bool})
    n = size(nonzero, 1)
    match = zeros(Int, n)
    return all(i -> augment_matching!(match, nonzero, i, falses(n)), 1:n)
end

"""
    diagnose(fc::FeasibleCommunities)

Necessary conditions for `fc` to have any point that fail, each as a sentence; empty when none
fails. An empty result does not prove the polytope nonempty, which needs a linear program (the
PolytopeSamplers.jl extension's `isempty`). The checks:

- with efficiency bands, the budget `Σ (λ - γ) x_i x_j + Σ σ_i x_i² = x ⋅ r` has a non-negative
  left side, so `x ⋅ r ≥ 0` must be possible;
- a species that eats nothing must be able to have a growth rate of at least 0, to pay for its
  losses;
- a species with no predators and no self-regulation must be able to have a growth rate of at most
  0, to balance its gains.
"""
function diagnose(fc::FeasibleCommunities)
    S = length(fc.x)
    problems = String[]
    r_lo, r_hi = fc.fixed_rates === nothing ? fc.rate_bounds : (fc.fixed_rates, fc.fixed_rates)
    # biomasses are positive, so x ⋅ r is largest at the upper rates
    if fc.efficiency !== nothing && sum(fc.x .* r_hi) < 0
        push!(problems, "with efficiency bands the budget needs x ⋅ r ≥ 0, but x ⋅ r is negative for every allowed rate")
    end
    for i in 1:S
        gains = any(fc.eats[:, i])
        losses = any(fc.eats[i, :]) || fc.self_regulated[i]
        # row i: gains - losses - self-regulation + r_i = 0
        !gains && r_hi[i] < 0 &&
            push!(problems, "species $i eats nothing, so its losses need a growth rate of at least 0, but its rate is below 0")
        !losses && r_lo[i] > 0 &&
            push!(problems, "species $i has no predators or self-regulation, so its gains need a growth rate of at most 0, but its rate is above 0")
    end
    return problems
end

"""
    CommunityPrior(; loss = nothing, gain = nothing, interaction = nothing,
                   self_regulation = nothing, growth = nothing)

A prior on the coordinates of a [`FeasibleCommunities`](@ref): a univariate distribution for each
kind of coordinate, or `nothing` for uniform. `interaction` sets both `loss` and `gain`. `growth`
is one distribution for every growth rate or a named tuple `(basal = d₁, consumer = d₂)` by role.
[`logdensity`](@ref) gives the density restricted to the polytope.

A density on the growth rates times a uniform distribution on the interaction strengths,
restricted to the polytope, conditions a joint prior on the equilibrium holding: each rate vector
is weighted by the volume of interaction strengths it allows. Drawing the rates first and then
sampling the interaction strengths uniformly does not weight by that volume, so the two differ.
"""
struct CommunityPrior{L, G, S, R}
    loss::L
    gain::G
    self_regulation::S
    growth::R
end
function CommunityPrior(; loss = nothing, gain = nothing, interaction = nothing, self_regulation = nothing, growth = nothing)
    loss = something(loss, interaction, Some(nothing))
    gain = something(gain, interaction, Some(nothing))
    return CommunityPrior(loss, gain, self_regulation, growth)
end

function coordinate_distribution(prior::CommunityPrior, fc::FeasibleCommunities, c)
    c.kind === :loss && return prior.loss
    c.kind === :gain && return prior.gain
    c.kind === :self && return prior.self_regulation
    prior.growth isa NamedTuple || return prior.growth
    return any(fc.eats[:, c.i]) ? prior.growth.consumer : prior.growth.basal
end

"""
    logdensity(prior::CommunityPrior, fc::FeasibleCommunities, θ; atol = 1e-8)

The log density of the point `θ` of `fc` under `prior`, up to a constant: the sum of the log
densities of the coordinates with a distribution, and `-Inf` outside the polytope (to within
`atol`).
"""
function logdensity(prior::CommunityPrior, fc::FeasibleCommunities, θ::AbstractVector{<:Real}; atol::Real = 1.0e-8)
    P = halfspaces(fc)
    inside = all(P.lower .- atol .<= θ .<= P.upper .+ atol) && all(P.A * θ .<= P.b .+ atol) &&
        all(abs.(P.C * θ .- P.d) .<= atol)
    inside || return -Inf
    ℓ = 0.0
    for (k, c) in enumerate(fc.coordinates)
        dist = coordinate_distribution(prior, fc, c)
        dist === nothing || (ℓ += logpdf(dist, θ[k]))
    end
    return ℓ
end

"""
    draw_growth_rates([rng], web; basal = Uniform(0, 5), consumer = Uniform(0, 1))

Growth rates drawn by role: `r_i ~ basal` for species that eat nothing and `r_i = -c`,
`c ~ consumer`, for the rest, as in the published constrained food web studies. Fixed rates
drawn this way and then a uniform distribution on the interaction strengths give a different
prior from rates as coordinates of [`FeasibleCommunities`](@ref) (see [`CommunityPrior`](@ref)).
"""
function draw_growth_rates(
        rng::AbstractRNG, web::AbstractMatrix; basal::Distribution = Uniform(0, 5),
        consumer::Distribution = Uniform(0, 1)
    )
    S = LinearAlgebra.checksquare(web)
    return [any(i != j && web[i, j] != 0 for i in 1:S) ? -rand(rng, consumer) : rand(rng, basal) for j in 1:S]
end
draw_growth_rates(web::AbstractMatrix; kwargs...) = draw_growth_rates(Random.default_rng(), web; kwargs...)

"""
    sample_communities([rng], fc::FeasibleCommunities[, prior::CommunityPrior]; target_ess = 1000,
                       sampler, kwargs...)

Sample the polytope of `fc`, uniformly or from the density of `prior`, with PolytopeSamplers.jl
(load it to enable this method), running chains until their smallest effective sample size
reaches `target_ess` with R̂ below 1.01. The default sampler is the billiard walk for the uniform
distribution and coordinate hit-and-run for a prior (`BilliardWalk` samples only uniform
targets); CRHMC uses the gradient of the prior's log density. Other keywords go to
`PolytopeSamplers.sample_until_converged`.

Returns a named tuple: `points`, a matrix with one sampled point per column after burn-in, pooled
over chains; `communities`, a named tuple `(A, r, J)` for each point; `diagnostics`; and
`converged`. Throws an `ArgumentError` if the polytope is empty or unbounded.
"""
function sample_communities end

"""
    feasible_growth_rates([rng], web, x; maxiter = 1000, basal, consumer, kwargs...)

Draw growth rates with [`draw_growth_rates`](@ref) until the [`FeasibleCommunities`](@ref) with
those fixed rates (built with the other keywords) has a point, checked with a linear program
through PolytopeSamplers.jl (load it to enable this method). Returns the rates `r`, the
`communities` and the number of `draws`, or throws an `ErrorException` after `maxiter` draws.
This reproduces the redraw loop of the published studies; rates as coordinates need no redraws.

With magnitude `bounds`, the scale of `x` sets whether a draw can be feasible: each row of
`r + A x = 0` divides a rate by biomasses, so with biomasses much smaller than the rates the
magnitudes exceed the bounds. `per_capita_bounds` bound the interactions on the scale of the
rates and do not depend on `x`.
"""
function feasible_growth_rates end
