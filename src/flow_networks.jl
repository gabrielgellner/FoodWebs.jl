# # Flow networks
# Ecological network analysis of a web of flows at steady state: each compartment's inputs and
# inflows equal its outflows, exports and respiration. The indices follow enaR (Borrett & Lau
# 2014), with trophic levels and the omnivory index as in NetIndices (Kones et al. 2009). Flows
# run from row to column, so the fluxes of `fluxes_from_metabolism` over `[prey, predator]` are
# flows.

"""
    FlowNetwork(flows; inputs, exports = 0, respiration = 0, living = true, names, tol = 0.05)

A network of compartments and the flows between them at steady state: `flows[i, j]` is the flow
from compartment `i` to compartment `j`, and each compartment has an input from outside, an export
and a respiration, numbers or one per compartment. `living` marks the compartments that are
organisms; the others, such as detritus, have trophic level 1 and do not prey on what flows into
them. Throws an `ArgumentError` when some compartment's inputs and inflows differ from its
outflows, exports and respiration by more than `tol` of its outflows, the check of enaR's
`ssCheck`; pass `tol = Inf` to skip it.

The indices are [`throughflow`](@ref), [`total_system_throughflow`](@ref),
[`total_system_throughput`](@ref), [`average_path_length`](@ref),
[`finn_cycling_index`](@ref), [`ascendency`](@ref), [`mixed_trophic_impact`](@ref), and the
flow-weighted [`trophic_levels`](@ref) and [`omnivory_index`](@ref).
"""
struct FlowNetwork
    flows::Matrix{Float64}
    inputs::Vector{Float64}
    exports::Vector{Float64}
    respiration::Vector{Float64}
    living::Vector{Bool}
    names::Vector{String}
end

function FlowNetwork(
        flows::AbstractMatrix{<:Real}; inputs, exports = 0.0, respiration = 0.0, living = true,
        names = string.(axes(flows, 1)), tol::Real = 0.05
    )
    n = LinearAlgebra.checksquare(flows)
    per_compartment(v, name) = v isa Real ? fill(float(v), n) :
        (length(v) == n ? float.(v) : throw(DimensionMismatch("$name must have length $n")))
    z, e, r = per_compartment(inputs, "inputs"), per_compartment(exports, "exports"), per_compartment(respiration, "respiration")
    alive = living isa Bool ? fill(living, n) :
        (length(living) == n ? Vector{Bool}(living) : throw(DimensionMismatch("living must have length $n")))
    length(names) == n || throw(DimensionMismatch("names must have length $n"))
    all(>=(0), flows) && all(>=(0), z) && all(>=(0), e) && all(>=(0), r) ||
        throw(ArgumentError("flows, inputs, exports and respiration must be non-negative"))
    net = FlowNetwork(Matrix{Float64}(flows), z, e, r, alive, string.(names))
    inflow = throughflow(net)
    outflow = vec(sum(net.flows; dims = 2)) .+ e .+ r
    unbalanced = [i for i in 1:n if abs(inflow[i] - outflow[i]) > tol * outflow[i]]
    isempty(unbalanced) || throw(
        ArgumentError("compartments $unbalanced are not at steady state: their inflows and outflows differ by more than $tol")
    )
    return net
end

Base.show(io::IO, net::FlowNetwork) = print(
    io, "FlowNetwork with ", length(net.inputs), " compartments and ", count(>(0), net.flows), " flows"
)

"""
    throughflow(net::FlowNetwork)

The throughflow of each compartment, its input plus its inflows from other compartments; at
steady state this equals its outflows, export and respiration.
"""
throughflow(net::FlowNetwork) = net.inputs .+ vec(sum(net.flows; dims = 1))

"""
    total_system_throughflow(net::FlowNetwork)

The sum of the compartments' [`throughflow`](@ref), `TST`.
"""
total_system_throughflow(net::FlowNetwork) = sum(throughflow(net))

"""
    total_system_throughput(net::FlowNetwork)

The sum of every flow, inputs, exports and respiration included, `TST_p` (Ulanowicz 1986), the
total of the flow matrix that [`ascendency`](@ref) extends with the boundary.
"""
total_system_throughput(net::FlowNetwork) = sum(net.flows) + sum(net.inputs) + sum(net.exports) + sum(net.respiration)

"""
    average_path_length(net::FlowNetwork)

The mean number of compartments a unit of input passes through before it leaves,
`TST / Σ inputs` (Finn 1976).
"""
average_path_length(net::FlowNetwork) = total_system_throughflow(net) / sum(net.inputs)

# The integral flow matrix N = (I - G)⁻¹, with G[i, j] the fraction of j's throughflow that flows to
# i; n[i, i] is the expected number of visits of a unit flowing through i
function integral_flow(net::FlowNetwork)
    T = throughflow(net)
    G = [T[j] > 0 ? net.flows[j, i] / T[j] : 0.0 for i in eachindex(T), j in eachindex(T)]
    return inv(I - G)
end

"""
    finn_cycling_index(net::FlowNetwork)

The fraction of the [`total_system_throughflow`](@ref) that cycles, `TST_c / TST` (Finn 1980),
following enaR's `enaFlow`: compartment `i` returns a fraction `1 - 1/n_ii` of its throughflow to
itself, where `n_ii` is the diagonal of the integral flow matrix `(I - G)⁻¹` and `G[i, j]` is the
fraction of `j`'s throughflow that flows to `i`, so `TST_c = Σ_i (1 - 1/n_ii) T_i`. NetIndices'
`FCIb` divides `TST_c` by the [`total_system_throughput`](@ref) instead, which gives a smaller
value.

Finn, J. T. Measures of ecosystem structure and function derived from analysis of flows. J. Theor.
Biol. 56, 363–380 (1976).

Finn, J. T. Flow analysis of models of the Hubbard Brook ecosystem. Ecology 61, 562–571 (1980).
"""
function finn_cycling_index(net::FlowNetwork)
    T = throughflow(net)
    n = LinearAlgebra.diag(integral_flow(net))
    return sum((1 .- 1 ./ n) .* T) / sum(T)
end

"""
    ascendency(net::FlowNetwork)

Ulanowicz's information measures of the flows, in bits (base-2 logarithms, as in enaR), as a named
tuple:

- `capacity`: the development capacity `C = -Σ T_ij log2(T_ij / T..)`;
- `ascendency`: `A = Σ T_ij log2(T_ij T.. / (T_i. T_.j))`, the part of `C` from the organization
  of the flows;
- `overhead`: `C - A`;
- `diversity`: the flow diversity `C / T..`;
- `mutual_information`: the average mutual information `A / T..`.

The sums run over the flow matrix extended with an input row and export and respiration columns,
`T..` is its total, the [`total_system_throughput`](@ref), and `T_i.` and `T_.j` are its row and
column sums.

Ulanowicz, R. E. Growth and Development: Ecosystems Phenomenology (Springer, 1986).
"""
function ascendency(net::FlowNetwork)
    n = length(net.inputs)
    # compartments, then export and respiration columns and an input row
    T = zeros(n + 3, n + 3)
    T[1:n, 1:n] = net.flows
    T[1:n, n + 1] = net.exports
    T[1:n, n + 2] = net.respiration
    T[n + 3, 1:n] = net.inputs
    total = sum(T)
    row, col = vec(sum(T; dims = 2)), vec(sum(T; dims = 1))
    capacity = -sum(T[i, j] * log2(T[i, j] / total) for i in 1:(n + 3), j in 1:(n + 3) if T[i, j] > 0)
    A = sum(T[i, j] * log2(T[i, j] * total / (row[i] * col[j])) for i in 1:(n + 3), j in 1:(n + 3) if T[i, j] > 0)
    return (;
        ascendency = A, capacity, overhead = capacity - A, diversity = capacity / total,
        mutual_information = A / total,
    )
end

# The diet matrix p[i, j], the share of compartment i's inflow from other compartments that comes
# from j, without cannibalism; zero for non-living compartments
function flow_diets(net::FlowNetwork)
    n = length(net.inputs)
    inflow = [sum(net.flows[j, i] for j in 1:n if j != i; init = 0.0) for i in 1:n]
    return [net.living[i] && inflow[i] > 0 && j != i ? net.flows[j, i] / inflow[i] : 0.0 for i in 1:n, j in 1:n]
end

"""
    trophic_levels(net::FlowNetwork)

The flow-weighted trophic level of each compartment: 1 for non-living compartments and those
with no inflow from other compartments, and otherwise 1 plus the mean level of the compartments it
takes inflow from, weighted by those inflows (NetIndices' `TrophInd`). Cannibalistic flows are
left out, as in [`trophic_levels`](@ref) of a web. NetIndices keeps a cannibal's flow to itself in
its diet but not in its level, which lowers the level by the share of that flow: for the blue crab
of Chesapeake Bay, 3.41 against 3.50 here.
"""
function trophic_levels(net::FlowNetwork)
    p = flow_diets(net)
    return (I - p) \ ones(size(p, 1))
end

"""
    omnivory_index(net::FlowNetwork)

The omnivory index of each compartment, the variance of the trophic levels of its prey weighted by
its diet, `Σ_j p_ij (TL_j - (TL_i - 1))²` (Pauly et al. 1993, NetIndices' `TrophInd`): 0 for a
compartment that feeds on one trophic level. The diets and levels are those of
[`trophic_levels`](@ref), without cannibalism.

Pauly, D., Soriano-Bartz, M. L. & Palomares, M. L. D. Improved construction, parametrization and
interpretation of steady-state ecosystem models. In Trophic Models of Aquatic Ecosystems, 1–13
(ICLARM, 1993).
"""
function omnivory_index(net::FlowNetwork)
    p = flow_diets(net)
    TL = (I - p) \ ones(size(p, 1))
    return [sum(p[i, j] * (TL[j] - (TL[i] - 1))^2 for j in eachindex(TL)) for i in eachindex(TL)]
end

"""
    mixed_trophic_impact(net::FlowNetwork)

The mixed trophic impact `M[i, j]` of compartment `i` on compartment `j`, direct and through every
path (Ulanowicz & Puccia 1990), following enaR's `enaMTI`: `M = (I - Q)⁻¹ - I` with the direct
impacts `Q[i, j] = d[i, j] - f[j, i]`. The diet share `d[i, j]` is the flow from `i` to `j` over
`j`'s [`throughflow`](@ref), the benefit of prey `i` to `j`; the host share `f[j, i]` is the flow
from `j` to `i` over `j`'s throughflow less its respiration, the cost to prey `j` of `i`, and is zero
when `i` is not living.

Ulanowicz, R. E. & Puccia, C. J. Mixed trophic impacts in ecosystems. Coenoses 5, 7–16 (1990).
"""
function mixed_trophic_impact(net::FlowNetwork)
    T = throughflow(net)
    n = length(T)
    d = [T[j] > 0 ? net.flows[i, j] / T[j] : 0.0 for i in 1:n, j in 1:n]
    f = [net.living[j] && T[i] > net.respiration[i] ? net.flows[i, j] / (T[i] - net.respiration[i]) : 0.0 for i in 1:n, j in 1:n]
    Q = d - f'
    return inv(I - Q) - I
end
