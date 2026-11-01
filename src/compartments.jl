# # Compartments
# Compartments are groups of species that interact more among themselves than with others. As
# elsewhere, `W[i, j]` is the link from prey `i` to predator `j`: 1 or 0 in a web, or a
# non-negative weight such as an energy flux. Cannibalistic links are left out.

"""
    modularity(W, groups; directed = true, γ = 1)

The modularity of the partition `groups` (one label per species) of the web or non-negative
weight matrix `W`, where `W[i, j]` is the link from prey `i` to predator `j`, such as an energy
flux. Cannibalistic links (the diagonal) are left out.

With `directed = true`, the directed modularity of Leicht & Newman (2008), weighted as in
Arenas et al. (2007),

    Q = (1/W) Σ_ij [W_ij − γ s_i^out s_j^in / W] δ(g_i, g_j),

where `s_i^out = Σ_j W_ij`, `s_j^in = Σ_i W_ij` and `W = Σ_ij W_ij`; `Q` is the same for `W` and its
transpose, so the direction convention does not matter. With `directed = false`, the modularity
of Newman & Girvan (2004) on the undirected web `W + Wᵀ`,

    Q = (1/2m) Σ_ij [A_ij − γ k_i k_j / 2m] δ(g_i, g_j),

in which each link counts once, so that two species that eat each other are joined by two links.
`γ` is the resolution of Reichardt & Bornholdt (2006): larger values favour smaller groups. `Q`
is 0 for a single group and for a web without links.

Leicht, E. A. & Newman, M. E. J. Community structure in directed networks. Phys. Rev. Lett. 100,
118703 (2008). Arenas, A., Duch, J., Fernández, A. & Gómez, S. Size reduction of complex networks
preserving modularity. New J. Phys. 9, 176 (2007). Newman, M. E. J. & Girvan, M. Finding and
evaluating community structure in networks. Phys. Rev. E 69, 026113 (2004). Reichardt, J. &
Bornholdt, S. Statistical mechanics of community detection. Phys. Rev. E 74, 016110 (2006).
"""
function modularity(W::AbstractMatrix, groups::AbstractVector; directed::Bool = true, γ::Real = 1)
    S = LinearAlgebra.checksquare(W)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S, got $(length(groups))"))
    γ >= 0 || throw(ArgumentError("the resolution γ must be non-negative, got $γ"))
    w = [i == j ? 0.0 : float(W[i, j]) for i in 1:S, j in 1:S]
    all(>=(0), w) || throw(ArgumentError("weights must be non-negative"))
    directed || (w = w + w')
    total = sum(w)
    total == 0 && return 0.0
    s_out = vec(sum(w; dims = 2))
    s_in = vec(sum(w; dims = 1))
    Q = 0.0
    for j in 1:S, i in 1:S
        groups[i] == groups[j] && (Q += w[i, j] - γ * s_out[i] * s_in[j] / total)
    end
    return Q / total
end
