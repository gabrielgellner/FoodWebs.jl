"""
    feeding_loops(J; min_length = 2, max_length = 5, both_directions = true, limit = 10^6)

Every feeding loop of the community matrix `J` with `min_length` to `max_length` species,
heaviest first. A loop is a directed cycle `i₁ → i₂ → ⋯ → iₖ → i₁` of nonzero off-diagonal
entries, and its weight is the geometric mean `|J[i₁, i₂] J[i₂, i₃] ⋯ J[iₖ, i₁]|^(1/k)` of its
entries (Neutel, Heesterbeek & de Ruiter 2002). A predator and its prey form a loop of length 2.
Longer loops are counted in both directions around the same species, since each direction takes
a different set of entries: around an omnivory triangle, one direction runs through two top-down
effects and one bottom-up effect, the other through one top-down and two bottom-up.
`both_directions = false` keeps only the heavier direction around each set of species.

Each loop is a named tuple with:
- `species`: the species in loop order, starting from the lowest index;
- `links`: the entries `(i₁, i₂), (i₂, i₃), …, (iₖ, i₁)` of `J` around the loop;
- `length`, `sign` (of the product of its entries) and `weight`.

The number of loops grows exponentially with their length, so `max_length` caps the search, and
more than `limit` loops throws an `ArgumentError`. [`loop_membership`](@ref) and
[`loop_mask`](@ref) give the species and entries of a set of loops, and `loopplot` draws them.
"""
function feeding_loops(
        J::AbstractMatrix; min_length::Integer = 2, max_length::Integer = 5, both_directions::Bool = true,
        limit::Integer = 10^6
    )
    S = LinearAlgebra.checksquare(J)
    2 <= min_length <= max_length ||
        throw(ArgumentError("need 2 ≤ min_length ≤ max_length, got min_length = $min_length and max_length = $max_length"))
    out = [[j for j in 1:S if j != i && J[i, j] != 0] for i in 1:S]
    loops = NamedTuple{
        (:species, :links, :length, :sign, :weight),
        Tuple{Vector{Int}, Vector{Tuple{Int, Int}}, Int, Int, Float64},
    }[]
    path = Int[]
    onpath = falses(S)
    # depth-first search from each start s through species above s only, so that each cycle is
    # found once, from its lowest species
    function extend!(s, v)
        for w in out[v]
            if w == s
                k = length(path)
                k >= min_length && push!(loops, make_loop(J, path))
                length(loops) > limit &&
                    throw(ArgumentError("more than $limit loops up to length $max_length; lower max_length or raise limit"))
            elseif w > s && !onpath[w] && length(path) < max_length
                push!(path, w)
                onpath[w] = true
                extend!(s, w)
                onpath[w] = false
                pop!(path)
            end
        end
        return nothing
    end
    for s in 1:S
        push!(path, s)
        onpath[s] = true
        extend!(s, s)
        onpath[s] = false
        pop!(path)
    end
    sort!(loops; by = l -> l.weight, rev = true)
    both_directions && return loops
    # the reverse of s → a → ⋯ → z is s → z → ⋯ → a; keep the first, heavier, of each pair
    seen = Set{Vector{Int}}()
    return filter(loops) do l
        p = l.species
        key = l.length > 2 && p[2] > p[end] ? [p[1]; reverse(p[2:end])] : p
        key in seen && return false
        push!(seen, key)
        return true
    end
end

function make_loop(J, path)
    k = length(path)
    links = [(path[m], path[mod1(m + 1, k)]) for m in 1:k]
    entries = [J[i, j] for (i, j) in links]
    weight = exp(sum(log ∘ abs, entries) / k)
    return (species = copy(path), links, length = k, sign = Int(prod(sign, entries)), weight)
end

"""
    heaviest_loop(J; min_length = 3, max_length = 8)

The heaviest of the [`feeding_loops`](@ref) of `J` with `min_length` to `max_length` species, or
`nothing` when `J` has none. Loops of two species, a predator and its prey, are left out by
default, as in Neutel, Heesterbeek & de Ruiter (2002). The search drops a partial loop when even
the largest entry of `J` on every remaining step would leave it lighter than the heaviest loop
found so far, so it reaches longer loops than [`feeding_loops`](@ref) can list.
"""
function heaviest_loop(J::AbstractMatrix; min_length::Integer = 3, max_length::Integer = 8)
    S = LinearAlgebra.checksquare(J)
    2 <= min_length <= max_length ||
        throw(ArgumentError("need 2 ≤ min_length ≤ max_length, got min_length = $min_length and max_length = $max_length"))
    out = [[j for j in 1:S if j != i && J[i, j] != 0] for i in 1:S]
    any(!isempty, out) || return nothing
    logw = [J[i, j] == 0 ? -Inf : log(abs(J[i, j])) for i in 1:S, j in 1:S]
    top = maximum(logw[i, j] for i in 1:S for j in out[i])
    best = -Inf
    heaviest = Int[]
    path = Int[]
    onpath = falses(S)
    # P is the sum of the log entries along the path so far
    function extend!(s, v, P)
        for w in out[v]
            if w == s
                k = length(path)
                if k >= min_length && (P + logw[v, s]) / k > best
                    best = (P + logw[v, s]) / k
                    heaviest = copy(path)
                end
            elseif w > s && !onpath[w] && length(path) < max_length
                Q = P + logw[v, w]
                # the best mean log weight of a loop through this path closing after k steps is
                # at most (Q + (k - q) top) / k, which is monotone in k
                q = length(path)
                kmin = max(q + 1, min_length)
                max((Q + (kmin - q) * top) / kmin, (Q + (max_length - q) * top) / max_length) <= best && continue
                push!(path, w)
                onpath[w] = true
                extend!(s, w, Q)
                onpath[w] = false
                pop!(path)
            end
        end
        return nothing
    end
    for s in 1:S
        push!(path, s)
        onpath[s] = true
        extend!(s, s, 0.0)
        onpath[s] = false
        pop!(path)
    end
    return isempty(heaviest) ? nothing : make_loop(J, heaviest)
end

"""
    loop_membership(loops, S)

A Boolean matrix with one row per species and one column per loop, true where the species is in
the loop. Its row sums count the loops each species is in.
"""
function loop_membership(loops::AbstractVector, S::Integer)
    member = falses(S, length(loops))
    for (k, loop) in enumerate(loops)
        member[loop.species, k] .= true
    end
    return member
end

"""
    loop_mask(loops, S)

A Boolean `S × S` matrix, true at the entries of a community matrix on any of `loops`. Entries
with the mask set lie inside the loops and the other nonzero off-diagonal entries outside them,
so `J[loop_mask(loops, S)] .*= c` scales the interaction strengths inside the loops alone.
"""
function loop_mask(loops::AbstractVector, S::Integer)
    mask = falses(S, S)
    for loop in loops, (i, j) in loop.links
        mask[i, j] = true
    end
    return mask
end

"""
    min_self_regulation(J; scale = ones(size(J, 1)), tol = 1e-10)

The smallest `s` for which the off-diagonal entries of `J` with diagonal `-s * scale` are stable:
the minimal intraspecific interaction strength of Neutel, Heesterbeek & de Ruiter (2002), relative
to `scale`, such as each species' specific death rate. The diagonal of `J` is ignored. With equal
`scale` this is the [`spectral_abscissa`](@ref) of the off-diagonal part divided by the scale. In
general it is found by bisection, to relative tolerance `tol`, between 0 and an `s` large enough
to be stable, and is the boundary of the stable range when stability persists as `s` grows. The
off-diagonal part has trace 0, so its spectral abscissa is at least 0; when it is 0, as for a
food chain, any positive `s` is stable and this returns 0.
"""
function min_self_regulation(J::AbstractMatrix; scale::AbstractVector = ones(size(J, 1)), tol::Real = 1.0e-10)
    S = LinearAlgebra.checksquare(J)
    length(scale) == S || throw(DimensionMismatch("scale must have length $S"))
    all(>(0), scale) || throw(ArgumentError("scale must be positive"))
    offdiag = Matrix{Float64}(J)
    offdiag[diagind(offdiag)] .= 0
    abscissa(s) = spectral_abscissa(offdiag - s * Diagonal(scale))
    abscissa(0.0) <= 0 && return 0.0
    hi = 1.0
    while abscissa(hi) >= 0
        hi *= 2
    end
    lo = 0.0
    while hi - lo > tol * hi
        mid = (lo + hi) / 2
        abscissa(mid) < 0 ? (hi = mid) : (lo = mid)
    end
    return hi
end

"""
    loopplot(web, loops; kwargs...)
    loopplot!(ax, web, loops; kwargs...)

Draw `loops`, such as the heaviest few from [`feeding_loops`](@ref), over a faded
[`foodwebplot`](@ref) of `web`. Each loop is a closed chain of arrows through its species in loop
order, in its own colour from `palette` (by default `Makie.wong_colors()`) and offset to the
right of its direction, so loops on the same species stay apart. The arrow from `i` to `j` stands
for the entry `J[i, j]` of the loop, the effect of `j` on `i`. Arrow widths scale with the loop
weights, the heaviest at `linewidth` (default 4), and the `k`th loop sits `spacing * k` pixels
(default 3) to the right. The web's links are faded to `edgealpha` (default 0.08). Requires
Makie: load a Makie backend such as CairoMakie to enable it. `loopplot` returns the figure and
axis, with a legend of the loops' species and weights; `loopplot!` draws into an existing axis.
Other keywords go to [`foodwebplot!`](@ref).
"""
function loopplot end

"""
    loopplot!(ax, web, loops; kwargs...)

Draw `loops` over `web` into the Makie axis `ax`; see [`loopplot`](@ref).
"""
function loopplot! end
