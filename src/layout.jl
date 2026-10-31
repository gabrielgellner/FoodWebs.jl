# Node positions for drawing a food web by trophic level: each species' height is its trophic
# level, and species of similar level share a band, ordered across the band to reduce edge
# crossings. Plotting itself is in the Makie extension (`foodwebplot`).

"""
    trophic_layout(web; levels = trophic_levels(web), band = 0.3, sweeps = 20)

Positions `(x, y)` for drawing `web`, where `web[i, j] != 0` means species `j` eats species `i`.
`y` is the trophic level `levels[i]`; prey-averaged levels by default, or any vector of levels.

Species are grouped into bands of levels at most `band` apart: going up the levels, a species
starts a new band when its level is more than `band` above the lowest level in the current band.
Measuring from the band's lowest species keeps bands short for continuous levels such as the
MacKay levels, where chaining each species to the next would merge most of the web into one
band. Each band is spread evenly across `x`, centred on zero with unit spacing, and shifted by
half a unit when it lies within `band` of the band below and would otherwise put nodes at the
same `x`. The order within each band comes from
the barycenter heuristic: `sweeps` passes, alternately up and down the bands, each ordering a
band by the mean `x` of every species' prey and predators in other bands, keeping the ordering
with the fewest edge crossings ([`FoodWebs.edge_crossings`](@ref)). The starting order is the
order of the species in `web`, so the layout is deterministic, and `sweeps = 0` keeps it.
"""
function trophic_layout(
        web::AbstractMatrix; levels::AbstractVector = trophic_levels(web), band::Real = 0.3,
        sweeps::Integer = 20,
    )
    S = LinearAlgebra.checksquare(web)
    length(levels) == S || throw(DimensionMismatch("levels must have length $S"))
    bands = trophic_bands(levels, band)
    neighbours = [[j for j in 1:S if j != i && (web[i, j] != 0 || web[j, i] != 0)] for i in 1:S]
    x = zeros(S)
    place!(x, bands, levels, band)
    positions() = [(x[i], float(levels[i])) for i in 1:S]
    # the heuristic can add crossings on a sweep, so keep the best ordering seen
    best, fewest = positions(), edge_crossings(web, positions())
    for sweep in 1:sweeps
        order = isodd(sweep) ? eachindex(bands) : reverse(eachindex(bands))
        for b in order
            members = bands[b]
            centre(i) = isempty(neighbours[i]) ? x[i] : mean(x[j] for j in neighbours[i])
            sort!(members; by = i -> (centre(i), x[i]))
            place!(x, bands, levels, band)
        end
        n = edge_crossings(web, positions())
        n < fewest && ((best, fewest) = (positions(), n))
    end
    return best
end

# species sorted by level, split where the level rises more than `band` above the band's lowest
function trophic_bands(levels, band)
    order = sortperm(levels)
    bands = [[first(order)]]
    for i in order[2:end]
        if levels[i] - levels[first(bands[end])] > band
            push!(bands, [i])
        else
            push!(bands[end], i)
        end
    end
    return bands
end

# Spread each band evenly across x, centred on zero, in its current order. A band within `band`
# of the band below would put nodes of nearly equal height at the same x when both have the same
# parity, so it is shifted by half a unit to interleave them.
function place!(x, bands, levels, band)
    previous_offset = 0.0
    for (b, members) in enumerate(bands)
        n = length(members)
        offset = mod(-(n + 1) / 2, 1)
        close = b > 1 && minimum(levels[members]) - maximum(levels[bands[b - 1]]) < band
        shift = close && offset == previous_offset ? 0.5 : 0.0
        for (k, i) in enumerate(members)
            x[i] = k - (n + 1) / 2 + shift
        end
        previous_offset = mod(offset + shift, 1)
    end
    return x
end

"""
    edge_crossings(web, positions)

The number of pairs of links of `web` whose straight segments between `positions` cross, not
counting pairs that share a species or cannibalistic links.
"""
function edge_crossings(web::AbstractMatrix, positions)
    S = LinearAlgebra.checksquare(web)
    links = [(i, j) for i in 1:S for j in 1:S if i != j && web[i, j] != 0]
    orientation(a, b, c) = sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
    crosses(p1, p2, q1, q2) = orientation(p1, p2, q1) * orientation(p1, p2, q2) < 0 &&
        orientation(q1, q2, p1) * orientation(q1, q2, p2) < 0
    n = 0
    for (k, (a, b)) in enumerate(links), (c, d) in links[(k + 1):end]
        length(unique((a, b, c, d))) == 4 || continue
        crosses(positions[a], positions[b], positions[c], positions[d]) && (n += 1)
    end
    return n
end

"""
    foodwebplot(web; kwargs...)
    foodwebplot!(ax, web; kwargs...)

Draw `web` with each species at the height of its trophic level and arrows from prey to predator
([`trophic_layout`](@ref)). Requires Makie: load a Makie backend such as CairoMakie to enable it.
`foodwebplot` returns the figure and axis; `foodwebplot!` draws into an existing axis.

Keywords:
- `levels` and `layout`: the trophic levels and node positions (by default from
  `trophic_levels` and `trophic_layout`);
- `color`: node colours, one per species or a single colour (by default by trophic level), with
  `colormap`;
- `markersize`: one per species or a single size (by default shrinking from 22 as `S` grows
  past 30);
- `edgecolor`, `edgealpha`, `edgewidth` and `tipsize`: for the arrows (by default the opacity
  and arrowhead size fall as the number of links grows, so dense webs stay legible). Arrows
  stop at the node edges in screen space, so they meet the nodes whatever the aspect ratio;
- `labels` and `fontsize`: text drawn on each node (`nothing` for none; by default the species
  indices for webs of up to 60 species and none for larger ones) and its size. Labels are
  white on dark nodes and black on light ones;
- `highlight_cannibals`: draw cannibals with a thick outline (default `true`).
"""
function foodwebplot end

"""
    foodwebplot!(ax, web; kwargs...)

Draw `web` into the Makie axis `ax`; see [`foodwebplot`](@ref) for the keywords.
"""
function foodwebplot! end
