# Node positions for drawing a food web by trophic level: each species' height is its trophic
# level, and species of similar level share a band, ordered across the band to reduce edge
# crossings. Plotting itself is in the Makie extension (`foodwebplot`).

"""
    trophic_layout(web; levels = trophic_levels(web), band = 0.3, sweeps = 20)

Positions `(x, y)` for drawing `web`, where `web[i, j] != 0` means species `j` eats species `i`.
`y` is the trophic level `levels[i]`; prey-averaged levels by default, or any vector of levels.

Species whose levels lie within `band` of the next lower species share a band, and each band is
spread evenly across `x`, centred on zero with unit spacing. The order within each band comes from
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
    place!(x, bands)
    positions() = [(x[i], float(levels[i])) for i in 1:S]
    # the heuristic can add crossings on a sweep, so keep the best ordering seen
    best, fewest = positions(), edge_crossings(web, positions())
    for sweep in 1:sweeps
        order = isodd(sweep) ? eachindex(bands) : reverse(eachindex(bands))
        for b in order
            members = bands[b]
            centre(i) = isempty(neighbours[i]) ? x[i] : mean(x[j] for j in neighbours[i])
            sort!(members; by = i -> (centre(i), x[i]))
            place!(x, bands)
        end
        n = edge_crossings(web, positions())
        n < fewest && ((best, fewest) = (positions(), n))
    end
    return best
end

# species sorted by level, split where the level rises by more than `band`
function trophic_bands(levels, band)
    order = sortperm(levels)
    bands = [[first(order)]]
    for i in order[2:end]
        if levels[i] - levels[last(bands[end])] > band
            push!(bands, [i])
        else
            push!(bands[end], i)
        end
    end
    return bands
end

# spread each band evenly across x, centred on zero, in its current order
function place!(x, bands)
    for members in bands
        n = length(members)
        for (k, i) in enumerate(members)
            x[i] = k - (n + 1) / 2
        end
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
- `markersize`: one per species or a single size;
- `edgecolor` and `edgewidth`: for the arrows;
- `labels`: text drawn on each node (`nothing` for none; by default the species indices);
- `highlight_cannibals`: draw cannibals with a thick outline (default `true`).
"""
function foodwebplot end

"""
    foodwebplot!(ax, web; kwargs...)

Draw `web` into the Makie axis `ax`; see [`foodwebplot`](@ref) for the keywords.
"""
function foodwebplot! end
