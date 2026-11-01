module FoodWebsMakieExt

using FoodWebs: FoodWebs, elliptic_law, trophic_layout, trophic_levels
using LinearAlgebra: LinearAlgebra, eigvals
using Makie: Makie, Axis, DataAspect, Figure, Point2f, Vec2f, arrows2d!, hidespines!, hidexdecorations!,
    limits!, lift, lines!, scatter!, text!, vlines!

function FoodWebs.foodwebplot(web::AbstractMatrix; size = (500, 500), title = "", kwargs...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title, ylabel = "trophic level")
    FoodWebs.foodwebplot!(ax, web; kwargs...)
    return fig, ax
end

# Defaults that thin out the drawing as webs grow: smaller nodes for more species, fainter
# edges and smaller arrowheads for more links
default_markersize(S) = clamp(22 * sqrt(30 / S), 9, 22)
default_edgealpha(L) = clamp(6 / sqrt(max(L, 1)), 0.12, 0.6)
default_tipsize(L) = clamp(10 * sqrt(100 / max(L, 1)), 4, 10)

# Origins and directions of arrows from `a[k]` to `b[k]`, shortened at each end by `ra[k]` and
# `rb[k]` pixels, given the pixels per data unit `scale` along each axis. The gap is computed in
# pixels so that arrows end at the node's edge whatever the axis aspect ratio. Arrows between
# overlapping nodes are dropped.
function arrow_geometry(a, b, ra, rb, scale)
    origins = Point2f[]
    directions = Vec2f[]
    for k in eachindex(a, b)
        A = a[k] .* scale
        B = b[k] .* scale
        d = B .- A
        n = sqrt(sum(abs2, d))
        n > ra[k] + rb[k] || continue
        u = d ./ n
        p = A .+ ra[k] .* u
        q = B .- rb[k] .* u
        push!(origins, Point2f(p ./ scale))
        push!(directions, Vec2f((q .- p) ./ scale))
    end
    return origins, directions
end

# Black or white, whichever reads better on a node of colour `c`
function label_color(c)
    rgb = Makie.to_color(c)
    luminance = 0.2126 * rgb.r + 0.7152 * rgb.g + 0.0722 * rgb.b
    return luminance < 0.5 ? Makie.RGBf(1, 1, 1) : Makie.RGBf(0, 0, 0)
end

# The colour of each node: values mapped through `colormap` over their range, or colours as given
function node_colors(color, colormap, S)
    color isa AbstractVector{<:Real} || return color isa AbstractVector ? color : fill(color, S)
    cmap = Makie.to_colormap(colormap)
    lo, hi = extrema(color)
    t = hi > lo ? (color .- lo) ./ (hi - lo) : fill(0.5, S)
    return [cmap[1 + round(Int, x * (length(cmap) - 1))] for x in t]
end

function FoodWebs.foodwebplot!(
        ax, web::AbstractMatrix; levels = trophic_levels(web), layout = trophic_layout(web; levels),
        color = levels, colormap = :viridis, markersize = nothing, edgecolor = :gray30,
        edgealpha = nothing, edgewidth = 1.2, tipsize = nothing, labels = :auto, fontsize = 10,
        highlight_cannibals = true,
    )
    S = size(web, 1)
    xs, ys = first.(layout), last.(layout)
    links = [(i, j) for i in 1:S, j in 1:S if i != j && web[i, j] != 0]
    L = length(links)
    markersize = something(markersize, default_markersize(S))
    edgealpha = something(edgealpha, default_edgealpha(L))
    tipsize = something(tipsize, default_tipsize(L))
    sizes = markersize isa AbstractVector ? markersize : fill(markersize, S)
    # pad the limits so that nodes at the edges are not clipped
    limits!(ax, minimum(xs) - 0.7, maximum(xs) + 0.7, minimum(ys) - 0.4, maximum(ys) + 0.4)
    # arrows stop at the node edges, recomputed whenever the axis is resized
    a = [Point2f(layout[i]) for (i, _) in links]
    b = [Point2f(layout[j]) for (_, j) in links]
    ra = [sizes[i] / 2 + 1 for (i, _) in links]
    rb = [sizes[j] / 2 + 2 for (_, j) in links]
    geometry = lift(ax.scene.viewport, ax.finallimits) do viewport, limits
        scale = Vec2f(viewport.widths ./ limits.widths)
        arrow_geometry(a, b, ra, rb, scale)
    end
    isempty(links) || arrows2d!(
        ax, lift(first, geometry), lift(last, geometry); color = (edgecolor, edgealpha),
        shaftwidth = edgewidth, tipwidth = tipsize, tiplength = tipsize,
    )
    cannibal = [web[i, i] != 0 for i in 1:S]
    strokewidth = highlight_cannibals ? [c ? 2.5 : 0.75 for c in cannibal] : 0.75
    scatter!(ax, xs, ys; color, colormap, markersize, strokewidth, strokecolor = :black)
    labels === :auto && (labels = S <= 60 ? string.(axes(web, 1)) : nothing)
    labels === nothing || text!(
        ax, xs, ys; text = labels, align = (:center, :center), fontsize,
        color = label_color.(node_colors(color, colormap, S)),
    )
    # only the trophic level axis carries information
    hidespines!(ax, :t, :r, :b)
    hidexdecorations!(ax)
    return ax
end


function FoodWebs.eigenplot(M::AbstractMatrix; size = (500, 500), title = "", kwargs...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title, xlabel = "Re λ", ylabel = "Im λ")
    FoodWebs.eigenplot!(ax, M; kwargs...)
    return fig, ax
end

function FoodWebs.eigenplot!(
        ax, M::AbstractMatrix; ellipse = true, outlier = true, stability = true,
        law = nothing, color = (:steelblue, 0.6), markersize = 5, ellipsecolor = :black,
        outliercolor = :firebrick, equal_aspect = true, expected = nothing, expectedcolor = :darkorange,
    )
    λ = eigvals(M)
    # points the window must contain: the eigenvalues and whatever else is drawn
    xs, ys = collect(real.(λ)), collect(imag.(λ))
    stability && (vlines!(ax, [0.0]; color = :gray60, linestyle = :dash); push!(xs, 0.0); push!(ys, 0.0))
    if ellipse || outlier
        law = something(law, elliptic_law(M))
        if ellipse
            t = range(0, 2π; length = 200)
            ex, ey = law.μ_diag .+ law.x_semi_axis .* cos.(t), law.y_semi_axis .* sin.(t)
            lines!(ax, ex, ey; color = ellipsecolor)
            append!(xs, ex)
            append!(ys, ey)
        end
        if outlier && !isnan(law.outlier)
            scatter!(
                ax, [law.outlier], [0.0]; markersize = 16, color = :transparent,
                strokecolor = outliercolor, strokewidth = 2,
            )
            push!(xs, law.outlier)
            push!(ys, 0.0)
        end
    end
    scatter!(ax, real.(λ), imag.(λ); color, markersize)
    if expected !== nothing
        μ = eigvals(expected)
        scatter!(ax, real.(μ), imag.(μ); marker = :xcross, markersize = 1.6 * markersize, color = expectedcolor)
        append!(xs, real.(μ))
        append!(ys, imag.(μ))
    end
    if equal_aspect
        # a square window with equal scales keeps the ellipse's shape without squeezing a narrow
        # spectrum into a thin axis
        cx, cy = (minimum(xs) + maximum(xs)) / 2, (minimum(ys) + maximum(ys)) / 2
        half = 0.55 * max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys), 1.0e-9)
        limits!(ax, cx - half, cx + half, cy - half, cy + half)
        ax.aspect = DataAspect()
    end
    return ax
end

end
