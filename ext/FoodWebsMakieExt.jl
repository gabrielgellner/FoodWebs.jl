module FoodWebsMakieExt

using FoodWebs: FoodWebs, FoodWeb, block_spectrum, elliptic_law, trophic_layout, trophic_levels
using LinearAlgebra: LinearAlgebra, eigvals
using Makie: Makie, Axis, Colorbar, DataAspect, Figure, Legend, MarkerElement, Point2f, Vec2f, arrows2d!,
    heatmap!, hidespines!, hidexdecorations!, hlines!, limits!, lift, lines!, scatter!, text!, vlines!

function FoodWebs.foodwebplot(
        web::AbstractMatrix; size = (500, 500), title = "", groups = nothing, roles = nothing,
        palette = Makie.wong_colors(), legend = groups !== nothing || roles !== nothing, kwargs...
    )
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title, ylabel = "trophic level")
    FoodWebs.foodwebplot!(ax, web; groups, roles, palette, kwargs...)
    if legend
        elements, names = Any[], String[]
        if groups !== nothing
            colors = group_colors(groups, palette)
            for g in unique(groups)
                count(==(g), groups) > 1 || continue
                push!(elements, MarkerElement(; marker = :circle, color = colors[g], strokecolor = :black, strokewidth = 0.75))
                push!(names, "group $g")
            end
            if any(g -> count(==(g), groups) == 1, groups)
                push!(elements, MarkerElement(; marker = :circle, color = SINGLE_COLOR, strokecolor = :black, strokewidth = 0.75))
                push!(names, "single species")
            end
        end
        if roles !== nothing
            for role in intersect(keys(ROLE_MARKERS), unique(roles))
                push!(elements, MarkerElement(; marker = ROLE_MARKERS[role], color = :white, strokecolor = :black, strokewidth = 0.75))
                push!(names, replace(string(role), "_" => " "))
            end
        end
        isempty(elements) || Legend(fig[1, 2], elements, names; framevisible = false)
    end
    return fig, ax
end

# markers for the species roles of Olesen et al. (2007)
const ROLE_MARKERS = Dict(:peripheral => :circle, :connector => :diamond, :module_hub => :rect, :network_hub => :star5)

# Group colours: groups of two or more species take the colours of `palette` in order of first
# appearance, cycling when there are more groups than colours; single-species groups are grey
const SINGLE_COLOR = Makie.RGBf(0.8, 0.8, 0.8)
function group_colors(groups, palette)
    colors = Dict{eltype(groups), Any}()
    k = 0
    for g in unique(groups)
        if count(==(g), groups) > 1
            k += 1
            colors[g] = palette[mod1(k, length(palette))]
        else
            colors[g] = SINGLE_COLOR
        end
    end
    return colors
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
        ax, web::AbstractMatrix; levels = trophic_levels(web), layout = nothing,
        color = nothing, colormap = :viridis, markersize = nothing, edgecolor = :gray30,
        edgealpha = nothing, edgewidth = 1.2, tipsize = nothing, labels = :auto, fontsize = 10,
        highlight_cannibals = true, groups = nothing, palette = Makie.wong_colors(),
        betweencolor = :firebrick, show_within = true, roles = nothing,
    )
    S = size(web, 1)
    groups === nothing || length(groups) == S || throw(DimensionMismatch("groups must have length $S"))
    roles === nothing || length(roles) == S || throw(DimensionMismatch("roles must have length $S"))
    if groups !== nothing && layout === nothing
        layout = trophic_layout(web; levels, groups)
    end
    layout = something(layout, trophic_layout(web; levels))
    if color === nothing
        color = if groups === nothing
            levels
        else
            colors = group_colors(groups, palette)
            [colors[g] for g in groups]
        end
    end
    xs, ys = first.(layout), last.(layout)
    links = [(i, j) for i in 1:S, j in 1:S if i != j && web[i, j] != 0]
    between = groups === nothing ? falses(length(links)) : [groups[i] != groups[j] for (i, j) in links]
    show_within || (links = links[between]; between = trues(length(links)))
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
    # links between groups are drawn in `betweencolor`, and more opaque
    edgecolors = [b ? Makie.RGBAf(Makie.to_color(betweencolor), min(1, 2edgealpha)) : Makie.RGBAf(Makie.to_color(edgecolor), edgealpha) for b in between]
    isempty(links) || arrows2d!(
        ax, lift(first, geometry), lift(last, geometry); color = edgecolors,
        shaftwidth = edgewidth, tipwidth = tipsize, tiplength = tipsize,
    )
    cannibal = [web[i, i] != 0 for i in 1:S]
    strokewidth = highlight_cannibals ? [c ? 2.5 : 0.75 for c in cannibal] : 0.75
    marker = roles === nothing ? :circle : [get(ROLE_MARKERS, r, :circle) for r in roles]
    scatter!(ax, xs, ys; color, colormap, markersize, strokewidth, strokecolor = :black, marker)
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


function FoodWebs.massplot(
        web::AbstractMatrix, M::AbstractVector, y::AbstractVector; size = (500, 500), title = "",
        ylabel = "", kwargs...
    )
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title, xlabel = "log₁₀ body mass", ylabel)
    FoodWebs.massplot!(ax, web, M, y; kwargs...)
    return fig, ax
end
function FoodWebs.massplot(fw::FoodWeb; y::Symbol = :abundance, kwargs...)
    values, ylabel = mass_plot_values(fw, y)
    return FoodWebs.massplot(fw, fw.traits["M"], values; ylabel, kwargs...)
end

function FoodWebs.massplot!(ax, web::AbstractMatrix, M::AbstractVector, y::AbstractVector; labels = :auto, kwargs...)
    S = LinearAlgebra.checksquare(web)
    length(M) == S && length(y) == S || throw(DimensionMismatch("M and y must have length $S"))
    keep = [i for i in 1:S if !ismissing(M[i]) && !ismissing(y[i])]
    labels === :auto && (labels = length(keep) <= 60 ? string.(keep) : nothing)
    levels = float.(y[keep])
    layout = [(log10(M[i]), float(y[i])) for i in keep]
    FoodWebs.foodwebplot!(ax, web[keep, keep]; levels, layout, labels, kwargs...)
    # foodwebplot! hides the horizontal axis, which here carries the masses
    ax.bottomspinevisible = true
    ax.xticksvisible = true
    ax.xticklabelsvisible = true
    ax.xlabelvisible = true
    return ax
end
function FoodWebs.massplot!(ax, fw::FoodWeb; y::Symbol = :abundance, kwargs...)
    values, _ = mass_plot_values(fw, y)
    return FoodWebs.massplot!(ax, fw, fw.traits["M"], values; kwargs...)
end

function mass_plot_values(fw::FoodWeb, y::Symbol)
    haskey(fw.traits, "M") || throw(ArgumentError("$(fw.name) has no body mass trait \"M\""))
    if y === :abundance
        haskey(fw.traits, "N") || throw(ArgumentError("$(fw.name) has no abundance trait \"N\""))
        return [ismissing(n) ? missing : log10(n) for n in fw.traits["N"]], "log₁₀ abundance"
    elseif y === :level
        return trophic_levels(fw), "trophic level"
    end
    throw(ArgumentError("y must be :abundance or :level, given :$y"))
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
        groups = nothing,
    )
    λ = eigvals(M)
    # points the window must contain: the eigenvalues and whatever else is drawn
    xs, ys = collect(real.(λ)), collect(imag.(λ))
    stability && (vlines!(ax, [0.0]; color = :gray60, linestyle = :dash); push!(xs, 0.0); push!(ys, 0.0))
    if groups !== nothing && (ellipse || outlier)
        # the predictions of block_spectrum for the partition `groups`
        p = block_spectrum(M, groups)
        if ellipse
            t = range(0, 2π; length = 200)
            ex, ey = p.center .+ p.horizontal .* cos.(t), p.vertical .* sin.(t)
            lines!(ax, ex, ey; color = ellipsecolor)
            append!(xs, ex)
            append!(ys, ey)
        end
        if outlier && !isempty(p.outliers)
            scatter!(
                ax, p.outliers, zeros(length(p.outliers)); markersize = 16, color = :transparent,
                strokecolor = outliercolor, strokewidth = 2,
            )
            append!(xs, p.outliers)
            append!(ys, zeros(length(p.outliers)))
        end
    elseif ellipse || outlier
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


function FoodWebs.blockplot(M::AbstractMatrix, groups::AbstractVector; size = (520, 450), title = "", kwargs...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title, aspect = DataAspect(), yreversed = true)
    hm = FoodWebs.blockplot!(ax, M, groups; kwargs...)
    Colorbar(fig[1, 2], hm)
    return fig, ax
end

function FoodWebs.blockplot!(
        ax, M::AbstractMatrix, groups::AbstractVector; colormap = nothing, linecolor = :black, linewidth = 1,
    )
    S = LinearAlgebra.checksquare(M)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S"))
    rank = Dict(g => k for (k, g) in enumerate(unique(groups)))
    order = sortperm(1:S; by = i -> (rank[groups[i]], i))
    A = float.(M[order, order])
    signed = any(<(0), A)
    m = maximum(abs, A; init = 0.0)
    colorrange = signed ? (-max(m, eps()), max(m, eps())) : (0.0, max(m, eps()))
    colormap = something(colormap, signed ? :RdBu : :Blues)
    # rows of M down the y axis and columns across x, as the matrix is written
    hm = heatmap!(ax, 1:S, 1:S, permutedims(A); colormap, colorrange)
    sizes = [count(==(g), groups) for g in unique(groups)]
    edges = cumsum(sizes)[1:(end - 1)] .+ 0.5
    isempty(edges) || (vlines!(ax, edges; color = linecolor, linewidth); hlines!(ax, edges; color = linecolor, linewidth))
    return hm
end

end
