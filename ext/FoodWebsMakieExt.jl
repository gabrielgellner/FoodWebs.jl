module FoodWebsMakieExt

using FoodWebs: FoodWebs, trophic_layout, trophic_levels
using Makie: Axis, Figure, arrows2d!, hidespines!, hidexdecorations!, limits!, scatter!, text!

function FoodWebs.foodwebplot(web::AbstractMatrix; size = (500, 500), title = "", kwargs...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title, ylabel = "trophic level")
    FoodWebs.foodwebplot!(ax, web; kwargs...)
    return fig, ax
end

function FoodWebs.foodwebplot!(
        ax, web::AbstractMatrix; levels = trophic_levels(web), layout = trophic_layout(web; levels),
        color = levels, colormap = :viridis, markersize = 22, edgecolor = (:gray30, 0.6),
        edgewidth = 1.2, labels = string.(axes(web, 1)), highlight_cannibals = true,
    )
    S = size(web, 1)
    xs, ys = first.(layout), last.(layout)
    # stop each arrow short of the nodes, in data units scaled to the layout's spacing
    gap = 0.18
    for i in 1:S, j in 1:S
        (i != j && web[i, j] != 0) || continue
        a, b = layout[i], layout[j]
        d = b .- a
        n = sqrt(sum(abs2, d))
        n > 2gap || continue
        arrows2d!(
            ax, [a[1] + gap * d[1] / n], [a[2] + gap * d[2] / n], [d[1] * (1 - 2gap / n)], [d[2] * (1 - 2gap / n)];
            color = edgecolor, shaftwidth = edgewidth, tipwidth = 4 + 3edgewidth, tiplength = 4 + 3edgewidth,
        )
    end
    cannibal = [web[i, i] != 0 for i in 1:S]
    strokewidth = highlight_cannibals ? [c ? 3.0 : 1.0 for c in cannibal] : 1.0
    scatter!(ax, xs, ys; color, colormap, markersize, strokewidth, strokecolor = :black)
    # a light glow keeps the labels readable on dark and light nodes alike
    labels === nothing || text!(
        ax, xs, ys; text = labels, align = (:center, :center), fontsize = 10, glowwidth = 2,
        glowcolor = (:white, 0.8),
    )
    # pad the limits so that nodes at the edges are not clipped
    limits!(ax, minimum(xs) - 0.7, maximum(xs) + 0.7, minimum(ys) - 0.4, maximum(ys) + 0.4)
    # only the trophic level axis carries information
    hidespines!(ax, :t, :r, :b)
    hidexdecorations!(ax)
    return ax
end

end
