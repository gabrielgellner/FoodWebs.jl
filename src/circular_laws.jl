
immutable EllipsisLawResult
    S::Int
    μ::Float64
    μ2::Float64
    μ_diag::Float64
    σ2::Float64
    ρ::Float64
    x_semi_axis::Float64
    y_semi_axis::Float64
end

function ellipsis_sample(mat::Matrix{Float64})
    @assert size(mat, 1) == size(mat, 2)
    S = size(mat, 1)
    # all calculations need to ignore the diagonal
    imat = copy(mat) #TODO: is this the best way to do this?
    imat[diagind(imat)] = 0.0
    #μ = sum(mat)/(S*(S - 1))
    # This requires that the upper and lower triangle are shuffled as far as being - or +
    #μ2 = 2*sum(triu(mat .* mat'))/(S*(S - 1))

    # try to calc everyting in a single pass of the matrix
    μ = 0.0
    μ2 = 0.0
    μ_diag = 0.0
    σ2 = 0.0
    for i = 1:S
        for j = i:S
            if i != j
                μ += mat[i, j] + mat[j, i]
                μ2 += mat[i, j]*mat[j, i]
                σ2 += mat[i, j]^2 + mat[j, i]^2
            else
                μ_diag += mat[i, i]
            end
        end
    end
    μ = μ/(S*(S - 1))
    μ2 = 2*μ2/(S*(S - 1))
    μ_diag = μ_diag/S
    σ2 = σ2/(S*(S - 1)) - μ^2

    #μ_diag = mean(diag(inmat))
    #σ2 = sum(mat.^2)/(S*(S - 1)) - μ^2
    ρ = (μ2 - μ^2)/σ2
    x_semi_axis = sqrt(S*σ2)*(1 + ρ)
    y_semi_axis = sqrt(S*σ2)*(1 - ρ)
    return EllipsisLawResult(S, μ, μ2, μ_diag, σ2, ρ, x_semi_axis, y_semi_axis)
end

#=
# Plotting Routines, depends on OdinSon/PyPlot
=#
using PyCall
using PyPlot
@pyimport matplotlib.patches as patches

function plot_eigs(mat)
    @assert size(mat, 1) == size(mat, 2)
    S = size(mat, 1)
    eigs = eigvals(mat)
    reigs = real(eigs)
    ieigs = imag(eigs)
    scatter(reigs, ieigs, s=5.0, color="black")
    #xlim(1.1*minimum(reigs), 1.1*maximum(reigs))
end

function eig_ellipse(mat)
    ellps = ellipsis_sample(mat)
    Θs = linspace(π/2.0, 0.0, 500)
    xbase = ellps.x_semi_axis*cos(Θs)
    ybase = ellps.y_semi_axis*sin(Θs)
    #TODO: add a commment as to why I am doing this
    x = vcat(xbase, reverse(xbase), -xbase, reverse(-xbase))
    y = vcat(ybase, reverse(-ybase), -ybase, reverse(ybase))
    pelip = patches.Polygon(hcat(x + ellps.μ + ellps.μ_diag, y), fill=false, edgecolor="#0072B2", linewidth=2)
    #TODO: should really pass this in
    ax = gca()
    ax[:add_patch](pelip)
end

function plot_tang(S, C, f, r)
    # calc the elipse points
    radiusx = hor_semi_an(S, C, f, r)
    radiusy = sqrt(S*V_an(C, f, r))*(1 - rho_an(C, f, r))
    centerx = -E_an(C, f, r)
    thetas = linspace(π/2.0, 0.0, 500)
    xbase = radiusx*cos(thetas)
    ybase = radiusy*sin(thetas)
    x = vcat(xbase, reverse(xbase), -xbase, reverse(-xbase))
    y = vcat(ybase, reverse(-ybase), -ybase, reverse(ybase))
    pelip = patches.Polygon(hcat(x + centerx, y), fill=false, edgecolor="#0072B2", linewidth=2)

    ax = gca()
    ax[:add_patch](pelip)
end

#=
function eig_ellipse(S, C, f, r; nsamp=50)
    max_eig = -100.0
    for i in 1:nsamp
        M = build_matrix(S, C, f, r)
        deig = real_dom_eig(M)
        if deig > max_eig
            max_eig = deig
        end
        plot_eigs(M)
    end
    @show tang_an(S, C, f, r)
    @show max_eig
    @show tau_naive(S, C, f, r)
    plot_tang(S, C, f, r)
    xlabel(L"Re(\lambda)")
    ylabel(L"Im(\lambda)")
end
=#
