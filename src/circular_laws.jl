
immutable EllipsisLawResult
    S::Int
    μ::Float64
    μ2::Float64
    μ_diag::Float64
    σ2::Float64
    ρ::Float64
    ρ2::Float64 # This is for the analytic formula (μ2 - μ^2)/σ2
    x_semi_axis::Float64
    y_semi_axis::Float64
end

"""
    offdiag_pairs(mat::Matrix)

return an array of the pairs of the offdiagonal pairs of the matrix `mat`

"""
function offdiag_pairs(mat::Matrix)
    # Allocate a matrix of columns [a_ij, a_ji] excluding the diagonals
    # the size is simply (S^2 - S)/2 = S(S - 1)/2
    S = size(mat, 1)
    @assert S == size(mat, 2) # must be square
    pairs = zeros(eltype(mat), Int(S*(S - 1)/2), 2)
    #TODO: likely there is a more efficient way of doing this? Can I skip over the diagonals
    # without an if statement?
    # use a simple counter to unwrap the matrix into the pairs matrix
    current_row = 1
    for i = 1:S
        for j = i:S
            if i != j
                pairs[current_row, :] = [mat[i, j], mat[j, i]]
                current_row += 1
            end
        end
    end
    return pairs
end

function ellipsis_sample(mat::Matrix{Float64})
    @assert size(mat, 1) == size(mat, 2)
    S = size(mat, 1)
    # all calculations need to ignore the diagonal
    #imat = copy(mat) #TODO: is this the best way to do this?
    #imat[diagind(imat)] = 0.0
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

    # this formula assumes that the marginal distributions of off diagonal pairs are equal
    # (which is why there is a single σ2)
    #ρ0 = (μ2 - μ^2)/σ2

    pairs = offdiag_pairs(mat)
    # Instead deal with the numerical covariance matrix of the off diagonal pairs
    σ11, ρ21, σ22 = cov(pairs, corrected=false)[[1, 2, 4]]
    #σ2 = σ2/(S*(S - 1)) - μ^2
    #σ2 = mean([σ1, σ2]) # take the mean if the marginals are not equal? Is this a good rule?
    σ2 = mean([σ11, σ22])
    ρ = ρ21/sqrt(σ11*σ22)
    ρ2 = (μ2 - μ^2)/σ2
    # So instead we use the actual pair distribution to calculate this
    #ρ2 = mean((pairs[:, 1] - mean(pairs[:, 1])).*(pairs[:, 2] - mean(pairs[:, 2])))/sqrt(σ11*σ22)
    #ρ3 = cor(pairs)[1, 2] #TODO: I am not sure if this is using a sample formula or population
    #@show ρ, ρ2, ρ3

    x_semi_axis = sqrt(S*σ2)*(1 + ρ)
    y_semi_axis = sqrt(S*σ2)*(1 - ρ)
    return EllipsisLawResult(S, μ, μ2, μ_diag, σ2, ρ, ρ2, x_semi_axis, y_semi_axis)
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

function eig_ellipse1(mat)
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

function eig_ellipse2(mat)
    ellps = ellipsis_sample(mat)
    Θs = linspace(π/2.0, 0.0, 500)
    x_semi_axis = sqrt(ellps.S*ellps.σ2)*(1 + ellps.ρ2)
    y_semi_axis = sqrt(ellps.S*ellps.σ2)*(1 - ellps.ρ2)
    xbase = x_semi_axis*cos(Θs)
    ybase = y_semi_axis*sin(Θs)
    #TODO: add a commment as to why I am doing this
    x = vcat(xbase, reverse(xbase), -xbase, reverse(-xbase))
    y = vcat(ybase, reverse(-ybase), -ybase, reverse(ybase))
    pelip = patches.Polygon(hcat(x + ellps.μ + ellps.μ_diag, y), fill=false, edgecolor="#A372B2", linewidth=2)
    #TODO: should really pass this in
    ax = gca()
    ax[:add_patch](pelip)
end

function plot_eigdist(mat)
    plot_eigs(mat)
    eig_ellipse1(mat)
    eig_ellipse2(mat)
end
