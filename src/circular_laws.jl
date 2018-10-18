import Random: shuffle
import Statistics: cov, mean

struct EllipsisLawResult
    S::Int
    μ::Float64
    μ2::Float64
    μ_diag::Float64
    σ2::Float64
    ρ::Float64
    x_semi_axis::Float64
    y_semi_axis::Float64
end

"""
    offdiag_pairs(mat::Matrix)

return an array of the pairs of the offdiagonal pairs of the square matrix `mat`

"""
function offdiag_pairs(mat::Matrix)
    S = size(mat, 1)
    @assert S == size(mat, 2) # must be square
    # Allocate a matrix of columns [a_ij, a_ji] excluding the diagonals
    # the size is simply (S^2 - S)/2 = S(S - 1)/2
    pairs = fill(zero(eltype(mat)), Int(S * (S - 1) / 2), 2)
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

"""
    shuffle_pairs(mat::Matrix)

Given a square matrix `mat` shuffle the offdiagonal pairs above and below the main diagonal
at random. This is often useful for looking at distributions of the off diagonal pairs that
that have equal marginal distributions.
"""
function shuffle_pairs(mat::Matrix)
    S = size(mat, 1)
    @assert S == size(mat, 2) # must be square
    shuff_mat = similar(mat)
    # leave the diagonal elements alone
    shuff_mat[diagind(mat)] = diag(mat)
    #TODO: it would be nice to have some kind of degree of shuffling in the future
    for i = 1:S
        for j = i:S
            if i != j
                shuff_pair = shuffle([mat[i, j], mat[j, i]])
                shuff_mat[i, j] = shuff_pair[1]
                shuff_mat[j, i] = shuff_pair[2]
            end
        end
    end
    return shuff_mat
end

function ellipsis_sample(mat::Matrix{Float64})
    S = size(mat, 1)
    @assert S == size(mat, 2)

    pairs = offdiag_pairs(mat)
    μx = mean(pairs[:, 1])
    μy = mean(pairs[:, 2])
    μ = mean([μx, μy]) # this is the same as mean(pairs[:]), that is, the mean of the pair distribution
    μxy = mean(pairs[:, 1] .* pairs[:, 2])
    μ_diag = mean(diag(mat))

    # Instead deal with the numerical covariance matrix of the off diagonal pairs
    # cov(mat, 1, false) uses the population formula for the variance, and goes over columns
    σ2x, ρ21, σ2y = cov(pairs, dims = 1, corrected = false)[[1, 2, 4]]
    # take the mean if the marginals are not equal? Is this a good rule?
    #σ2 = mean([σ2x, σ2y])
    # I for some reason feel that it should be the geometric mean, as I see the variance
    # being multiplicative for some reason, note this is also the denominator of the correllation
    σ2 = sqrt(σ2x * σ2y)
    ρ = ρ21 / sqrt(σ2x * σ2y)
    #ρ3 = (μxy - μx * μy) / sqrt(σ2x * σ2y) same as ρ from scaled cov function, like above

    #TODO: now that I am using the geometric mean of var, I am assumign that this makes sense
    # in the sqrt(S*V) part
    x_semi_axis = sqrt(S * σ2) * (1 + ρ)
    y_semi_axis = sqrt(S * σ2) * (1 - ρ)

    #return EllipsisLawResult(S, μx, μy, μ, μxy, μ_diag, σ2x, σ2y, σ2, ρ, ρ2, x_semi_axis, y_semi_axis)
    return Dict(
        :S => S,
        :μx => μx, :μy => μy, :μxy => μxy, :μ => μ, :μ_diag => μ_diag,
        :σ2x => σ2x, :σ2y => σ2y, :σ2 => σ2,
        :ρ => ρ
    )
end
