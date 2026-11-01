"""
    random_predprey([rng], web, dist, f; shuffle_prey = false)

A random community (Jacobian) matrix with the structure of the food web `web`, using the
standard convention that `J[i, j]` is the effect of species `j` on the growth of species `i`.

For each link where predator `j` eats prey `i` (`web[i, j] == 1`), the effect of the predator
on the prey is `J[i, j] = -|a|` and the effect of the prey on the predator is
`J[j, i] = f |b|`, with `a, b ~ dist` independent and `f` the conversion efficiency. The
diagonal is `-1`. Cannibalistic links are ignored.

When two species eat each other, each entry of the pair sums the two links' effects: species
`i` loses `|a|` to `j` and gains `f |b'|` from eating `j`, so `J[i, j] = f |b'| - |a|`, of
either sign, with all four magnitudes drawn independently.

With `shuffle_prey = true` the two entries of each pair are randomly swapped, so the
predator–prey sign structure is kept but not which species is the predator.
"""
function random_predprey(
        rng::AbstractRNG, web::PredationMatrix, dist::Distribution, f::Real;
        shuffle_prey = false
    )
    S = size(web, 1)
    J = -Matrix{Float64}(I, S, S)
    # links can be on either side of the diagonal (e.g. niche webs where a prey has a higher
    # niche value than its predator)
    for i in 1:S, j in 1:S
        if i != j && web[i, j]
            neg = -abs(rand(rng, dist))
            pos = f * abs(rand(rng, dist))
            if shuffle_prey && rand(rng, Bool)
                neg, pos = pos, neg
            end
            # accumulate, so that a pair of species that eat each other gets both links' effects
            J[i, j] += neg
            J[j, i] += pos
        end
    end
    return J
end
function random_predprey(web::PredationMatrix, dist::Distribution, f::Real; kwargs...)
    return random_predprey(Random.default_rng(), web, dist, f; kwargs...)
end

"""
    expected_community_matrix(web, dist, f; shuffle_prey = false)

The mean of [`random_predprey`](@ref)`(web, dist, f; shuffle_prey)` given the web: `-E|a|` on
each link's entry for the predator's effect on its prey, `f E|a|` on the prey's effect on its
predator, both summed for species that eat each other, and `-1` on the diagonal, with `E|a|`
the mean magnitude of `dist`, computed by the midpoint rule over 20 000 quantiles (relative error
about `1e-5` for normal and exponential magnitudes, below `1e-3` for heavy tails such as a t
distribution with 3 degrees of freedom). With `shuffle_prey = true` each entry of a linked pair is the mean
of the two, `(f - 1) E|a| / 2`.

For structured webs such as cascade and niche webs, where the elliptic law's exchangeable pairs
do not hold, the eigenvalues of this matrix give the largest imaginary part of `random_predprey`'s
spectrum to about 1%, and the spread of imaginary parts to within about 15% (the fluctuations
widen it). The real parts, and so stability, come from the fluctuations around the mean: the
expected matrix is stable, and with `f = 1` it is the identity's negative plus a skew-symmetric
matrix, with every eigenvalue at real part `-1`.
"""
function expected_community_matrix(web::PredationMatrix, dist::Distribution, f::Real; shuffle_prey = false)
    S = size(web, 1)
    m = mean_magnitude(dist)
    neg, pos = shuffle_prey ? ((f - 1) * m / 2, (f - 1) * m / 2) : (-m, f * m)
    J = -Matrix{Float64}(I, S, S)
    for i in 1:S, j in 1:S
        if i != j && web[i, j]
            J[i, j] += neg
            J[j, i] += pos
        end
    end
    return J
end

# E|X| for a univariate distribution, as the mean of |quantile| over a midpoint grid of
# probabilities
mean_magnitude(dist::Distribution; n::Int = 20_000) = sum(abs(quantile(dist, (k - 0.5) / n)) for k in 1:n) / n
