"""
    random_predprey([rng], web, dist, f; shuffle_prey = false)

A random community (Jacobian) matrix with the structure of the food web `web`, using the
standard convention that `J[i, j]` is the effect of species `j` on the growth of species `i`.

For each link where predator `j` eats prey `i` (`web[i, j] == 1`), the effect of the predator
on the prey is `J[i, j] = -|a|` and the effect of the prey on the predator is
`J[j, i] = f |b|`, with `a, b ~ dist` independent and `f` the conversion efficiency. The
diagonal is `-1`. Cannibalistic links are ignored.

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
            J[i, j] = neg
            J[j, i] = pos
        end
    end
    return J
end
function random_predprey(web::PredationMatrix, dist::Distribution, f::Real; kwargs...)
    return random_predprey(Random.default_rng(), web, dist, f; kwargs...)
end
