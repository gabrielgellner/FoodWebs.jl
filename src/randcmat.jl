#TODO: I need to think of a better name for this.
function random_predprey(pmat::PredationMatrix, dist::Distribution, f::Float64; shuffle_prey = false)
    #TODO: allow the diagonals to be set from an passed in argument
    cmat = -Matrix{Float64}(I, pmat.S, pmat.S)
    for i in 1:pmat.S, j in (i + 1):pmat.S
        if pmat[i, j] == 1
            neg = -abs(rand(dist))
            pos = f * abs(rand(dist))
            if shuffle_prey
                # Don't assume a given species is a predator or a prey, randomly choose for
                # each pair
                cmat[i, j], cmat[j, i] = shuffle([neg, pos])
            else
                cmat[j, i] = neg
                cmat[i, j] = pos
            end
        end
    end
    return cmat
end
