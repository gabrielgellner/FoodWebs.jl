#TODO: I need to think of a better name for this.
function random_predprey(adj::Matrix{Int}, dist::Distribution, f::Float64, shuffle=false)
    S = size(adj, 1)
    @assert S == size(adj, 2) # only allow square matrices
    cmat = diagm(-ones(S)) # todo allow the diagonals to be set from an passed in argument
    for i in 1:S, j in (i + 1):S
        if adj[i, j] == 1
            neg = -abs(rand(dist))
            pos = f*abs(rand(dist))
            if shuffle
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
