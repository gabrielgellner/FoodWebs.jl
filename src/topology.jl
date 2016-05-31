"""
    may_network(S::Int, C::Float64)

Returns an `S`x`S` Erdos-Renlyi adjacency matrix with connectance `C`.
"""
function may_network(S::Int, C::Float64)
    adj = zeros(Int, S, S)
    for i in 1:S, j in (i + 1):S
        if rand() < C
            adj[i, j] = 1
        end
    end
    return adj
end

"""
    cascade_network(S::Int, C::Float64)

Returns an `S`x`S` adjacency matrix with network structure described by Cohen, Briand & Newman, 1989.
"""
function cascade_network(S::Int, C::Float64)
    adj = zeros(Int, S, S)
    for i in 1:S, j in (i + 1):S
        # I need to get rid of the 2 from the orignal model to get the correct
        # connectance, need to think about his clearly this is related to the
        # fact that we double the number of entries when we move from an
        # adjancency matrix to a community matrix. Maybe I just need to fix my
        # connectance calculator to only look at the lower triangle of the
        # community matrix
        if rand() < S*C/(S - 1.0)
            adj[i, j] = 1
        end
    end
    adj
end

"""
    niche_network(S::Int, C::Float64)

Returns an `S`x`S` adjacency matrix with network structure described by Williams & Martinez 2000.
"""
function niche_network(S::Int, C::Float64)
    eta = sort(rand(S))
    beta = 1.0/C - 1.0;
    r = eta .* rand(Beta(1, beta), S)
    # set r[position of min[eta]] = 0, so that every web has at least one basal species
    r[1] = 0.0

    adj = zeros(Int, S, S)
    for i in 1:S, j in 1:S
        center = rand(Uniform(r[i]/2.0, eta[i]))
        #TODO I currently don't allow things on the diagonal, but really
        # this is not correct as this would simply be canabals
        if ((center - r[i]/2) <= eta[j] <= (center + r[i]/2)) && (i != j)
            adj[i, j] = 1
        end
    end
    return adj
end

#TODO: I need to think of a better name for this
"""
    random_predprey(adj::Array{Int, 2}, dist::Distribution, f::Float64)
"""
function random_predprey(adj::Array{Int, 2}, dist::Distribution, f::Float64)
    S = size(adj, 1) # I assume it is square
    cmat = diagm(-ones(S))
    for i in 1:S, j in (i + 1):S
        if adj[i, j] == 1
            cmat[i, j] = -abs(rand(dist))
            cmat[j, i] = f*abs(rand(dist))
        end
    end
    return cmat
end
