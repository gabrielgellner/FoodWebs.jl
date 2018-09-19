struct PredationMatrix <: AbstractArray{Int, 2}
    links::Matrix{Int}
    S::Int
end

function PredationMatrix(links::Matrix{Int})
    S1, S2 = size(links)
    if S1 != S2
        error("PredationMatrix must be square, given dimensions ($S1, S2)")
    end
    return PredationMatrix(links, S1)
end

Base.size(pmat::PredationMatrix) = size(pmat.links)
Base.IndexStyle(::Type{<:PredationMatrix}) = IndexCartesian()
Base.getindex(pmat::PredationMatrix, I...) = getindex(pmat.links, I...)

"""
    may_network(S, C)

Generate a `S`×`S` adjancency matrix, with connectance `C`. Following the random network
described by:
May, R. M. Will a large complex system be stable? Nature 238, 413–414 (1972).
"""
function may_network(S::Int, C::Float64)
    adj = fill(0, S, S)
    for i in 1:S
        for j in 1:S
            if rand() < C / 2
                adj[i, j] = 1
            end
        end
    end
    return PredationMatrix(adj)
end

"""
cascade_network(S, C)::PredationMatrix

Generate a `S`×`S` adjancency matrix, with connectance `C`, following the network structure
described by:
Cohen, J. E. 1978. Food Webs and Niche Space. Princeton University Press, Princeton, N.J.
"""
function cascade_network(S, C)
    adj = fill(0, S, S)
    #NOTE: we only use C * S / (S - 1) (instead of 2CS/(S - 1) ) since we return a
    #      PredationMatrix and we want our connectance to be over both the predator and
    #      prey link, that is implicit in this formulation
    p = C * S / (S - 1)
    for i in 1:(S - 1)
        for j in (i + 1):S
            if rand() < p
                adj[i, j] = 1
            end
        end
    end
    return PredationMatrix(adj)
end

"""
    generalized_cascade_network(S, C)::PredationMatrix

Generate a `S`×`S` adjancency matrix, with connectance `C`, following the network structure
described by:
Stouffer, D. B., Camacho, J., Guimerà, R. & Ng, C. Quantitative patterns in the structure of
model and empirical food webs. Ecology 86, 1301–1311 (2005).
"""
function generalized_cascade_network(S::Int, C::Float64)
    β = 1 / C - 1
    adj = fill(0, S, S)
    for i = 1:S
        for j = i:S
            # I need to get rid of the 2 from the orignal model to get the correct connectance
            # this is related to the fact that we double the number of entries when we move
            # from an adjancency matrix to a  # community matrix. Maybe I just need to fix the
            # connectance calculator to only look at the lower triangle of the community matrix
            if rand() < rand(Beta(1, β))
                adj[i, j] = 1
            end
        end
    end
    return PredationMatrix(adj)
end

"""
    niche_network(S, C)::PredationMatrix

Generate a `S`×`S` adjancency matrix, with connectance `C`, following the network
structure described by:
Williams, R. J. & Martinez, N. D. Simple rules yield complex food webs. Nature 404, 180–3 (2000).
"""
function niche_network(S, C)
    adj = fill(0, S, S) #not sure why this is needed outside
    cond = false
    while !cond
        n_i = sort(rand(S))
        r_i = rand(Beta(1, ((1 / (2 * C)) - 1)), S) .* n_i
        c_i = [rand(Uniform(r_i[i] / 2, n_i[i])) for i in 1:S]

        adj = fill(0, S, S)

        for i in 2:S
            for j in 1:S
                if n_i[j] > (c_i[i] - (0.5 * r_i[i])) && n_i[j] < (c_i[i] + 0.5 * r_i[i])
                    adj[j, i] = 1
                end
            end
        end

        cond = is_connected(DiGraph(adj))
    end

    return PredationMatrix(adj)
end

function niche_network(S::Int, C::Float64)
    adj = fill(0, S, S) #not sure why this is needed outside
    cond = false
    while !cond
        η = sort(rand(S))
        β = 1 / C - 1
        r = η .* rand(Beta(1, β), S)
        # set r[position of min[η]] = 0, so that every web has at least one basal species
        r[1] = 0.0

        adj = fill(0, S, S)
        for i in 1:S
            for j in 1:S
                center = rand(Uniform(r[i]/2, η[i]))
                if ((center - r[i] / 2) <= η[j] <= (center + r[i] / 2))
                    adj[j, i] = 1
                end
            end
        end

        cond = is_connected(DiGraph(adj))
    end
    return PredationMatrix(adj)
 end

"""
    generalized_niche_network(S, C, diet)

Generate a `S`×`S` adjancency matrix with connectance `C`, and diet continuity `diet`
following the network structure described by:
Guimerà, R. et al. Origin of compartmentalization in food webs. Ecology 91, 2941–51 (2010).
"""
function generalized_niche_network(S::Int, C::Float64, diet::Float64)
    η = sort(rand(S))
    β = 1/C - 1
    r = η .* rand(Beta(1, β), S)
    # set r[position of min[η]] = 0, so that every web has at least one basal species
    r[1] = 0.0
    reduced_r = diet * r
    centers = [rand(Uniform(reduced_r[i] / 2, η[i])) for i = 1:S]
    adj = fill(0, S, S)
    for i = 1:S
        for j = 1:S
            # add the reduced Niche style links
            if i != j && (centers[i] - reduced_r[i] / 2) <= η[j] && η[j] <= (centers[i] + reduced_r[i] / 2)
                adj[i, j] = 1
            end
        end
        # now we pick the cascade style links to make up for the missing links from the
        # reduced nich ranges
        n_missing_prey = round(Int, (r - reduced_r) * S)
        if n_missing_prey[i] > 0
            pos = findall(x -> x == 0, adj[i, 1:i]) #TODO: this is likely not efficient
            # If missing prey is larger than available: eat them all ... need to check if
            # this is a bug in my code or normal behavior
            for j in sample(pos, min(length(pos), n_missing_prey[i]))
                adj[i, j] = 1
            end
        end
    end
    #TODO: check for being connected
    return PredationMatrix(adj)
end
