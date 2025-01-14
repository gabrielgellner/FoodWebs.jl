struct PredationMatrix <: AbstractArray{Int, 2}
    links::Matrix{Int}
    S::Int
    C::Float64
end

function PredationMatrix(links::Matrix{Int})
    S1, S2 = size(links)
    if S1 != S2
        error("PredationMatrix must be square, given dimensions ($S1, S2)")
    end
    # if there is canibals the connectance measure might be off
    return PredationMatrix(links, S1, 2 * sum(links) / S1 ^ 2)
end

Base.size(pmat::PredationMatrix) = size(pmat.links)
Base.IndexStyle(::Type{<:PredationMatrix}) = IndexCartesian()
Base.getindex(pmat::PredationMatrix, I...) = getindex(pmat.links, I...)
Base.strides(pmat::PredationMatrix) = strides(pmat.links)

"""
    food_chain_network(S)

Generate a linear food chain of length `S`
"""
function food_chain_network(S::Int)
    @assert S > 0
    return PredationMatrix(diagm(1 => fill(1, S - 1)))
end

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
function cascade_network(S, C; maxiter = 1000)
    for iter in 1:maxiter
        adj = fill(0, S, S)
        for i in 1:(S - 1)
            for j in (i + 1):S
                # Since we are only looping over the top triangle we are only checking
                # S * (S - 1) / 2 possible links and we know that there are C * (S - 1) / 2 in
                # total so the probability is just C that we have a link at the specific location
                if rand() < C
                    adj[i, j] = 1
                end
            end
        end
        if is_connected(DiGraph(adj))
            return PredationMatrix(adj)
        end
    end
    error("Could not find a connected network for ($S, $C) in $maxiter iteractions")
end

"""
    generalized_cascade_network(S, C)::PredationMatrix

Generate a `S`×`S` adjancency matrix, with connectance `C`, following the network structure
described by:
Stouffer, D. B., Camacho, J., Guimerà, R. & Ng, C. Quantitative patterns in the structure of
model and empirical food webs. Ecology 86, 1301–1311 (2005).
"""
function generalized_cascade_network(S::Int, C::Float64; maxiter = 1000)
    for iter in 1:maxiter
        β = 1 / C - 1
        adj = fill(0, S, S)
        for i = 1:S
            for j = (i + 1):S
                # I need to get rid of the 2 from the orignal model to get the correct connectance
                # this is related to the fact that we double the number of entries when we move
                # from an adjancency matrix to a  # community matrix. Maybe I just need to fix the
                # connectance calculator to only look at the lower triangle of the community matrix
                if rand() < rand(Beta(1, β))
                    adj[i, j] = 1
                end
            end
        end
        if is_connected(DiGraph(adj))
            return PredationMatrix(adj)
        end
    end
    error("Could not find a connected network for ($S, $C) in $maxiter iteractions")
end

"""
    niche_network(S, C)::PredationMatrix

Generate a `S`×`S` adjancency matrix, with connectance `C`, following the network
structure described by:
Williams, R. J. & Martinez, N. D. Simple rules yield complex food webs. Nature 404, 180–3 (2000).
"""
function niche_network(S::Int, C::Float64; maxiter = 1000)
    adj = fill(0, S, S) #not sure why this is needed outside
    cond = false
    for i in 1:maxiter
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

        if is_connected(DiGraph(adj))
            return PredationMatrix(adj)
        end
    end
    error("Could not find a connected network for ($S, $C) in $maxiter iteractions")
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
