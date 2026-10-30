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
        # each predator j draws a single probability x_j of eating each lower ranked species.
        # Drawing per pair instead would marginalize to Bernoulli(C), i.e. the plain cascade.
        x = rand(Beta(1, β), S)
        adj = fill(0, S, S)
        for i = 1:S
            for j = (i + 1):S
                if rand() < x[j]
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

        # each species gets a single feeding center, so its diet is an interval of niche space
        centers = [rand(Uniform(r[i] / 2, η[i])) for i in 1:S]
        adj = fill(0, S, S)
        for i in 1:S
            for j in 1:S
                if ((centers[i] - r[i] / 2) <= η[j] <= (centers[i] + r[i] / 2))
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
function generalized_niche_network(S::Int, C::Float64, diet::Float64; maxiter = 1000)
    β = 1 / C - 1
    for iter in 1:maxiter
        η = sort(rand(S))
        r = η .* rand(Beta(1, β), S)
        # set r[position of min[η]] = 0, so that every web has at least one basal species
        r[1] = 0.0
        reduced_r = diet * r
        centers = [rand(Uniform(reduced_r[i] / 2, η[i])) for i = 1:S]
        n_missing_prey = round.(Int, (r - reduced_r) * S)
        adj = fill(0, S, S)
        for i = 1:S
            # add the reduced niche style links, prey j of predator i
            for j = 1:S
                if i != j && (centers[i] - reduced_r[i] / 2) <= η[j] <= (centers[i] + reduced_r[i] / 2)
                    adj[j, i] = 1
                end
            end
            # now we pick cascade style links among lower ranked species to make up for the
            # missing links from the reduced niche ranges
            pos = [j for j in 1:(i - 1) if adj[j, i] == 0]
            k = min(length(pos), n_missing_prey[i])
            if k > 0
                for j in sample(pos, k; replace = false)
                    adj[j, i] = 1
                end
            end
        end
        if is_connected(DiGraph(adj))
            return PredationMatrix(adj)
        end
    end
    error("Could not find a connected network for ($S, $C, $diet) in $maxiter iteractions")
end
