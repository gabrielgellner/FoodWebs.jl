"""
    connectance(web)

The connectance `L / S^2` of the `S` species food web `web`, where `L` is the number of links.
"""
function connectance(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    return count(!iszero, web) / S^2
end

"""
    basal_species(web)

The species that eat nothing, where `web[i, j] != 0` means species `j` eats species `i`.
"""
basal_species(web::AbstractMatrix) = [j for j in axes(web, 2) if all(iszero, view(web, :, j))]

"""
    top_species(web)

The species that nothing eats, where `web[i, j] != 0` means species `j` eats species `i`.
"""
top_species(web::AbstractMatrix) = [i for i in axes(web, 1) if all(iszero, view(web, i, :))]

"""
    disconnected_from_basal(web)

The species in `web` with no path of feeding links up from a basal species, where
`web[i, j] != 0` means species `j` eats species `i`. These are species in closed loops that
only eat each other, so no energy can flow to them. Cannibalistic links are ignored, so a
species that only eats itself counts as basal.
"""
function disconnected_from_basal(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    eats(pred, prey) = pred != prey && web[prey, pred] != 0
    basal = [j for j in 1:S if !any(eats(j, i) for i in 1:S)]
    # search up the web from the basal species along prey -> predator links
    reached = falses(S)
    reached[basal] .= true
    frontier = copy(basal)
    while !isempty(frontier)
        prey = pop!(frontier)
        for pred in 1:S
            if eats(pred, prey) && !reached[pred]
                reached[pred] = true
                push!(frontier, pred)
            end
        end
    end
    return findall(!, reached)
end

"""
    trophic_levels(web)

The prey-averaged trophic levels of the species in `web`, where `web[i, j] != 0` means
species `j` eats species `i`. Basal species have trophic level 1 and every other species has
1 plus the mean trophic level of its prey (Levine 1980), that is, the solution of
`(I - P) t = 1` where `P` is the diet matrix with predators as rows normalized to sum to 1.

Cannibalistic links are ignored. Throws an `ArgumentError` if some species have no path of
feeding links up from a basal species (see [`disconnected_from_basal`](@ref)), as their
trophic levels are then undefined.
"""
function trophic_levels(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    unfed = disconnected_from_basal(web)
    if !isempty(unfed)
        throw(ArgumentError("species $unfed have no feeding path up from a basal species, so " *
                            "their trophic levels are undefined"))
    end
    # diet matrix with predators as rows and without cannibalism
    P = Matrix{Float64}(transpose(web) .!= 0)
    P[diagind(P)] .= 0
    P ./= max.(sum(P, dims = 2), 1)
    return (I - P) \ ones(S)
end
