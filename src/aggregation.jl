# # Trophic species and aggregation
# Species with identical prey and identical predators form a trophic species (Briand & Cohen
# 1984; Williams & Martinez 2000). Lumping species into groups, trophic species or any other
# partition, gives a web on the groups. Aggregation changes connectance, motif counts and
# stability, so it is a separate, explicit step.

"""
    trophic_species(web)

The trophic species of `web`: a group label for each species, where species with identical prey
and identical predators share a label. Labels are `1, 2, …` in order of each group's first
species. A cannibalistic link counts among a species' prey and predators like any other, so two
cannibals with otherwise identical links, each eating only itself, are different trophic
species. This follows `TrophicSpecies` of cheddar (Hudson, Reuman & Emerson 2013), and isolated
species with no links form one trophic species. Lump the groups with [`lump_species`](@ref).
"""
function trophic_species(web::AbstractMatrix)
    S = LinearAlgebra.checksquare(web)
    A = web .!= 0
    labels = zeros(Int, S)
    signatures = Dict{Tuple{BitVector, BitVector}, Int}()
    for j in 1:S
        key = (BitVector(view(A, :, j)), BitVector(view(A, j, :)))
        labels[j] = get!(signatures, key, length(signatures) + 1)
    end
    return labels
end

"""
    lump_species(web, groups)
    lump_species(fw::FoodWeb, groups; names, combine = Dict())

The web on the groups of `groups`, one label per species (such as [`trophic_species`](@ref) or
taxonomic or functional groups): group `b` eats group `a` when some species of `b` eats some
species of `a`, and a group eats itself when some species in it eats another or itself. Groups
are ordered by their first species. A web returns a [`PredationMatrix`](@ref).

A `FoodWeb` returns a `FoodWeb` with one species per group, named by `names` (by default the
member names joined by `", "`), the web's name and metadata, and its traits combined over each
group's members:

- `"N"` (abundance) is summed;
- `"M"` (body mass) is the mean weighted by `"N"` when the web has abundances, so that the
  group's biomass `M N` is the sum of its members', and the plain mean otherwise;
- other numeric traits take the mean;
- other traits join their distinct non-empty values with `", "`.

A missing value among a group's members makes its numeric trait `missing`. `combine` maps a
trait name to a function of the members' values, overriding these rules, for example
`Dict("M" => maximum)`. cheddar's `LumpNodes` takes the same weighted mean of `"M"` but the
plain mean of `"N"`.
"""
function lump_species(web::AbstractMatrix, groups::AbstractVector)
    S = LinearAlgebra.checksquare(web)
    length(groups) == S || throw(DimensionMismatch("groups must have length $S"))
    index, members = group_members(groups)
    A = zeros(Bool, length(members), length(members))
    for i in 1:S, j in 1:S
        web[i, j] != 0 && (A[index[groups[i]], index[groups[j]]] = true)
    end
    return PredationMatrix(A)
end
function lump_species(
        fw::FoodWeb, groups::AbstractVector; names = nothing, combine::AbstractDict = Dict{String, Function}()
    )
    web = lump_species(fw.web, groups)
    _, members = group_members(groups)
    names = something(names, [join(fw.species[m], ", ") for m in members])
    length(names) == length(members) ||
        throw(DimensionMismatch("$(length(names)) names for $(length(members)) groups"))
    traits = Dict{String, Vector}(
        trait => [combine_trait(fw, trait, values, m, combine) for m in members] for (trait, values) in fw.traits
    )
    return FoodWeb(web, names, fw.name, fw.metadata, traits)
end

# The index of each group label, in order of first appearance, and the members of each group
function group_members(groups)
    index = Dict{eltype(groups), Int}()
    members = Vector{Int}[]
    for (k, g) in enumerate(groups)
        haskey(index, g) || (index[g] = length(members) + 1; push!(members, Int[]))
        push!(members[index[g]], k)
    end
    return index, members
end

function combine_trait(fw::FoodWeb, trait, values, members, combine)
    v = values[members]
    haskey(combine, trait) && return combine[trait](v)
    if all(x -> x isa Union{Missing, Real}, v)
        trait == "N" && return sum(v)
        if trait == "M" && haskey(fw.traits, "N")
            N = fw.traits["N"][members]
            return sum(v .* N) / sum(N)
        end
        return mean(v)
    end
    return join(unique(string(x) for x in v if !ismissing(x) && !isempty(string(x))), ", ")
end
