# # Empirical food webs
# A web on disk is a directory with `nodes.csv`, whose first column names the species, and
# `trophic.links.csv` with columns `resource,consumer`. A collection is a directory of webs with a
# `catalog.toml` giving each web's title, citation, source and license, and a
# `datapackage.json` (Frictionless Data Package) declaring the files' columns. Bundled collections live
# in the package's `data/` directory; downloaded ones live in `foodweb_cache()`.

"""
    FoodWeb(web; species = string.(axes(web, 1)), name = "", metadata = Dict())

A food web with named species: an `S`×`S` binary matrix where `web[i, j] == 1` means species
`j` eats species `i`, together with the species names, a name for the web and a metadata
dictionary (for webs read from a catalog, its title, citation, source and license). A
`FoodWeb` is an `AbstractMatrix{Bool}`, so every function taking a web accepts it.
"""
struct FoodWeb <: AbstractMatrix{Bool}
    web::PredationMatrix
    species::Vector{String}
    name::String
    metadata::Dict{String, Any}
    function FoodWeb(web::PredationMatrix, species::AbstractVector, name::AbstractString, metadata::AbstractDict)
        length(species) == size(web, 1) || throw(
            DimensionMismatch("$(length(species)) species names for a web of $(size(web, 1)) species")
        )
        allunique(species) || throw(ArgumentError("species names must be unique"))
        return new(web, string.(species), name, Dict{String, Any}(metadata))
    end
end
function FoodWeb(
        web::AbstractMatrix; species = string.(axes(web, 1)), name::AbstractString = "",
        metadata::AbstractDict = Dict{String, Any}()
    )
    return FoodWeb(PredationMatrix(web), species, name, metadata)
end

Base.size(fw::FoodWeb) = size(fw.web)
Base.getindex(fw::FoodWeb, i::Int, j::Int) = fw.web[i, j]

function Base.show(io::IO, ::MIME"text/plain", fw::FoodWeb)
    title = isempty(fw.name) ? "FoodWeb" : "FoodWeb \"$(fw.name)\""
    return print(io, title, " with ", size(fw, 1), " species and ", count(fw.web.links), " links")
end

# Split one CSV line into fields, following RFC 4180 quoting: a quoted field may hold commas,
# and a doubled quote inside it is a literal quote.
function split_csv_line(line::AbstractString)
    fields = String[]
    field = IOBuffer()
    quoted = false
    chars = collect(line)
    k = 1
    while k <= length(chars)
        c = chars[k]
        if quoted
            if c == '"' && k < length(chars) && chars[k + 1] == '"'
                write(field, '"')
                k += 1
            elseif c == '"'
                quoted = false
            else
                write(field, c)
            end
        elseif c == '"'
            quoted = true
        elseif c == ','
            push!(fields, String(take!(field)))
        else
            write(field, c)
        end
        k += 1
    end
    push!(fields, String(take!(field)))
    return fields
end

# The header and rows of a CSV file without quoted line breaks
function read_csv(path::AbstractString)
    lines = filter(!isempty, readlines(path))
    isempty(lines) && throw(ArgumentError("$path is empty"))
    return split_csv_line(lines[1]), [split_csv_line(l) for l in lines[2:end]]
end

"""
    read_foodweb(dir; name = basename(dir), metadata = Dict())

Read the food web stored in the directory `dir`, which holds `nodes.csv`, whose first column
names the species, and `trophic.links.csv` with columns `resource,consumer`, one row per feeding link.
Returns a [`FoodWeb`](@ref) with the species in the order of `nodes.csv`.
"""
function read_foodweb(
        dir::AbstractString; name::AbstractString = basename(normpath(dir)),
        metadata::AbstractDict = Dict{String, Any}()
    )
    _, nodes = read_csv(joinpath(dir, "nodes.csv"))
    species = [first(row) for row in nodes]
    header, links = read_csv(joinpath(dir, "trophic.links.csv"))
    r = findfirst(==("resource"), header)
    c = findfirst(==("consumer"), header)
    (r === nothing || c === nothing) &&
        throw(ArgumentError("$(joinpath(dir, "trophic.links.csv")) needs resource and consumer columns"))
    index = Dict(s => k for (k, s) in enumerate(species))
    web = zeros(Bool, length(species), length(species))
    for row in links
        for s in (row[r], row[c])
            haskey(index, s) || throw(ArgumentError("species \"$s\" in trophic.links.csv is missing from nodes.csv"))
        end
        web[index[row[r]], index[row[c]]] = true
    end
    return FoodWeb(web; species, name, metadata)
end

"""
    foodweb_cache()

The directory holding downloaded collections of food webs: the environment variable
`FOODWEBS_CACHE` if set, otherwise a scratch space in the first Julia depot. The scripts in the
repository's `data/scripts` directory download webs here.
"""
function foodweb_cache()
    uuid = "6d07a039-1ac1-577b-a10e-abe3f3d3fc7d"
    return get(ENV, "FOODWEBS_CACHE", joinpath(first(DEPOT_PATH), "scratchspaces", uuid, "webs"))
end

# The directories holding collections, bundled ones first
collection_roots() = filter(isdir, [joinpath(dirname(@__DIR__), "data"), foodweb_cache()])

# id => (directory, catalog entry) for every web in every collection
function foodweb_catalog()
    catalog = Dict{String, Tuple{String, Dict{String, Any}}}()
    for root in collection_roots(), collection in readdir(root)
        file = joinpath(root, collection, "catalog.toml")
        isfile(file) || continue
        for (id, entry) in TOML.parsefile(file)
            catalog["$collection/$id"] = (joinpath(root, collection, id), entry)
        end
    end
    return catalog
end

"""
    foodweb_ids()

The identifiers `"collection/web"` of the empirical food webs available to [`foodweb`](@ref):
the collections bundled with the package and those downloaded to [`foodweb_cache`](@ref).
"""
foodweb_ids() = sort!(collect(keys(foodweb_catalog())))

"""
    foodweb(id)

Read the empirical food web `id`, of the form `"collection/web"` (see [`foodweb_ids`](@ref)),
as a [`FoodWeb`](@ref) whose metadata is the web's catalog entry: its title, citation, source
and license. For example `foodweb("williams_martinez_2000/skipwith_pond")`.
"""
function foodweb(id::AbstractString)
    catalog = foodweb_catalog()
    haskey(catalog, id) || throw(ArgumentError("no food web \"$id\"; see foodweb_ids()"))
    dir, entry = catalog[id]
    return read_foodweb(dir; name = get(entry, "title", id), metadata = entry)
end
