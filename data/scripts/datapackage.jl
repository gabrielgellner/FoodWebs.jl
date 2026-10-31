# Write a Frictionless Data Package (https://specs.frictionlessdata.io/data-package/) describing a
# collection: one `nodes` and one `trophic.links` resource per web, with column types inferred
# from the files and foreign keys from both link columns to the node names.

using CSV: CSV
using DataFrames: DataFrame
using JSON: JSON

# "number" if every non-empty value parses as a number, otherwise "string"
function column_type(values)
    present = filter(v -> !ismissing(v) && !isempty(strip(string(v))), values)
    isempty(present) && return "string"
    return all(v -> v isa Real || tryparse(Float64, string(v)) !== nothing, present) ? "number" : "string"
end

function table_schema(path)
    df = CSV.read(path, DataFrame; stringtype = String)
    return [Dict("name" => n, "type" => column_type(df[!, n])) for n in names(df)]
end

"""
    write_datapackage(dir; name, title, description, licenses, sources)

Describe every web in the collection directory `dir` (each subdirectory with `nodes.csv` and
`trophic.links.csv`) in `dir/datapackage.json`.
"""
function write_datapackage(dir; name, title, description, licenses, sources)
    resources = Any[]
    for web in sort(filter(d -> isdir(joinpath(dir, d)), readdir(dir)))
        nodes = joinpath(web, "nodes.csv")
        links = joinpath(web, "trophic.links.csv")
        isfile(joinpath(dir, nodes)) || continue
        node_schema = table_schema(joinpath(dir, nodes))
        # species names are text even when they are numbers
        node_schema[1]["type"] = "string"
        push!(
            resources, Dict(
                "name" => "$web-nodes", "path" => nodes, "format" => "csv", "mediatype" => "text/csv",
                "schema" => Dict("fields" => node_schema, "primaryKey" => node_schema[1]["name"]),
            )
        )
        push!(
            resources, Dict(
                "name" => "$web-trophic-links", "path" => links, "format" => "csv", "mediatype" => "text/csv",
                "schema" => Dict(
                    "fields" => [
                        Dict("name" => "resource", "type" => "string", "description" => "species eaten"),
                        Dict("name" => "consumer", "type" => "string", "description" => "species eating"),
                    ],
                    "primaryKey" => ["resource", "consumer"],
                    "foreignKeys" => [
                        Dict("fields" => f, "reference" => Dict("resource" => "$web-nodes", "fields" => node_schema[1]["name"]))
                            for f in ("resource", "consumer")
                    ],
                ),
            )
        )
    end
    package = Dict(
        "name" => name, "title" => title, "description" => description,
        "licenses" => licenses, "sources" => sources, "resources" => resources,
    )
    open(joinpath(dir, "datapackage.json"), "w") do io
        JSON.json(io, package; pretty = true)
        println(io)
    end
    return package
end
