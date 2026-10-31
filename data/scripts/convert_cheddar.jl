# Convert the food webs shipped with the R package cheddar (BSD-2-Clause) to the package's
# data format: one directory per web with `nodes.csv` (the node table, species name first) and
# `trophic.links.csv` (`resource,consumer`), a `catalog.toml` describing each web and a
# `datapackage.json` schema.
#
# Run from the repository root:
#     julia --project=data/scripts data/scripts/convert_cheddar.jl

using CSV: CSV
using DataFrames: DataFrame, select
using Downloads: download
using RData: RData
using TOML: TOML

include("datapackage.jl")

const VERSION = "0.1-640"
const URL = "https://cran.r-project.org/src/contrib/cheddar_$VERSION.tar.gz"
const OUT = joinpath(@__DIR__, "..", "cheddar")

# Citations of the original data, from cheddar's documentation. Two are corrected: Yodzis's
# Benguela paper is 1998 (J. Anim. Ecol. 67), and Bersier et al. 2002 is Ecology 83.
const CITATIONS = Dict(
    "TL84" => "Carpenter, S. R. & Kitchell, J. F. (eds) The Trophic Cascade in Lakes (Cambridge University Press, 1996); Cohen, J. E., Jonsson, T. & Carpenter, S. R. Ecological community description using the food web, species abundance, and body size. PNAS 100, 1781–1786 (2003); Jonsson, T., Cohen, J. E. & Carpenter, S. R. Food webs, body size, and species abundance in ecological community description. Adv. Ecol. Res. 36, 1–84 (2005).",
    "YthanEstuary" => "Hall, S. J. & Raffaelli, D. Food-web patterns: lessons from a species-rich web. J. Anim. Ecol. 60, 823–841 (1991); Emmerson, M. C. & Raffaelli, D. Predator–prey body size, interaction strength and the stability of a real food web. J. Anim. Ecol. 73, 399–409 (2004).",
    "Benguela" => "Yodzis, P. Local trophodynamics and the interaction of marine mammals and fisheries in the Benguela ecosystem. J. Anim. Ecol. 67, 635–658 (1998).",
    "SkipwithPond" => "Warren, P. H. Spatial and temporal variation in the structure of a freshwater food web. Oikos 55, 299–311 (1989).",
    "BroadstoneStream" => "Woodward, G., Speirs, D. C. & Hildrew, A. G. Quantification and resolution of a complex, size-structured food web. Adv. Ecol. Res. 36, 85–135 (2005).",
    "ChesapeakeBay" => "Baird, D. & Ulanowicz, R. E. The seasonal dynamics of the Chesapeake Bay ecosystem. Ecol. Monogr. 59, 329–364 (1989); Bersier, L.-F., Banašek-Richter, C. & Cattin, M.-F. Quantitative descriptors of food-web matrices. Ecology 83, 2394–2407 (2002).",
    "Millstream" => "Ledger, M. E., Brown, L. E., Edwards, F. K., Milner, A. M. & Woodward, G. Drought alters the structure and functioning of complex food webs. Nat. Clim. Change 3, 223–227 (2013); Woodward, G. et al. Climate change impacts in multispecies systems: drought alters food web size-structure in a field experiment. Phil. Trans. R. Soc. B 367, 2990–2997 (2012).",
    "pHWebs" => "Layer, K., Riede, J. O., Hildrew, A. G. & Woodward, G. Food web structure and stability in 20 streams across a wide pH gradient. Adv. Ecol. Res. 42, 265–299 (2010).",
)
CITATIONS["TL86"] = CITATIONS["TL84"]

const CHEDDAR = "Hudson, L. N., Reuman, D. C. & Emerson, R. Cheddar: analysis and visualisation of ecological communities in R. Methods Ecol. Evol. 4, 99–104 (2013)."

slug(s) = replace(lowercase(s), r"[^a-z0-9]+" => "_") |> x -> strip(x, '_')

# A cheddar community is a dictionary with "nodes", "trophic.links" and "properties"
function write_web(dir, community)
    mkpath(dir)
    nodes = community["nodes"]
    links = community["trophic.links"]
    CSV.write(joinpath(dir, "nodes.csv"), nodes)
    CSV.write(joinpath(dir, "trophic.links.csv"), select(links, "resource", "consumer"))
    props = community["properties"]
    # a DictoVec iterates over its values, so look the keys up
    props = Dict(String(k) => (v isa AbstractString ? String(v) : v) for k in keys(props) for v in (props[k],))
    return props, size(nodes, 1), size(links, 1)
end

function main()
    work = mktempdir()
    tarball = download(URL, joinpath(work, "cheddar.tar.gz"))
    run(`tar -xzf $tarball -C $work`)
    datadir = joinpath(work, "cheddar", "data")
    rm(OUT; force = true, recursive = true)
    mkpath(OUT)
    cp(joinpath(work, "cheddar", "LICENSE"), joinpath(OUT, "LICENSE"))
    catalog = Dict{String, Any}()
    for file in sort(readdir(datadir))
        name = splitext(file)[1]
        data = RData.load(joinpath(datadir, file))[name]
        # a CommunityCollection is a dictionary of communities; a Community has "nodes"
        communities = haskey(data, "nodes") ? [(name, data)] : [(String(k), data[k]) for k in keys(data)]
        for (title, community) in communities
            id = length(communities) == 1 ? slug(name) : slug(name) * "_" * slug(title)
            props, S, L = write_web(joinpath(OUT, id), community)
            catalog[id] = Dict(
                "title" => get(props, "title", title),
                "dataset" => name,
                "species" => S,
                "links" => L,
                "citation" => CITATIONS[name],
                "source" => "cheddar $VERSION ($URL), $CHEDDAR",
                "license" => "BSD-2-Clause (see LICENSE)",
                "properties" => props,
            )
        end
    end
    open(joinpath(OUT, "catalog.toml"), "w") do io
        println(io, "# Generated by data/scripts/convert_cheddar.jl; do not edit by hand.")
        TOML.print(io, catalog; sorted = true)
    end
    write_datapackage(
        OUT; name = "cheddar", title = "Food webs of the R package cheddar",
        description = "The food webs shipped with cheddar $VERSION, converted by data/scripts/convert_cheddar.jl.",
        licenses = [Dict("name" => "BSD-2-Clause", "path" => "LICENSE", "title" => "BSD 2-Clause License")],
        sources = [Dict("title" => CHEDDAR, "path" => URL)],
    )
    return println("wrote $(length(catalog)) webs to $(normpath(OUT))")
end

main()
