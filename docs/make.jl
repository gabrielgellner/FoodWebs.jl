using Documenter
using FoodWebs
using CairoMakie: CairoMakie

CairoMakie.activate!(type = "png")

DocMeta.setdocmeta!(FoodWebs, :DocTestSetup, :(using FoodWebs); recursive = true)

makedocs(;
    modules = [FoodWebs, Base.get_extension(FoodWebs, :FoodWebsMakieExt)],
    sitename = "FoodWebs.jl",
    authors = "Gabriel Gellner",
    repo = Remotes.GitHub("gabrielgellner", "FoodWebs.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://gabrielgellner.github.io/FoodWebs.jl",
        edit_link = "main",
        size_threshold = 2^20,
    ),
    pages = [
        "Home" => "index.md",
        "Food web models" => "models.md",
        "Empirical food webs" => "data.md",
        "Trophic levels and structure" => "properties.md",
        "Equilibrium biomass" => "biomass.md",
        "Feasible communities" => "feasible.md",
        "Motifs" => "motifs.md",
        "Community matrices and stability" => "matrices.md",
        "Plotting" => "plotting.md",
        "API index" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/gabrielgellner/FoodWebs.jl.git", devbranch = "main")
