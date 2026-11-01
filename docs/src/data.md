# Empirical food webs

FoodWebs.jl bundles 26 empirical food webs:

- `cheddar`: the 19 webs of the R package cheddar (BSD-2-Clause): Tuesday Lake 1984 and 1986,
  Ythan Estuary, Benguela, Skipwith Pond, Broadstone Stream, Chesapeake Bay, two Millstream
  mesocosms and ten pH-gradient streams, with body masses, abundances and taxonomy.
- `williams_martinez_2000`: the seven trophic-species webs of Williams & Martinez (2000), the
  reference webs for the structural properties: Skipwith Pond, Little Rock Lake, Bridge Brook
  Lake, Chesapeake Bay, Ythan Estuary, Coachella Valley and St. Martin Island.

```@example data
using FoodWebs

fw = foodweb("williams_martinez_2000/skipwith_pond")
```

```@example data
fw.metadata["citation"]
```

```@example data
first(foodweb_ids(), 5)
```

A [`FoodWeb`](@ref) is an `AbstractMatrix{Bool}`, so every function that takes a web takes it.

## Format

Each web is a directory with `nodes.csv`, whose first column names the species, and
`trophic.links.csv` with columns `resource,consumer`, one row per feeding link, the layout of
cheddar. Each collection has a `catalog.toml` with every web's title, species and link counts,
citation, source and license, and a `datapackage.json`
([Frictionless Data Package](https://specs.frictionlessdata.io/data-package/)) declaring the
files' columns, types and foreign keys. [`read_foodweb`](@ref) reads any directory in this
format.

## Downloaded webs

The Web of Life and Mangal databases allow use with attribution but do not clearly allow
redistribution, so scripts in the repository's `data/scripts` directory download their food
webs to [`FoodWebs.foodweb_cache`](@ref), where [`foodweb_ids`](@ref) and [`foodweb`](@ref)
find them:

```sh
julia --project=data/scripts data/scripts/fetch_web_of_life.jl   # 16 webs
julia --project=data/scripts data/scripts/fetch_mangal.jl        # 218 webs
```

See `data/README.md` for sources, licenses and how the bundled data were converted.

```@docs
FoodWeb
foodweb
foodweb_ids
read_foodweb
FoodWebs.foodweb_cache
```
