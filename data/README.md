# Empirical food webs

Each directory below `data/` is a collection of webs. A web is a directory with `nodes.csv`,
whose first column names the species, and `trophic.links.csv` with columns `resource,consumer`,
one row per feeding link, the layout of the R package cheddar. Each collection has a
`catalog.toml` giving every web's title, species and link counts, citation, source and license,
which the package reads, and a `datapackage.json`
([Frictionless Data Package](https://specs.frictionlessdata.io/data-package/)) declaring each
file's columns and types, with foreign keys from both link columns to the node names. Load webs
with `foodweb("collection/web")`; list them with `foodweb_ids()`.

## Bundled collections

| Collection | Webs | Source | License |
|---|---|---|---|
| `cheddar` | 19: Tuesday Lake 1984 and 1986, Ythan Estuary, Benguela, Skipwith Pond, Broadstone Stream, Chesapeake Bay, two Millstream mesocosms, ten pH-gradient streams | R package cheddar 0.1-640 (Hudson, Reuman & Emerson 2013), converted by `scripts/convert_cheddar.jl` | BSD-2-Clause, `cheddar/LICENSE` |
| `williams_martinez_2000` | 7 trophic-species webs: Skipwith Pond, Coachella Valley, Chesapeake Bay, St. Martin Island, Bridge Brook Lake, Little Rock Lake, Ythan Estuary | The webs of Williams & Martinez (2000) and Dunne, Williams & Martinez (2002), converted by `scripts/convert_williams_martinez.jl` from the original `.dat` files | None given; redistributed with citations of the original data |

The nodes of the cheddar webs keep cheddar's columns: category, body mass `M`, numerical
abundance `N` and taxonomy where available. Two cheddar citations are corrected in the catalog:
Yodzis's Benguela paper is from 1998, and Bersier et al. (2002) is Ecology 83. The
Williams & Martinez webs name no species, so their species are numbered.

## Downloaded collections

These databases allow use with attribution but do not clearly allow redistribution, so their webs
are downloaded to `FoodWebs.foodweb_cache()` (the `FOODWEBS_CACHE` environment variable, or a
scratch space in the Julia depot) instead of being committed:

```sh
julia --project=data/scripts data/scripts/fetch_web_of_life.jl
julia --project=data/scripts data/scripts/fetch_mangal.jl [dataset ...]
```

- `web_of_life`: the 16 binary food webs of the Web of Life database (Thompson & Townsend stream
  webs). Its terms ask users to acknowledge "This work has used the Web of Life database
  (www.web-of-life.es)" and to cite the original authors.
- `mangal`: every Mangal network whose interactions are all trophic (predation, herbivory,
  detritivory or scavenging). Mangal's older documentation releases its data under the CC0
  waiver; the current site states no license.

KONECT's Little Rock Lake web is left out: it has no species names and adds a node feeding every
producer. Its Florida Bay webs are unnamed carbon-flow networks with non-living compartments.

## References

- Hudson, L. N., Reuman, D. C. & Emerson, R. Cheddar: analysis and visualisation of ecological
  communities in R. Methods Ecol. Evol. 4, 99–104 (2013).
- Williams, R. J. & Martinez, N. D. Simple rules yield complex food webs. Nature 404, 180–183
  (2000).
- Dunne, J. A., Williams, R. J. & Martinez, N. D. Food-web structure and network theory: the
  role of connectance and size. PNAS 99, 12917–12922 (2002).
