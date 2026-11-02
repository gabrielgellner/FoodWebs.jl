# Structural properties of the bundled cheddar webs, for the tests against cheddar in
# test/properties.jl. Computes with the R package cheddar (Hudson, Reuman & Emerson 2013,
# BSD-2-Clause), version 0.1-640 as bundled in data/cheddar, and writes:
#
# - structural_properties.csv: per web, the numbers of basal, top, intermediate and isolated
#   nodes (cheddar ignores cannibalistic links in these), of cannibals, the SD of normalized
#   generality and vulnerability, and the mean maximum trophic similarity (MeanMaximumTrophicSimilarity);
# - trophic_levels.csv: per node, the prey-averaged, shortest and short-weighted trophic levels
#   with cannibalistic links removed, which FoodWebs.jl ignores (NA where cheddar has none).
#
# cheddar's trophic similarity pools each node's resources and consumers into one set, where
# Williams & Martinez (2000) count shared prey and shared predators separately.
#
# Run from the repository root:
#
#     Rscript test/fixtures/make_structural_properties.R
suppressPackageStartupMessages(library(cheddar))
stopifnot(as.character(packageVersion("cheddar")) == "0.1.640")
data(Benguela, BroadstoneStream, ChesapeakeBay, Millstream, SkipwithPond, TL84, TL86, YthanEstuary, pHWebs)

webs <- list(
  benguela = Benguela, broadstonestream = BroadstoneStream, chesapeakebay = ChesapeakeBay,
  millstream_c4 = Millstream[["c4"]], millstream_d4 = Millstream[["d4"]], skipwithpond = SkipwithPond,
  tl84 = TL84, tl86 = TL86, ythanestuary = YthanEstuary,
  phwebs_old_lodge = pHWebs[["Old Lodge"]], phwebs_afon_hafren = pHWebs[["Afon Hafren"]],
  phwebs_broadstone = pHWebs[["Broadstone"]], phwebs_dargall_lane = pHWebs[["Dargall Lane"]],
  phwebs_mosedal_beck = pHWebs[["Mosedal Beck"]], phwebs_duddon_pike_beck = pHWebs[["Duddon Pike Beck"]],
  phwebs_allt_a_mharcaidh = pHWebs[["Allt a'Mharcaidh"]], phwebs_hardknott_gill = pHWebs[["Hardknott Gill"]],
  phwebs_bere_stream = pHWebs[["Bere Stream"]], phwebs_mill_stream = pHWebs[["Mill Stream"]]
)

number <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

properties <- do.call(rbind, lapply(names(webs), function(id) {
  community <- webs[[id]]
  data.frame(
    web = paste0("cheddar/", id), basal = length(BasalNodes(community)),
    top = length(TopLevelNodes(community)), intermediate = length(IntermediateNodes(community)),
    isolated = length(IsolatedNodes(community)), cannibals = length(Cannibals(community)),
    GenSD = number(sd(NormalisedTrophicGenerality(community))),
    VulSD = number(sd(NormalisedTrophicVulnerability(community))),
    MxSim = number(MeanMaximumTrophicSimilarity(community))
  )
}))

levels <- do.call(rbind, lapply(names(webs), function(id) {
  community <- RemoveCannibalisticLinks(webs[[id]])
  tl <- TrophicLevels(community)
  data.frame(
    web = paste0("cheddar/", id), node = seq_len(NumberOfNodes(community)),
    prey_averaged = number(PreyAveragedTrophicLevel(community)),
    shortest = number(tl[, "ShortestTL"]), short_weighted = number(tl[, "ShortWeightedTL"])
  )
}))

write.csv(properties, file.path("test", "fixtures", "structural_properties.csv"), row.names = FALSE, quote = FALSE)
write.csv(levels, file.path("test", "fixtures", "trophic_levels.csv"), row.names = FALSE, quote = FALSE)
