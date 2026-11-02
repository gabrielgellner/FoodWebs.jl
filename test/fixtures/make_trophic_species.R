# Trophic species of the bundled cheddar webs and their lumped webs, for the trophic species
# tests in test/aggregation.jl. Computes with the R package cheddar (Hudson, Reuman & Emerson
# 2013, BSD-2-Clause), version 0.1-640 as bundled in data/cheddar, and writes:
#
# - trophic_species.csv: the trophic species of every node (TrophicSpecies), with nodes in
#   cheddar's order, which data/cheddar keeps;
# - trophic_species_lumped.csv: every lumped node (LumpTrophicSpecies) with its mass M and
#   abundance N where the web has them, and its resources as lumped node numbers separated by
#   spaces. cheddar weights M by N (the plain mean in webs without N) and takes the mean of N.
#
# Run from the repository root:
#
#     Rscript test/fixtures/make_trophic_species.R
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

number <- function(x) ifelse(is.na(x), "missing", sprintf("%.17g", x))

# TrophicSpecies compares the resources and consumers of two nodes as lists in the order of the
# trophic links, so it separates nodes with the same resources and consumers listed in different
# orders (two nodes of Dargall Lane and of Hardknott Gill). Sorting the links by resource and then
# consumer lists both in node order.
sorted_links <- function(community) {
  nodes <- NP(community, "node")
  links <- TLPS(community)
  links <- links[order(match(links$resource, nodes), match(links$consumer, nodes)), , drop = FALSE]
  Community(nodes = NPS(community), trophic.links = links, properties = CPS(community))
}
webs <- lapply(webs, sorted_links)

species <- do.call(rbind, lapply(names(webs), function(id) {
  community <- webs[[id]]
  data.frame(web = paste0("cheddar/", id), node = seq_len(NumberOfNodes(community)),
             trophic_species = unname(TrophicSpecies(community)))
}))

lumped <- do.call(rbind, lapply(names(webs), function(id) {
  # webs without abundances take the plain mean of M
  weight.by <- if ("N" %in% NodePropertyNames(webs[[id]])) "N" else NULL
  community <- LumpTrophicSpecies(webs[[id]], weight.by = weight.by)
  nodes <- NPS(community)
  resources <- ResourcesByNode(community)
  index <- setNames(seq_len(nrow(nodes)), nodes$node)
  data.frame(
    web = paste0("cheddar/", id), group = seq_len(nrow(nodes)),
    M = if ("M" %in% colnames(nodes)) number(nodes$M) else "missing",
    N = if ("N" %in% colnames(nodes)) number(nodes$N) else "missing",
    resources = sapply(nodes$node, function(n) paste(sort(index[resources[[n]]]), collapse = " "))
  )
}))

write.csv(species, file.path("test", "fixtures", "trophic_species.csv"), row.names = FALSE, quote = FALSE)
write.csv(lumped, file.path("test", "fixtures", "trophic_species_lumped.csv"), row.names = FALSE, quote = FALSE)
