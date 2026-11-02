# Body mass summaries of the cheddar webs with body mass and abundance, for the body mass tests
# in test/biomass.jl. Computes with the R package cheddar (Hudson, Reuman & Emerson 2013,
# BSD-2-Clause), version 0.1-640 as bundled in data/cheddar, and writes body_mass.csv:
#
# - slope, intercept: the regression of log10(N) on log10(M) over species (NvMSlopeAndIntercept);
# - ratios, sum_log10_ratio: the number of predator-to-prey mass ratios without cannibalism and
#   the sum of their log10 (CRMRatio);
# - links, mean_length, mean_angle: the links in the log M-log N plane, without cannibalism and
#   species missing M or N, and their mean length and angle (NvMTriTrophicStatistics, after
#   Cohen, Schittler, Raffaelli & Reuman 2009).
#
# Run from the repository root:
#
#     Rscript test/fixtures/make_body_mass.R
suppressPackageStartupMessages(library(cheddar))
stopifnot(as.character(packageVersion("cheddar")) == "0.1.640")
data(TL84, TL86, YthanEstuary, BroadstoneStream, Millstream, pHWebs)

webs <- list(
  tl84 = TL84, tl86 = TL86, ythanestuary = YthanEstuary, broadstonestream = BroadstoneStream,
  millstream_c4 = Millstream[["c4"]], millstream_d4 = Millstream[["d4"]],
  phwebs_old_lodge = pHWebs[["Old Lodge"]], phwebs_afon_hafren = pHWebs[["Afon Hafren"]],
  phwebs_broadstone = pHWebs[["Broadstone"]], phwebs_dargall_lane = pHWebs[["Dargall Lane"]],
  phwebs_mosedal_beck = pHWebs[["Mosedal Beck"]], phwebs_duddon_pike_beck = pHWebs[["Duddon Pike Beck"]],
  phwebs_allt_a_mharcaidh = pHWebs[["Allt a'Mharcaidh"]], phwebs_hardknott_gill = pHWebs[["Hardknott Gill"]],
  phwebs_bere_stream = pHWebs[["Bere Stream"]], phwebs_mill_stream = pHWebs[["Mill Stream"]]
)

rows <- do.call(rbind, lapply(names(webs), function(id) {
  community <- webs[[id]]
  fit <- NvMSlopeAndIntercept(community)
  ratios <- CRMRatio(RemoveCannibalisticLinks(community))
  ratios <- ratios[!is.na(ratios)]
  links <- NvMTriTrophicStatistics(community)$links
  data.frame(
    web = paste0("cheddar/", id), slope = fit[["slope"]], intercept = fit[["intercept"]],
    ratios = length(ratios), sum_log10_ratio = sum(log10(ratios)), links = nrow(links),
    mean_length = mean(links$length), mean_angle = mean(links$angle)
  )
}))
out <- file.path("test", "fixtures", "body_mass.csv")
numeric <- sapply(rows, is.numeric)
rows[numeric] <- lapply(rows[numeric], function(x) sprintf("%.17g", x))
write.csv(rows, out, row.names = FALSE, quote = FALSE)
