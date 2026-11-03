# Fluxes of the three example webs of the R package fluxweb (Gauzens et al. 2019, GPL >= 2),
# version 2.0.1, for the flux tests in test/energetics.jl. Writes:
#
# - fluxweb_species.csv: web, species, biomass, body mass, loss per unit biomass and assimilation
#   efficiency of each species, as passed to fluxing;
# - fluxweb_links.csv: web, prey, predator, the preference in fluxweb's matrix and the flux that
#   fluxing returns, one row per link;
# - fluxweb_jacobian.csv: the nonzero entries of create.jacob for species.level, and its
#   stability.value as the row with row = col = 0.
#
# simple.case gives its losses per unit biomass directly (met.rate, with body mass 1);
# species.level and groups.level use the average allometric losses 0.71 M^-0.25 of the fluxweb
# vignette. groups.level has integer preferences, which fluxing weights by prey biomass.
#
# Run from the repository root:
#
#     Rscript test/fixtures/make_fluxweb.R
suppressPackageStartupMessages(library(fluxweb))
stopifnot(as.character(packageVersion("fluxweb")) == "2.0.1")
data(simple.case, species.level, groups.level)

webs <- list(
  simple_case = list(
    mat = simple.case$mat, biomasses = simple.case$biomasses, M = rep(1, 4),
    losses = simple.case$met.rate, efficiencies = simple.case$efficiencies
  ),
  species_level = list(
    mat = species.level$mat, biomasses = species.level$biomasses, M = species.level$bodymasses,
    losses = 0.71 * species.level$bodymasses^-0.25, efficiencies = species.level$efficiencies
  ),
  groups_level = list(
    mat = groups.level$mat, biomasses = groups.level$biomasses, M = groups.level$bodymasses,
    losses = 0.71 * groups.level$bodymasses^-0.25, efficiencies = groups.level$efficiencies
  )
)

fmt <- function(df) {
  numeric <- sapply(df, is.double)
  df[numeric] <- lapply(df[numeric], function(x) sprintf("%.17g", x))
  df
}

species <- list()
links <- list()
for (id in names(webs)) {
  w <- webs[[id]]
  mat <- matrix(as.numeric(w$mat), nrow(w$mat))
  F <- fluxing(mat, w$biomasses, w$losses, w$efficiencies)
  S <- nrow(mat)
  species[[id]] <- data.frame(
    web = id, species = 1:S, biomass = w$biomasses, M = w$M, loss = w$losses,
    efficiency = w$efficiencies
  )
  idx <- which(mat != 0, arr.ind = TRUE)
  links[[id]] <- data.frame(
    web = id, prey = idx[, 1], predator = idx[, 2], preference = mat[idx], flux = F[idx]
  )
}
write.csv(fmt(do.call(rbind, species)), file.path("test", "fixtures", "fluxweb_species.csv"), row.names = FALSE, quote = FALSE)
write.csv(fmt(do.call(rbind, links)), file.path("test", "fixtures", "fluxweb_links.csv"), row.names = FALSE, quote = FALSE)

# create.jacob and stability.value need metabolic types; with none marked plant or detritus the
# diagonal holds only cannibalism, (e - 1) F[i, i] / x[i] for the 26 cannibals of species.level
w <- webs$species_level
F <- fluxing(w$mat, w$biomasses, w$losses, w$efficiencies)
types <- rep("animal", nrow(F))
J <- create.jacob(F, w$biomasses, w$efficiencies, types)
idx <- which(J != 0, arr.ind = TRUE)
jac <- data.frame(row = c(idx[, 1], 0), col = c(idx[, 2], 0), value = c(J[idx], stability.value(F, w$biomasses, w$efficiencies, types)))
write.csv(fmt(jac), file.path("test", "fixtures", "fluxweb_jacobian.csv"), row.names = FALSE, quote = FALSE)
