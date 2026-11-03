# Flow-based network indices of the Chesapeake Bay carbon flow network (Baird & Ulanowicz 1989),
# for the tests in test/flow_indices.jl. Computes with the R packages enaR (Borrett & Lau 2014,
# GPL-3), version 3.1.1, and NetIndices (Kones et al. 2009, GPL >= 2), version 1.4.4.1, and writes:
#
# - flow_indices_nodes.csv: each compartment's name, whether it is living, its input, export and
#   respiration, and its trophic level and omnivory index (NetIndices TrophInd, with the non-living
#   compartments as Dead);
# - flow_indices_flows.csv: the flows between compartments, from compartment `from` to `to`;
# - flow_indices.csv: total system throughflow and throughput, average path length and Finn's
#   cycling index (enaFlow), and flow diversity, average mutual information, development capacity,
#   ascendency and overhead in bits (enaAscendency);
# - flow_indices_mti.csv: the mixed trophic impact of compartment `from` on `to` (enaMTI).
#
# Run from the repository root:
#
#     Rscript test/fixtures/make_flow_indices.R
suppressPackageStartupMessages({
  library(enaR)
  library(NetIndices)
  library(network)
})
stopifnot(as.character(packageVersion("enaR")) == "3.1.1")
stopifnot(as.character(packageVersion("NetIndices")) == "1.4.4.1")
data(troModels)
x <- troModels[["Chesapeake Bay "]]
stopifnot(ssCheck(x))

fmt <- function(df) {
  numeric <- sapply(df, is.double)
  df[numeric] <- lapply(df[numeric], function(v) sprintf("%.17g", v))
  df
}
out <- function(df, name) {
  write.csv(fmt(df), file.path("test", "fixtures", name), row.names = FALSE, quote = FALSE)
}

flow <- as.matrix(x, attrname = "flow")
living <- x %v% "living"
n <- nrow(flow)
# NetIndices takes flows with the boundary as extra rows and columns: Tij[i, j] is the flow from
# j to i
tij <- rbind(
  cbind(t(flow), matrix(0, n, 1)),
  c(x %v% "export" + x %v% "respiration", 0)
)
tij[1:n, n + 1] <- x %v% "input"
rownames(tij) <- c(1:n, "out")
colnames(tij) <- c(1:n, "in")
troph <- TrophInd(Tij = tij, Import = "in", Export = "out", Dead = which(!living))
names <- gsub(",", ";", trimws(x %v% "vertex.names"))
out(data.frame(
  index = 1:n, name = names, living = living, input = x %v% "input", export = x %v% "export",
  respiration = x %v% "respiration", trophic_level = troph$TL, omnivory = troph$OI
), "flow_indices_nodes.csv")

idx <- which(flow != 0, arr.ind = TRUE)
out(data.frame(from = idx[, 1], to = idx[, 2], flow = flow[idx]), "flow_indices_flows.csv")

f <- enaFlow(x)$ns
a <- enaAscendency(x)
out(data.frame(
  index = c("TST", "TSTp", "APL", "FCI", "H", "AMI", "CAP", "ASC", "OH"),
  value = c(f[, "TST"], f[, "TSTp"], f[, "APL"], f[, "FCI"], a[, "H"], a[, "AMI"], a[, "CAP"], a[, "ASC"], a[, "OH"])
), "flow_indices.csv")

M <- enaMTI(x)$M
out(data.frame(from = rep(1:n, n), to = rep(1:n, each = n), value = as.vector(M)), "flow_indices_mti.csv")
