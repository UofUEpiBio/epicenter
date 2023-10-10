library(ergm)
library(sna)
library(data.table)

nets <- readRDS("models/2022-04-25-bipartite-ergms.rds")[["Model 5"]]

degs <- degree(nets$network, cmode = "outdegree")

degs <- data.table(
  degree = degs,
  type   = fifelse(nets$network %v% "is_actor", "HCW", "Resident")
)

degs[, mean(degree), by = "type"]
