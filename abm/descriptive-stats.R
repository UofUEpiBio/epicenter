library(ergm)
library(sna)
library(data.table)

load("LTCF_final_report/res10_new.RData")

nets <- res10_new

degs <- degree(nets$network, cmode = "outdegree")

degs <- data.table(
  degree = degs,
  type   = fifelse(nets$network %v% "is_actor", "HCW", "Resident")
)

degs[, mean(degree), by = "type"]
