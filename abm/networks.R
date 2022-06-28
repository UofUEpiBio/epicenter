library(ergm)
library(sna)
library(data.table)

nets <- readRDS("models/2022-04-25-bipartite-ergms.rds")[["Model 5"]]
log_n <- log(summary(nets$network ~ degrange(0, by = "net_id")))

# Saving attributes
X <- as.data.frame(nets$network, unit = "vertices") |>
  as.data.table()

X[, .(as.integer(is_actor), as.integer(res_age <= 65), ventilator)] |>
  unlist() |> unname() |> cat(file = "abm/actor_attributes.txt", sep = "\n")



# Simulating ERGM networks
set.seed(123)
z <- simulate(nets, 10)

for (i in seq_along(z)) {

  fwrite(
    as.edgelist(z[[i]]) - 1,
    sprintf("abm/networks/ergm-%04i.txt", i),
    sep = " ",
    col.names = FALSE
    )
    
}
