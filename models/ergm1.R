library(data.table)
library(igraph)
library(netplot)
library(ergm)


# Reading the data in
dat <- fread("data/nelsons-modeling/2020-11-24/matched MI data 20201124.csv")

# HCW data
cnames <- colnames(dat)
cnames[grepl("hcw", cnames, ignore.case = TRUE)]
  

# Getting nodes
edges <- dat[, .(state, date, Fac_Name, Fac_Unit, hcw_id, hcw_type, res_id)]

namemilize <- function(...) {
  do.call(
    paste,
    c(lapply(list(...), \(x) gsub("[^[:alnum:]]", "_", trimws(x))),
    list(sep = "-"))
  ) |>
    tolower() |>
    gsub(pattern = "_+", replacement = "_")
    
}

edges[, loc_id := namemilize(state, date, Fac_Name, Fac_Unit)]

nodes <- rbind(
  edges[, .(id = hcw_id, loc_id, state, Fac_Name, Fac_Unit, type = hcw_type, group = "hwc")],
  edges[, .(id = res_id, loc_id, state, Fac_Name, Fac_Unit, type = "resident", group = "residents")]
) |> unique()

networks <- vector("list", length(unique(nodes$loc_id)))
names(networks) <- unique(nodes$loc_id)
networks_proj <- networks

nerrors <- 0L
for (n in names(networks)) {
  
  edges_n <- edges[loc_id == n]
  nodes_n <- nodes[loc_id == n]
    
  
  # Trying to read the data
  networks[[n]] <- tryCatch(graph_from_data_frame(
    d        = edges_n,
    vertices = nodes_n,
    directed = FALSE
  ), error = function(e) e)
  
  # Notification!
  if (inherits(networks[[n]], what = "error")) {
    message("Network ", n, " failed to be built")
    nerrors <- nerrors + 1L
  }
  
  # Generating projection matrix
  networks_proj[[n]] <- merge(
    x = edges_n[,.(ego   = res_id, hcw_id)] |> unique(),
    y = edges_n[,.(alter = res_id, hcw_id)] |> unique(),
    by = "hcw_id", allow.cartesian = TRUE
  )[ego < alter]
  
}

message("Done! ", nerrors, "/", length(networks), " failed.")

# Example of a graph
nplot(networks$`ga-6_6_2019-ag_rhodes_health_rehab-rehab_unit_3rd_floor`)

# Projection networks ----------------------------------------------------------

