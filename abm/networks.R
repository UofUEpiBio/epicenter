library(ergm)
library(sna)
library(data.table)

nets <- readRDS("../models/2022-04-25-bipartite-ergms.rds")[["Model 5"]]
log_n <- log(summary(nets$network ~ degrange(0, by = "net_id")))

# Saving attributes
X <- as.data.frame(nets$network, unit = "vertices") |>
  as.data.table()

X[, .(as.integer(is_actor), res_age, ventilator)] |>
  unlist() |> unname() |> cat(file = "actor_attributes.txt", sep = "\n")



# Simulating ERGM networks
set.seed(123)

fnames <- list.files(
  "networks",
  pattern = "ergm-[0-9]+\\.txt",
  full.names = TRUE
  )

for (i in 1:1000) {

  fname <- sprintf("networks/ergm-%04i.txt", i)
  if (fname %in% fnames)
    next

  n <- simulate(nets, 1)

  fwrite(
    as.edgelist(n) - 1,
    fname,
    sep = " ",
    col.names = FALSE
    )

  if (!i %% 50)
    message("Network ", sprintf("% 5i", i), " done...")

}

# Simulating random from the baseline preserving degree
# sequence

# Splitting the networks
nets_individual <- nets$network
vnames          <- nets_individual %v% "vertex.names"
nets_names      <- unique(nets_individual %v% "net_id")

nets_individual <- parallel::mclapply(nets_names, \(id) {
  get.inducedSubgraph(
    nets_individual,
    which(nets_individual %v% "net_id" == id)
    )
}, mc.cores = 4L)

# We will use netdiffuseR's rewiring algorithm which is
# fast
library(netdiffuseR)

fnames <- list.files(
  "networks",
  pattern = "degseq-[0-9]+\\.txt",
  full.names = TRUE
  )

for (i in 1:1000) {

  fname <- sprintf("networks/degseq-%04i.txt", i)
  if (fname %in% fnames)
    next

  # Empty network
  empty <- nets$network
  empty[,] <- 0L

  # Getting the new ties
  count <- 1L
  for (n in nets_individual) {

    # Retrieving the edgelist
    el_tmp <- as.edgelist(n)
    n_sparse <- matrix(0L, nrow = network.size(n), ncol = network.size(n))

    dimnames(n_sparse) <- list(
      attr(el_tmp, "vnames"),
      attr(el_tmp, "vnames")
    )

    n_sparse[el_tmp] <- 1L

    # Rewiring with degree-sequence
    n_sparse <- rewire_graph(
      n_sparse,
      p = network.edgecount(n) * 15,
      algorithm = "swap"
      ) |> as.matrix()

    ids <- match(n %v% "vertex.names", vnames)

    e  <- which(n_sparse != 0, arr.ind = TRUE)
    e[] <- ids[as.vector(e)]

    add.edges(empty, e[, 1], e[, 2])

    # if (!count %% 10)
    #   message("Network ", count, " in ", i, " done.")

    count <- count + 1

  }

  # Must have the same sequence
  all(degree(empty) == degree(nets$network))

  fwrite(
    as.edgelist(empty) - 1,
    fname,
    sep = " ",
    col.names = FALSE
    )

  if (!i %% 50)
    message("Network ", sprintf("% 5i", i), " done...")

}

?rewire_graph
