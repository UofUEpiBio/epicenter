library(ergm)
library(sna)
library(data.table)

models <- "ergm" # c("degseq", "permute", "ergm")

nets <- readRDS("../models/2022-04-25-bipartite-ergms.rds")[["Model 5"]]
log_n <- log(summary(nets$network ~ degrange(0, by = "net_id")))

# Saving attributes
X <- as.data.frame(nets$network, unit = "vertices") |>
  as.data.table()

# Networks of individuals who are in ventilators
X[, net_has_vent := sum(ventilator), by = "net_id"]
X[, net_nothas_vent := .N - net_has_vent, by = "net_id"]

X[net_has_vent >= 4 & net_nothas_vent >= 4, unique(net_id)]

# Only those who have ventilators
X <- X[net_has_vent >= 4 & net_nothas_vent >= 4]

X_int <- X[, .(as.integer(is_actor), res_age, ventilator, as.integer(as.factor(net_id)) - 1)]

X_int |>
  unlist() |> unname() |> cat(file = "actor_attributes.txt", sep = "\n")

# Generating permuted version of actors attributes
set.seed(331)
for (i in 1:1000) {

  # Generating permutation
  Z <- X[, .(is_actor, net_id, o = runif(.N))]
  
  Z <- cbind(X_int, Z)[, {
    data.table(V1,res_age, ventilator, V4)[order(o), ]
  }, by = .(is_actor, net_id)][, -c(1,2)]

  # Saving the data
  Z |> unlist() |> unname() |> cat(
    file = sprintf(
      "actor_attributes/permuted_attributes-%04i.txt",
      i),
      sep = "\n"
      )
}

# We will use netdiffuseR's rewiring algorithm which is
# fast
library(netdiffuseR)

# ------------------------------------------------------------------------------
# Original network
# ------------------------------------------------------------------------------
empty <- get.inducedSubgraph(
  nets$network,
  which((nets$network %v% "vertex.names") %in% X$vertex.names)
  )

fwrite(
  as.edgelist(empty) - 1,
  "networks/original.txt",
  sep = " ",
  col.names = FALSE
  )

# ------------------------------------------------------------------------------
# Simulating ERGM networks 
# ------------------------------------------------------------------------------
set.seed(123)

if ("ergm" %in% models) {

  for (i in 1:1000) {

    fname <- sprintf("networks/ergm-%04i.txt", i)
    n <- simulate(nets, 1)

    # Subsetting networks with ventilator
    n <- get.inducedSubgraph(
      n,
      which((n %v% "vertex.names") %in% X$vertex.names)
      )

    fwrite(
      as.edgelist(n) - 1,
      fname,
      sep = " ",
      col.names = FALSE
      )

    
    # Splitting the networks
    nets_individual <- n
    vnames          <- nets_individual %v% "vertex.names"
    nets_names      <- unique(nets_individual %v% "net_id")

    nets_individual <- parallel::mclapply(nets_names, \(id) {
      get.inducedSubgraph(
        nets_individual,
        which(nets_individual %v% "net_id" == id)
        )
    }, mc.cores = 4L)
    
    # Generating rewired version -------------------------------------------------
    fname <- sprintf("networks/ergm+degseq-%04i.txt", i)

    # Empty network
    empty    <- n
    empty[,] <- 0L

    # Getting the new ties
    for (n_i in nets_individual) {

      # Retrieving the edgelist
      el_tmp <- as.edgelist(n_i)
      n_sparse <- matrix(0L, nrow = network.size(n_i), ncol = network.size(n_i))

      dimnames(n_sparse) <- list(
        attr(el_tmp, "vnames"),
        attr(el_tmp, "vnames")
      )

      n_sparse[el_tmp] <- 1L

      # Rewiring with degree-sequence
      n_sparse <- rewire_graph(
        n_sparse,
        p = network.edgecount(n_i) * 15,
        algorithm = "swap"
        ) |> as.matrix()

      ids <- match(n_i %v% "vertex.names", vnames)

      e  <- which(n_sparse != 0, arr.ind = TRUE)
      e[] <- ids[as.vector(e)]

      add.edges(empty, e[, 1], e[, 2])

    }

    # Must have the same sequence
    stopifnot(all(degree(empty) == degree(n)))

    empty <- get.inducedSubgraph(
      empty,
      which((empty %v% "vertex.names") %in% X$vertex.names)
      )

    fwrite(
      as.edgelist(empty) - 1,
      fname,
      sep = " ",
      col.names = FALSE
      )

    # Generating the fully permuted version ------------------------------------
    fname <- sprintf("networks/ergm+bernoulli-%04i.txt", i)
    
    empty    <- n
    empty[,] <- 0L

    # Getting the new ties
    for (n_i in nets_individual) {

      # Actor count
      nactors <- sum(n_i %v% "is_actor")
      netsize <- network.size(n_i)
      indices_n_i <- lapply(1:nactors, \(j) {
        cbind(j, nactors:netsize)
      }) |> do.call(what=rbind)

      # Retrieving the edgelist
      el_tmp <- as.edgelist(n_i)
      n_sparse <- matrix(0L, nrow = network.size(n_i), ncol = network.size(n_i))

      dimnames(n_sparse) <- list(
        attr(el_tmp, "vnames"),
        attr(el_tmp, "vnames")
      )

      n_sparse[el_tmp] <- 1L

      # Rewiring free
      n_sparse[indices_n_i] <- sample(n_sparse[indices_n_i])

      ids <- match(n_i %v% "vertex.names", vnames)

      e  <- which(n_sparse != 0, arr.ind = TRUE)
      e[] <- ids[as.vector(e)]

      add.edges(empty, e[, 1], e[, 2])

    }

    if (!i %% 50)
      message("Network ", sprintf("% 5i", i), " done...")

    
    fwrite(
      as.edgelist(empty) - 1,
      fname,
      sep = " ",
      col.names = FALSE
      )

  }

}

# ------------------------------------------------------------------------------
# Simulating random from the baseline preserving degree 
# sequence
# ------------------------------------------------------------------------------

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

if ("degseq" %in% models) {

  for (i in 1:1000) {

    fname <- sprintf("networks/degseq-%04i.txt", i)

    # Empty network
    empty <- nets$network
    empty[,] <- 0L

    # Getting the new ties
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

    }

    # Must have the same sequence
    stopifnot(all(degree(empty) == degree(nets$network)))

    empty <- get.inducedSubgraph(
      empty,
      which((empty %v% "vertex.names") %in% X$vertex.names)
      )

    fwrite(
      as.edgelist(empty) - 1,
      fname,
      sep = " ",
      col.names = FALSE
      )

    if (!i %% 50)
      message("Network ", sprintf("% 5i", i), " done...")

  }

}

# ------------------------------------------------------------------------------
# Simulating permute
# ------------------------------------------------------------------------------
net_induced <- get.inducedSubgraph(
  nets$network, 
  which((nets$network %v% "vertex.names") %in% X$vertex.names)
  )

net_degrees <- degree(net_induced)

set.seed(881)
if ("permute" %in% models) {

  for (i in 1:1000) {

    fname <- sprintf("networks/permute-%04i.txt", i)

    # Empty network
    empty <- network.copy(nets$network)
    empty[,] <- 0L

    # Getting the new ties
    for (n in nets_individual) {

      # Skipping this network
      if (!all((n %v% "vertex.names") %in% X$vertex.names))
        next

      nn <- as.edgelist(n)

      # Permuting residents and hcw
      nactors <- sum(n %v% "is_actor")
      p <- c(
        sample(1:nactors),
        sample((nactors + 1):network.size(n))
        )
      
      nn[] <- p[as.vector(nn)]

      ids <- match((n %v% "vertex.names"), vnames)

      nn[] <- ids[as.vector(nn)]

      add.edges(empty, nn[, 1], nn[, 2])

    }

    # Not shared degree sequence
    empty <- get.inducedSubgraph(
      empty,
      which((empty %v% "vertex.names") %in% X$vertex.names)
      )

    fwrite(
      as.edgelist(empty) - 1,
      fname,
      sep = " ",
      col.names = FALSE
      )

    if (!i %% 50)
      message("Network ", sprintf("% 5i", i), " done...")

  }

}