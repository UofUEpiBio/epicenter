library(data.table)
library(ergm.multi)
library(ggplot2)
library(ggridges)
library(ggrepel)
library(texreg)
library(parallel)

#' Function to configure the ERGM model for the paper
#' @param nets A list of [network] objects.
#' @param boot Logical. Whether to return a boostrap sample or not
#' @return A list with the following three elements:
#' - `ids`: Integer vector of the selected networks.
#' - `nets`: A list of the selected networks.
#' - `filters`: A list of index vectors with the needed filtering for
#' the ERGMs
#' - `n_hcp`: Number of healthcare providers per network.
#' - `n_patients`: Number of patients per network.
#' - `facility_state`: Character vector with facilities' states. 
set_ergm <- function(
  nets,
  idx = NULL
  ) {

  if (is.null(idx))
    idx <- seq_along(nets)
  
  nets <- nets[idx] |> unname()

  # Updating category of PA/NP in type to "Nurse" (only 1)
  # and Physician with "Other" (only 3)
  for (i in seq_along(nets)) {

    types <- nets[[i]] %v% "type"
    types[types == "PA/NP"] <- "Nurse"
    types[types == "Physician"] <- "Other"
    nets[[i]] %v% "type" <- types

  }

  # Re-building the dat to be sure it is right
  dat <- lapply(seq_along(nets), \(i) {
    net <- nets[[i]]
    data.table(
      NetID = i,
      Nurse = net %v% "Nurse",
      CNA = net %v% "CNA",
      PTOTRT = net %v% "PTOTRT",
      Other = net %v% "Other",
      dialysis = net %v% "dialysis",
      wound_care = net %v% "wound_care",
      mdro = net %v% "mdro",
      bedridden = net %v% "bedridden",
      ventilator = net %v% "ventilator",
      wheelchair = net %v% "wheelchair",
      state = net %v% "state",
      facility = net %v% "Fac_Name",
      diabetes = net %v% "diabetes"
    )
  }) |> rbindlist()

  ans <- list()

  ans$Nurse  = dat[, n:=sum(Nurse==1), NetID][n>=1, unique(NetID)]
  ans$CNA    = dat[, n:=sum(CNA==1), NetID][n>=1, unique(NetID)]
  ans$PTOTRT = dat[, n:=sum(PTOTRT==1), NetID][n>=1, unique(NetID)]
  ans$Other  = dat[, n:=sum(Other==1), NetID][n>=1, unique(NetID)] 

  # For b2factor (resident type), the networks need to have at least 1 resident of that type. 
  ans$dialysis   <- dat[, n:=sum(dialysis=='Yes'), NetID][n>0, unique(NetID)]
  ans$wound_care <- dat[, n:=sum(wound_care=='Yes'), NetID][n>0, unique(NetID)] 
  ans$mdro       <- dat[, n:=sum(mdro=='Yes'), NetID][n>0, unique(NetID)]
  ans$bedridden  <- dat[, n:=sum(bedridden=='Yes'), NetID][n>0, unique(NetID)]  
  ans$diabetes   <- dat[, n:=sum(diabetes=='Yes'), NetID][n>0, unique(NetID)]
  ans$ventilator <- dat[, n:=sum(ventilator=='Yes'), NetID][n>0, unique(NetID)]
  ans$wheelchair <- dat[, n:=sum(wheelchair=='Yes'), NetID][n>0, unique(NetID)]

  # For b1startmix(2), at least two residents of such type need to be present in the network.
  ans$dialysis2   <- dat[, n:=sum(dialysis=='Yes'), NetID][n>1, unique(NetID)]
  ans$diabetes2   <- dat[, n:=sum(diabetes=='Yes'), NetID][n>1, unique(NetID)]
  ans$wound_care2 <- dat[, n:=sum(wound_care=='Yes'), NetID][n>1, unique(NetID)]
  ans$ventilator2 <- dat[, n:=sum(ventilator=='Yes'), NetID][n>1, unique(NetID)]
  ans$mdro2       <- dat[, n:=sum(mdro=='Yes'), NetID][n>1, unique(NetID)]
  ans$wheelchair2 <- dat[, n:=sum(wheelchair=='Yes'), NetID][n>1, unique(NetID)]

  n_hcp <- sapply(nets, \(net) summary(net~b1degrange(0)))
  n_patients <- sapply(nets, \(net) summary(net ~ b2degrange(0)))

  c(
    list(ids = idx, nets = nets, filters = ans),
    ans,
    list(
      n_hcp = n_hcp |> unname(),
      n_patients = n_patients |> unname(),
      n_tot = (n_hcp + n_patients) |> unname(),
      facility_state = sapply(nets, \(net) unique(net %v% "state")) |>
        as.factor() |> unname(),
      facility = sapply(nets, \(net) unique(net %v% "Fac_Name")) |>
        as.factor() |> unname()
    )
  )

}

# Wrapper function of texreg::* to remove the
# offset terms and make the table easier to read.
tabulator <- function(x, ..., tmethod = texreg::screenreg) {

  .mapper <- function(x) {
    if (!inherits(x, "ergm") && length(x) > 1) {
      ans <- unique(unlist(sapply(x, .mapper)))
      return(
        setNames(
          as.list(ans),
          ans
        )
      )
    }

    x_coefs <- coef(x) |> names()
    x_coefs <- x_coefs[!grepl("^offset", x_coefs)]
    setNames(
      as.list(x_coefs),
      x_coefs
    )
  }

  maps <- .mapper(x)

  nnets <- list("N nets" = sapply(x, \(n) {
    length(n$newnetwork$gal$.subnetcache$.NetworkID)
  }))

  # Check if method was MCMLE
  mcmc_results <- lapply(x, \(n) {
    if (n$control$main.method == "MCMLE") {
      list(
        "N steps" = prettyNum(summary(n$sample)$end, big.mark = ","),
        "Joint Geweke" = suppressWarnings(ergm::geweke.diag.mv(n$sample)$p.value[1]) |> unname()
      )
    } else
      list("N steps" = 0, "Joint Geweke" = 0)
  })

  nnets$"N steps" <- sapply(mcmc_results, \(y) y[["N steps"]])
  nnets$"Joint Geweke" <- sapply(mcmc_results, \(y) y[["Joint Geweke"]])

  tmethod(
    x,
    custom.coef.map = maps,
    ...,
    custom.gof.rows = nnets
    )

}

# Session-level registry: maps output name -> hash of its call expression.
# Created once when helpers.R is sourced; persists for the session.
.with_restore_registry <- new.env(hash = TRUE, parent = emptyenv())

#' Conditionally run and cache an ERGM model
#'
#' Wraps \code{with()} so that the result is saved to
#' \code{path/<name>.rds} the first time and simply read back on subsequent
#' runs.  Within a session, calling \code{with_restore} twice with the same
#' \code{name} but different code raises an error to prevent accidental
#' name re-use.
#'
#' @param name  Character string used as the variable name in the calling
#'   frame and as the file stem of the cached \code{.rds}.
#' @param data  The data environment forwarded to \code{with()}.
#' @param ...   Expression to evaluate inside \code{with(data, ...)} —
#'   typically an \code{ergm()} call.
#' @param path  Directory where \code{.rds} files are stored.
with_restore <- function(name, data, ..., path = "models/03-pooled-ergms-rds") {

  # Build a string from the unevaluated call expressions, then MD5-hash it
  # using tools::md5sum() (always available, no extra dependencies).
  expr_str <- paste(
    deparse(substitute(data)),
    deparse(substitute(list(...))),
    collapse = "\n"
  )
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(expr_str, tmp)
  call_hash <- unname(tools::md5sum(tmp))

  fpath <- file.path(path, paste0(name, ".rds"))

  # Registry check: error if the name appears again with different code.
  if (file.exists(fpath) && exists(name, envir = .with_restore_registry, inherits = FALSE)) {
    prev_hash <- get(name, envir = .with_restore_registry, inherits = FALSE)
    if (!identical(prev_hash, call_hash)) {
      stop(
        "with_restore: name '", name, "' was already used in this session ",
        "with different code. Use a different name or restart the session ",
        "to clear the registry."
      )
    }
  } else {
    assign(name, call_hash, envir = .with_restore_registry)
  }

  
  if (file.exists(fpath)) {
    message("Loading '", name, "' from cache.")
    res <- readRDS(fpath)
  } else {
    res <- with(data, ...)
    saveRDS(res, fpath)
  }
  assign(name, res, envir = parent.frame())
  invisible(res)
}

#' Restore a cached result
#' @param i The index of the cached result to restore.
#' @param cache The directory where the cached results are stored.
#' @return The cached result, or NULL if it doesn't exist.
restore_cache <- function(i, cache = NULL) {

  if (!length(cache))
    return(NULL)

  # Checking the filepath exists
  if (!dir.exists(cache)) {
    stop("Cache directory does not exist. To make things cleaner, we require the user create the directory first.")
  }

  fn <- file.path(cache, sprintf("%04i.rds", i))
  
  if (file.exists(fn)) {
    message("Loading cached result ", fn, ".")
    return(readRDS(fn))
  }
  NULL
}

#' Save a result to the cache
#' @param i The index of the cached result to save.
#' @param res The result to save.
#' @param cache The directory where the cached results are stored.
save_cache <- function(i, res, cache = NULL) {

  if (!length(cache))
    return(invisible(NULL))

  fn <- file.path(cache, sprintf("%04i.rds", i))
  saveRDS(res, fn)
}

#' Run a function with caching
#' @param i The index of the cached result to use.
#' @param cache The directory where the cached results are stored.
#' @param ... The arguments to pass to the function.
#' @return The result of the function, either from the cache or by running it.
run_with_cache <- function(i, cache, ...) {

  res <- restore_cache(i, cache)

  if (!length(res)) {
    res <- eval(...)
    save_cache(i, res, cache)
  }

  res

}

# Bootstrapping function of the ERGM model
boot_ergm <- function(
  networks,
  ergm_formula,
  ...,
  nboot = 100,
  cl = NULL,
  cache = NULL
) {

  # Preparing the boot
  nnets <- length(networks)

  # Generating the indices for the bootstrapping
  boot_indices <- sample.int(nnets, nboot * nnets, replace = TRUE) |>
    matrix(ncol = nboot, nrow = nnets)

  # Sometimes, the execution order can make things off,
  # so we will shuffle the execution order
  boot_order <- sample(seq_len(nboot), nboot)

  fitfun <- function(idx, networks,...) {

    # Preparing the networks
    dat <- set_ergm(nets = networks, idx = idx)

    # Unpack all dat elements into the local frame so the formula
    # environment can resolve Networks(nets) and any other objects.
    list2env(dat, envir = environment())

    # Ensuring it runs locally
    environment(ergm_formula) <- environment()

    res <- tryCatch(
      with(dat, ergm::ergm(ergm_formula, ...)),
      error = function(e) e)

    if (inherits(res, "error"))
      return(res)

    coef(res)
  }

  # Passing the data to the cl
  if (!is.null(cl)) {

    # Exporting the needed data for the cluster
    parallel::clusterExport(
      cl, varlist = c(
        "boot_indices", "fitfun", "networks", "ergm_formula",
        "restore_cache", "save_cache", "cache", "run_with_cache"
        ),
      envir = environment()
    )

    parallel::parLapplyLB(
      cl,
      boot_order,
      fun = \(i, ...) {
        run_with_cache(
          i, cache,
          fitfun(boot_indices[, i], networks, ...)
          )
      },
      ...
    )
    
  } else {
    lapply(boot_order, \(i, ...) 
      run_with_cache(
        i, cache,
        fitfun(boot_indices[, i], networks, ...)
      ),
    ...
    )
  }

}
