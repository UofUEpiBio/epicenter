# convert NA to 0
naton <- function(x){y=x;y[is.na(x)]=0; return(y)} 

# count unique values
nuni <- function(x){length(unique(x))}

# count missings in all columns
findmiss <- function(dat){M=apply(dat,2,function(x){sum(is.na(x))}); M[M>0]}

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
  boot = FALSE
  ) {

  # Re-building the dat to be sure it is right
  dat <- lapply(nets, \(net) {
    data.table(
      NetID = net %v% "net_id",
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
      diabetes = net %v% "diabetes"
    )
  }) |> rbindlist()

  dat[, NetID := as.integer(NetID)]

  idx <- seq_along(nets)
  if (boot) {
    # Sampling the networks
    idx <- sample(idx, replace = TRUE) |> sort()
    nets <- nets[idx]

    # Redoing indices
    dat[, id. := 1:.N, by = "NetID"]

    # Merging in
    dat_idx <- data.table(
      NetID = idx,
      NetID_new = 1:length(idx)
    )

    dat_new <- merge(
      x = dat_idx,
      y = dat,
      by = "NetID",
      all.x = TRUE,
      all.y = FALSE
    )

    setorder(dat_new, "NetID", "id.")

    # Saving the new position id
    dat_new[, NetID := NetID_new]
    dat_new[, NetID_new := NULL]
    dat_new[, id. := NULL]
    
    # Saving
    dat <- dat_new
  }

  ans <- list()

  ans$Nurse  = dat[, n:=sum(Nurse==1), NetID][n>1, unique(NetID)]
  ans$CNA    = dat[, n:=sum(CNA==1), NetID][n>1, unique(NetID)]
  ans$PTOTRT = dat[, n:=sum(PTOTRT==1), NetID][n>1, unique(NetID)]
  ans$Other  = dat[, n:=sum(Other==1), NetID][n>1, unique(NetID)] 

  # For b2factor (resident type), the networks need to have at least 1 resident of that type. 
  ans$dialysis   <- dat[, n:=sum(dialysis=='Yes'), NetID][n>0, unique(NetID)]
  ans$wound_care <- dat[, n:=sum(wound_care=='Yes'), NetID][n>0, unique(NetID)] 
  ans$mdro       <- dat[, n:=sum(mdro=='Yes'), NetID][n>0, unique(NetID)]
  ans$bedridden  <- dat[, n:=sum(bedridden=='Yes'), NetID][n>0, unique(NetID)]  

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
      n_hcp = n_hcp,
      n_patients = n_patients,
      n_tot = n_hcp + n_patients,
      facility_state = sapply(nets, \(net) unique(net %v% "state")) |>
        as.factor()
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

  tmethod(
    x,
    custom.coef.map = maps,
    ...,
    custom.gof.rows = nnets
    )

}
