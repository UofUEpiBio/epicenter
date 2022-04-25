library(data.table)

dat <- fread("data/Matched Interaction and Resident.csv")

# Trimming all character
classes <- sapply(dat, class)

# Looking at multiclass (mostly dates)
classes[which(sapply(classes, length) > 1)]

# HCW data
cnames <- colnames(dat)
cnames[grepl("hcw", cnames, ignore.case = TRUE)]

# Function to collapse multiple fields into a single string
# keeping only alphanumeric characters.
namemilize <- function(...) {
  do.call(
    paste,
    c(lapply(list(...), \(x) gsub("[^[:alnum:]]", "_", trimws(x))),
      list(sep = "-"))
  ) |>
    tolower() |>
    gsub(pattern = "_+", replacement = "_")
  
}

# Getting nodes
edges <- dat[, .(state.x, date, Fac_Name, Fac_Unit, hcw_id, hcw_type, res_id)]
dat[, loc_id := namemilize(state, date, Fac_Name, Fac_Unit)]

edges[, loc_id := namemilize(state, date, Fac_Name, Fac_Unit)]
