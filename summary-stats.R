library(ergm.multi)
library(ergm)
library(data.table)

load("LTCF_final_report/res10_new.RData")

# Counting observed sufficient statistics
counts <- summary_formula(
  res10_new$network ~
    edges + 
    N(~edges, ~log(n)) +
    b1dsp(1:2) +
    b1starmix(2, "wound_care", diff = FALSE) +
    b1starmix(2, "dialysis", diff = FALSE) +
    b1starmix(2, "wheelchair", diff = FALSE) +
    b1starmix(2, "bedridden", diff = FALSE) +
    b1starmix(2, "ventilator", diff = FALSE) +
    b1starmix(2, "mdro", diff = FALSE) +
    b2cov("sqrt_age") +
    b2star(2) +
    b2starmix(2, "CNA", diff = FALSE) +
    b2starmix(2, "Nurse", diff = FALSE) +
    b2starmix(2, "Other", diff = FALSE) +
    b2starmix(2, "PTOTRT", diff = FALSE) 
  )

n <- uncombine_network(res10_new$network) |>
  length()

mytab <- data.table(
  Statistic = paste("\\hspace{.5cm}", names(counts)),
  Value     = as.integer(counts),
  Avg       = counts/n
)

mytab[, Statistic := gsub("b.starmix.{3}", "", Statistic)]

# Creating the latex table
mytab <- xtable::xtable(mytab)

ans <- list(
  mytab[1:5,],
  mytab[6:11,],
  mytab[11:nrow(mytab),]
  )

lapply(ans, print, 
  include.colnames = FALSE,
  include.rownames = FALSE,
  only.contents = TRUE,
  hline.after = NULL,
  sanitize.text.function	= \(x) x
  ) |> invisible()



