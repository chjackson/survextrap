library(survextrap)
library(DT)

mspline_defaults <- list(
  df = 10,
  degree = 3,
  bsmooth = TRUE,
  bknot = 10
)

selection_defaults <- list(
  degree_unsmooth = 2:4,
  df = 5:12
)
