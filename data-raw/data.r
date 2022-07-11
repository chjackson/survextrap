library(dplyr)
library(survival)
set.seed(1)
colons <- colon |>
    group_by(rx) |>
    slice_sample(prop=0.2) |>
    mutate(years = time/365.25) |>
    filter(etype == 1) |>
    mutate(status = ifelse(years > 3, 0, status),
           years = ifelse(years > 3, 3, years))
use_data(colons, overwrite=TRUE)
