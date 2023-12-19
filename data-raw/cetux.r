library(tidyverse)

cetux <- read.table("data-raw/ipd_recon.dat",header = TRUE) %>%
  mutate(treat = fct_recode(factor(treat), "Control"="0", "Cetuximab"="1"),
         years = t / 12) %>%
  rename(months = t) %>%
  select(months, years, d, treat)

agemed <- 57
agesd <- (83 - 34)/6 # interpret range as 6 SDs. not used
p_male <- 0.8

# Use mortality rate for median age at randomisation + years of follow up
# averaged over gender balance at rand
# current cetux_bh obtained with data from
# https://mortality.org/Country/Country?cntr=USA downloaded in early 2023
cetux_bh <- read.table("data-raw/mort_usa.txt",header=TRUE,skip=2) %>%
  filter(Year==1999) %>%
  mutate(Age = as.numeric(ifelse(Age=="110+", "110", Age)),
         Agenext = ifelse(Age==110, Inf, Age+1),
         Trialyear = Age - agemed,
         hazard = p_male*Male + (1 - p_male)*Female) %>%
  filter(Trialyear >= 0) %>%
  select(time=Trialyear, hazard)
use_data(cetux_bh, overwrite=TRUE)

cetux_seer <- read.table("data-raw/seer.dat",skip=10,
                         col.names=c("t","num","denom"),nrows = 21) %>%
  mutate(start = t-1) %>%
  select(start,stop=t,r=num,n=denom) %>%
  mutate(treat = factor("Control"),
         haz = -log(r/n),
         haz_upper = -log(qbeta(0.025, r, n-r)),
         haz_lower = -log(qbeta(0.975, r, n-r)))

use_data(cetux, overwrite=TRUE)
use_data(cetux_seer, overwrite=TRUE)
