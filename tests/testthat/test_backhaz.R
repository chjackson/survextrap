
test_that("get_cum_backhaz",{
  backhaz <- data.frame(hazard = c(0.01, 0.02, 0.03), time = c(0, 5, 10))
  expect_equal(get_cum_backhaz(2, backhaz), 0.02)
  expect_equal(get_cum_backhaz(6, backhaz), 0.07)
  expect_equal(get_cum_backhaz(12, backhaz), 0.21)
  expect_equal(get_cum_backhaz(c(2,6,12), backhaz), c(0.02, 0.07, 0.21))
  expect_error(get_cum_backhaz(2, backhaz, strata="Male"),
               "`strata` provided, but no column named `stratum` in `backhaz`")
  backhaz <- data.frame(hazard = c(0.01, 0.02, 0.03), time = c(0, 5, 10),
                        stratum = c("Male","Female","Male"))
  # no requirement to sort backhaz within strata
  expect_equal(get_cum_backhaz(6, backhaz, strata = "Male"), 0.06)
  expect_error(get_cum_backhaz(6, backhaz, strata = "Foo"),
               "No stratum with value")
})


alpha <- -0.35
coefs <- c(0.0012, 0.029, 0.1, 0.11, 0.23, 0.0074, 0.052, 0.1)
knots <- c(0.36, 0.62, 0.91, 1.2, 1.54, 2.02, 3)

test_that("distribution functions with backhaz",{
  backhaz <- data.frame(hazard = c(0.01, 0.02, 0.03), time = c(0, 5, 10))
  h1 <-  hsurvmspline(5, alpha, coefs, knots)
  h2 <- hsurvmspline(5, alpha, coefs, knots, backhaz=backhaz)
  expect_lt(h1, h2)
})

### Set up data for testing stratified background hazards

## Create age and sex strata in colon data
colons$sex <- factor(colons$sex, labels=c("Female","Male"))
colons$agegroup <- cut(colons$age, breaks=c(0,70,Inf),
                       right=FALSE, labels=c("Under 70","Over 70"))

## Simulate stratified background hazard data with same risk for each stratum
## so can check if get same results as unstratified analysis
bh <- data.frame(hazard = c(0.01, 0.02, 0.03), time=c(0, 5, 10))
age4 <- rep(c("Under 70","Over 70"), 2)
sex4 <- rep(c("Male","Male","Female","Female"))
hrnull <- rep(1, 4)
bh_strata <-
  bh[rep(1:nrow(bh), 4),] %>%
  mutate(agegroup = factor(rep(age4, each=nrow(bh))),
         sex = factor(rep(sex4, each=nrow(bh))),
         hazard = hazard*rep(hrnull, each=nrow(bh)),
         stratum = paste(agegroup, sex, sep=";"))
rownames(bh_strata) <- NULL

## unstratified
set.seed(1)   #why does posterior mode depend on the seed??
mod1 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                   backhaz=bh)
## stratified but with no strata effect: expect same results
set.seed(1)
mod1st <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                     backhaz=bh_strata, backhaz_strata=c("agegroup","sex"))

bh_stratad <- bh_strata
bh_stratad$hazard <- rep(c(0.8, 1.3, 0.8, 1.2), each=nrow(bh)) * bh_strata$hazard
## stratified, with no strata effect: expect different results
set.seed(1)
mod1std <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      backhaz=bh_stratad, backhaz_strata=c("agegroup","sex"))

test_that("Stratified background hazards: just IPD",{
  expect_equal(
    get_cum_backhaz(c(2,5,8), bh_strata,
                    strata=c("Under 70;Male","Under 70;Female","Over 70;Male")),
    get_cum_backhaz(c(2,5,8), bh))
  expect_equal(summary(mod1)$mode[1], summary(mod1st)$mode[1])
  expect_true(summary(mod1)$mode[1] != summary(mod1std)$mode[1])

  nd <- data.frame(agegroup="Over 70", sex="Female")
  nd2 <- data.frame(agegroup="Under 70", sex="Female")

  set.seed(1); r1 <- rmst(mod1, newdata=nd, t=5, niter=50)
  set.seed(1); r1st <- rmst(mod1st, newdata=nd, t=5, niter=50)
  set.seed(1); r1std <- rmst(mod1std, newdata=nd, t=5, niter=50)
  expect_equal(r1, r1st)
  set.seed(1); s1 <- survival(mod1, newdata=nd, t=5, niter=50)

  ## TODO dim of coefs.  its 4d array on exit from prepars
  ## we at least need it to be subsetted by nvars in rmst
  ## TODO document vectorisation of dpqr.  coefs ends up as a matrix with
  ## nvars cols,  so nvar should be the final dimension
  ## what was it before: a matrix with nvars cols
  ## suspect it's going wrong in _dist_setup
  ##    coefs <- matrix(rep(as.numeric(t(coefs)), length.out = ncol(coefs) * nret),
 ##             ncol = ncol(coefs), byrow = TRUE)
  ## input is nt x niter x nvars x nvals
  ## output is niter x (niter x nvals),  weird
  ## it doesn't handle 3+d arrays
  ## (a) if supply vector, interpret as 1-row matrix
  ## (b) if supply matrix, why is this line done?  It's doing the vectorisation.
  ## Output should end up with nret rows, nvals cols
  ## eg if coefs has 2 rows and alpha is length 8, rep 4 times
  ## combining t() and byrow() is opaque.   Better with an index

  set.seed(1); s1st <- survival(mod1st, newdata=nd, t=5, niter=50)
  set.seed(1); s1std <- survival(mod1std, newdata=nd, t=5, niter=50)
  expect_equal(s1, s1st)
  set.seed(1); h1 <- hazard(mod1, newdata=nd, t=5, niter=50)
  set.seed(1); h1st <- hazard(mod1st, newdata=nd, t=5, niter=50)
  set.seed(1); h1std <- hazard(mod1std, newdata=nd, t=5, niter=50)
  expect_equal(h1, h1st)

  set.seed(1);  r1st <- rmst(mod1std, newdata=nd, t=5, niter=50)
  set.seed(1);  r1st2 <- rmst(mod1std, newdata=nd2, t=5, niter=50)
  expect_gt(r1st2$median, r1st$median)

  set.seed(1); s1std <- survival(mod1std, newdata=nd, t=5, niter=50)
  set.seed(1); s1std2 <- survival(mod1std, newdata=nd2, t=5, niter=50)
  expect_gt(s1std2$median, s1std$median)

  set.seed(1); s1std <- survival(mod1std, newdata=nd, t=1:3, niter=50)

  ## default strata in newdata
  set.seed(1); s1 <- survival(mod1std, t=1:3, niter=50)
  set.seed(1); s2 <- survival(mod1std, t=1:3, niter=50,
                              newdata = list(agegroup="Under 70",sex="Male"))
  expect_equal(s1$median, s2$median)
})

test_that("Stratified background hazards with external data",{
  ext <- data.frame(start=5, stop=10, n=30, r=5,
                    backsurv_start = 0.4, backsurv_stop = 0.3,
                    agegroup = "Under 70", sex = "Female") |>
    mutate(stratum = paste(agegroup, sex, sep=";"))
  set.seed(1)
  mod1 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                     external = ext, backhaz=bh)
  set.seed(1)
  mod1st <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                       external = ext,
                       backhaz=bh_strata, backhaz_strata=c("agegroup","sex"))
  expect_equal(summary(mod1)$mode[1], summary(mod1st)$mode[1])

  set.seed(1)
  mod1st <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                       external = ext,
                       backhaz=bh_stratad, backhaz_strata=c("agegroup","sex"))
  expect_true(summary(mod1)$mode[1] != summary(mod1st)$mode[1])
})

test_that("Stratified background hazards with no individual data",{
  extdat <- data.frame(start = c(5, 10, 15, 20),
                       stop =  c(10, 15, 20, 25),
                       n = c(100, 100, 100, 100),
                       r = c(50, 40, 30, 20),
                       agegroup = rep("Under 70",4),
                       sex = rep("Male",4))
  set.seed(1)
  nde_mod1 <- survextrap(~1, external = extdat, backhaz=bh,
                         mspline = list(df=5), fit_method="opt")
  set.seed(1)
  nde_mod2 <- survextrap(~1, external = extdat,
                         backhaz=bh_strata, backhaz_strata = c("agegroup","sex"),
                         mspline = list(df=5), fit_method="opt")
  expect_equal(summary(nde_mod1)$mode[1], summary(nde_mod2)$mode[1])
  set.seed(1)
  nde_mod3 <- survextrap(~1, external = extdat,
                         backhaz=bh_stratad, backhaz_strata = c("agegroup","sex"),
                         mspline = list(df=5), fit_method="opt")
  expect_true(summary(nde_mod3)$mode[1] != summary(nde_mod2)$mode[1])
  nd <- data.frame(agegroup="Under 70", sex="Male")
  plot_hazard(nde_mod3, niter=100, newdata = nd)
})

test_that("backhaz error handling",{
  bh <- data.frame(hazard = c(0.01, 0.02, 0.03), time=c(1, 5, 10))
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons, backhaz=bh),
               "The first element of the \"time\" column of `backhaz`")

  bh <- data.frame(hazard = c(0.01, 0.02, 0.03, 0.04), time=c(0, 5, 1, 10))
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons, backhaz=bh),
               "backhaz.time should be sorted in increasing order")

  bhe <- bh_strata
  bhe$time[4] <- 1
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons,
                          fit_method = "opt",
                          backhaz=bhe, backhaz_strata=c("agegroup","sex")),
               "The first element of the \"time\" column of `backhaz`  within stratum \"Over 70;Male\" should be 0")

  bhe <- bh_strata
  bhe$time[5:6] <- 6:5
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons,
                          fit_method = "opt",
                          backhaz=bhe, backhaz_strata=c("agegroup","sex")),
               "backhaz.time should be sorted in increasing order within stratum \"Over 70.Male\"")

  expect_error(survextrap(Surv(years, status) ~ 1, data=colons,
                          fit_method = "opt", backhaz=bhe, backhaz_strata=3),
               "`strata` should be a character vector of variables indicating background hazard strata")

  expect_error(survextrap(Surv(years, status) ~ 1, data=colons,
                          fit_method = "opt", backhaz=bhe, backhaz_strata=c("age","sex")),
               "strata variable \"age\" not found in")

  expect_error(survextrap(Surv(years, status) ~ 1, data=colons,
                          fit_method = "opt", backhaz=bhe, backhaz_strata=c("age","boo")),
               "strata variables \"age\",\"boo\" not found in")

  bhe <- bh_strata
  bhe$agegroup <- 1
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons,
                          fit_method = "opt", backhaz=bhe, backhaz_strata=c("agegroup","sex")),
               "stratifying variable \"agegroup\" has class \"factor\" in `data`, but class \"numeric\" in `backhaz`, so cannot match")

  bhe <- bh_strata
  bhe$sex[bhe$sex=="Female"] <- "Male"
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons,
                          fit_method = "opt", backhaz=bhe, backhaz_strata=c("agegroup","sex")),
               "rows 1,3,4 and others in `data` have no `sex` value in `backhaz` matching their value for `sex` in `data`")
})
