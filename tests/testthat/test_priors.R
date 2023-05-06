test_that("Error handling for priors",{
    expect_error(
        survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                   prior_loghr = list(p_normal(1,2))),
        "1 component, but there are 2 coefficients")
    expect_error(survextrap(Surv(t, status) ~ 1, data=curedata, fit_method="opt", cure=~x,
                            prior_hsd = p_normal(0.2, 10)), "should be p_gamma")
    expect_error(
        survextrap(Surv(t, status) ~ 1, data=curedata, fit_method="opt", cure=~x,
                   prior_logor_cure = p_gamma(0.2, 10)), "should be one of")
    expect_error(
        survextrap(Surv(t, status) ~ 1, data=curedata, fit_method="opt", cure=~x,
                   prior_logor_cure = "gamma"), "should be a call")
})

test_that("Priors on hscale behave",{
    nd1 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                      prior_hscale = p_normal(0, 2))
    m1 <- summary(nd1) %>% filter(variable=="alpha") %>% pull(median)
    nd2 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                      prior_hscale = p_normal(4, 0.1))
    m2 <- summary(nd2) %>% filter(variable=="alpha") %>% pull(median)
    expect_gt(m2, m1)
})

test_that("Priors on loghr behave",{
    nd1 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                      prior_loghr = p_normal(1, 2))
    nd2 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                      prior_loghr = p_normal(1, 0.1))
    m1 <- summary(nd1) %>% filter(variable=="loghr", term=="rxLev") %>% pull(median)
    m2 <- summary(nd2) %>% filter(variable=="loghr", term=="rxLev") %>% pull(median)
    expect_gt(m2, m1)
    nd3 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                      prior_loghr = p_normal(2, 0.1))
    nd4 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                      prior_loghr = list("rxLev"=p_normal(2, 0.1),
                                         "rxLev+5FU"=p_normal(0, 2.5)))
    m3 <- summary(nd3) %>% filter(variable=="loghr", term=="rxLev+5FU") %>% pull(median)
    m4 <- summary(nd4) %>% filter(variable=="loghr", term=="rxLev+5FU") %>% pull(median)
    expect_gt(m3, m4)
})

test_that("Priors on cure prob behave",{
    nd1 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt", cure=TRUE,
                      prior_cure = p_beta(10, 10))
    nd2 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt", cure=TRUE,
                      prior_cure = p_beta(1, 10))
    nd3 <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt", cure=TRUE,
                      prior_cure = p_beta(100, 1)) # dominated by data
    m1 <- summary(nd1) %>% filter(variable=="pcure") %>% pull(mode)
    m2 <- summary(nd2) %>% filter(variable=="pcure") %>% pull(mode)
    m3 <- summary(nd3) %>% filter(variable=="pcure") %>% pull(mode)
    expect_gt(m1, m2)
    expect_gt(m3, m1)
})

test_that("Priors on cure prob covariates behave",{
    nd1 <- survextrap(Surv(t, status) ~ 1, data=curedata, fit_method="opt", cure=~x,
                      prior_logor_cure = p_normal(0, 0.1))
    nd2 <- survextrap(Surv(t, status) ~ 1, data=curedata, fit_method="opt", cure=~x,
                      prior_logor_cure = p_normal(4, 0.1))
    m1 <- summary(nd1) %>% filter(variable=="logor_cure") %>% pull(median)
    m2 <- summary(nd2) %>% filter(variable=="logor_cure") %>% pull(median)
    expect_gt(m2, m1)
})

test_that("Priors on smoothing SD behave",{
    nd1 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      prior_hsd = p_gamma(1, 2))
    nd2 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      prior_hsd = p_gamma(0.8, 2))
    m1 <- summary(nd1) %>% filter(variable=="hsd") %>% pull(median)
    m2 <- summary(nd2) %>% filter(variable=="hsd") %>% pull(median)
    expect_lt(m2, m1)
})

test_that("Errors in priors on nonprop",{
  expect_error(survextrap(Surv(years,status)~rx, data=colons, prior_hrsd = p_normal(4, 2), nonprop=TRUE),
               "p_gamma")
  expect_error(survextrap(Surv(years,status)~rx, data=colons, prior_hrsd = list(x=p_gamma(4, 2)), nonprop=TRUE),
               "1 component")
  expect_error(survextrap(Surv(years,status)~rx, data=colons,
                          prior_hrsd = list(xwrong1=p_gamma(4, 2), xwrong2=p_gamma(4, 2)), nonprop=TRUE),
               "names of prior_hrsd do not match")
})

test_that("Priors on non-proportional hazards smoothing SD behave",{
    nd1 <- survextrap(Surv(years, status) ~ rx, data=colons, nonprop=TRUE, fit_method="opt",
                      prior_hrsd = p_gamma(2, 1))
    nd2 <- survextrap(Surv(years, status) ~ rx, data=colons, nonprop=TRUE, fit_method="opt",
                      prior_hrsd = p_gamma(3, 1))
    m1 <- summary(nd1) %>% filter(variable=="hrsd", term=="rxLev") %>% pull(median)
    m2 <- summary(nd2) %>% filter(variable=="hrsd", term=="rxLev") %>% pull(median)
    expect_lt(m1, m2)
})

