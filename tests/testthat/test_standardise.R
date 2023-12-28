
test_that("Basic spline model, standardised output",{
    mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
    ref_pop <- data.frame(rx = c("Obs","Lev+5FU"))
    nostd <- survival(mod, t=c(5,10), newdata=ref_pop)
    std <- survival(mod, t = c(5,10), newdata = standardise_to(ref_pop))
    expect_lt(nostd$median[nostd$t==5 & nostd$rx=="Obs"],
              std$median[std$t==5])
    expect_lt(std$median[std$t==5],
              nostd$median[nostd$t==5 & nostd$rx=="Lev+5FU"])

    set.seed(1)
    stdr <- survival(mod, t = c(5,10), newdata = standardise_to(ref_pop, random=TRUE))
    expect_lt(nostd$median[nostd$t==5 & nostd$rx=="Obs"],
              stdr$median[stdr$t==5])
    expect_lt(stdr$median[stdr$t==5],
              nostd$median[nostd$t==5 & nostd$rx=="Lev+5FU"])
})
