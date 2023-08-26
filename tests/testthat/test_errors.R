test_that("Errors in formula",{
    expect_error(survextrap(formula="boo"), "must be a formula")
    expect_error(survextrap(Surv(years, status, type="left") ~ rx, data=colons), "Cannot handle \'left\'")
})

test_that("Errors in external data",{
    external <- data.frame(x=1:2)
    expect_error(survextrap(Surv(years, status) ~ rx, data=colons, external=external),
                 "not found in `external`")
    external <- data.frame(start=10, stop=15)
    expect_error(survextrap(Surv(years, status) ~ rx, data=colons, external=external),
                 "not found in `external`")
    external <- data.frame(start=10, stop=15, n=20, r=5)
    expect_error(survextrap(Surv(years, status) ~ rx, data=colons, external=external),
                 "Covariate \"rx\" not found")
    expect_error(survextrap(Surv(years, status) ~ 1, data=colons, external=external, cure=~rx),
                 "Covariate \"rx\" not found")

    colonse <- colons
    colonse$bh <- rep(0.01, nrow(colons))
    ext <- data.frame(start=5, stop=10, n=30, r=5)
    expect_error(
        survextrap(Surv(years, status) ~ 1, data=colonse, external=ext, backhaz="bh"),
        "Columns \"backsurv_start\"")

})
