
test_that("Random walk priors",{
  skip_on_cran()
  expect_no_error({
    mod <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                      smooth_model = "exchangeable")
    hm <- hazard(mod, niter=100) |> mutate(model="exchangeable")
    modr <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                       smooth_model = "random_walk")
    hmr <- hazard(modr, niter=1000) |> mutate(model="random_walk")
    haz <- rbind(hm, hmr)
    ggplot(haz,aes(x=t,col=model)) +
      geom_ribbon(aes(ymin=lower,ymax=upper,fill=model), alpha=0.2) +
      geom_line(aes(y=median))
    plot_survival(modr)

    ## Compare different gamma priors
    modr <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                       smooth_model = "random_walk")
    hmr <- hazard(modr, niter=1000) |> mutate(model="gamma")
    modr5 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                        smooth_model = "random_walk", prior_hsd = p_gamma(2,5))
    hmr5 <- hazard(modr5, niter=1000) |> mutate(model="gamma5")
    modr5 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                        smooth_model = "random_walk", prior_hsd = p_gamma(2,0.1))
    hmr01 <- hazard(modr5, niter=1000) |> mutate(model="gamma01")
    haz <- rbind(hmr, hmr5)
    ggplot(haz,aes(x=t,col=model)) +
      geom_ribbon(aes(ymin=lower,ymax=upper,fill=model), alpha=0.2) +
      geom_line(aes(y=median)) + ylim(0,1)

    ## Non-proportional hazards models
    rxnph_mod <- survextrap(Surv(years, status) ~ rx, data=colons,
                            nonprop=TRUE, fit_method = "opt",
                            smooth_model="exchangeable")
    rxnphr_mod <- survextrap(Surv(years, status) ~ rx, data=colons,
                            nonprop=TRUE, fit_method = "opt",
                            smooth_model="random_walk")
    nd <- data.frame(rx = c("Lev+5FU","Lev"))
    hr1 <- hazard_ratio(rxnph_mod,newdata=nd, niter=1000) |> mutate(model="exchangeable")
    hr2 <- hazard_ratio(rxnphr_mod,newdata=nd,niter=1000) |> mutate(model="random_walk")
    hrs <- rbind(hr1, hr2)
    ggplot(hrs,aes(x=t,col=model)) +
      geom_ribbon(aes(ymin=lower,ymax=upper,fill=model), alpha=0.2) +
      geom_line(aes(y=median)) + coord_cartesian(ylim=c(0,10))

  })
})
