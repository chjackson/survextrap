mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
nd <- data.frame(rx = c("Lev+5FU","Lev"))

test_that("hazard ratio plot", {
  expect_error({
    plot_hazard_ratio(mod, newdata=nd, niter=20, t=1:5)
  }, NA)
})

test_that("survival and hazard plots", {
  expect_error({
    p <- plot_survival(mod, niter=20)
    p <- plot_survival(mod, niter=20, km=FALSE)
    p <- plot_hazard(mod, niter=20)
    p <- plot(mod, niter=20)
  }, NA)
})

test_that("deconstruct fitted spline",{
  expect_error({
    hazdf <-  deconstruct_mspline(mod)
    knots <- unname(attr(hazdf,"knots"))
    p <- ggplot(hazdf, aes(x=time, y=value, group=term)) +
      geom_line(alpha=0.5) +
      scale_x_continuous(breaks=c(knots)) +
      theme_minimal() +
      theme(panel.grid.minor = element_blank()) +
      geom_vline(xintercept = max(knots), col="gray50") +
      xlab("Time") + ylab("Hazard")
  }, NA)
})
