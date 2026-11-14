test_that("equivalent trim pairs share the existing pair-owner contract", {
  withr::local_options(np.messages = FALSE, np.tree = FALSE)

  norm <- get(".np_plot_normalize_public_dots", asNamespace("np"))
  pair <- c(.1, .9)
  for(context in c("plot.singleindex", "plot.sibandwidth", "plot.conmode")) {
    both <- list(xtrim = pair, grid.control = np_grid_control(xtrim = pair))
    single <- list(grid.control = np_grid_control(xtrim = pair))
    expect_identical(norm(both, context), norm(single, context))
    both$xtrim <- c(.2, .8)
    expect_error(norm(both, context), "conflicting plot arguments", fixed = TRUE)
  }
  for(context in c("plot.npdensity", "plot.npregression", "plot.condistribution"))
    expect_error(norm(list(xtrim = pair, grid.control = np_grid_control(xtrim = pair)),
      context), "conflicting plot arguments", fixed = TRUE)
  d <- data.frame(x = seq(.02,.98,length.out=24), z = cos(seq(.1,4,length.out=24)))
  d$y <- sin(3*d$x)
  fit <- npindex(y ~ x + z, data = d, bws = c(1,.5,.3), bandwidth.compute = FALSE)
  extract <- function(x) lapply(x$panels, function(p) p[c("grid","values","intervals")])
  direct <- plot(fit, neval = 7L, output = "data", xtrim = pair)
  both <- plot(fit, neval = 7L, output = "data", xtrim = pair,
    grid.control = np_grid_control(xtrim = pair))
  expect_identical(extract(direct), extract(both))
})
