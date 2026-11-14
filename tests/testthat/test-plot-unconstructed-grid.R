test_that("grid controls cannot silently alter stored copula or conmode samples", {
  if (identical(environmentName(environment(.np_plot_reject_unconstructed_grid)), "npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  withr::local_options(np.messages = FALSE)
  withr::local_seed(75)
  d <- data.frame(x = seq(.02, .98, length.out = 24),
                  z = cos(seq(.1, 4, length.out = 24)),
                  g = factor(rep(c("a", "b"), 12)))
  cm <- npconmode(g ~ x, data = d, bws = c(.2, .3),
                 bandwidth.compute = FALSE, probabilities = TRUE)
  co <- npcopula(~x+z, data = d, bws = c(.3, .4), neval = 4L)
  snapshot <- function(x) serialize(x, NULL, refhook = function(e) format(e))
  before <- snapshot(list(cm, co))
  seed <- .Random.seed
  controls <- list(list(xtrim = c(.1, .9)), list(xq = .3),
                   list(grid.control = np_grid_control(xtrim = c(.1, .9))),
                   list(grid.control = np_grid_control(xq = .3)))
  copula <- plot(co, output = "data")
  sample <- plot(cm, view = "sample", output = "data")
  for (name in c("xtrim", "xq")) {
    a <- setNames(list(NULL), name)
    expect_identical(do.call(plot, c(list(co, output="data"), a))$panels, copula$panels)
    expect_identical(do.call(plot, c(list(cm, view="sample", output="data"), a))$panels, sample$panels)
  }
  expect_identical(plot(cm, view="sample", neval=NULL, output="data")$panels, sample$panels)
  for (a in controls) {
    expect_error(do.call(plot, c(list(x = co, output = "data"), a)), "does not construct a grid")
    expect_error(do.call(plot, c(list(x = cm, view = "sample", output = "data"), a)), "does not construct a grid")
    expect_error(do.call(plot, c(list(x = copula), a)), "does not construct a grid")
    expect_error(do.call(plot, c(list(x = sample), a)), "does not construct a grid")
  }
  expect_error(plot(cm, view = "sample", neval = 4, output = "data"), "does not construct a grid")
  expect_error(plot(co, grid.control = list(xq = .3), output = "data"), "np_grid_control")
  expect_error(plot(co, grid.control = np_grid_control(slices = 1), output = "data"), "does not construct a grid")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add=TRUE)
  for (saved in list(copula, sample)) {
    expect_error(plot(saved, neval=5), "does not construct a grid")
    suppressWarnings(expect_no_error(plot(saved, xtrim=NULL, xq=NULL, neval=NULL)))
  }
  # This repair deliberately leaves other replay families' existing policy alone.
  reg <- npreg(z~x, data=d, bws=.3, bandwidth.compute=FALSE)
  saved.reg <- plot(reg, neval=4, output="data")
  suppressWarnings(expect_no_error(plot(saved.reg, xq=.3)))
  suppressWarnings(expect_no_error(plot(saved.reg, neval=5)))
  for (ctrl in list(NULL, np_grid_control())) {
    expect_identical(plot(co, grid.control = ctrl, output = "data")$panels, copula$panels)
    if (!is.null(ctrl))
      expect_identical(plot(cm, view = "sample", grid.control = ctrl, output = "data")$panels, sample$panels)
    else
      expect_error(plot(cm, view = "sample", grid.control = ctrl, output = "data"), "np_grid_control")
  }
  fixed <- plot(cm, view = "fixed", neval = 7, output = "data")
  trimmed <- plot(cm, view = "fixed", neval = 7,
                  grid.control = np_grid_control(xtrim = c(.1, .9)), output = "data")
  expect_false(identical(fixed$panels[[1]]$grid, trimmed$panels[[1]]$grid))
  expect_equal(nrow(trimmed$panels[[1]]$grid), 7L)
  expect_identical(snapshot(list(cm, co)), before)
  expect_identical(.Random.seed, seed)
})
