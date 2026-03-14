fake_unconditional_density_bws <- function(ndim,
                                           ncon = ndim,
                                           nord = 0L,
                                           nuno = 0L,
                                           iord = rep(FALSE, ndim),
                                           inumord = rep(FALSE, ndim),
                                           xnames = paste0("x", seq_len(ndim))) {
  structure(
    list(
      ndim = as.integer(ndim),
      ncon = as.integer(ncon),
      nord = as.integer(nord),
      nuno = as.integer(nuno),
      xnames = xnames,
      xdati = list(
        iord = as.logical(iord),
        inumord = as.logical(inumord)
      )
    ),
    class = "bandwidth"
  )
}

test_that("unconditional density bootstrap planner emits one slice target per dimension", {
  planner <- getFromNamespace(".np_plot_bootstrap_plan_unconditional_density", "np")

  xdat <- data.frame(x1 = runif(5), x2 = runif(5), x3 = runif(5))
  bws <- fake_unconditional_density_bws(ndim = 3L)

  plan <- planner(
    bws = bws,
    xdat = xdat,
    neval = 25L,
    perspective = FALSE,
    plot.errors.boot.method = "inid",
    plot.errors.type = "simultaneous"
  )

  expect_identical(plan$family, "unconditional_density")
  expect_identical(plan$method_label, "inid")
  expect_identical(plan$band_type, "simultaneous")
  expect_identical(plan$target_total, 3L)
  expect_identical(vapply(plan$targets, `[[`, character(1L), "kind"), rep("slice", 3L))
  expect_identical(vapply(plan$targets, `[[`, character(1L), "label"), c("x1", "x2", "x3"))
  expect_identical(vapply(plan$targets, `[[`, integer(1L), "eval_size"), rep(25L, 3L))
  expect_identical(unname(plan$phase_labels[["constructing"]]), "constructing simultaneous bands")
})

test_that("unconditional density bootstrap planner collapses eligible 2D perspective work to one surface target", {
  planner <- getFromNamespace(".np_plot_bootstrap_plan_unconditional_density", "np")

  xdat <- data.frame(x1 = runif(6), x2 = runif(6))
  bws <- fake_unconditional_density_bws(ndim = 2L, ncon = 2L, nord = 0L, nuno = 0L)

  plan <- planner(
    bws = bws,
    xdat = xdat,
    neval = 40L,
    perspective = TRUE,
    plot.errors.boot.method = "geom",
    plot.errors.type = "all"
  )

  expect_identical(plan$target_total, 1L)
  expect_identical(plan$targets[[1L]]$kind, "surface")
  expect_identical(plan$targets[[1L]]$label, "x1,x2 surface")
  expect_identical(plan$targets[[1L]]$eval_size, 1600L)
  expect_identical(unname(plan$phase_labels[["constructing"]]), "constructing all bands")
})

test_that("unconditional density bootstrap planner uses factor cardinality for slice eval sizes", {
  planner <- getFromNamespace(".np_plot_bootstrap_plan_unconditional_density", "np")

  xdat <- data.frame(
    x1 = factor(c("a", "b", "c", "a")),
    x2 = runif(4)
  )
  bws <- fake_unconditional_density_bws(ndim = 2L)

  plan <- planner(
    bws = bws,
    xdat = xdat,
    neval = 30L,
    perspective = FALSE,
    plot.errors.boot.method = "fixed",
    plot.errors.type = "pointwise"
  )

  expect_identical(vapply(plan$targets, `[[`, integer(1L), "eval_size"), c(3L, 30L))
  expect_identical(vapply(plan$targets, `[[`, character(1L), "label"), c("x1", "x2"))
})

test_that("unconditional density bootstrap planner stays agnostic to execution-path details", {
  planner <- getFromNamespace(".np_plot_bootstrap_plan_unconditional_density", "np")

  planner.formals <- names(formals(planner))

  expect_false("plot.errors.boot.nonfixed" %in% planner.formals)
  expect_false("counts.drawer" %in% planner.formals)
  expect_false("H" %in% planner.formals)
  expect_false("helper" %in% planner.formals)
  expect_false("use.frozen.nonfixed" %in% planner.formals)
})
