library(npRmpi)

make_proper_cdist_test_bws <- function() {
  list(
    xbw = 0.25,
    ybw = 0.25,
    xnames = "x",
    ynames = "y",
    xndim = 1L,
    yndim = 1L,
    xnord = 0L,
    xnuno = 0L,
    xncon = 1L,
    ynord = 0L,
    ynuno = 0L,
    yncon = 1L,
    pscaling = "Bandwidth(s)",
    ptype = "Fixed",
    pcxkertype = "Gaussian",
    puxkertype = "Aitchison-Aitken",
    poxkertype = "Wang-Van Ryzin",
    pcykertype = "Gaussian",
    puykertype = "Aitchison-Aitken",
    poykertype = "Wang-Van Ryzin",
    regtype = "lp",
    degree = 3L,
    regtype.engine = "lp",
    degree.engine = 3L
  )
}

make_proper_cdist_test_object <- function(condist,
                                          x.grid = c(0.2, 0.8),
                                          y.grid = seq(-1, 1, length.out = length(condist) / length(x.grid)),
                                          trainiseval = FALSE,
                                          gradients = FALSE) {
  stopifnot(length(condist) == length(x.grid) * length(y.grid))
  eval.grid <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))
  getFromNamespace("condistribution", "npRmpi")(
    bws = make_proper_cdist_test_bws(),
    xeval = eval.grid["x"],
    yeval = eval.grid["y"],
    condist = condist,
    conderr = rep(NA_real_, nrow(eval.grid)),
    ntrain = 10L,
    trainiseval = trainiseval,
    gradients = gradients
  )
}

test_that("proper helper projection enforces monotone bounded slices", {
  p.fun <- getFromNamespace(".np_condist_project_bounded_isotonic", "npRmpi")
  w.fun <- getFromNamespace(".np_condens_trapezoid_weights", "npRmpi")

  w <- w.fun(c(0, 1, 2, 3))
  f <- c(-0.2, 0.7, 0.4, 1.2)
  g <- p.fun(f = f, w = w, tol = 1e-12)

  expect_true(all(g >= -1e-10))
  expect_true(all(g <= 1 + 1e-10))
  expect_true(all(diff(g) >= -1e-10))
})

test_that("proper finalizer repairs supported synthetic grids and preserves raw values", {
  raw <- c(-0.3, 0.7, 0.4, -0.1, 0.8, 1.2)
  obj <- make_proper_cdist_test_object(condist = raw)

  out <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list()
  )

  expect_true(isTRUE(out$proper.requested))
  expect_true(isTRUE(out$proper.applied))
  expect_identical(out$proper.method, "isotonic")
  expect_equal(out$condist.raw, raw, tolerance = 1e-12)
  expect_true(all(out$condist >= -1e-10))
  expect_true(all(out$condist <= 1 + 1e-10))

  split.idx <- split(seq_along(out$condist), out$xeval[[1L]])
  for (idx in split.idx)
    expect_true(all(diff(out$condist[idx]) >= -1e-10))
})

test_that("proper finalizer is a no-op for already-proper synthetic lc cdf objects", {
  raw <- c(-0.3, 0.7, 0.4, -0.1, 0.8, 1.2)
  obj <- make_proper_cdist_test_object(condist = raw)
  obj$bws$regtype <- "lc"
  obj$bws$regtype.engine <- "lc"
  obj$bws$degree <- NULL
  obj$bws$degree.engine <- 0L

  out <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list()
  )

  expect_true(isTRUE(out$proper.requested))
  expect_false(isTRUE(out$proper.applied))
  expect_identical(out$proper.info$reason, "already_proper")
  expect_equal(out$condist, raw, tolerance = 1e-12)
  expect_null(out$condist.raw)
})

test_that("proper finalizer is a no-op for already-proper synthetic lp-degree-zero cdf objects", {
  raw <- c(-0.3, 0.7, 0.4, -0.1, 0.8, 1.2)
  obj <- make_proper_cdist_test_object(condist = raw)
  obj$bws$regtype <- "lp"
  obj$bws$regtype.engine <- "lp"
  obj$bws$degree <- 0L
  obj$bws$degree.engine <- 0L

  out <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list()
  )

  expect_true(isTRUE(out$proper.requested))
  expect_false(isTRUE(out$proper.applied))
  expect_identical(out$proper.info$reason, "already_proper")
  expect_equal(out$condist, raw, tolerance = 1e-12)
  expect_null(out$condist.raw)
})

test_that("proper plan projects bootstrap-style distribution matrices slice-wise", {
  raw <- c(-0.3, 0.7, 0.4, -0.1, 0.8, 1.2)
  obj <- make_proper_cdist_test_object(condist = raw)

  plan <- getFromNamespace(".np_condist_prepare_proper_plan", "npRmpi")(
    object = obj,
    proper.control = list()
  )
  expect_true(isTRUE(plan$supported))

  projected <- getFromNamespace(".np_condist_project_values_with_plan", "npRmpi")(
    values = rbind(raw, raw + 0.05),
    plan = plan
  )

  expect_equal(dim(projected), c(2L, length(raw)))
  expect_true(all(projected >= -1e-10))
  expect_true(all(projected <= 1 + 1e-10))

  split.idx <- split(seq_len(length(raw)), obj$xeval[[1L]])
  for (row in seq_len(nrow(projected))) {
    for (idx in split.idx)
      expect_true(all(diff(projected[row, idx]) >= -1e-10))
  }
})

test_that("mode='slice' still defers to exact-grid cdf repair on supported synthetic grids", {
  raw <- c(-0.3, 0.7, 0.4, -0.1, 0.8, 1.2)
  obj <- make_proper_cdist_test_object(condist = raw)

  out.grid <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list()
  )
  out.slice <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list(mode = "slice")
  )

  expect_true(isTRUE(out.slice$proper.applied))
  expect_equal(out.slice$condist, out.grid$condist, tolerance = 1e-12)
})

test_that("proper finalizer stores request-only metadata on unsupported geometry", {
  obj.train <- make_proper_cdist_test_object(condist = rep(0.2, 6L), trainiseval = TRUE)
  out.train <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj.train,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list()
  )
  expect_true(isTRUE(out.train$proper.requested))
  expect_false(isTRUE(out.train$proper.applied))
  expect_identical(out.train$proper.info$reason, "no_eval_grid")

  obj.grad <- make_proper_cdist_test_object(condist = rep(0.2, 6L), gradients = TRUE)
  out.grad <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj.grad,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list()
  )
  expect_true(isTRUE(out.grad$proper.requested))
  expect_false(isTRUE(out.grad$proper.applied))
  expect_identical(out.grad$proper.info$reason, "gradients_unsupported")
})

test_that("apply='fitted' leaves explicit-evaluation synthetic distribution objects unchanged", {
  raw <- c(-0.3, 0.7, 0.4, -0.1, 0.8, 1.2)
  obj <- make_proper_cdist_test_object(condist = raw)

  out <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list(mode = "slice", apply = "fitted"),
    slice.context = list(
      txdat = data.frame(x = c(0.2, 0.8)),
      tydat = data.frame(y = seq(-1, 1, length.out = 3L)),
      exdat = obj$xeval,
      eydat = obj$yeval
    )
  )

  expect_true(isTRUE(out$proper.requested))
  expect_false(isTRUE(out$proper.applied))
  expect_identical(out$proper.info$reason, "scope_not_selected")
  expect_true(isTRUE(out$proper.info$supported))
  expect_equal(out$condist, raw, tolerance = 1e-12)
})

test_that("apply='fitted' also leaves exact-grid synthetic distribution objects unchanged", {
  raw <- c(-0.3, 0.7, 0.4, -0.1, 0.8, 1.2)
  obj <- make_proper_cdist_test_object(condist = raw)

  out <- getFromNamespace(".np_condist_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list(mode = "slice", apply = "fitted"),
    slice.context = list(
      txdat = data.frame(x = c(0.2, 0.8)),
      tydat = data.frame(y = seq(-1, 1, length.out = 3L)),
      exdat = obj$xeval,
      eydat = obj$yeval
    )
  )

  expect_true(isTRUE(out$proper.requested))
  expect_false(isTRUE(out$proper.applied))
  expect_identical(out$proper.info$reason, "scope_not_selected")
  expect_true(isTRUE(out$proper.info$supported))
  expect_equal(out$condist, raw, tolerance = 1e-12)
})
