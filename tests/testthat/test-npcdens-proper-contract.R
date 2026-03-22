library(npRmpi)

make_proper_test_bws <- function() {
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

make_proper_test_object <- function(condens,
                                    x.grid = c(0.2, 0.8),
                                    y.grid = seq(-1, 1, length.out = length(condens) / length(x.grid)),
                                    trainiseval = FALSE,
                                    gradients = FALSE) {
  stopifnot(length(condens) == length(x.grid) * length(y.grid))
  eval.grid <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))
  getFromNamespace("condensity", "npRmpi")(
    bws = make_proper_test_bws(),
    xeval = eval.grid["x"],
    yeval = eval.grid["y"],
    condens = condens,
    conderr = rep(NA_real_, nrow(eval.grid)),
    ntrain = 10L,
    trainiseval = trainiseval,
    gradients = gradients
  )
}

test_that("proper helper weights and projection satisfy core invariants", {
  w.fun <- getFromNamespace(".np_condens_trapezoid_weights", "npRmpi")
  p.fun <- getFromNamespace(".np_condens_project_weighted_simplex", "npRmpi")

  w <- w.fun(c(0, 1, 2))
  expect_equal(w, c(0.5, 1, 0.5), tolerance = 1e-12)

  f <- c(-0.2, 0.4, 1.8)
  g <- p.fun(f = f, w = w, mass = 1, tol = 1e-12)
  expect_true(all(g >= -1e-10))
  expect_equal(sum(w * g), 1, tolerance = 1e-8)
})

test_that("proper finalizer repairs supported synthetic grids and preserves raw values", {
  raw <- c(-0.3, 0.2, 0.9, -0.1, 0.4, 1.1)
  obj <- make_proper_test_object(condens = raw)

  out <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "project",
    proper.control = list()
  )

  expect_true(isTRUE(out$proper.requested))
  expect_true(isTRUE(out$proper.applied))
  expect_identical(out$proper.method, "project")
  expect_equal(out$condens.raw, raw, tolerance = 1e-12)
  expect_true(all(out$condens >= -1e-10))

  w.fun <- getFromNamespace(".np_condens_trapezoid_weights", "npRmpi")
  w <- w.fun(unique(out$yeval[[1L]]))
  split.idx <- split(seq_len(nrow(out$xeval)), out$xeval[[1L]])
  mass <- vapply(split.idx, function(idx) sum(w * out$condens[idx]), numeric(1))
  expect_equal(unname(mass), rep(1, length(mass)), tolerance = 1e-8)
})

test_that("proper finalizer is a no-op for already-proper synthetic lc objects", {
  raw <- c(-0.3, 0.2, 0.9, -0.1, 0.4, 1.1)
  obj <- make_proper_test_object(condens = raw)
  obj$bws$regtype <- "lc"
  obj$bws$regtype.engine <- "lc"
  obj$bws$degree <- NULL
  obj$bws$degree.engine <- 0L

  out <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "project",
    proper.control = list()
  )

  expect_true(isTRUE(out$proper.requested))
  expect_false(isTRUE(out$proper.applied))
  expect_identical(out$proper.info$reason, "already_proper")
  expect_equal(out$condens, raw, tolerance = 1e-12)
  expect_null(out$condens.raw)
})

test_that("proper finalizer is a no-op for already-proper synthetic lp-degree-zero objects", {
  raw <- c(-0.3, 0.2, 0.9, -0.1, 0.4, 1.1)
  obj <- make_proper_test_object(condens = raw)
  obj$bws$regtype <- "lp"
  obj$bws$regtype.engine <- "lp"
  obj$bws$degree <- 0L
  obj$bws$degree.engine <- 0L

  out <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "project",
    proper.control = list()
  )

  expect_true(isTRUE(out$proper.requested))
  expect_false(isTRUE(out$proper.applied))
  expect_identical(out$proper.info$reason, "already_proper")
  expect_equal(out$condens, raw, tolerance = 1e-12)
  expect_null(out$condens.raw)
})

test_that("proper plan projects bootstrap-style density matrices slice-wise", {
  raw <- c(-0.3, 0.2, 0.9, -0.1, 0.4, 1.1)
  obj <- make_proper_test_object(condens = raw)

  plan <- getFromNamespace(".np_condens_prepare_proper_plan", "npRmpi")(
    object = obj,
    proper.control = list()
  )
  expect_true(isTRUE(plan$supported))

  projected <- getFromNamespace(".np_condens_project_values_with_plan", "npRmpi")(
    values = rbind(raw, raw + 0.05),
    plan = plan
  )

  expect_equal(dim(projected), c(2L, length(raw)))
  expect_true(all(projected >= -1e-10))

  w <- getFromNamespace(".np_condens_trapezoid_weights", "npRmpi")(unique(obj$yeval[[1L]]))
  split.idx <- split(seq_len(length(raw)), obj$xeval[[1L]])
  for (row in seq_len(nrow(projected))) {
    mass <- vapply(split.idx, function(idx) sum(w * projected[row, idx]), numeric(1))
    expect_equal(unname(mass), rep(1, length(mass)), tolerance = 1e-8)
  }
})

test_that("mode='slice' still defers to exact-grid repair on supported synthetic grids", {
  raw <- c(-0.3, 0.2, 0.9, -0.1, 0.4, 1.1)
  obj <- make_proper_test_object(condens = raw)

  out.grid <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "project",
    proper.control = list()
  )
  out.slice <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "project",
    proper.control = list(mode = "slice")
  )

  expect_true(isTRUE(out.slice$proper.applied))
  expect_equal(out.slice$condens, out.grid$condens, tolerance = 1e-12)
})

test_that("proper finalizer stores request-only metadata on unsupported geometry", {
  obj.train <- make_proper_test_object(condens = rep(0.2, 6L), trainiseval = TRUE)
  out.train <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj.train,
    proper = TRUE,
    proper.method = "project",
    proper.control = list()
  )
  expect_true(isTRUE(out.train$proper.requested))
  expect_false(isTRUE(out.train$proper.applied))
  expect_identical(out.train$proper.info$reason, "no_eval_grid")

  obj.grad <- make_proper_test_object(condens = rep(0.2, 6L), gradients = TRUE)
  out.grad <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj.grad,
    proper = TRUE,
    proper.method = "project",
    proper.control = list()
  )
  expect_true(isTRUE(out.grad$proper.requested))
  expect_false(isTRUE(out.grad$proper.applied))
  expect_identical(out.grad$proper.info$reason, "gradients_unsupported")
})

test_that("apply='fitted' leaves explicit-evaluation synthetic density objects unchanged", {
  raw <- c(-0.3, 0.2, 0.9, -0.1, 0.4, 1.1)
  obj <- make_proper_test_object(condens = raw)

  out <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "project",
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
  expect_equal(out$condens, raw, tolerance = 1e-12)
})

test_that("apply='fitted' also leaves exact-grid synthetic density objects unchanged", {
  raw <- c(-0.3, 0.2, 0.9, -0.1, 0.4, 1.1)
  obj <- make_proper_test_object(condens = raw)

  out <- getFromNamespace(".np_condens_finalize_proper_object", "npRmpi")(
    object = obj,
    proper = TRUE,
    proper.method = "project",
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
  expect_equal(out$condens, raw, tolerance = 1e-12)
})
