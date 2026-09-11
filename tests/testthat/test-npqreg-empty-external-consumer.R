test_that("selected-CDF inversion preserves healthy work and typed empty-row identity", {
  ns <- asNamespace("npRmpi")
  inv <- get(".npqreg_invert_selected_cdf", ns)
  cache.new <- get(".npqreg_selected_cdf_cache_new", ns)
  keys <- get(".npqreg_selected_cdf_cache_row_keys", ns)
  b <- list(regtype.engine = "lc", basis.engine = "glp", degree.engine = 0L,
            bernstein.basis.engine = FALSE, xncon = 1L, cykerorder = 2L)
  x <- data.frame(x = c(0, .5, 1))
  y <- data.frame(y = c(-1, 1))
  ex <- data.frame(x = c(.2, .8, .2, .5, .8))
  invoke <- function(fun, exdat = ex, cache = NULL, allow = FALSE, tau = .5, itmax = 100L) {
    inv(b, x, y, exdat, tau, tol = 1e-4, small = 1e-5, itmax = itmax,
        cdf.cache = cache, cdf.row.keys = keys(exdat),
        cdf.values = fun, allow.external = allow)
  }
  make.cdf <- function(contrast = FALSE, conflict = FALSE) {
    function(bws, xdat, ydat, exdat, ycand) {
      value <- pmax(0, pmin(1, (as.double(ycand) + 1)/2 + .2*(exdat$x - .5)))
      flags <- as.integer(exdat$x > 1)
      if (conflict && all(ycand == 1)) flags[] <- 0L
      if (any(flags == 1L)) {
        if (!contrast) {
          value[flags == 1L] <- NA_real_
          attr(value, ".np.empty.base.rows") <- flags
        }
        attr(value, ".np.empty.rows") <- flags
      }
      value
    }
  }
  controls <- list()
  for (enabled in c(FALSE, TRUE)) {
    cache <- cache.new(enabled)
    calls <- 0L
    rows <- 0L
    healthy <- function(bws, xdat, ydat, exdat, ycand) {
      calls <<- calls + 1L
      rows <<- rows + nrow(exdat)
      make.cdf()(bws, xdat, ydat, exdat, ycand)
    }
    controls[[as.character(enabled)]] <- lapply(c(.25, .5, .75), function(tau)
      invoke(healthy, cache = cache, tau = tau))
    expect_identical(calls, if (enabled) 45L else 51L)
    expect_identical(rows, if (enabled) 132L else 255L)
    expect_identical(cache$hits, if (enabled) 35L else 0L)
  }
  expect_identical(controls[["TRUE"]], controls[["FALSE"]])
  mixed <- ex
  mixed$x[c(2, 5)] <- 10
  for (enabled in c(FALSE, TRUE)) {
    cache <- cache.new(enabled)
    for (j in seq_along(c(.25, .5, .75))) {
      value <- invoke(make.cdf(), mixed, cache, TRUE, c(.25, .5, .75)[[j]])
      expect_identical(which(is.na(value)), c(2L, 5L))
      expect_identical(attr(value, ".np.empty.base.rows"), c(0L, 1L, 0L, 0L, 1L))
      expect_identical(as.double(value[c(1, 3, 4)]),
                       as.double(controls[["TRUE"]][[j]][c(1, 3, 4)]))
    }
    expect_true(all(is.na(invoke(make.cdf(), mixed[c(2, 5), , drop = FALSE], cache, TRUE))))
    contrast <- invoke(make.cdf(contrast = TRUE), mixed, cache.new(enabled), TRUE)
    expect_true(all(is.finite(contrast)))
    expect_null(attr(contrast, ".np.empty.base.rows"))
  }
  expect_error(invoke(make.cdf(), mixed), "required CDF evaluation", fixed = TRUE)
  expect_error(invoke(make.cdf(conflict = TRUE), mixed, allow = TRUE),
               "inconsistent base-support", fixed = TRUE)
  expect_error(invoke(function(bws, xdat, ydat, exdat, ycand)
    rep(NA_real_, nrow(exdat)), allow = TRUE), "non-finite bracket", fixed = TRUE)
  expect_error(invoke(function(bws, xdat, ydat, exdat, ycand)
    ifelse(abs(ycand) == 1, (ycand + 1)/2, NaN), allow = TRUE),
    "non-finite refinement", fixed = TRUE)
  expect_error(invoke(function(bws, xdat, ydat, exdat, ycand)
    structure(rep(NA_real_, nrow(exdat)), .np.empty.base.rows = rep(1, nrow(exdat))),
    allow = TRUE), "malformed empty-row", fixed = TRUE)
  expect_error(invoke(make.cdf(), itmax = 1L), "failed to converge", fixed = TRUE)
  for (bound in c(0, 1))
    expect_true(all(invoke(function(...) rep(bound, nrow(ex))) == if (bound == 0) 1 else -1))
  constant <- inv(b, x, data.frame(y = rep(1, 3)), ex, .5,
    tol = 1e-4, small = 1e-5, itmax = 100L,
    cdf.values = function(...) stop("constant response invoked CDF"),
    allow.external = TRUE)
  expect_true(all(constant == 1))
})

test_that("quantile delta evaluates only typed-supported rows and restores public arrays", {
  ns <- asNamespace("npRmpi")
  env <- new.env(parent = ns)
  fn <- get(".npqreg_quantile_delta_from_conditional", ns)
  environment(fn) <- env
  env$.npqreg_quantile_delta_from_conditional <- fn
  calls <- list()
  env$.np_conditional_eval_selected <- function(bws, xdat, ydat, exdat, eydat,
      cdf, gradients, se, allow.external, .np.defer.empty.rows, ...) {
    expect_true(allow.external && .np.defer.empty.rows)
    expect_true(all(is.finite(eydat[[1L]])))
    calls[[length(calls) + 1L]] <<- exdat$x
    if (cdf) list(condist = rep(.5, nrow(exdat)), conderr = rep(.2, nrow(exdat)),
                  congrad = matrix(.3, nrow(exdat), 1L),
                  congerr = matrix(.1, nrow(exdat), 1L))
    else list(condens = rep(2, nrow(exdat)))
  }
  b <- list(regtype.engine = "lc", basis.engine = "glp", degree.engine = 0L,
            bernstein.basis.engine = FALSE, xncon = 1L, cykerorder = 2L,
            xnuno = 0L, xnord = 0L, xndim = 1L, ixuno = FALSE, ixord = FALSE)
  x <- data.frame(x = c(.2, 10, .2))
  y <- data.frame(y = c(-1, 1))
  q <- structure(c(0, NA_real_, 0), .np.empty.base.rows = c(0L, 1L, 0L),
                 .np.empty.rows = c(0L, 1L, 0L))
  out <- fn(b, x, y, x, q, gradients = TRUE, se = TRUE, allow.external = TRUE)
  expect_identical(calls, list(c(.2, .2), c(.2, .2)))
  expect_identical(out$evaluated.rows, c(1L, 3L))
  expect_identical(out$quanterr, c(.1, NA_real_, .1))
  expect_identical(as.double(out$quantgrad), c(-.15, NA_real_, -.15))
  expect_identical(dim(out$quantgrad), c(3L, 1L))
  expect_error(fn(b, x, y, x, q, se = TRUE), "required delta evaluation", fixed = TRUE)
  all.empty <- structure(rep(NA_real_, 3), .np.empty.base.rows = rep(1L, 3))
  empty <- fn(b, x, y, x, all.empty, gradients = TRUE, se = TRUE, allow.external = TRUE)
  expect_length(calls, 2L)
  expect_null(empty$cdf)
  expect_null(empty$dens)
  expect_identical(empty$evaluated.rows, integer(0L))
  expect_true(all(is.na(empty$quantgrad)))
  expect_identical(dim(empty$quantgrad), c(3L, 1L))
})

test_that("unsupported categorical counterfactuals do not erase the base quantile", {
  ns <- asNamespace("npRmpi")
  env <- new.env(parent = ns)
  fn <- get(".npqreg_categorical_first_differences", ns)
  environment(fn) <- env
  env$.npqreg_invert_selected_cdf <- function(bws, xdat, ydat, exdat, tau,
      tol, small, itmax, cdf.cache, cdf.row.keys, allow.external) {
    expect_true(allow.external)
    value <- as.double(exdat$z == "b")
    empty <- as.integer(exdat$z == "a" & exdat$x > .4)
    if (any(empty == 1L)) {
      value[empty == 1L] <- NA_real_
      attr(value, ".np.empty.base.rows") <- empty
      attr(value, ".np.empty.rows") <- empty
    }
    value
  }
  ex <- data.frame(x = c(.2, .8), z = factor(c("b", "b"), levels = c("a", "b")))
  b <- list(ixuno = c(FALSE, TRUE), ixord = c(FALSE, FALSE), xndim = 2L)
  out <- fn(b, ex, data.frame(y = c(-1, 1)), ex, .5,
            tol = 1e-4, small = 1e-5, itmax = 100L, allow.external = TRUE)
  expect_identical(out[, 2L], c(1, NA_real_))
  expect_identical(attr(out, ".np.empty.rows"), c(0L, 1L))
  expect_null(attr(out, ".np.empty.base.rows"))
})

test_that("quantile fanout metadata follows query rows without changing numeric layouts", {
  ns <- asNamespace("npRmpi")
  fanout <- get(".npRmpi_bootstrap_run_fanout", ns)
  expect_null(formals(fanout)$metadata.reducer)
  expect_gt(match("metadata.reducer", names(formals(fanout))),
            match("...", names(formals(fanout))))
  worker.called <- FALSE
  expect_error(fanout(list(list(start = 1L, bsz = 1L)),
    function(...) { worker.called <<- TRUE }, ncol.out = 1L,
    metadata.reducer = "invalid"), "invalid internal fan-out metadata reducer", fixed = TRUE)
  expect_false(worker.called)
  collect <- get(".npqreg_collect_empty_rows", ns)
  decode <- get(".npqreg_fit_tau_vector_from_parallel_matrix", ns)
  tasks <- list(list(start = 1L, bsz = 2L), list(start = 3L, bsz = 3L))
  parts <- list(matrix(as.double(1:4), 2, 2), matrix(as.double(5:10), 3, 2))
  values <- do.call(rbind, parts)
  expect_identical(collect(values, parts, tasks, 2L), values)
  parts[[1L]][2L, ] <- NA_real_
  values <- do.call(rbind, parts)
  attr(parts[[1L]], ".npqreg.empty.tau.base.rows") <- matrix(c(0L, 1L, 0L, 1L), 2, 2)
  attr(parts[[1L]], ".npqreg.empty.tau.rows") <- matrix(c(0L, 1L, 0L, 1L), 2, 2)
  attr(parts[[2L]], ".npqreg.empty.tau.rows") <- matrix(c(0L, 0L, 1L, 1L, 0L, 0L), 3, 2)
  parts <- unserialize(serialize(parts, NULL))
  out <- collect(values, parts, tasks, 2L)
  expect_identical(as.vector(out), as.vector(values))
  expect_identical(dim(out), dim(values))
  expect_identical(attr(out, ".np.empty.base.rows"), c(0L, 1L, 0L, 0L, 0L))
  expect_identical(attr(out, ".np.empty.rows"), c(0L, 1L, 1L, 0L, 1L))
  expect_identical(dim(attr(out, ".npqreg.empty.tau.rows")), c(5L, 2L))
  decoded <- decode(out, c(.25, .75), se = FALSE)
  expect_identical(as.vector(decoded$yq), as.vector(values))
  expect_identical(attr(decoded, ".npqreg.empty.tau.rows"), attr(out, ".npqreg.empty.tau.rows"))
  bad.tasks <- tasks
  bad.tasks[[2L]]$start <- 2L
  expect_error(collect(values, parts, bad.tasks, 2L), "chunk row identity", fixed = TRUE)
  bad.tasks[[2L]]$start <- 4L
  expect_error(collect(values, parts, bad.tasks, 2L), "chunk row identity", fixed = TRUE)
  bad.parts <- parts
  attr(bad.parts[[2L]], ".npqreg.empty.tau.rows") <- matrix(1L, 2, 2)
  expect_error(collect(values, bad.parts, tasks, 2L), "tau empty-row status", fixed = TRUE)
  bad.parts <- parts
  attr(bad.parts[[1L]], ".npqreg.empty.tau.base.rows")[2, 2] <- 0L
  expect_error(collect(values, bad.parts, tasks, 2L), "differs across tau", fixed = TRUE)
  bad.values <- values
  bad.values[2, 1] <- 3
  expect_error(collect(bad.values, parts, tasks, 2L), "finite quantiles", fixed = TRUE)
  cdf.parts <- list(c(0, NA_real_), c(1, .5, 0))
  attr(cdf.parts[[1L]], ".np.empty.base.rows") <- c(0L, 1L)
  attr(cdf.parts[[1L]], ".np.empty.rows") <- c(0L, 1L)
  cdf.values <- matrix(c(0, NA_real_, 1, .5, 0), ncol = 1L)
  cdf.out <- collect(cdf.values, cdf.parts, tasks)
  expect_identical(as.vector(cdf.out), as.vector(cdf.values))
  expect_identical(attr(cdf.out, ".np.empty.base.rows"), c(0L, 1L, 0L, 0L, 0L))
  bad.parts <- cdf.parts
  attr(bad.parts[[1L]], ".np.empty.rows") <- c(0, 1)
  expect_error(collect(cdf.values, bad.parts, tasks), "malformed empty-row status", fixed = TRUE)
})
