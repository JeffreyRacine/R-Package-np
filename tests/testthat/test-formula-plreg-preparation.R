plreg_preparation_internal <- function(name) {
  get(name, envir = environment(getS3method("npplreg", "formula")), inherits = TRUE)
}

test_that("partially linear residuals retain training rows independently of evaluation rows", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  set.seed(1789)
  d <- data.frame(y = rnorm(28), x = runif(28), z = runif(28))
  d$x[c(3L, 8L)] <- NA
  keep <- complete.cases(d)
  bw <- npplregbw(y ~ x | z, data = d, na.action = na.exclude,
                  bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  nd <- data.frame(x = c(.2, NA, .5, .7, .8), z = c(.3, .4, .6, NA, .5))
  native <- npplreg(bw, txdat = d[keep, "x", drop = FALSE], tydat = d$y[keep],
    tzdat = d[keep, "z", drop = FALSE], exdat = nd["x"], ezdat = nd["z"],
    residuals = TRUE)
  expected.residuals <- rep(NA_real_, nrow(d))
  expected.residuals[keep] <- native$resid
  training <- npplreg(bw, residuals = TRUE)
  expect_identical(training$resid, expected.residuals)
  for (extra in list(list(newdata = nd),
                    list(exdat = nd["x"], ezdat = nd["z"]),
                    list(newdata = data.frame(wrong = 1),
                         exdat = nd["x"], ezdat = nd["z"]))) {
    value <- do.call(npplreg, c(list(bws = bw, residuals = TRUE), extra))
    expect_identical(fitted(value), fitted(native))
    expect_identical(which(is.na(fitted(value))), c(2L, 4L))
    expect_identical(value$resid, expected.residuals)
    expect_identical(which(is.na(value$resid)), c(3L, 8L))
    expect_length(value$resid, nrow(d))
  }
})

test_that("partially linear automatic preparation contexts clear after errors", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1790)
  d <- data.frame(y = rnorm(28), x1 = runif(28), x2 = runif(28), z = runif(28))
  h <- matrix(.6, 3, 1)
  states <- list()
  fail.preparation <- FALSE
  capture <- function(x) {
    candidates <- lapply(sys.frames(), function(frame) {
      if (exists("frame.state", frame, inherits = FALSE))
        get("frame.state", frame, inherits = FALSE) else NULL
    })
    candidates <- Filter(function(state)
      is.environment(state) && identical(parent.env(state), emptyenv()), candidates)
    if (!length(candidates))
      stop("test could not observe formula preparation context")
    states[[length(states) + 1L]] <<- candidates[[1L]]
    if (fail.preparation) stop("deterministic plreg preparation failure")
    x
  }
  f <- y ~ capture(x1) + x2 | z
  run <- function(dat) npplreg(formula = f, data = dat, bws = h,
                                bandwidth.compute = FALSE, residuals = TRUE)
  oracle <- npplreg(txdat = d[c("x1", "x2")], tydat = d$y, tzdat = d["z"],
                    bws = h, bandwidth.compute = FALSE, residuals = TRUE)
  for (stage in c("preparation", "numeric")) {
    fail.preparation <- identical(stage, "preparation")
    bad <- d
    if (!fail.preparation) bad$x2 <- bad$x1
    previous <- length(states)
    expect_error(run(bad), if (fail.preparation) "deterministic plreg preparation failure"
                           else "rank deficient after smoothing")
    expect_length(states, previous + 1L)
    failed.state <- states[[length(states)]]
    expect_identical(ls(failed.state, all.names = TRUE), character())
    expect_error(plreg_preparation_internal(".np_formula_frame_take")(failed.state),
                 "no training frame")
    fail.preparation <- FALSE
    value <- run(d)
    expect_length(states, previous + 2L)
    successful.state <- states[[length(states)]]
    expect_false(identical(failed.state, successful.state))
    expect_identical(ls(successful.state, all.names = TRUE), character())
    expect_identical(fitted(value), fitted(oracle))
    expect_identical(unname(coef(value)), unname(coef(oracle)))
    expect_identical(value$resid, oracle$resid)
    expect_false(grepl("formula.state", paste(deparse(value$bws$call), collapse = "")))
    expect_false(grepl("formula.state", paste(deparse(value$call), collapse = "")))
  }
})

test_that("partially linear trained scalar calls retain their training parameters", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1787)
  d <- data.frame(y = rnorm(32), x = runif(32), z = runif(32))
  centered <- function(x, center = mean(x))
    structure(x - center, class = c("np_plreg_test_center", "numeric"), center = center)
  make.call <- function(var, call) {
    if (is.call(call) && identical(call[[1L]], quote(centered)))
      call$center <- attr(var, "center")
    call
  }
  methods <- get(".__S3MethodsTable__.", envir = asNamespace("stats"))
  method.name <- "makepredictcall.np_plreg_test_center"
  previous <- if (exists(method.name, methods, inherits = FALSE))
    get(method.name, methods) else NULL
  on.exit(if (is.null(previous)) rm(list = method.name, envir = methods)
          else assign(method.name, previous, methods), add = TRUE)
  registerS3method("makepredictcall", "np_plreg_test_center", make.call,
                   envir = asNamespace("stats"))
  bw <- npplregbw(y ~ x | centered(z), data = d,
                  bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  expect_identical(attr(bw$terms, "predvars")[[3L]]$center, mean(d$z))
  nd <- transform(d[1:6, ], z = z + 1)
  value <- npplreg(bw, newdata = nd)
  tz <- setNames(data.frame(centered(d$z, mean(d$z))), bw$znames)
  ez <- setNames(data.frame(centered(nd$z, mean(d$z))), bw$znames)
  oracle <- npplreg(bw, txdat = d["x"], tydat = d$y, tzdat = tz,
                    exdat = nd["x"], ezdat = ez)
  expect_identical(fitted(value), fitted(oracle))
  expect_identical(value$evalz[[1]], ez[[1]])
  expect_identical(plreg_preparation_internal(".np_plot_plreg_training_data")(bw)$zdat[[1]],
                   tz[[1]])
})

test_that("partially linear native evaluation owns omissions and response precedence", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1788)
  d <- data.frame(y = rnorm(28), x = runif(28), z = runif(28))
  d$x[c(3L, 8L)] <- NA
  bw <- npplregbw(y ~ x | z, data = d, na.action = na.exclude,
                  bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  keep <- complete.cases(d)
  ex <- data.frame(x = c(.3, NA, .8, .5, .2))
  ez <- data.frame(z = c(.2, .4, .7, .5, .6))
  native <- npplreg(bw, txdat = d[keep, "x", drop = FALSE], tydat = d$y[keep],
                    tzdat = d[keep, "z", drop = FALSE], exdat = ex, ezdat = ez)
  value <- npplreg(bw, newdata = data.frame(wrong = 1), exdat = ex, ezdat = ez)
  expect_identical(fitted(value), fitted(native))
  expect_length(fitted(value), 5L)
  expect_identical(value$eval.rows.omit, native$eval.rows.omit)
  fit <- npplreg(bw)
  expect_length(fitted(fit), nrow(d))
  expect_identical(predict(fit, newdata = data.frame(wrong = 1), exdat = ex, ezdat = ez),
                   fitted(native))
  nd <- d[20:25, c("x", "z")]
  ey <- seq_len(nrow(nd)) / 10
  value <- npplreg(bw, newdata = nd, y.eval = TRUE, eydat = ey)
  oracle <- npplreg(bw, txdat = d[keep, "x", drop = FALSE], tydat = d$y[keep],
    tzdat = d[keep, "z", drop = FALSE], exdat = nd["x"], ezdat = nd["z"], eydat = ey)
  expect_identical(fitted(value), fitted(oracle))
  expect_identical(value$MSE, oracle$MSE)
})

test_that("partially linear formula transactions prepare every expression once", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1781)
  d <- data.frame(y = rnorm(36), x = runif(36), z = runif(36))
  counts <- setNames(integer(3), c("y", "x", "z"))
  counted <- function(value, name) { counts[[name]] <<- counts[[name]] + 1L; value }
  f <- counted(y, "y") ~ counted(x, "x") | counted(z, "z")
  h <- matrix(.6, 2L, 1L)
  bw <- npplregbw(f, data = d, bws = h, bandwidth.compute = FALSE)
  expect_identical(unname(counts), rep(1L, 3L))
  counts[] <- 0L
  ref <- npplreg(bw, se = TRUE, residuals = TRUE)
  expect_identical(unname(counts), rep(1L, 3L))
  for (args in list(list(f, data = d, bws = h),
                   list(formula = f, data = d, bws = h),
                   list(data = d, bws = h, formula = f))) {
    counts[] <- 0L; rng <- .Random.seed
    actual <- do.call(npplreg, c(args,
      list(bandwidth.compute = FALSE, se = TRUE, residuals = TRUE)))
    expect_identical(unname(counts), rep(1L, 3L))
    expect_identical(.Random.seed, rng)
    for (field in c("mean", "xcoef", "xcoeferr", "xcoefvcov", "resid", "R2", "MSE"))
      expect_identical(actual[[field]], ref[[field]], info = field)
    expect_false(grepl("formula.state", paste(deparse(actual$bws$call), collapse = "")))
  }
  counts[] <- 0L
  evaluated <- npplreg(bw, newdata = d[1:7, ], se = TRUE)
  expect_identical(unname(counts), c(1L, 2L, 2L))
  counts[] <- 0L
  expect_identical(predict(ref, newdata = d[1:7, ]), fitted(evaluated))
  expect_identical(unname(counts), c(1L, 2L, 2L))
  expect_error(predict(ref, newdata = data.frame(wrong = 1:7)), "columns.*x")
  for (helper in c(".np_plreg_predict_train_data", ".np_plot_plreg_training_data")) {
    counts[] <- 0L
    plreg_preparation_internal(helper)(bw)
    expect_identical(unname(counts), rep(1L, 3L), info = helper)
  }
  counts[] <- 0L
  npplreg(bw, newdata = d[1:7, ], y.eval = TRUE)
  expect_identical(unname(counts), rep(2L, 3L))
})

test_that("partially linear time alignment agrees with integer observation oracles", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1782)
  y <- ts(rnorm(30), frequency = 4)
  nd <- data.frame(y = rnorm(18)); nd$y <- ts(nd$y, frequency = 4)
  for (shift in c(-2L, 1L)) {
    f <- if (shift == -2L) y ~ lag(y, -1) | lag(y, -2) else y ~ lag(y, -1) | lag(y, 1)
    yi <- if (shift == -2L) 3:30 else 2:29
    xi <- if (shift == -2L) 2:29 else 1:28
    zi <- if (shift == -2L) 1:28 else 3:30
    exi <- if (shift == -2L) 2:18 else 1:16
    ezi <- if (shift == -2L) 1:17 else 3:18
    bw <- npplregbw(f, bws = matrix(.8, 2, 1), bandwidth.compute = FALSE)
    x <- setNames(data.frame(as.numeric(y)[xi]), bw$xnames)
    z <- setNames(data.frame(as.numeric(y)[zi]), bw$znames)
    ex <- setNames(data.frame(as.numeric(nd$y)[exi]), bw$xnames)
    ez <- setNames(data.frame(as.numeric(nd$y)[ezi]), bw$znames)
    fit <- npplreg(bw, se = TRUE, residuals = TRUE)
    ref <- npplreg(bw, txdat = x, tydat = as.numeric(y)[yi], tzdat = z,
                  se = TRUE, residuals = TRUE)
    val <- npplreg(bw, newdata = nd, se = TRUE)
    oracle <- npplreg(bw, txdat = x, tydat = as.numeric(y)[yi], tzdat = z,
                     exdat = ex, ezdat = ez, se = TRUE)
    for (field in c("mean", "xcoef", "xcoeferr", "xcoefvcov", "R2", "MSE")) {
      expect_identical(fit[[field]], ref[[field]], info = field)
      expect_identical(val[[field]], oracle[[field]], info = field)
    }
    expect_identical(predict(fit, newdata = data.frame(wrong = 1),
                             exdat = ex, ezdat = ez), fitted(oracle))
    ev.yi <- if (shift == -2L) 3:18 else 2:17
    ev.xi <- if (shift == -2L) 2:17 else 1:16
    ev.zi <- if (shift == -2L) 1:16 else 3:18
    ey <- npplreg(bw, newdata = nd, y.eval = TRUE)
    expected <- npplreg(bw, txdat = x, tydat = as.numeric(y)[yi], tzdat = z,
      exdat = setNames(data.frame(as.numeric(nd$y)[ev.xi]), bw$xnames),
      ezdat = setNames(data.frame(as.numeric(nd$y)[ev.zi]), bw$znames),
      eydat = as.numeric(nd$y)[ev.yi])
    expect_identical(fitted(ey), fitted(expected))
  }
})

test_that("partially linear preparation uses joint subset missingness and saved data ownership", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1783)
  d <- data.frame(y = rnorm(38), x = runif(38), z = runif(38),
                  g = factor(rep(letters[1:2], 19)))
  d$y[5] <- NA; d$x[9] <- NA; d$z[12] <- NA
  f <- y ~ x | z + g
  wrap <- function(input, ff) {
    dd <- input
    npplregbw(ff, data = dd, subset = seq_along(y) > 2L,
              bws = matrix(c(.6, .2), 2, 2, byrow = TRUE), bandwidth.compute = FALSE)
  }
  bw <- wrap(d, f)
  train <- plreg_preparation_internal(".np_plreg_formula_training")(bw)
  keep <- seq_len(nrow(d)) > 2L & complete.cases(d)
  expect_identical(as.numeric(train$tydat), d$y[keep])
  expect_identical(train$txdat[[1]], d$x[keep])
  expect_identical(train$tzdat$z, d$z[keep])
  expect_identical(train$tzdat$g, d$g[keep])
  fitted.formula <- npplreg(bw, se = TRUE)
  fitted.native <- npplreg(bw, txdat = d[keep, "x", drop = FALSE],
    tydat = d$y[keep], tzdat = d[keep, c("z", "g")], se = TRUE)
  expect_identical(fitted(fitted.formula), fitted(fitted.native))
  expect_identical(coef(fitted.formula), coef(fitted.native))
  override <- transform(d, x = x + .05)
  override.fit <- npplreg(bw, data = override)
  expect_identical(override.fit$evalx[[1]], override$x[keep])
  expect_identical(plreg_preparation_internal(".np_plot_plreg_training_data")(bw)$xdat,
                   train$txdat)
  ps <- predict(fitted.formula, newdata = d[20:24, ], se.fit = TRUE)
  expect_length(ps$fit, 5L)
  expect_true(all(is.finite(ps$se.fit)))
})

test_that("partially linear one-call stochastic and repeated-role expressions are single use", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1784)
  d <- data.frame(y = rnorm(30), x = runif(30), z = runif(30))
  jittered <- function(x) x + rnorm(length(x), sd = .01)
  set.seed(1785)
  x <- data.frame(jittered(d$x)); names(x) <- "jittered(x)"
  rng <- .Random.seed
  ref <- npplreg(txdat = x, tydat = d$y, tzdat = d["z"],
                 bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  set.seed(1785)
  fit <- npplreg(y ~ jittered(x) | z, data = d,
                 bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  expect_identical(.Random.seed, rng)
  expect_identical(fitted(fit), fitted(ref))
  counts <- 0L
  counted <- function(x) { counts <<- counts + 1L; x }
  bw <- npplregbw(y ~ counted(x) | counted(x) + z, data = d,
                  bws = matrix(.6, 2, 2), bandwidth.compute = FALSE)
  expect_identical(counts, 1L)
  counts <- 0L
  npplreg(bw)
  expect_identical(counts, 1L)
  bad <- bw
  attr(bad$xterms, "predvars")[[2L]] <- quote(x + 100)
  expect_error(plreg_preparation_internal(".np_plreg_formula_terms")(bad),
               "inconsistent.*trained")
})

test_that("partially linear legacy scaffolding is recognized exactly", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1786)
  y <- ts(rnorm(30), frequency = 4)
  bw <- npplregbw(y ~ lag(y, -1) | lag(y, -2),
                  bws = matrix(.8, 2, 1), bandwidth.compute = FALSE)
  original <- npplreg(bw)
  legacy <- bw
  joint <- plreg_preparation_internal(".np_plreg_formula_spec")(bw$formula)$joint
  arguments <- as.list(attr(joint, "variables"))[-1L]
  inner <- as.call(list(quote(as.data.frame),
    as.call(c(list(quote(ts.intersect)), arguments))))
  attr(legacy$terms, "predvars") <- substitute(INNER[, INDEX, drop = FALSE],
                                               list(INNER = inner, INDEX = c(1L, 3L)))
  attr(legacy$xterms, "predvars") <- substitute(INNER[, 2L, drop = FALSE],
                                                list(INNER = inner))
  expect_identical(fitted(npplreg(legacy)), fitted(original))
  bad <- legacy
  attr(bad$xterms, "predvars")[[4L]] <- 1L
  before <- attr(bad$xterms, "predvars")
  expect_error(plreg_preparation_internal(".np_plreg_formula_terms")(bad),
               "unsupported.*prediction")
  expect_identical(attr(bad$xterms, "predvars"), before)
  arbitrary <- bw
  attr(arbitrary$xterms, "predvars") <- quote(identity(list(lag(y, -1))))
  before <- attr(arbitrary$xterms, "predvars")
  expect_error(plreg_preparation_internal(".np_plreg_formula_terms")(arbitrary),
               "unsupported.*prediction")
  expect_identical(attr(arbitrary$xterms, "predvars"), before)
  x <- runif(28)
  mixed <- npplregbw(y ~ x | lag(y, -2),
                     bws = matrix(.8, 2, 1), bandwidth.compute = FALSE)
  old.mixed <- mixed
  inner <- substitute((cbind(as.data.frame(ts.intersect(y, lag(y, -2))), x,
                        check.rows = TRUE)[, INDEX]), list(INDEX = c(1L, 3L, 2L)))
  attr(old.mixed$terms, "predvars") <- substitute(INNER[, INDEX, drop = FALSE],
                                                  list(INNER = inner, INDEX = c(1L, 3L)))
  attr(old.mixed$xterms, "predvars") <- substitute(INNER[, 2L, drop = FALSE],
                                                   list(INNER = inner))
  expect_identical(fitted(npplreg(old.mixed)), fitted(npplreg(mixed)))
  far <- ts(rnorm(8), start = 100, frequency = 4)
  expect_error(suppressWarnings(npplregbw(y ~ far | lag(y, -2),
    bws = matrix(.8, 2, 1), bandwidth.compute = FALSE)))
  expect_error(npplregbw(y ~ x, bws = matrix(.8, 2, 1), bandwidth.compute = FALSE),
               "improper formula")
})
