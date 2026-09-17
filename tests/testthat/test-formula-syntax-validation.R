b_formula_fixture <- function() {
  set.seed(92017)
  d <- data.frame(y = exp(rnorm(40, sd = .3)), x = runif(40), z = rnorm(40),
                  g = factor(rep(1:2, 20)), keep = rep(c(TRUE, FALSE), 20),
                  offset = seq_len(40)/40)
  d[["response value"]] <- d$y
  d
}

test_that("conditional response syntax is rejected before expression evaluation", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(list(np.messages = FALSE))
  d <- b_formula_fixture()
  count <- 0L
  touched <- function(x) { count <<- count + 1L; x }
  forms <- list(log(y) ~ x, I(y^2) ~ x, (y * z) ~ x, touched(y) ~ x)
  for (family in c("npcdensbw", "npcdistbw", "npcdens", "npcdist", "npqreg", "npconmode")) {
    fun <- get(family, mode = "function")
    for (f in forms) {
      expect_error(fun(formula = f, data = d, bws = c(.4, .4), bandwidth.compute = FALSE),
                   "conditional formula responses", fixed = TRUE)
    }
  }
  expect_identical(count, 0L)
  # A constructor after failure still owns exactly its valid sample.
  good <- npcdens(y ~ x, data = d, bws = c(.4, .4), bandwidth.compute = FALSE)
  native <- npcdens(txdat = d["x"], tydat = d$y, bws = c(.4, .4), bandwidth.compute = FALSE)
  expect_equal(fitted(good), fitted(native), tolerance = 1e-12)
})

test_that("actual offset specials fail before values while ordinary names remain valid", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(list(np.messages = FALSE))
  d <- b_formula_fixture()
  count <- 0L
  touched <- function(x) { count <<- count + 1L; x }
  formulas <- list(
    npudensbw = ~ x + offset(touched(z)), npudistbw = ~ x + offset(touched(z)),
    npregbw = y ~ x + offset(touched(z)), npcdensbw = y ~ x + offset(touched(z)),
    npcdistbw = y ~ x + offset(touched(z)), npindexbw = y ~ x + offset(touched(z)),
    npplregbw = y ~ x + offset(touched(z)) | z,
    npscoefbw = y ~ x | z + offset(touched(x)),
    npudens = ~ x + offset(touched(z)), npudist = ~ x + offset(touched(z)),
    npreg = y ~ x + offset(touched(z)), npcdens = y ~ x + offset(touched(z)),
    npcdist = y ~ x + offset(touched(z)), npindex = y ~ x + offset(touched(z)),
    npqreg = y ~ x + offset(touched(z)), npconmode = g ~ x + offset(touched(z)),
    npplreg = y ~ x | z + offset(touched(x)),
    npscoef = y ~ x + offset(touched(z)) | z,
    npcopula = ~ x + offset(touched(z)), npksum = y ~ x + offset(touched(z)),
    nplsqregbw = y ~ x + offset(touched(z)), nplsqreg = y ~ x + offset(touched(z)),
    npregiv = y ~ x + offset(touched(z)) | z,
    npregivderiv = y ~ x | z + offset(touched(x)))
  for (family in names(formulas)) {
    fun <- get(family, mode = "function")
    expect_error(fun(formulas[[family]], data = d), "offset() terms are not supported", fixed = TRUE)
  }
  expect_identical(count, 0L)
  validate <- get(".np_formula_validate_syntax", envir = environment(npregbw))
  expect_error(validate(y ~ x + (offset(z))), "offset() terms", fixed = TRUE)
  for (f in list(offset(y) ~ x, y ~ x + stats::offset(z), y ~ x + offset,
                 y ~ x + I(offset(z)), y ~ x + `offset(z)`))
    expect_identical(validate(f), f)
  for (f in list(offset(y) ~ x, y ~ x + stats::offset(z), y ~ x + offset)) {
    frame <- model.frame(f, d)
    tx <- frame[attr(terms(f), "term.labels")]
    widths <- rep(.4, ncol(tx))
    actual <- npreg(f, data = d, bws = widths, bandwidth.compute = FALSE)
    native <- npreg(txdat = tx, tydat = model.response(frame), bws = widths,
                    bandwidth.compute = FALSE)
    expect_equal(fitted(actual), fitted(native), tolerance = 1e-12)
  }
})

test_that("conditional literal stored dynamic dot and grouped formulas retain native targets", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(list(np.messages = FALSE))
  d <- b_formula_fixture()[c("y", "x", "z")]
  for (family in c("npcdens", "npcdist")) {
    fun <- get(family, mode = "function")
    bwfun <- get(paste0(family, "bw"), mode = "function")
    native <- fun(txdat = d[c("x", "z")], tydat = d$y,
                  bws = c(.4, .4, .4), bandwidth.compute = FALSE)
    dynamic <- as.formula(paste("y ~", paste(c("x", "z"), collapse = "+")))
    for (f in list(y ~ x + z, dynamic, y ~ ., (y) ~ x + z, +y ~ x + z)) {
      bw <- bwfun(f, data = d, bws = c(.4, .4, .4), bandwidth.compute = FALSE)
      actual <- fun(bws = bw)
      expect_equal(fitted(actual), fitted(native), tolerance = 1e-12)
      expect_equal(fitted(fun(formula = f, data = d, bws = c(.4, .4, .4),
                             bandwidth.compute = FALSE)), fitted(native), tolerance = 1e-12)
      nd <- d[1:5, , drop = FALSE]
      expected <- fun(bws = native$bws, txdat = d[c("x", "z")], tydat = d$y,
                      exdat = nd[c("x", "z")], eydat = nd$y)
      expect_equal(as.numeric(predict(actual, newdata = nd)), as.numeric(fitted(expected)), tolerance = 1e-12)
    }
    multivariate <- bwfun((y + z) ~ x, data = d, bws = c(.4, .4, .4), bandwidth.compute = FALSE)
    native.bw <- bwfun(xdat = d["x"], ydat = d[c("y", "z")],
                      bws = c(.4, .4, .4), bandwidth.compute = FALSE)
    expect_equal(multivariate$bw, native.bw$bw, tolerance = 1e-12)
    expect_identical(multivariate$ynames, c("y", "z"))
  }
})

test_that("formula validation preserves lexical terms subset promises and reentry", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(list(np.messages = FALSE))
  d <- b_formula_fixture()
  count <- 0L
  counted <- function(x) { count <<- count + 1L; x }
  f <- y ~ counted(x)
  actual <- npcdens(formula = f, data = d, subset = keep, bws = c(.4, .4), bandwidth.compute = FALSE)
  expect_identical(count, 1L)
  native <- npcdens(txdat = d[d$keep, "x", drop = FALSE], tydat = d$y[d$keep],
                    bws = c(.4, .4), bandwidth.compute = FALSE)
  expect_equal(fitted(actual), fitted(native), tolerance = 1e-12)
  wrapper <- function(data) {
    local.x <- data$x
    formula <- as.formula("y ~ local.x", env = environment())
    npcdens(formula = formula, data = data, bws = c(.4, .4), bandwidth.compute = FALSE)
  }
  local.fit <- wrapper(d)
  expect_equal(fitted(npcdens(bws = local.fit$bws)), fitted(local.fit), tolerance = 1e-12)
  reentrant <- function(x) {
    out <- npreg(txdat = data.frame(x = x), tydat = x, bws = .4, bandwidth.compute = FALSE)
    expect_length(fitted(out), length(x))
    x
  }
  nested <- npcdens(y ~ reentrant(x), data = d, bws = c(.4, .4), bandwidth.compute = FALSE)
  plain <- npcdens(y ~ x, data = d, bws = c(.4, .4), bandwidth.compute = FALSE)
  expect_equal(fitted(nested), fitted(plain), tolerance = 1e-12)
  constructed <- 0L
  make.formula <- function(text) {
    constructed <<- constructed + 1L
    as.formula(text, env = environment())
  }
  resolved <- make.formula("y ~ x")
  built <- npcdens(formula = resolved, data = d,
                   bws = c(.4, .4), bandwidth.compute = FALSE)
  expect_identical(constructed, 1L)
  expect_equal(fitted(built), fitted(plain), tolerance = 1e-12)
  constructed <- 0L
  resolved.pipe <- make.formula("y ~ x | z")
  built.pipe <- npplregbw(formula = resolved.pipe, data = d,
                         bws = matrix(.4, 2L, 1L), bandwidth.compute = FALSE)
  expect_identical(constructed, 1L)
  expect_s3_class(built.pipe, "plbandwidth")
  d$ly <- log(d$y)
  transformed.data <- npcdens(ly ~ x, data = d, bws = c(.4, .4), bandwidth.compute = FALSE)
  explicit <- npcdens(txdat = d["x"], tydat = d$ly, bws = c(.4, .4), bandwidth.compute = FALSE)
  expect_equal(fitted(transformed.data), fitted(explicit), tolerance = 1e-12)
  regression <- npreg(log(y) ~ x, data = d, bws = .4, bandwidth.compute = FALSE)
  reg.native <- npreg(txdat = d["x"], tydat = log(d$y), bws = .4, bandwidth.compute = FALSE)
  expect_equal(fitted(regression), fitted(reg.native), tolerance = 1e-12)
  density <- npudens(~ log(y), data = d, bws = .4, bandwidth.compute = FALSE)
  den.native <- npudens(tdat = data.frame(y = log(d$y)), bws = .4, bandwidth.compute = FALSE)
  expect_equal(fitted(density), fitted(den.native), tolerance = 1e-12)
  # Syntax-only validation does not evaluate backticked response names.
  validate <- get(".np_formula_validate_syntax", envir = environment(npregbw))
  backtick <- `response value` ~ x
  expect_identical(validate(backtick, conditional.response = TRUE), backtick)
})

test_that("saved conditional formulas are checked only when formula metadata is used", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(list(np.messages = FALSE))
  d <- b_formula_fixture()
  for (family in c("npcdens", "npcdist")) {
    bwfun <- get(paste0(family, "bw"), mode = "function")
    fun <- get(family, mode = "function")
    bw <- bwfun(y ~ x, data = d, bws = c(.4, .4), bandwidth.compute = FALSE)
    plain <- fun(bws = bw, txdat = d["x"], tydat = d$y)
    # The old implementation stored bare-y terms even when formula said log(y).
    legacy <- bw
    legacy$formula <- log(y) ~ x
    legacy$call$formula <- legacy$formula
    original <- legacy
    expect_error(fun(bws = legacy), "conditional formula responses", fixed = TRUE)
    expect_error(plot(legacy, plot.behavior = "data", neval = 5L),
                 "conditional formula responses", fixed = TRUE)
    expect_identical(legacy, original)
    # Explicit native data bypasses formula replay, retaining native precedence.
    native <- fun(bws = legacy, txdat = d["x"], tydat = d$y)
    expect_equal(fitted(native), fitted(plain), tolerance = 1e-12)
  }
})
