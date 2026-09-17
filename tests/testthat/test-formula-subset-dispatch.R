p2_subset_fixture <- function() {
  set.seed(91621)
  d <- data.frame(y = rnorm(36), x = runif(36), z = runif(36),
    f = factor(rep(letters[1:3], 12)), keep = rep(c(TRUE, FALSE, TRUE), 12))
  rownames(d) <- paste0("case", seq_len(nrow(d)))
  d
}

test_that("MPI unconditional distribution raw bandwidth-formula branch owns subset", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(np.messages = FALSE)
  d <- data.frame(f = ordered(rep(1:4, 9)), keep = rep(c(TRUE, FALSE, TRUE), 12))
  # For ordered data, the default unsearched zero width is an exact
  # categorical limit. Supplying numeric bws would bypass this MPI branch.
  bw <- npudistbw(~ f, data = d, subset = keep, bandwidth.compute = FALSE)
  first <- npudist(bw)
  native <- npudist(bw, tdat = d[d$keep, "f", drop = FALSE])
  probe <- new.env(parent = emptyenv())
  probe$n <- 0L
  actual <- npudist(bws = ~ f, data = d,
    subset = { probe$n <- probe$n + 1L; keep }, bandwidth.compute = FALSE)
  expect_identical(probe$n, 1L)
  expect_identical(fitted(actual), fitted(first))
  expect_identical(fitted(actual), fitted(native))
  expect_identical(actual$bws$bw, bw$bw)
  expect_true(all(is.finite(fitted(actual))))
  expect_identical(actual$ntrain, native$ntrain)
})

p2_subset_specs <- function() {
  list(
    npudens = list(formula = quote(~ x), constructor = "npudensbw", widths = .5, kind = "unconditional"),
    npudist = list(formula = quote(~ x), constructor = "npudistbw", widths = .5, kind = "unconditional"),
    npreg = list(formula = quote(y ~ x), constructor = "npregbw", widths = .5, kind = "regression"),
    npcdens = list(formula = quote(y ~ x), constructor = "npcdensbw", widths = c(.6, .5), kind = "regression"),
    npcdist = list(formula = quote(y ~ x), constructor = "npcdistbw", widths = c(.6, .5), kind = "regression"),
    npplreg = list(formula = quote(y ~ x | z), constructor = "npplregbw", widths = matrix(.5, 2, 1), kind = "multipart"),
    npscoef = list(formula = quote(y ~ x | z), constructor = "npscoefbw", widths = .5, kind = "multipart"),
    npindex = list(formula = quote(y ~ x + z), constructor = "npindexbw", widths = c(1, .25, .5), kind = "index"),
    npqreg = list(formula = quote(y ~ x), constructor = "npcdistbw", widths = c(.6, .5), kind = "regression"),
    npconmode = list(formula = quote(f ~ x), constructor = "npcdensbw", widths = c(.3, .5), kind = "classification"),
    npcopula = list(formula = quote(~ x + z), constructor = "npudensbw", widths = c(.5, .5), kind = "copula"),
    npksum = list(formula = quote(y ~ x), constructor = NULL, widths = .5, kind = "kernel"),
    nplsqreg = list(formula = quote(y ~ x), constructor = "nplsqregbw", widths = NULL, kind = "lsq")
  )
}

p2_subset_eval <- function(fun, args, envir) {
  eval(as.call(c(list(as.name(fun)), args)), envir = envir)
}

p2_subset_case <- function(family, route) {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_preserve_seed()
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  d <- p2_subset_fixture()
  selected <- d[d$keep, , drop = FALSE]
  spec <- p2_subset_specs()[[family]]
  controls <- if (family == "npksum") list() else list(bandwidth.compute = FALSE)
  fit.controls <- if (family == "npcopula") list(target = "density", evaluation = "sample") else
    if (family == "npindex") list(se = FALSE) else list()
  if (family == "nplsqreg")
    controls <- c(controls, list(delta = .5, scale = rep(1, nrow(selected)),
                                regtype = "ll", nomad = FALSE))
  # LSQ's named formula is bws=; it does not take a second numeric bws.
  formula.args <- setNames(list(spec$formula),
    if (family == "nplsqreg") "bws" else "formula")
  width.args <- if (is.null(spec$widths)) list() else list(bws = spec$widths)
  if (!is.null(spec$constructor)) {
    bw <- p2_subset_eval(spec$constructor,
      c(formula.args, list(data = quote(d), subset = quote(keep)), width.args, controls),
      environment())
    first.data <- if (family == "npcopula") list(data = quote(selected[c("x", "z")])) else list()
    first <- p2_subset_eval(family, c(list(bws = quote(bw)), first.data, fit.controls), environment())
  } else {
    # npksum has no bandwidth constructor: an explicit model frame is its
    # constructor-first control, separate from the row-filtered native oracle.
    frame <- model.frame(y ~ x, data = d, subset = keep)
    first <- npksum(txdat = frame["x"], tydat = frame$y, bws = spec$widths)
  }
  native <- switch(spec$kind,
    unconditional = list(tdat = selected["x"]),
    copula = list(data = selected[c("x", "z")]),
    multipart = list(txdat = selected["x"], tydat = selected$y, tzdat = selected["z"]),
    classification = list(txdat = selected["x"], tydat = selected$f),
    index = list(txdat = selected[c("x", "z")], tydat = selected$y),
    list(txdat = selected["x"], tydat = selected$y))
  native$bws <- if (is.null(spec$constructor)) spec$widths else bw
  oracle <- do.call(get(family, mode = "function"), c(native, fit.controls))
  payload <- function(value) if (family == "npksum") value$ksum else fitted(value)
  expect_equal(payload(first), payload(oracle), tolerance = 1e-10)
  if (route == "positional") names(formula.args) <- ""
  actual <- p2_subset_eval(family,
    c(formula.args, list(data = quote(d), subset = quote(keep)), width.args,
      controls, fit.controls), environment())
  expect_equal(payload(actual), payload(oracle), tolerance = 1e-10)
  expect_equal(payload(actual), payload(first), tolerance = 1e-10)
}

for (family in names(p2_subset_specs())) {
  for (route in c("named", "positional")) {
    local({
      family <- family
      route <- route
      test_that(paste("data-only subset dispatch preserves", family, route), {
        p2_subset_case(family, route)
      })
    })
  }
}

test_that("subset expressions preserve data masks and forwarded promise owners once", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_preserve_seed()
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  d <- p2_subset_fixture()
  keep <- !d$keep
  probe <- new.env(parent = emptyenv())
  probe$n <- 0L
  counted <- function(value) { probe$n <- probe$n + 1L; value }
  forward <- function(data, rows) npreg(formula = y ~ x, data = data,
    subset = rows, bws = .5, bandwidth.compute = FALSE)
  expected <- npreg(txdat = d[d$keep, "x", drop = FALSE],
    tydat = d$y[d$keep], bws = .5, bandwidth.compute = FALSE)
  actual <- npreg(formula = y ~ x, data = d, subset = counted(keep),
    bws = .5, bandwidth.compute = FALSE)
  expect_identical(probe$n, 1L)
  expect_equal(fitted(actual), fitted(expected), tolerance = 1e-12)
  probe$n <- 0L
  nested <- forward(d, counted(keep))
  expect_identical(probe$n, 1L)
  # A forwarded promise is evaluated in its original caller, not recursively
  # reinterpreted in the data mask. Freeze the working constructor contract.
  constructor <- function(data, rows) npreg(npregbw(y ~ x, data = data,
    subset = rows, bws = .5, bandwidth.compute = FALSE))
  probe$n <- 0L
  nested.first <- constructor(d, counted(keep))
  expect_identical(probe$n, 1L)
  nested.oracle <- npreg(txdat = d[keep, "x", drop = FALSE], tydat = d$y[keep],
    bws = .5, bandwidth.compute = FALSE)
  expect_equal(fitted(nested), fitted(nested.first), tolerance = 1e-12)
  expect_equal(fitted(nested), fitted(nested.oracle), tolerance = 1e-12)
  lexical <- which(d$keep)
  probe$n <- 0L
  fallback <- forward(d, counted(lexical))
  expect_identical(probe$n, 1L)
  expect_equal(fitted(fallback), fitted(expected), tolerance = 1e-12)
  probe$n <- 0L
  expect_error(forward(d, { probe$n <- probe$n + 1L; stop("subset sentinel") }),
               "subset sentinel")
  expect_identical(probe$n, 1L)
  recovered <- forward(d, lexical)
  expect_equal(fitted(recovered), fitted(expected), tolerance = 1e-12)
  bound <- new.env(parent = environment())
  probe$n <- 0L
  makeActiveBinding("active.rows", function(value) {
    if (!missing(value)) stop("read-only subset binding")
    probe$n <- probe$n + 1L
    lexical
  }, env = bound)
  active <- eval(quote(npreg(formula = y ~ x, data = d, subset = active.rows,
    bws = .5, bandwidth.compute = FALSE)), envir = bound)
  expect_identical(probe$n, 1L)
  expect_equal(fitted(active), fitted(expected), tolerance = 1e-12)
})

test_that("unconditional subset expressions retain index and no-subset contracts", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_preserve_seed()
  withr::local_options(np.messages = FALSE)
  d <- p2_subset_fixture()
  choices <- list(which(d$keep), -which(!d$keep), d$keep)
  for (selected in choices) {
    actual <- npudens(formula = ~ x, data = d, subset = selected,
      bws = .5, bandwidth.compute = FALSE)
    oracle <- npudens(tdat = d[selected, "x", drop = FALSE],
      bws = .5, bandwidth.compute = FALSE)
    expect_equal(fitted(actual), fitted(oracle), tolerance = 1e-12)
  }
  actual <- npudens(~ x, data = d, subset = keep & x > .25,
    bws = .5, bandwidth.compute = FALSE)
  selected <- d$keep & d$x > .25
  oracle <- npudens(tdat = d[selected, "x", drop = FALSE],
    bws = .5, bandwidth.compute = FALSE)
  expect_equal(fitted(actual), fitted(oracle), tolerance = 1e-12)
  expect_equal(fitted(npudens(~ x, data = d, bws = .5, bandwidth.compute = FALSE)),
    fitted(npudens(tdat = d["x"], bws = .5, bandwidth.compute = FALSE)), tolerance = 1e-12)
})

test_that("subset row ownership survives NA actions prediction and native evaluation", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_preserve_seed()
  withr::local_options(np.messages = FALSE)
  d <- p2_subset_fixture()
  d$y[1L] <- NA_real_
  d$x[7L] <- NA_real_
  selected <- d[d$keep, c("y", "x"), drop = FALSE]
  complete <- selected[complete.cases(selected), , drop = FALSE]
  expected <- npreg(txdat = complete["x"], tydat = complete$y,
    bws = .5, bandwidth.compute = FALSE)
  for (action in list(na.omit, na.exclude)) {
    actual <- npreg(formula = y ~ x, data = d, subset = keep, na.action = action,
      bws = .5, bandwidth.compute = FALSE)
    value <- fitted(actual)
    expect_equal(as.numeric(value[!is.na(value)]), as.numeric(fitted(expected)), tolerance = 1e-12)
    expect_identical(as.integer(actual$bws$rows.omit), which(!complete.cases(selected)))
    if (identical(action, na.exclude))
      expect_identical(which(is.na(value)), which(!complete.cases(selected)))
    nd <- data.frame(x = c(.2, .5, .8))
    expect_equal(predict(actual, newdata = nd), predict(expected, exdat = nd), tolerance = 1e-12)
    expect_equal(predict(actual, newdata = data.frame(wrong = 1:3), exdat = nd),
                 predict(actual, exdat = nd), tolerance = 1e-12)
  }
  expect_error(npreg(formula = y ~ x, data = d, subset = keep, na.action = na.fail,
    bws = .5, bandwidth.compute = FALSE), "missing")
})

test_that("time-series subset indices apply after formula alignment", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(np.messages = FALSE)
  d <- list(x = ts(sin(seq_len(20L)), start = 3))
  selected <- c(1L, 4L, 8L, 12L)
  aligned <- as.data.frame(ts.intersect(x = d$x, previous = lag(d$x, -1)))
  names(aligned) <- c("x", "lag(x, -1)")
  actual <- npudens(formula = ~ x + lag(x, -1), data = d, subset = selected,
    bws = c(.5, .5), bandwidth.compute = FALSE)
  oracle <- npudens(tdat = aligned[selected, , drop = FALSE],
    bws = c(.5, .5), bandwidth.compute = FALSE)
  expect_equal(fitted(actual), fitted(oracle), tolerance = 1e-12)
})

test_that("significance formula subset uses the constructor-selected sample", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_preserve_seed()
  withr::local_options(np.messages = FALSE)
  d <- p2_subset_fixture()
  bw <- npregbw(y ~ x, data = d, subset = keep, bws = .5, bandwidth.compute = FALSE)
  expected <- npsigtest(bw, index = 1L, B = 9L, random.seed = 943L)
  native.bw <- npregbw(xdat = d[d$keep, "x", drop = FALSE], ydat = d$y[d$keep],
    bws = .5, bandwidth.compute = FALSE)
  oracle <- npsigtest(native.bw, index = 1L, B = 9L, random.seed = 943L)
  probe <- new.env(parent = emptyenv()); probe$n <- 0L
  actual <- npsigtest(formula = y ~ x, data = d, subset = { probe$n <- probe$n + 1L; keep },
    bws = .5, bandwidth.compute = FALSE, index = 1L, B = 9L, random.seed = 943L)
  expect_identical(probe$n, 1L)
  expect_equal(actual[c("In", "P", "In.bootstrap")],
               expected[c("In", "P", "In.bootstrap")], tolerance = 1e-12)
  expect_equal(actual[c("In", "P", "In.bootstrap")],
               oracle[c("In", "P", "In.bootstrap")], tolerance = 1e-12)
})
