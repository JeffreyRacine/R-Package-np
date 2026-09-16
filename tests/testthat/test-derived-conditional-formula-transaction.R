test_that("derived conditional formulas share one constructor and fit sample", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(887L)
  d <- data.frame(y = rnorm(30), cy = factor(rep(1:2, 15)), x = runif(30))
  counter <- new.env(parent = emptyenv()); counter$n <- 0L
  ns <- environment(npqreg)
  if (exists(".npRmpi_autodispatch_materialize_call", ns, inherits = FALSE)) {
    dispatch <- new.env(parent = emptyenv()); dispatch$calls <- list()
    old.dispatch <- options(np.test.derived.dispatch = dispatch)
    on.exit(options(old.dispatch), add = TRUE)
    trace(".npRmpi_autodispatch_materialize_call", where = ns, tracer = quote({
      .record <- getOption("np.test.derived.dispatch")
      .record$calls <- c(.record$calls, list(mc))
    }), print = FALSE)
    on.exit(untrace(".npRmpi_autodispatch_materialize_call", where = ns), add = TRUE)
  }
  counted <- function(x) { counter$n <- counter$n + 1L; x }
  for (family in c("npqreg", "npconmode")) {
    fun <- get(family)
    f <- if (family == "npqreg") y ~ counted(x) else cy ~ counted(x)
    h <- if (family == "npqreg") c(.7, .6) else c(.2, .6)
    response <- if (family == "npqreg") d$y else d$cy
    for (route in c("positional", "named")) {
      counter$n <- 0L
      args <- if (route == "positional") list(f, data = d, bws = h) else
        list(formula = f, data = d, bws = h)
      fit <- do.call(fun, args)
      expect_identical(counter$n, 1L)
      if (family == "npqreg") {
        expect_identical(unname(fit$bws$ybw), h[1L])
        expect_identical(unname(fit$bws$xbw), h[2L])
        expect_true(is.na(fit$bws$fval))
        expect_true(is.na(fit$bws$num.feval))
      }
      native <- fun(bws = fit$bws, txdat = d["x"], tydat = response)
      expect_equal(fitted(fit), fitted(native), tolerance = 1e-12)
      counter$n <- 0L
      refit <- fun(bws = fit$bws)
      expect_identical(counter$n, 1L)
      expect_identical(fitted(refit), fitted(fit))
      expect_false(".np.formula.state" %in% names(fit$bws$call$...))
    }
    controls <- if (family == "npqreg") list(bwmethod = "normal-reference") else
      list(bandwidth.compute = FALSE)
    for (route in c("positional", "named", "bws")) {
      counter$n <- 0L
      args <- switch(route, positional = list(f, data = d),
        named = list(formula = f, data = d), bws = list(bws = f, data = d))
      fit <- do.call(fun, c(args, controls))
      expect_identical(counter$n, 1L)
      native <- fun(bws = fit$bws, txdat = d["x"], tydat = response)
      expect_equal(fitted(fit), fitted(native), tolerance = 1e-12)
    }
  }
  if (exists(".npRmpi_autodispatch_materialize_call", ns, inherits = FALSE)) {
    expect_false(any(vapply(dispatch$calls, function(x)
      is.symbol(x[[1L]]) && startsWith(as.character(x[[1L]]), "npqreg") &&
        "formula" %in% names(x), logical(1L))))
    expect_false(any(vapply(dispatch$calls, function(x)
      ".np.formula.state" %in% names(x), logical(1L))))
  }
})

test_that("mode validates its prepared response before bandwidth computation", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  counter <- new.env(parent = emptyenv()); counter$n <- 0L
  old.counter <- options(np.test.derived.counter = counter)
  on.exit(options(old.counter), add = TRUE)
  ns <- environment(npcdensbw)
  trace("npcdensbw.default", where = ns, tracer = quote({
    .counter <- getOption("np.test.derived.counter")
    .counter$n <- .counter$n + 1L
  }), print = FALSE)
  on.exit(untrace("npcdensbw.default", where = ns), add = TRUE)
  d <- data.frame(y = seq_len(24), x = seq(-1, 1, length.out = 24))
  for (route in c("positional", "named", "bws")) {
    args <- switch(route, positional = list(y ~ x, data = d),
      named = list(formula = y ~ x, data = d), bws = list(bws = y ~ x, data = d))
    expect_error(do.call(npconmode, args), "categorical response")
    expect_identical(counter$n, 0L)
  }
  d$y <- ordered(rep(1:3, 8))
  fit <- npconmode(y ~ x, data = d, bws = c(.2, .6))
  expect_length(fit$conmode, 24L)
  expect_gt(counter$n, 0L)
})

test_that("training validation is one-use and cannot leave a stored failed frame", {
  ns <- environment(npregbw)
  store <- get(".np_formula_frame_store", ns)
  take <- get(".np_formula_frame_take", ns)
  state <- new.env(parent = emptyenv())
  counter <- new.env(parent = emptyenv()); counter$n <- 0L
  state$validate.training <- function(frame) {
    counter$n <- counter$n + 1L
    stopifnot(identical(frame, data.frame(x = 1:3)))
  }
  store(state, data.frame(x = 1:3))
  expect_identical(counter$n, 1L)
  expect_identical(ls(state), "frame")
  expect_identical(take(state), data.frame(x = 1:3))
  expect_length(ls(state), 0L)
  state$validate.training <- function(frame) stop("invalid prepared response")
  expect_error(store(state, data.frame(x = 1:3)), "invalid prepared response")
  expect_length(ls(state), 0L)
  state$validate.training <- TRUE
  expect_error(store(state, data.frame(x = 1:3)), "invalid formula training validator")
  expect_false(exists("frame", state, inherits = FALSE))
})
