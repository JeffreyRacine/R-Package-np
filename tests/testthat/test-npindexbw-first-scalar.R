test_that("first-scalar restoration has exact bounded work and no rejected stencils", {
  for (method in c("Nelder-Mead", "BFGS", "CG")) {
    for (nreject in c(0L, 1L, 3L, 8L)) {
      calls <- gradients <- numeric()
      controls <- list(maxit = 17L, reltol = 1e-7)
      factory <- .npindexbw_first_scalar_guard
      scope <- new.env(parent = environment(factory))
      environment(factory) <- scope
      scope$.np_progress_bandwidth_notice <- function(labels) invisible(NULL)
      scope$optim <- function(par, fn, gr, method, control) {
        expect_identical(control, controls)
        expect_identical(fn, args$fn)
        value <- fn(par)
        if (method != "Nelder-Mead") {
          gradients <<- c(gradients, par[2L])
          fn(par + c(0, .001))
          fn(par - c(0, .001))
        }
        list(par = par, value = value, convergence = 0L)
      }
      guard <- factory("ichimura")
      fn <- function(par) {
        calls <<- c(calls, par[2L])
        raw <- if (par[2L] < 2^nreject) .Machine$double.xmax else 1
        guard$observe(raw, par)
        if (identical(raw, .Machine$double.xmax)) 10 else raw
      }
      args <- list(par = c(.75, 1), fn = fn, gr = NULL, method = method,
                   control = controls)
      set.seed(606)
      before <- .Random.seed
      got <- guard$run(args, TRUE, FALSE, 2, .1, 1, 1L, 0L)
      expect_identical(got$par, c(.75, 2^nreject))
      expect_identical(calls[seq_len(nreject + 1L)], 2^(0:nreject))
      expect_length(calls, nreject + 1L + if (method == "Nelder-Mead") 0L else 2L)
      expect_identical(gradients, if (method == "Nelder-Mead") numeric() else 2^nreject)
      expect_identical(.Random.seed, before)
      expect_false(guard$active)
    }
  }
})

test_that("valid starts stay in the original real optim invocation", {
  for (method in c("Nelder-Mead", "BFGS", "CG")) {
    guard <- .npindexbw_first_scalar_guard("kleinspady")
    old.calls <- new.calls <- list()
    ordinary <- function(par) sum((par - c(.3, 4))^2)
    old.fn <- function(par) { old.calls[[length(old.calls)+1L]] <<- par; ordinary(par) }
    new.fn <- function(par) {
      new.calls[[length(new.calls)+1L]] <<- par
      raw <- ordinary(par)
      guard$observe(raw, par)
      raw
    }
    args <- list(par = c(.7, 2), fn = old.fn, gr = NULL, method = method,
                 control = list(maxit = 50L, reltol = 1e-7))
    old <- do.call(optim, args)
    args$fn <- new.fn
    got <- guard$run(args, TRUE, FALSE, 1, .1, 2, 1L, 0L)
    expect_identical(got, old)
    expect_identical(new.calls, old.calls)
    expect_false(guard$active)
  }
})

test_that("first-scalar failures never become a survivor or generic rescue", {
  factory <- .npindexbw_first_scalar_guard
  scope <- new.env(parent = environment(factory))
  environment(factory) <- scope
  invocations <- 0L
  scope$optim <- function(par, fn, gr, method, control, ...) {
    invocations <<- invocations + 1L
    list(par = par, value = fn(par, ...), convergence = 0L)
  }
  scope$.np_progress_bandwidth_notice <- function(labels) invisible(NULL)
  guard <- factory("ichimura")
  fn <- function(par, h = NULL) {
    guard$observe(.Machine$double.xmax, if (is.null(h)) par else c(par, h))
    stop("unreachable mapped first value")
  }
  args <- list(par = c(.75, 1), fn = fn, gr = NULL, method = "BFGS", control = list())
  run <- function(...) guard$run(args, scale = 1, lower = .1, h = 1,
                                 start = 2L, retry = 1L, ...)
  expect_error(run(automatic = TRUE, held = FALSE), "9 starting scalar.*eight-doubling")
  expect_identical(invocations, 9L)
  expect_false(guard$active)
  invocations <- 0L
  expect_error(run(automatic = FALSE, held = FALSE), "1 starting scalar.*explicit initial")
  expect_identical(invocations, 1L)
  invocations <- 0L
  expect_error(run(automatic = TRUE, held = FALSE, upper = 1.5), "existing search bounds")
  expect_identical(invocations, 1L)
  for (scale in c(0, NA_real_, Inf)) {
    expect_error(guard$run(args, TRUE, FALSE, scale, .1, 1, 1L, 0L), "positive finite index scale")
    expect_false(guard$active)
  }
  huge <- .Machine$double.xmax / 1.5
  big.args <- args; big.args$par[2L] <- huge
  expect_error(guard$run(big.args, TRUE, FALSE, 1, .1, huge, 1L, 0L), "finite strictly larger")
  held.args <- args; held.args$par <- .75; held.args$h <- 1
  expect_error(guard$run(held.args, FALSE, TRUE, 1, .1, 1, 1L, 0L), "raw-invalid held bandwidth")
  missing <- args; missing$fn <- function(par) 10
  expect_error(guard$run(missing, TRUE, FALSE, 1, .1, 1, 1L, 0L), "without its current ordinary raw")
  wrong.beta <- args; wrong.beta$fn <- function(par) guard$observe(1, par + c(1, 0))
  expect_error(guard$run(wrong.beta, TRUE, FALSE, 1, .1, 1, 1L, 0L), "stale or mismatched raw witness")
  converted <- factory("ichimura", function(beta) beta * 2)
  converted.args <- args; converted.args$fn <- function(par) {
    converted$observe(1, c(par[1L] * 2, par[2L])); 1
  }
  expect_equal(converted$run(converted.args, TRUE, FALSE, 1, .1, 1, 1L, 0L)$value, 1)
  malformed <- args; malformed$fn <- function(par) guard$observe(NULL, par)
  expect_error(guard$run(malformed, TRUE, FALSE, 1, .1, 1, 1L, 0L), "normal finite raw")
  unexpected <- args; unexpected$fn <- function(par) stop("ordinary unrelated failure")
  expect_error(guard$run(unexpected, TRUE, FALSE, 1, .1, 1, 1L, 0L), "ordinary unrelated failure")
  forged <- args; forged$fn <- function(par) stop(structure(
    list(message = "forged same-class failure", call = NULL, token = new.env()),
    class = c("np_index_first_scalar_invalid", "error", "condition")))
  expect_error(guard$run(forged, TRUE, FALSE, 1, .1, 1, 1L, 0L), "forged same-class")
  stale <- NULL
  stale.args <- args
  stale.args$fn <- function(par) {
    if (!is.null(stale)) stop(stale)
    tryCatch(guard$observe(.Machine$double.xmax, par), error = function(e) {
      stale <<- e; stop(e)
    })
  }
  expect_error(guard$run(stale.args, TRUE, FALSE, 1, .1, 1, 1L, 0L), "private invalid first scalar")
  expect_false(guard$active)
  valid <- args; valid$fn <- function(par) {guard$observe(1, par); 1}
  expect_equal(guard$run(valid, TRUE, FALSE, 1, .1, 1, 3L, 0L)$value, 1)
  expect_false(guard$active)
})

test_that("restoration admission depends only on fixed bandwidth and ownership", {
  b <- list(type = "fixed", ckerbound = "none", ckertype = "epanechnikov", ckerorder = 2L)
  eligible <- function(bw = b, owner = TRUE)
    .npindexbw_restore_start_eligible(bw, owner)
  expect_true(eligible())
  expect_false(eligible(owner = FALSE))
  for (type in c("generalized_nn", "adaptive_nn")) {
    alt <- b; alt$type <- type; expect_false(eligible(alt))
  }
  for (bound in c("range", "fixed")) {
    alt <- b; alt$ckerbound <- bound; expect_true(eligible(alt))
  }
  for (kernel in c("gaussian", "beta")) {
    alt <- b; alt$ckertype <- kernel; expect_true(eligible(alt))
  }
  for (order in c(4L, 6L, 8L)) {
    alt <- b; alt$ckerorder <- order; expect_true(eligible(alt))
  }
  alt <- b; alt$ckertype <- "uniform"; alt$ckerorder <- 8L
  expect_true(eligible(alt))
})

test_that("retained restart scale uses the existing projection and U draw", {
  x <- cbind(seq_len(40), sin(seq_len(40)))
  beta <- .75
  fit <- x[,1L]
  controls <- .npindexbw_h_start_controls(.7, .9, .5, 0)
  helper <- .npindex_random_restart_bandwidth
  scope <- new.env(parent = environment(helper))
  environment(helper) <- scope
  project <- .npindex_index_from_beta_tail
  random <- .npindex_random_start_bandwidth
  scale <- .npindex_start_bandwidth_scale
  nproject <- nscale <- 0L
  scope$.npindex_index_from_beta_tail <- function(...) {nproject <<- nproject+1L; project(...)}
  environment(random) <- scope
  scope$.npindex_random_start_bandwidth <- random
  scope$.npindex_start_bandwidth_scale <- function(...) {nscale <<- nscale+1L; scale(...)}
  set.seed(606)
  expected <- runif(1, .7, .9) * scale(as.double(x %*% c(1,beta)), nrow(x))
  after <- .Random.seed
  set.seed(606)
  details <- helper(x, beta, fit, "fixed", nrow(x), controls, 0, retain.scale = TRUE)
  expect_identical(details$h, expected)
  expect_identical(details$scale, scale(as.double(x %*% c(1,beta)), nrow(x)))
  expect_identical(.Random.seed, after)
  expect_identical(nproject, 1L)
  expect_identical(nscale, 1L)
  set.seed(606)
  expect_identical(helper(x, beta, fit, "fixed", nrow(x), controls, 0), expected)
})
