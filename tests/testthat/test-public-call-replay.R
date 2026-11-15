# Broad producer and MPI mode matrices are retained with the R84 evidence.
test_that("public replay calls preserve syntax and originating package", {
  pkg <- environmentName(environment(npreg))
  normalize <- getFromNamespace(".np_public_replay_call", pkg)
  capture <- function(x, ..., z = NULL)
    normalize(match.call(expand.dots = FALSE), "npreg")
  call <- capture(stop("not forced"), q = NULL, na.action = , z = quote(value()))
  expect_identical(call[[1L]], as.call(list(quote(`::`), as.name(pkg), quote(npreg))))
  expect_false("..." %in% names(call))
  expect_false("na.action" %in% names(call))
  expect_true("q" %in% names(call))
  expect_null(call$q)
  expect_identical(call$x, quote(stop("not forced")))
  expect_identical(call$z, quote(quote(value())))
})

test_that("updates replay public calls without evaluating inert metadata", {
  pkg <- environmentName(environment(npreg))
  if (pkg == "npRmpi") {
    skip_on_cran()
    if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  withr::local_options(np.messages = FALSE)
  d <- data.frame(x = seq(-1, 1, length.out = 30L))
  d$y <- sin(3 * d$x)
  d$z <- cos(seq_len(nrow(d)) / 4)
  count <- new.env(parent = emptyenv())
  count$n <- 0L
  data.value <- function() { count$n <- count$n + 1L; d }
  b <- npregbw(y ~ x, data = data.value(), bws = .4, bandwidth.compute = FALSE)
  set.seed(8421)
  rng <- .Random.seed
  call <- update(b, evaluate = FALSE)
  expect_identical(count$n, 1L)
  expect_identical(.Random.seed, rng)
  expect_identical(call[[1L]], as.call(list(quote(`::`), as.name(pkg), quote(npregbw))))
  again <- update(b)
  expect_identical(count$n, 2L)
  expect_identical(again$bw, b$bw)
  local({
    npregbw <- function(...) stop("masked generic")
    expect_identical(update(b)$bw, b$bw)
  })
  for (name in c("npreg", "npplreg")) {
    fun <- get(name)
    formula <- if (name == "npreg") y ~ x else y ~ x | z
    h <- if (name == "npreg") .4 else matrix(.4, 2L, 1L)
    fit <- fun(formula, data = d, bws = h, bandwidth.compute = FALSE)
    expect_identical(fitted(update(fit)), fitted(fit), info = name)
    expect_identical(getCall(fit)[[1L]], as.call(list(quote(`::`), as.name(pkg), as.name(name))))
  }
  lsq <- nplsqregbw(bws = y ~ x, data = d, bw = .4, scale = rep(1, nrow(d)),
                    delta = .5, regtype = "ll", bandwidth.compute = FALSE)
  fit <- nplsqreg(lsq)
  expect_identical(fitted(update(fit)), fitted(fit))
})
