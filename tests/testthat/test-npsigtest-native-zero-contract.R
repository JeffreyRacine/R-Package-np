npsig_native_zero_contract <- function(package, expect.zero = TRUE, smoke = FALSE) {
  ns <- asNamespace(package)
  bwfun <- get("npregbw", ns)
  fitfun <- get("npreg", ns)
  tile <- get(".np_npsig_streamed_iid_tile", ns)
  statistic <- get(".np_npsig_statistic", ns)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  attempt <- function(expr) tryCatch(force(expr), error = function(e)
    structure(conditionMessage(e), class = "native_zero_error"))
  failed <- function(value) inherits(value, "native_zero_error")
  call.tile <- function(bw, x, response, pivotal = TRUE)
    tile(bw, x, 1L, response.matrix = response,
         null.mean = response[, 1L], residual.pool = response[, 1L],
         pivotal = pivotal)
  if (smoke) {
    n <- 1024L
    x <- data.frame(x = seq(-1, 1, length.out = n))
    response <- vapply(1:8, function(k)
      sin((k + 1) * x$x) + .2 * cos(seq_len(n) * k), numeric(n))
    bw <- bwfun(xdat = x, ydat = response[, 1L], bws = .45,
                bandwidth.compute = FALSE, regtype = "ll")
    invisible(call.tile(bw, x, response))
    elapsed <- system.time(value <- call.tile(bw, x, response))[["elapsed"]]
    stopifnot(all(is.finite(value)), all(value > 0))
    return(list(input = list(x = x, response = response, bw = .45),
                ordinary = value, elapsed = elapsed))
  }
  n <- 40L
  x <- data.frame(x = seq(-1, 1, length.out = n))
  y <- sin(2 * x$x) + .2 * cos(seq_len(n))
  response <- cbind(y, rev(y), y + .1 * sin(seq_len(n) * 2))
  ordinary <- inputs <- outcomes <- list()
  for (regtype in c("lc", "ll", "lp")) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      key <- paste(regtype, bwtype, sep = "/")
      h <- if (bwtype == "fixed") .45 else 9
      bw <- bwfun(xdat = x, ydat = y, bws = h, bandwidth.compute = FALSE,
        regtype = regtype, degree = if (regtype == "lp") 2L else NULL,
        bwtype = bwtype, ckertype = "gaussian")
      value <- call.tile(bw, x, response)
      scalar <- vapply(seq_len(ncol(response)), function(j) {
        fit <- fitfun(bws = bw, txdat = x, tydat = response[, j],
                      gradients = TRUE, se = TRUE)
        statistic(fit, 1L, TRUE)
      }, numeric(1L))
      tol <- if (bwtype == "adaptive_nn") 5e-7 else 2e-10
      stopifnot(all(is.finite(value)), all(value > 0),
                isTRUE(all.equal(value, scalar, tolerance = tol)))
      raw <- call.tile(bw, x, response, FALSE)
      zero <- attempt(call.tile(bw, x, matrix(0, n, 1L)))
      mixed <- cbind(0, response, 0, response[, 3:1, drop = FALSE])
      mixed.value <- attempt(call.tile(bw, x, mixed))
      permutation <- c(8L, 1L, 7L, 2L, 6L, 3L, 5L, 4L)
      permuted <- attempt(call.tile(bw, x, mixed[, permutation]))
      if (expect.zero) {
        stopifnot(identical(zero, 0), !failed(mixed.value), !failed(permuted),
          identical(mixed.value, c(0, value, 0, rev(value))),
          identical(permuted, mixed.value[permutation]),
          identical(call.tile(bw, x, response[, 1L, drop = FALSE]), value[[1L]]))
      } else {
        stopifnot(failed(zero), failed(mixed.value), failed(permuted))
      }
      tiny.response <- matrix(y * 1e-200, ncol = 1L)
      tiny <- call.tile(bw, x, tiny.response)
      tiny.fit <- fitfun(bws = bw, txdat = x, tydat = tiny.response[, 1L],
                         gradients = TRUE, se = FALSE)
      stopifnot(all(is.finite(tiny.fit$grad)), any(tiny.fit$grad != 0),
                is.finite(tiny), tiny > 0)
      ordinary[[key]] <- list(pivotal = value, raw = raw, tiny = tiny)
      inputs[[key]] <- list(regtype = regtype, bwtype = bwtype, bw = h,
                            x = x, response = response)
      outcomes[[key]] <- list(zero = zero, mixed = mixed.value, permuted = permuted)
    }
  }
  bw <- bwfun(xdat = x, ydat = y, bws = .45,
              bandwidth.compute = FALSE, regtype = "ll")
  bad.payload <- matrix(y, ncol = 1L)
  bad.payload[1L, 1L] <- Inf
  nonfinite <- attempt(call.tile(bw, x, bad.payload))
  oversized <- attempt(call.tile(bw, x, matrix(rep(y, 9L), ncol = 9L)))
  stopifnot(failed(nonfinite), failed(oversized))

  # A genuine zero derivative/SE at an isolated evaluation row does not
  # forgive a nonzero gradient elsewhere in the same response column.
  isolated.x <- data.frame(x = c(seq(-1, 1, length.out = 39L), 10))
  isolated.y <- sin(2 * isolated.x$x) + .2 * cos(seq_len(n))
  isolated.bw <- bwfun(xdat = isolated.x, ydat = isolated.y, bws = .35,
    bandwidth.compute = FALSE, regtype = "lc", ckertype = "epanechnikov")
  isolated.fit <- fitfun(bws = isolated.bw, txdat = isolated.x,
    tydat = isolated.y, gradients = TRUE, se = FALSE)
  stopifnot(all(is.finite(isolated.fit$grad)), any(isolated.fit$grad != 0),
            isolated.fit$grad[n, 1L] == 0)
  isolated <- attempt(call.tile(isolated.bw, isolated.x,
                                matrix(isolated.y, ncol = 1L)))
  isolated.mixed <- attempt(call.tile(isolated.bw, isolated.x,
    cbind(0, isolated.y, 0, isolated.y, 0, isolated.y, 0, isolated.y)))
  stopifnot(failed(isolated), failed(isolated.mixed))

  # Zero responses do not convert an invalid nearest-neighbor radius to zero.
  tied.x <- data.frame(x = c(0, 0, 0, 1:5))
  tied.bw <- bwfun(xdat = tied.x, ydat = numeric(8L), bws = 2,
    bandwidth.compute = FALSE, regtype = "ll", bwtype = "generalized_nn")
  geometry <- attempt(call.tile(tied.bw, tied.x, matrix(0, 8L, 1L)))
  stopifnot(failed(geometry))
  list(input = inputs, ordinary = ordinary, outcomes = outcomes,
       invalid = list(nonfinite = nonfinite, oversized = oversized,
                      isolated = isolated, isolated.mixed = isolated.mixed,
                      geometry = geometry))
}

test_that("native pivotal response columns use complete finite zero effects", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  result <- npsig_native_zero_contract("npRmpi")
  expect_length(result$ordinary, 9L)
  expect_length(result$invalid, 5L)
})
