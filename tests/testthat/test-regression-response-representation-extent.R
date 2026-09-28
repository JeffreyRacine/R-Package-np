test_that("response coding refuses missing factor metadata and preserves fitted codes", {
  ns <- asNamespace("np")
  adjust <- get(".np_regression_response_levels", ns)
  numeric.response <- get(".np_plreg_numeric_response", ns)
  hat.response <- get(".np_hat_response", ns)
  ut <- get("untangle", ns)
  numeric.dati <- ut(data.frame(y = c(10, 30, 50)))
  for (ordered in c(FALSE, TRUE)) for (lev in list(c("10", "30", "50"), c("a", "b", "c"))) {
    y <- factor(rep(lev, 8L), levels = lev, ordered = ordered)
    dati <- ut(data.frame(y))
    expected <- rep(if (ordered && lev[1L] == "10") c(10, 30, 50) else 1:3, 8L)
    expect_equal(numeric.response(y, dati), expected, tolerance = 0)
    expect_equal(hat.response(y, dati), expected, tolerance = 0)
    expect_error(adjust(data.frame(y), numeric.dati), "factor response requires")
    expect_error(numeric.response(y, numeric.dati), "factor response requires")
    expect_identical(hat.response(y, numeric.dati), y)
  }
  y <- ordered(c(10, 30, 50))
  dati <- ut(data.frame(y))
  expect_error(adjust(data.frame(y = ordered(c(10, 20))), dati), "unknown factors")
  ev <- adjust(data.frame(y = ordered(c(10, 20))), dati, allowNewCells = TRUE)
  expect_equal(as.vector(get("toMatrix", ns)(ev)), c(10, 20), tolerance = 0)
  payload <- cbind(1:24, (1:24)^2)
  expect_identical(hat.response(payload, numeric.dati), payload)
})

test_that("regression rejects incompatible coding while explicit hat payloads retain their contract", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.05, .95, length.out = 24L))
  y <- rep(c(10, 30, 50), 8L)
  f <- ordered(y)
  b <- npregbw(xdat = x, ydat = y, bws = .25, bandwidth.compute = FALSE)
  expect_error(npreg(b, txdat = x, tydat = f, se = FALSE), "factor response requires")
  expect_error(npreg(b, txdat = x, tydat = y, eydat = f, se = FALSE), "factor response requires")
  H <- npreghat(b, txdat = x)
  expect_equal(as.numeric(npreghat(b, txdat = x, y = f, output = "apply")),
               as.numeric(H %*% as.numeric(f)), tolerance = 1e-11)
  # Cached numeric-fitted hats already refuse factor matrix multiplication;
  # preserve that existing R condition rather than inventing a coding policy.
  expect_error(predict(H, y = f, output = "apply"), "numeric/complex")
  valid <- npreg(b, txdat = x, tydat = y, se = FALSE)
  expect_equal(as.numeric(fitted(valid)), as.numeric(H %*% y), tolerance = 1e-11)
  external <- npreg(b, txdat = x, tydat = y, exdat = x[c(3, 12, 20), , drop = FALSE], se = FALSE)
  expect_length(fitted(external), 3L)
  expect_true(all(is.finite(fitted(external))))
})

test_that("native regression response extents are checked independently of R coding", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  ns <- asNamespace("np")
  x <- data.frame(x = seq(.05, .95, length.out = 24L))
  y <- sin(x$x)
  b <- npregbw(xdat = x, ydat = y, bws = .25, bandwidth.compute = FALSE)
  # Intercept a valid argument builder before native entry. The MPI overrides
  # are capture-only: no computational claim is made from this local seam.
  fn <- get("npreg.rbandwidth", ns)
  env <- new.env(parent = environment(fn))
  environment(fn) <- env
  captured <- NULL
  env$.Call <- function(.NAME, ...) {
    stopifnot(identical(.NAME, "C_np_regression"))
    captured <<- c(list(.NAME), list(...))
    stop(structure(list(message = "captured", call = NULL),
                   class = c("response_capture", "error", "condition")))
  }
  if (identical("np", "npRmpi")) {
    env$.npRmpi_require_active_slave_pool <- function(...) invisible(NULL)
    env$.npRmpi_master_local_entry_needed <- function(...) FALSE
    env$.npRmpi_autodispatch_active <- function(...) FALSE
    env$.npRmpi_guard_no_auto_object_in_manual_bcast <- function(...) invisible(NULL)
    env$.npRmpi_rank_local_regression_context <- function(...) FALSE
    env$mpi.comm.size <- function(...) 1L
  }
  expect_error(fn(b, txdat = x, tydat = y, eydat = y, se = FALSE), class = "response_capture")
  expect_length(captured, 26L)
  short.train <- captured; short.train[[5L]] <- double()
  expect_error(do.call(base::.Call, short.train), "training-response buffer is too short")
  short.eval <- captured; short.eval[[9L]] <- double()
  expect_error(do.call(base::.Call, short.eval), "evaluation-response buffer is too short")
  # Public valid calls test native recovery and legitimate absent ey.
  expect_true(all(is.finite(fitted(npreg(b, txdat = x, tydat = y, se = FALSE)))))
})
