test_that("ordinary LL index starts recover at fixed beta in both objectives", {
  skip_on_cran()
  spawn_mpi_slaves(1L)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(619L)
  n <- 64L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  x$x1[n] <- 5
  index <- x$x1 + .6*x$x2
  state <- new.env(parent = emptyenv())
  factory <- .npindexbw_first_scalar_guard
  local_mocked_bindings(.npindexbw_first_scalar_guard = function(...) {
    guard <- factory(...)
    observe <- guard$observe
    guard$observe <- function(raw, point) {
      state$raw <- c(state$raw, raw)
      state$points[[length(state$points) + 1L]] <- point
      observe(raw, point)
    }
    guard
  }, .package = "npRmpi")
  for (method in c("ichimura", "kleinspady")) {
    y <- if (method == "ichimura") sin(index) + .1*cos(seq_len(n)) else
      as.double(index + sin(seq_len(n)*2.3) > 0)
    state$raw <- numeric()
    state$points <- list()
    bw <- npindexbw(xdat = x, ydat = y, method = method,
                    regtype = "ll", ckertype = "epanechnikov",
                    optim.method = "BFGS", nmulti = 1L)
    expect_true(is.finite(bw$fval) && bw$fval < .Machine$double.xmax)
    expect_identical(state$raw[1L], .Machine$double.xmax)
    first.valid <- which(state$raw < .Machine$double.xmax)[1L]
    expect_true(!is.na(first.valid) && first.valid > 1L)
    if (!is.na(first.valid) && first.valid > 1L) for (i in 2:first.valid) {
      expect_identical(head(state$points[[i]], -1L), head(state$points[[1L]], -1L))
      expect_identical(tail(state$points[[i]], 1L),
                       2*tail(state$points[[i-1L]], 1L))
    }
  }
})
