suppressPackageStartupMessages(library(np))
options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)

objective <- function(xdat, y, state, accelerate) {
  options(np.macMseries.accelerate = accelerate)
  as.numeric(np:::.npregbw_eval_only(xdat, y, state)$objective)
}

make_state <- function(n, kernel, order) {
  set.seed(2026072821L + n + order)
  x <- replicate(3L, runif(n, -0.8, 0.8))
  colnames(x) <- paste0("x", seq_len(ncol(x)))
  xdat <- as.data.frame(x)
  y <- sin(1.8 * x[, 1L]) + rowSums(x^2) / 6 +
    rnorm(n, sd = 0.2)
  k <- vapply(
    seq_len(ncol(x)),
    function(j) as.integer(max(12L, round(n / (2.8 + 0.3 * j)))),
    integer(1L)
  )
  state <- npregbw(
    xdat = xdat,
    ydat = y,
    bws = k,
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "adaptive_nn",
    regtype = "lp",
    basis = "glp",
    degree = rep.int(1L, 3L),
    bernstein.basis = FALSE,
    ckertype = kernel,
    ckerorder = order
  )
  list(xdat = xdat, y = y, state = state)
}

active <- make_state(2048L, "gaussian", 2L)
active_scalar <- objective(active$xdat, active$y, active$state, FALSE)
active_blas <- objective(active$xdat, active$y, active$state, TRUE)
stopifnot(
  abs(active_scalar - active_blas) <=
    64 * .Machine$double.eps * max(1, abs(active_scalar), abs(active_blas))
)

threshold <- make_state(2047L, "gaussian", 2L)
stopifnot(identical(
  objective(threshold$xdat, threshold$y, threshold$state, FALSE),
  objective(threshold$xdat, threshold$y, threshold$state, TRUE)
))

signed <- make_state(2048L, "gaussian", 6L)
stopifnot(identical(
  objective(signed$xdat, signed$y, signed$state, FALSE),
  objective(signed$xdat, signed$y, signed$state, TRUE)
))

cat(
  "ADAPTIVE_REGRESSION_FULLROW_BLAS_INSTALLED_PASS ",
  "active_abs=", format(abs(active_scalar - active_blas), digits = 17L),
  " threshold_exact=TRUE signed_exact=TRUE\n",
  sep = ""
)
