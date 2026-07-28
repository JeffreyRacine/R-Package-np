suppressPackageStartupMessages(library(npRmpi))
options(
  np.messages = FALSE,
  np.tree = FALSE,
  np.largeh = FALSE,
  npRmpi.autodispatch = TRUE,
  npRmpi.reuse.slaves = FALSE
)

nslaves <- as.integer(Sys.getenv(
  "NP_RMPI_SENTINEL_NSLAVES",
  unset = "1"
))
npRmpi.init(nslaves = nslaves, quiet = TRUE)
closed <- FALSE
on.exit({
  if (!closed)
    try(npRmpi.quit(force = TRUE), silent = TRUE)
}, add = TRUE)

set_accelerate <- function(value) {
  expr <- substitute(
    options(np.macMseries.accelerate = VALUE),
    list(VALUE = value)
  )
  invisible(npRmpi:::.npRmpi_bcast_cmd_expr(
    expr,
    comm = 1L,
    caller.execute = TRUE
  ))
}

make_case <- function(n, kernel, order, binary = FALSE) {
  set.seed(2026072821L + n + order)
  x <- replicate(3L, runif(n, -0.8, 0.8))
  colnames(x) <- paste0("x", seq_len(ncol(x)))
  xdat <- as.data.frame(x)
  y <- sin(1.8 * x[, 1L]) + rowSums(x^2) / 6 +
    rnorm(n, sd = 0.2)
  if (isTRUE(binary))
    y <- as.double(y > median(y))
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

objective_pair <- function(case) {
  shadow <- npRmpi:::.npregbw_nomad_shadow_begin(
    xdat = case$xdat,
    ydat = case$y,
    bws = case$state,
    broadcast = TRUE
  )
  on.exit(npRmpi:::.npregbw_nomad_shadow_end(
    shadow,
    broadcast = TRUE
  ), add = TRUE)

  set_accelerate(FALSE)
  scalar <- npRmpi:::.npregbw_nomad_shadow_eval(
    shadow,
    bw = case$state$bw,
    degree = case$state$degree,
    broadcast = TRUE
  )[1L]
  set_accelerate(TRUE)
  blas <- npRmpi:::.npregbw_nomad_shadow_eval(
    shadow,
    bw = case$state$bw,
    degree = case$state$degree,
    broadcast = TRUE
  )[1L]
  c(scalar = scalar, blas = blas)
}

active <- objective_pair(make_case(2048L, "gaussian", 2L))
stopifnot(
  abs(active[["scalar"]] - active[["blas"]]) <=
    64 * .Machine$double.eps *
      max(1, abs(active[["scalar"]]), abs(active[["blas"]]))
)

threshold <- objective_pair(make_case(2047L, "gaussian", 2L))
stopifnot(identical(
  unname(threshold[["scalar"]]),
  unname(threshold[["blas"]])
))

signed <- objective_pair(make_case(2048L, "gaussian", 6L))
stopifnot(identical(
  unname(signed[["scalar"]]),
  unname(signed[["blas"]])
))

ks_pair <- function(case) {
  expr <- substitute(
    npRmpi:::.npregbw_eval_only(
      XDAT, YDAT, STATE, localize = FALSE, objective = "ks"
    )$objective,
    list(XDAT = case$xdat, YDAT = case$y, STATE = case$state)
  )
  set_accelerate(FALSE)
  scalar <- npRmpi:::.npRmpi_bcast_cmd_expr(
    expr, comm = 1L, caller.execute = TRUE
  )
  set_accelerate(TRUE)
  blas <- npRmpi:::.npRmpi_bcast_cmd_expr(
    expr, comm = 1L, caller.execute = TRUE
  )
  c(scalar = scalar, blas = blas)
}
ks <- ks_pair(make_case(2048L, "gaussian", 2L, binary = TRUE))
stopifnot(
  abs(ks[["scalar"]] - ks[["blas"]]) <=
    64 * .Machine$double.eps *
      max(1, abs(ks[["scalar"]]), abs(ks[["blas"]]))
)

check_pair <- function(case) {
  evaluate <- function(accelerate) {
    set_accelerate(accelerate)
    nplsqregbw(
      xdat = case$xdat,
      ydat = case$y,
      bws = case$state$bw,
      bandwidth.compute = FALSE,
      scale = rep.int(1, nrow(case$xdat)),
      tau = 0.35,
      delta = 0.4,
      bwtype = "adaptive_nn",
      regtype = "lp",
      basis = "glp",
      degree = rep.int(1L, ncol(case$xdat)),
      bernstein.basis = FALSE
    )$fval
  }
  c(scalar = evaluate(FALSE), blas = evaluate(TRUE))
}
check <- check_pair(make_case(2048L, "gaussian", 2L))
stopifnot(
  abs(check[["scalar"]] - check[["blas"]]) <=
    64 * .Machine$double.eps *
      max(1, abs(check[["scalar"]]), abs(check[["blas"]]))
)

info <- npRmpi.session.info()
cat(
  "ADAPTIVE_REGRESSION_FULLROW_BLAS_MPI_INSTALLED_PASS ",
  "comm_size=", info$comm_size,
  " active_abs=",
  format(abs(active[["scalar"]] - active[["blas"]]), digits = 17L),
  " ks_abs=",
  format(abs(ks[["scalar"]] - ks[["blas"]]), digits = 17L),
  " check_abs=",
  format(abs(check[["scalar"]] - check[["blas"]]), digits = 17L),
  " threshold_exact=TRUE signed_exact=TRUE\n",
  sep = ""
)
npRmpi.quit(force = TRUE)
closed <- TRUE
