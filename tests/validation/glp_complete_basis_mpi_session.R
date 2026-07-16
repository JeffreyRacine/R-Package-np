args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(paste(
    "usage:",
    "glp_complete_basis_mpi_session.R <library> serial <result.rds>",
    "or <library> mpi <result.rds> <nslaves>",
    "or <library> compare <summary.rds> <serial.rds> <mpi1.rds> <mpi3.rds>"
  ))
}

lib <- args[1L]
mode <- args[2L]
output <- args[3L]
.libPaths(c(lib, .libPaths()))

make_fixture <- function() {
  set.seed(2026071608L)
  n <- as.integer(Sys.getenv("NP_GLP_MPI_N", "600"))
  tx <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- 0.5 + tx$x1 - 0.4 * tx$x1^2 + 0.6 * tx$x2 +
    0.8 * tx$x1 * tx$x2 + 0.03 * sin(seq_len(n) * 0.9)
  ex <- data.frame(
    x1 = seq(-0.8, 0.8, length.out = 36L),
    x2 = seq(-0.7, 0.7, length.out = 36L)
  )
  list(n = n, tx = tx, y = y, ex = ex)
}

fit_fixture <- function(package) {
  fixture <- make_fixture()
  bw <- npregbw(
    xdat = fixture$tx,
    ydat = fixture$y,
    regtype = "lp",
    degree = c(2L, 1L),
    degree.select = "manual",
    basis = "glp",
    bernstein.basis = TRUE,
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE
  )
  fit <- npreg(
    txdat = fixture$tx,
    tydat = fixture$y,
    exdat = fixture$ex,
    bws = bw,
    gradients = TRUE,
    gradient.order = c(2L, 1L)
  )
  list(
    package = package,
    version = as.character(packageVersion(package)),
    n = fixture$n,
    degree = c(2L, 1L),
    dimension = get("npLpBasisNcol", asNamespace(package))("glp", c(2L, 1L)),
    mean = fit$mean,
    grad = fit$grad,
    merr = fit$merr,
    gerr = fit$gerr
  )
}

if (identical(mode, "serial")) {
  suppressPackageStartupMessages(library(np))
  result <- fit_fixture("np")
  saveRDS(result, output, version = 3L)
  cat("PASS: serial np complete-GLP reference written.\n")
} else if (identical(mode, "mpi")) {
  if (length(args) < 4L)
    stop("mpi mode requires nslaves")
  nslaves <- as.integer(args[4L])
  if (is.na(nslaves) || !(nslaves %in% c(1L, 3L)))
    stop("nslaves must be 1 or 3")

  suppressPackageStartupMessages(library(npRmpi))
  options(
    np.messages = FALSE,
    np.tree = FALSE,
    npRmpi.autodispatch = TRUE,
    npRmpi.reuse.slaves = FALSE
  )
  npRmpi.init(nslaves = nslaves, quiet = TRUE)
  closed <- FALSE
  on.exit({
    if (!closed)
      try(npRmpi.quit(force = TRUE), silent = TRUE)
  }, add = TRUE)

  info <- npRmpi.session.info()
  result <- fit_fixture("npRmpi")
  result$nslaves <- nslaves
  result$comm_size <- info$comm_size

  # Write substantive results before teardown. A known teardown exit 137 is
  # classifiable only when this file is complete and no route-scoped workers
  # remain after the subprocess exits.
  saveRDS(result, output, version = 3L)
  npRmpi.quit(force = TRUE)
  closed <- TRUE
  cat(sprintf("PASS: npRmpi complete-GLP session result written (nslaves=%d).\n",
              nslaves))
} else if (identical(mode, "compare")) {
  if (length(args) < 6L)
    stop("compare mode requires serial, mpi1, and mpi3 result paths")
  serial <- readRDS(args[4L])
  mpi1 <- readRDS(args[5L])
  mpi3 <- readRDS(args[6L])

  compare_one <- function(candidate) {
    stopifnot(
      identical(as.integer(candidate$degree), as.integer(serial$degree)),
      as.numeric(candidate$dimension) == as.numeric(serial$dimension),
      candidate$n == serial$n
    )
    c(
      mean = max(abs(candidate$mean - serial$mean)),
      grad = max(abs(candidate$grad - serial$grad)),
      merr = max(abs(candidate$merr - serial$merr)),
      gerr = max(abs(candidate$gerr - serial$gerr))
    )
  }

  errors1 <- compare_one(mpi1)
  errors3 <- compare_one(mpi3)
  stopifnot(
    errors1[["mean"]] < 1e-9,
    errors1[["grad"]] < 1e-8,
    errors1[["merr"]] < 1e-9,
    errors1[["gerr"]] < 1e-8,
    errors3[["mean"]] < 1e-9,
    errors3[["grad"]] < 1e-8,
    errors3[["merr"]] < 1e-9,
    errors3[["gerr"]] < 1e-8
  )

  summary <- list(serial = serial[c("package", "version", "n", "degree", "dimension")],
                  mpi1 = errors1, mpi3 = errors3)
  saveRDS(summary, output, version = 3L)
  cat("PASS: installed np/npRmpi complete-GLP session parity.\n")
} else {
  stop("mode must be serial, mpi, or compare")
}
