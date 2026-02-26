## the idea is that you have done bandwidth selection
## you just need to supply training data and the bandwidth
## this tool will help you visualize the result

.np_seed_enter <- function(random.seed = 42L) {
  save.seed <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    save.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    exists.seed <- TRUE
  } else {
    exists.seed <- FALSE
  }

  set.seed(random.seed)
  list(exists.seed = exists.seed, save.seed = save.seed)
}

.np_seed_exit <- function(state, remove_if_absent = FALSE) {
  if (isTRUE(state$exists.seed)) {
    assign(".Random.seed", state$save.seed, envir = .GlobalEnv)
  } else if (isTRUE(remove_if_absent) &&
             exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

.np_with_seed <- function(random.seed = 42L, code) {
  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state), add = TRUE)
  force(code)
}

.np_mammen_draws <- function(n, B) {
  a <- (1 - sqrt(5)) / 2
  p.a <- (sqrt(5) + 1) / (2 * sqrt(5))
  u <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  out <- matrix(1 - a, nrow = n, ncol = B)
  out[u <= p.a] <- a
  out
}

.np_rademacher_draws <- function(n, B) {
  u <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  out <- matrix(1.0, nrow = n, ncol = B)
  out[u <= 0.5] <- -1.0
  out
}

.np_wild_draws <- function(n, B, wild = c("mammen", "rademacher")) {
  if (length(wild) > 1L)
    wild <- wild[1L]
  wild <- match.arg(wild, c("mammen", "rademacher"))
  if (identical(wild, "mammen")) {
    return(.np_mammen_draws(n = n, B = B))
  }
  .np_rademacher_draws(n = n, B = B)
}

.np_wild_chunk_size <- function(n, B) {
  chunk.opt <- getOption("np.plot.wild.chunk.size")
  if (!is.null(chunk.opt)) {
    chunk.opt <- as.integer(chunk.opt)
    if (length(chunk.opt) != 1L || is.na(chunk.opt) || chunk.opt < 1L)
      stop("option 'np.plot.wild.chunk.size' must be a positive integer")
    return(min(B, chunk.opt))
  }

  if (n < 1L || B < 1L)
    return(1L)

  # Keep the temporary n x chunk response matrix in a moderate memory range.
  target.bytes <- 64 * 1024 * 1024
  chunk <- as.integer(floor(target.bytes / (8 * n)))
  if (!is.finite(chunk) || is.na(chunk) || chunk < 1L)
    chunk <- 1L
  min(B, chunk)
}

.np_wild_boot_t <- function(H, fit.mean, residuals, B, wild = c("mammen", "rademacher")) {
  B <- as.integer(B)
  n <- length(residuals)
  if (length(fit.mean) != n)
    stop("length mismatch between fitted means and residuals for wild bootstrap")
  if (B < 1L)
    stop("argument 'plot.errors.boot.num' must be a positive integer")

  t.mpi <- .npRmpi_wild_boot_t_parallel(
    H = H,
    fit.mean = fit.mean,
    residuals = residuals,
    B = B,
    wild = wild,
    comm = 1L
  )
  if (is.matrix(t.mpi))
    return(t.mpi)

  chunk.size <- .np_wild_chunk_size(n = n, B = B)
  out <- matrix(NA_real_, nrow = B, ncol = nrow(H))
  fit.mean <- as.double(fit.mean)
  residuals <- as.double(residuals)

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    draws <- .np_wild_draws(n = n, B = bsz, wild = wild)
    ystar <- matrix(fit.mean, nrow = n, ncol = bsz) +
      matrix(residuals, nrow = n, ncol = bsz) * draws
    out[start:stopi, ] <- t(H %*% ystar)
    start <- stopi + 1L
  }

  out
}

.np_plot_is_wild_method <- function(method) {
  isTRUE(length(method) == 1L && !is.na(method) && method == "wild")
}

.np_plot_reject_wild_unsupervised <- function(method, where) {
  if (.np_plot_is_wild_method(method)) {
    stop(sprintf("plot.errors.boot.method='wild' is not supported for %s; use one of 'inid', 'fixed', or 'geom'", where))
  }
  invisible(NULL)
}

.np_plot_inid_fastpath_enabled <- function() {
  !isTRUE(getOption("np.plot.inid.fastpath.disable", FALSE))
}

.np_plot_require_bws <- function(bws, where) {
  if (is.null(bws))
    stop(sprintf("required argument 'bws' is missing or NULL in %s", where))
  invisible(TRUE)
}

.npRmpi_plot_inid_ksum_fastpath_enabled <- function() {
  if (isFALSE(getOption("np.plot.inid.ksum.fastpath.nprmpi", FALSE)))
    return(FALSE)
  TRUE
}

.npRmpi_bootstrap_worker_count <- function(comm = 1L) {
  size <- tryCatch(as.integer(mpi.comm.size(comm = comm)), error = function(e) NA_integer_)
  if (is.na(size) || size <= 1L)
    return(0L)
  size - 1L
}

.npRmpi_bootstrap_fanout_enabled <- function(comm = 1L) {
  # Experimental gate: session-mode daemon execution has unresolved hangs
  # for apply-style fan-out in this environment.
  if (!isTRUE(getOption("np.plot.bootstrap.mpi.experimental", FALSE)))
    return(FALSE)
  if (isTRUE(getOption("np.plot.bootstrap.mpi.disable", FALSE)))
    return(FALSE)
  if (isTRUE(.npRmpi_autodispatch_called_from_bcast()))
    return(FALSE)
  if (!isTRUE(.npRmpi_has_active_slave_pool(comm = comm)))
    return(FALSE)
  .npRmpi_bootstrap_worker_count(comm = comm) >= 1L
}

.npRmpi_bootstrap_chunk_tasks <- function(B, chunk.size) {
  B <- as.integer(B)
  chunk.size <- as.integer(chunk.size)
  if (B < 1L || chunk.size < 1L)
    stop("invalid chunk configuration")

  starts <- seq.int(1L, B, by = chunk.size)
  lens <- pmin(chunk.size, B - starts + 1L)
  seeds <- sample.int(.Machine$integer.max, length(starts))

  lapply(seq_along(starts), function(i) {
    list(
      start = as.integer(starts[i]),
      bsz = as.integer(lens[i]),
      seed = as.integer(seeds[i])
    )
  })
}

.npRmpi_bootstrap_collect_chunks <- function(parts, tasks, ncol.out, what = "bootstrap") {
  if (!is.list(parts) || length(parts) != length(tasks)) {
    warning(sprintf("MPI %s fan-out returned malformed chunk results; using local path", what))
    return(NULL)
  }

  has.try.error <- vapply(parts, function(x) inherits(x, "try-error"), logical(1))
  if (any(has.try.error)) {
    warning(sprintf("MPI %s fan-out worker error detected; using local path", what))
    return(NULL)
  }

  total.rows <- sum(vapply(tasks, function(tt) as.integer(tt$bsz), integer(1)))
  out <- matrix(NA_real_, nrow = total.rows, ncol = as.integer(ncol.out))
  rowi <- 1L

  for (i in seq_along(parts)) {
    bsz <- as.integer(tasks[[i]]$bsz)
    chunk <- parts[[i]]
    if (!is.matrix(chunk))
      chunk <- as.matrix(chunk)

    if (!identical(dim(chunk), c(bsz, as.integer(ncol.out)))) {
      if (length(chunk) != (bsz * as.integer(ncol.out))) {
        warning(sprintf("MPI %s fan-out chunk dimension mismatch; using local path", what))
        return(NULL)
      }
      chunk <- matrix(as.numeric(chunk), nrow = bsz, ncol = as.integer(ncol.out))
    }

    out[rowi:(rowi + bsz - 1L), ] <- chunk
    rowi <- rowi + bsz
  }

  out
}

.npRmpi_wild_boot_t_parallel <- function(H, fit.mean, residuals, B, wild, comm = 1L) {
  if (isTRUE(getOption("np.plot.wild.mpi.parallel.disable", FALSE)))
    return(NULL)
  if (!.npRmpi_bootstrap_fanout_enabled(comm = comm))
    return(NULL)

  n <- length(residuals)
  p <- nrow(H)
  chunk.size <- .np_wild_chunk_size(n = n, B = B)
  tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
  if (!length(tasks))
    return(NULL)

  H <- as.matrix(H)
  fit.mean <- as.double(fit.mean)
  residuals <- as.double(residuals)
  wild <- match.arg(if (length(wild) > 1L) wild[1L] else wild,
                    c("mammen", "rademacher"))

  worker <- function(task) {
    set.seed(as.integer(task$seed))
    bsz <- as.integer(task$bsz)
    draws <- .np_wild_draws(n = n, B = bsz, wild = wild)
    ystar <- matrix(fit.mean, nrow = n, ncol = bsz) +
      matrix(residuals, nrow = n, ncol = bsz) * draws
    t(H %*% ystar)
  }

  t.comm <- proc.time()
  parts <- tryCatch(mpi.applyLB(tasks, worker, comm = comm), error = function(e) e)
  .npRmpi_profile_add_comm_elapsed(
    elapsed_sec = unname(as.double((proc.time() - t.comm)[["elapsed"]])),
    where = "mpi.applyLB:wild"
  )
  if (inherits(parts, "error")) {
    warning(sprintf("MPI wild bootstrap fan-out failed (%s); using local path",
                    conditionMessage(parts)))
    return(NULL)
  }

  .npRmpi_bootstrap_collect_chunks(parts = parts, tasks = tasks,
                                   ncol.out = p, what = "wild")
}

.np_inid_chunk_size <- function(n, B) {
  chunk.opt <- getOption("np.plot.inid.chunk.size")
  if (!is.null(chunk.opt)) {
    chunk.opt <- as.integer(chunk.opt)
    if (length(chunk.opt) != 1L || is.na(chunk.opt) || chunk.opt < 1L)
      stop("option 'np.plot.inid.chunk.size' must be a positive integer")
    return(min(B, chunk.opt))
  }

  if (n < 1L || B < 1L)
    return(1L)

  target.bytes <- 64 * 1024 * 1024
  chunk <- as.integer(floor(target.bytes / (8 * n)))
  if (!is.finite(chunk) || is.na(chunk) || chunk < 1L)
    chunk <- 1L
  min(B, chunk)
}

.np_inid_counts_matrix <- function(n, B, counts = NULL) {
  n <- as.integer(n)
  B <- as.integer(B)
  if (n < 1L || B < 1L)
    stop("invalid inid bootstrap dimensions")

  if (!is.null(counts)) {
    counts <- as.matrix(counts)
    if (!is.numeric(counts) ||
        nrow(counts) != n ||
        ncol(counts) != B)
      stop("counts must be an n x B numeric matrix")
    return(counts)
  }

  stats::rmultinom(n = B, size = n, prob = rep.int(1 / n, n))
}

.np_inid_lc_boot_from_hat <- function(H, ydat, B, counts = NULL) {
  H <- as.matrix(H)
  ydat <- as.double(ydat)
  B <- as.integer(B)
  n <- length(ydat)
  if (B < 1L)
    stop("argument 'plot.errors.boot.num' must be a positive integer")
  if (ncol(H) != n)
    stop("hat matrix columns must match length of ydat")

  W <- t(H)
  Wy <- W * ydat
  t0 <- as.vector(H %*% ydat)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    den <- crossprod(counts.mat, W)
    num <- crossprod(counts.mat, Wy)
    return(list(
      t = num / pmax(den, .Machine$double.eps),
      t0 = t0
    ))
  }

  if (!isTRUE(getOption("np.plot.inid.mpi.parallel.disable", FALSE)) &&
      .npRmpi_bootstrap_fanout_enabled(comm = 1L)) {
    chunk.size <- .np_inid_chunk_size(n = n, B = B)
    tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
    prob <- rep.int(1 / n, n)

    worker <- function(task) {
      set.seed(as.integer(task$seed))
      bsz <- as.integer(task$bsz)
      counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
      den <- crossprod(counts.chunk, W)
      num <- crossprod(counts.chunk, Wy)
      num / pmax(den, .Machine$double.eps)
    }

    t.comm <- proc.time()
    parts <- tryCatch(mpi.applyLB(tasks, worker, comm = 1L), error = function(e) e)
    .npRmpi_profile_add_comm_elapsed(
      elapsed_sec = unname(as.double((proc.time() - t.comm)[["elapsed"]])),
      where = "mpi.applyLB:inid"
    )
    if (!inherits(parts, "error")) {
      t.mpi <- .npRmpi_bootstrap_collect_chunks(parts = parts,
                                                tasks = tasks,
                                                ncol.out = nrow(H),
                                                what = "inid")
      if (is.matrix(t.mpi))
        return(list(t = t.mpi, t0 = t0))
    } else {
      warning(sprintf("MPI inid bootstrap fan-out failed (%s); using local path",
                      conditionMessage(parts)))
    }
  }

  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  prob <- rep.int(1 / n, n)
  tmat <- matrix(NA_real_, nrow = B, ncol = nrow(H))

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
    den <- crossprod(counts.chunk, W)
    num <- crossprod(counts.chunk, Wy)
    tmat[start:stopi, ] <- num / pmax(den, .Machine$double.eps)
    start <- stopi + 1L
  }

  list(t = tmat, t0 = t0)
}

.np_inid_lp_unpack_sym_row <- function(mrow, p) {
  A <- matrix(0.0, nrow = p, ncol = p)
  idx <- 1L
  for (a in seq_len(p)) {
    for (b in a:p) {
      A[a, b] <- mrow[idx]
      A[b, a] <- mrow[idx]
      idx <- idx + 1L
    }
  }
  A
}

.np_inid_lp_solver_backend <- function() {
  backend <- getOption("np.plot.inid.lp.solver", "solve")
  backend <- as.character(backend)[1L]
  if (is.na(backend) || !(backend %in% c("auto", "chol", "solve", "qr")))
    backend <- "auto"
  backend
}

.np_inid_lp_solve_once <- function(A, z, backend) {
  if (identical(backend, "chol")) {
    R <- tryCatch(chol(A), error = function(e) NULL)
    if (!is.null(R))
      return(backsolve(R, forwardsolve(t(R), z)))
    return(NULL)
  }
  if (identical(backend, "solve"))
    return(tryCatch(solve(A, z), error = function(e) NULL))
  if (identical(backend, "qr"))
    return(tryCatch(qr.solve(A, z, tol = .Machine$double.eps), error = function(e) NULL))

  R <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(R)) {
    beta <- tryCatch(backsolve(R, forwardsolve(t(R), z)), error = function(e) NULL)
    if (!is.null(beta) && all(is.finite(beta)))
      return(beta)
  }

  beta <- tryCatch(solve(A, z), error = function(e) NULL)
  if (!is.null(beta) && all(is.finite(beta)))
    return(beta)

  tryCatch(qr.solve(A, z, tol = .Machine$double.eps), error = function(e) NULL)
}

.np_inid_lp_predict_row <- function(A, z, rhs, ridge.base = 1.0e-12) {
  ridge <- max(0, as.double(ridge.base))
  backend <- .np_inid_lp_solver_backend()

  for (attempt in 0:8) {
    Ar <- A
    if (ridge > 0)
      diag(Ar) <- diag(Ar) + ridge

    beta <- .np_inid_lp_solve_once(A = Ar, z = z, backend = backend)
    if (!is.null(beta) && all(is.finite(beta)))
      return(sum(rhs * beta))

    ridge <- if (ridge > 0) ridge * 10 else 1.0e-12
  }

  NA_real_
}

.np_inid_lp_predict_chunk_general <- function(Mvals, Zvals, rhs, ridge.base = 1.0e-12) {
  Mvals <- as.matrix(Mvals)
  Zvals <- as.matrix(Zvals)
  rhs <- as.double(rhs)

  bsz <- nrow(Mvals)
  p <- ncol(Zvals)
  out <- numeric(bsz)

  for (ii in seq_len(bsz)) {
    A <- .np_inid_lp_unpack_sym_row(mrow = Mvals[ii, ], p = p)
    out[ii] <- .np_inid_lp_predict_row(
      A = A,
      z = as.double(Zvals[ii, ]),
      rhs = rhs,
      ridge.base = ridge.base
    )
  }

  out
}

.np_inid_lp_predict_chunk <- function(Mvals, Zvals, rhs, ridge.base = 1.0e-12) {
  Mvals <- as.matrix(Mvals)
  Zvals <- as.matrix(Zvals)
  rhs <- as.double(rhs)

  bsz <- nrow(Mvals)
  p <- ncol(Zvals)
  out <- rep(NA_real_, bsz)

  if (p == 1L) {
    den <- as.double(Mvals[, 1L])
    out <- rhs[1L] * as.double(Zvals[, 1L]) / pmax(den, .Machine$double.eps)
    return(out)
  }

  if (p == 2L) {
    a <- as.double(Mvals[, 1L])
    b <- as.double(Mvals[, 2L])
    c <- as.double(Mvals[, 3L])
    u <- as.double(Zvals[, 1L])
    v <- as.double(Zvals[, 2L])

    det <- a * c - b * b
    good <- is.finite(det) & (abs(det) > .Machine$double.eps)
    if (any(good)) {
      invdet <- 1 / det[good]
      beta1 <- (c[good] * u[good] - b[good] * v[good]) * invdet
      beta2 <- (a[good] * v[good] - b[good] * u[good]) * invdet
      out[good] <- rhs[1L] * beta1 + rhs[2L] * beta2
    }
  } else if (p == 3L) {
    a <- as.double(Mvals[, 1L])
    b <- as.double(Mvals[, 2L])
    c <- as.double(Mvals[, 3L])
    d <- as.double(Mvals[, 4L])
    e <- as.double(Mvals[, 5L])
    f <- as.double(Mvals[, 6L])
    u <- as.double(Zvals[, 1L])
    v <- as.double(Zvals[, 2L])
    w <- as.double(Zvals[, 3L])

    det <- a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d)
    good <- is.finite(det) & (abs(det) > .Machine$double.eps)
    if (any(good)) {
      c11 <- d[good] * f[good] - e[good] * e[good]
      c12 <- c[good] * e[good] - b[good] * f[good]
      c13 <- b[good] * e[good] - c[good] * d[good]
      c22 <- a[good] * f[good] - c[good] * c[good]
      c23 <- b[good] * c[good] - a[good] * e[good]
      c33 <- a[good] * d[good] - b[good] * b[good]
      invdet <- 1 / det[good]

      beta1 <- (c11 * u[good] + c12 * v[good] + c13 * w[good]) * invdet
      beta2 <- (c12 * u[good] + c22 * v[good] + c23 * w[good]) * invdet
      beta3 <- (c13 * u[good] + c23 * v[good] + c33 * w[good]) * invdet
      out[good] <- rhs[1L] * beta1 + rhs[2L] * beta2 + rhs[3L] * beta3
    }
  }

  bad <- which(!is.finite(out))
  if (length(bad)) {
    p <- ncol(Zvals)
    for (ii in bad) {
      A <- .np_inid_lp_unpack_sym_row(mrow = Mvals[ii, ], p = p)
      out[ii] <- .np_inid_lp_predict_row(
        A = A,
        z = as.double(Zvals[ii, ]),
        rhs = rhs,
        ridge.base = ridge.base
      )
    }
  }

  out
}

.np_inid_boot_from_regression <- function(xdat,
                                          exdat,
                                          bws,
                                          ydat,
                                          B,
                                          counts = NULL,
                                          ridge = 1.0e-12) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  ydat <- as.double(ydat)
  B <- as.integer(B)

  n <- nrow(xdat)
  neval <- nrow(exdat)
  if (length(ydat) != n)
    stop("length of ydat must match training rows")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid inid regression bootstrap dimensions")

  regtype <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  ncon <- bws$ncon

  degree <- if (identical(regtype, "lc")) {
    rep.int(0L, ncon)
  } else if (identical(regtype, "ll")) {
    rep.int(1L, ncon)
  } else {
    npValidateGlpDegree(
      regtype = "lp",
      degree = bws$degree,
      ncon = ncon
    )
  }

  basis <- npValidateLpBasis(
    regtype = "lp",
    basis = if (is.null(bws$basis)) "glp" else bws$basis
  )
  bernstein.basis <- npValidateGlpBernstein(
    regtype = "lp",
    bernstein.basis = isTRUE(bws$bernstein.basis)
  )

  npksum.fun <- .npRmpi_bootstrap_estimator("npksum.default")
  kw <- npksum.fun(
    txdat = xdat,
    exdat = exdat,
    bws = bws,
    return.kernel.weights = TRUE,
    bandwidth.divide = TRUE,
    leave.one.out = FALSE
  )$kw
  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = n)
  if (nrow(kw) != n || ncol(kw) != neval)
    stop("kernel-weight matrix shape mismatch")

  W <- W.lp(
    xdat = xdat,
    degree = degree,
    basis = basis,
    bernstein.basis = bernstein.basis
  )
  W.eval <- W.lp(
    xdat = xdat,
    exdat = exdat,
    degree = degree,
    basis = basis,
    bernstein.basis = bernstein.basis
  )
  W <- as.matrix(W)
  W.eval <- as.matrix(W.eval)

  if (nrow(W) != n || nrow(W.eval) != neval || ncol(W.eval) != ncol(W))
    stop("regression moment design matrix shape mismatch")

  p <- ncol(W)
  mcols <- p * (p + 1L) / 2L
  rhs <- W.eval
  ones <- matrix(1.0, nrow = n, ncol = 1L)

  Mfeat <- vector("list", neval)
  Zfeat <- vector("list", neval)
  t0 <- numeric(neval)

  for (i in seq_len(neval)) {
    k <- as.double(kw[, i])
    WK <- W * k
    Zfeat[[i]] <- WK * ydat

    mf <- matrix(0.0, nrow = n, ncol = mcols)
    idx <- 1L
    for (a in seq_len(p)) {
      for (b in a:p) {
        mf[, idx] <- WK[, a] * W[, b]
        idx <- idx + 1L
      }
    }
    Mfeat[[i]] <- mf

    M0 <- crossprod(ones, mf)
    Z0 <- crossprod(ones, Zfeat[[i]])
    t0[i] <- if (p > 3L) {
      .np_inid_lp_predict_chunk_general(
        Mvals = M0,
        Zvals = Z0,
        rhs = rhs[i, ],
        ridge.base = ridge
      )[1L]
    } else {
      .np_inid_lp_predict_chunk(
        Mvals = M0,
        Zvals = Z0,
        rhs = rhs[i, ],
        ridge.base = ridge
      )[1L]
    }
  }

  tmat <- matrix(NA_real_, nrow = B, ncol = neval)

  fill_chunk <- function(counts.chunk, start, stopi) {
    for (i in seq_len(neval)) {
      Mvals <- crossprod(counts.chunk, Mfeat[[i]])
      Zvals <- crossprod(counts.chunk, Zfeat[[i]])
      tmat[start:stopi, i] <<- if (p > 3L) {
        .np_inid_lp_predict_chunk_general(
          Mvals = Mvals,
          Zvals = Zvals,
          rhs = rhs[i, ],
          ridge.base = ridge
        )
      } else {
        .np_inid_lp_predict_chunk(
          Mvals = Mvals,
          Zvals = Zvals,
          rhs = rhs[i, ],
          ridge.base = ridge
        )
      }
    }
  }

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    fill_chunk(counts.chunk = counts.mat, start = 1L, stopi = B)
  } else {
    chunk.size <- .np_inid_chunk_size(n = n, B = B)
    prob <- rep.int(1 / n, n)
    start <- 1L
    while (start <= B) {
      stopi <- min(B, start + chunk.size - 1L)
      bsz <- stopi - start + 1L
      counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
      fill_chunk(counts.chunk = counts.chunk, start = start, stopi = stopi)
      start <- stopi + 1L
    }
  }

  if (any(!is.finite(t0)) || any(!is.finite(tmat)))
    stop("inid regression fast path produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_inid_scoef_numeric_y <- function(ydat, bws) {
  if (is.factor(ydat)) {
    if (is.null(bws$ydati))
      stop("factor response requires bws$ydati for smooth coefficient inid helper")
    yadj <- adjustLevels(data.frame(ydat), bws$ydati)
    return((bws$ydati$all.dlev[[1L]])[as.integer(yadj[, 1L])])
  }
  as.double(ydat)
}

.np_inid_scoef_predict_row <- function(mrow, zrow, rhs, epsilon) {
  A <- .np_inid_lp_unpack_sym_row(mrow = mrow, p = length(rhs))
  tyw <- as.double(zrow)
  nc <- ncol(A)

  maxPenalty <- sqrt(.Machine$double.xmax)
  coef.ii <- rep(maxPenalty, nc)
  ridge <- 0.0
  doridge <- TRUE

  while (doridge) {
    doridge <- FALSE
    ridge.val <- ridge * tyw[1L] / NZD(A[1L, 1L])
    coef.try <- tryCatch(
      solve(
        A + diag(rep(ridge, nc)),
        tyw + c(ridge.val, rep(0, nc - 1L))
      ),
      error = function(e) e
    )
    if (inherits(coef.try, "error") || any(!is.finite(coef.try))) {
      ridge <- ridge + epsilon
      doridge <- TRUE
      coef.try <- rep(maxPenalty, nc)
    }
    coef.ii <- as.double(coef.try)
  }

  sum(as.double(rhs) * coef.ii)
}

.np_inid_scoef_predict_chunk <- function(Mvals, Zvals, rhs) {
  Mvals <- as.matrix(Mvals)
  Zvals <- as.matrix(Zvals)
  rhs <- as.double(rhs)

  bsz <- nrow(Mvals)
  p <- ncol(Zvals)
  out <- rep(NA_real_, bsz)

  if (p == 1L) {
    den <- as.double(Mvals[, 1L])
    good <- is.finite(den) & (abs(den) > .Machine$double.eps)
    out[good] <- rhs[1L] * as.double(Zvals[good, 1L]) / den[good]
    return(out)
  }

  if (p == 2L) {
    a <- as.double(Mvals[, 1L])
    b <- as.double(Mvals[, 2L])
    c <- as.double(Mvals[, 3L])
    u <- as.double(Zvals[, 1L])
    v <- as.double(Zvals[, 2L])
    det <- a * c - b * b
    good <- is.finite(det) & (abs(det) > .Machine$double.eps)
    if (any(good)) {
      invdet <- 1 / det[good]
      beta1 <- (c[good] * u[good] - b[good] * v[good]) * invdet
      beta2 <- (a[good] * v[good] - b[good] * u[good]) * invdet
      out[good] <- rhs[1L] * beta1 + rhs[2L] * beta2
    }
    return(out)
  }

  if (p == 3L) {
    a <- as.double(Mvals[, 1L])
    b <- as.double(Mvals[, 2L])
    c <- as.double(Mvals[, 3L])
    d <- as.double(Mvals[, 4L])
    e <- as.double(Mvals[, 5L])
    f <- as.double(Mvals[, 6L])
    u <- as.double(Zvals[, 1L])
    v <- as.double(Zvals[, 2L])
    w <- as.double(Zvals[, 3L])

    det <- a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d)
    good <- is.finite(det) & (abs(det) > .Machine$double.eps)
    if (any(good)) {
      c11 <- d[good] * f[good] - e[good] * e[good]
      c12 <- c[good] * e[good] - b[good] * f[good]
      c13 <- b[good] * e[good] - c[good] * d[good]
      c22 <- a[good] * f[good] - c[good] * c[good]
      c23 <- b[good] * c[good] - a[good] * e[good]
      c33 <- a[good] * d[good] - b[good] * b[good]
      invdet <- 1 / det[good]

      beta1 <- (c11 * u[good] + c12 * v[good] + c13 * w[good]) * invdet
      beta2 <- (c12 * u[good] + c22 * v[good] + c23 * w[good]) * invdet
      beta3 <- (c13 * u[good] + c23 * v[good] + c33 * w[good]) * invdet
      out[good] <- rhs[1L] * beta1 + rhs[2L] * beta2 + rhs[3L] * beta3
    }
    return(out)
  }

  out
}

.np_inid_boot_from_scoef <- function(txdat,
                                     ydat,
                                     tzdat,
                                     exdat,
                                     ezdat,
                                     bws,
                                     B,
                                     counts = NULL,
                                     leave.one.out = FALSE) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  B <- as.integer(B)

  miss.z <- missing(tzdat) || is.null(tzdat)
  if (miss.z) {
    tzdat <- txdat
    ezdat <- exdat
  } else {
    tzdat <- toFrame(tzdat)
    ezdat <- toFrame(ezdat)
  }

  if (nrow(txdat) != nrow(tzdat))
    stop("smooth coefficient inid helper requires aligned txdat/tzdat rows")
  if (nrow(exdat) != nrow(ezdat))
    stop("smooth coefficient inid helper requires aligned exdat/ezdat rows")
  if (ncol(txdat) != ncol(exdat))
    stop("smooth coefficient inid helper requires matching txdat/exdat columns")
  if (nrow(txdat) < 1L || nrow(exdat) < 1L || B < 1L)
    stop("invalid smooth coefficient inid helper dimensions")

  if (length(ydat) != nrow(txdat))
    stop("length of ydat must match training rows in smooth coefficient inid helper")

  txdat <- adjustLevels(txdat, bws$xdati)
  exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
  if (!miss.z) {
    tzdat <- adjustLevels(tzdat, bws$zdati)
    ezdat <- adjustLevels(ezdat, bws$zdati, allowNewCells = TRUE)
  }

  y.num <- .np_inid_scoef_numeric_y(ydat = ydat, bws = bws)
  X.train <- toMatrix(txdat)
  X.eval <- toMatrix(exdat)
  W.train <- as.matrix(data.frame(1, X.train))
  W.eval <- as.matrix(data.frame(1, X.eval))

  npksum.fun <- .npRmpi_bootstrap_estimator("npksum.default")
  kw <- npksum.fun(
    txdat = tzdat,
    exdat = ezdat,
    bws = bws,
    return.kernel.weights = TRUE,
    bandwidth.divide = TRUE,
    leave.one.out = leave.one.out
  )$kw
  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(txdat))

  n <- nrow(W.train)
  neval <- nrow(W.eval)
  if (nrow(kw) != n || ncol(kw) != neval)
    stop("smooth coefficient inid helper kernel-weight matrix shape mismatch")

  p <- ncol(W.train)
  mcols <- p * (p + 1L) / 2L
  ones <- matrix(1.0, nrow = n, ncol = 1L)
  epsilon <- 1.0 / neval

  Mfeat <- vector("list", neval)
  Zfeat <- vector("list", neval)
  t0 <- numeric(neval)

  for (i in seq_len(neval)) {
    k <- as.double(kw[, i])
    WK <- W.train * k
    zf <- WK * y.num

    mf <- matrix(0.0, nrow = n, ncol = mcols)
    idx <- 1L
    for (a in seq_len(p)) {
      for (b in a:p) {
        mf[, idx] <- WK[, a] * W.train[, b]
        idx <- idx + 1L
      }
    }

    Mfeat[[i]] <- mf
    Zfeat[[i]] <- zf

    M0 <- crossprod(ones, mf)
    Z0 <- crossprod(ones, zf)
    t0i <- .np_inid_scoef_predict_chunk(Mvals = M0, Zvals = Z0, rhs = W.eval[i, ])[1L]
    if (!is.finite(t0i)) {
      t0i <- .np_inid_scoef_predict_row(
        mrow = M0[1L, ],
        zrow = Z0[1L, ],
        rhs = W.eval[i, ],
        epsilon = epsilon
      )
    }
    t0[i] <- t0i
  }

  tmat <- matrix(NA_real_, nrow = B, ncol = neval)

  fill_chunk <- function(counts.chunk, start, stopi) {
    bsz <- ncol(counts.chunk)
    for (i in seq_len(neval)) {
      Mvals <- crossprod(counts.chunk, Mfeat[[i]])
      Zvals <- crossprod(counts.chunk, Zfeat[[i]])
      if (bsz == 1L) {
        Mvals <- matrix(Mvals, nrow = 1L)
        Zvals <- matrix(Zvals, nrow = 1L)
      }
      out <- .np_inid_scoef_predict_chunk(Mvals = Mvals, Zvals = Zvals, rhs = W.eval[i, ])
      bad <- which(!is.finite(out))
      if (length(bad)) {
        for (bb in bad) {
          out[bb] <- .np_inid_scoef_predict_row(
            mrow = Mvals[bb, ],
            zrow = Zvals[bb, ],
            rhs = W.eval[i, ],
            epsilon = epsilon
          )
        }
      }
      tmat[start:stopi, i] <<- out
    }
  }

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    fill_chunk(counts.chunk = counts.mat, start = 1L, stopi = B)
  } else {
    chunk.size <- .np_inid_chunk_size(n = n, B = B)
    prob <- rep.int(1 / n, n)
    start <- 1L
    while (start <= B) {
      stopi <- min(B, start + chunk.size - 1L)
      bsz <- stopi - start + 1L
      counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
      fill_chunk(counts.chunk = counts.chunk, start = start, stopi = stopi)
      start <- stopi + 1L
    }
  }

  if (any(!is.finite(t0)) || any(!is.finite(tmat)))
    stop("inid smooth coefficient helper produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_plreg_numeric_x_matrix <- function(txdat, exdat, bws) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)

  p <- ncol(txdat)
  if (p < 1L)
    stop("plreg inid fast path requires at least one linear regressor")
  if (ncol(exdat) != p)
    stop("training/evaluation linear regressor dimensions do not match")

  x.train.num <- matrix(0.0, nrow = nrow(txdat), ncol = p)
  x.eval.num <- matrix(0.0, nrow = nrow(exdat), ncol = p)

  for (j in seq_len(p)) {
    if (is.factor(txdat[[j]])) {
      trj <- adjustLevels(txdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati)
      evj <- adjustLevels(exdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati, allowNewCells = TRUE)
      lev <- bws$bw[[j + 1L]]$ydati$all.dlev[[1L]]
      x.train.num[, j] <- lev[as.integer(trj[, 1L])]
      x.eval.num[, j] <- lev[as.integer(evj[, 1L])]
    } else {
      x.train.num[, j] <- as.double(txdat[[j]])
      x.eval.num[, j] <- as.double(exdat[[j]])
    }
  }

  list(train = x.train.num, eval = x.eval.num)
}

.np_plreg_weighted_coef <- function(X, y, w, ridge = 1.0e-12) {
  X <- as.matrix(X)
  y <- as.double(y)
  w <- as.double(w)
  if (nrow(X) != length(y) || length(w) != length(y))
    stop("weighted plreg solve dimension mismatch")

  w <- pmax(w, 0.0)
  sw <- sqrt(w)
  Xw <- X * sw
  yw <- y * sw

  beta <- tryCatch(
    qr.solve(Xw, yw, tol = .Machine$double.eps),
    error = function(e) NULL
  )
  if (!is.null(beta) && all(is.finite(beta)))
    return(as.double(beta))

  XtWX <- crossprod(X, X * w)
  XtWy <- crossprod(X, y * w)
  ridge <- max(0.0, as.double(ridge))

  for (attempt in 0:8) {
    A <- XtWX
    if (ridge > 0)
      diag(A) <- diag(A) + ridge
    beta <- tryCatch(
      qr.solve(A, XtWy, tol = .Machine$double.eps),
      error = function(e) NULL
    )
    if (!is.null(beta) && all(is.finite(beta)))
      return(as.double(beta))
    ridge <- if (ridge > 0) ridge * 10 else 1.0e-12
  }

  stop("plreg weighted solve failed")
}

.np_inid_boot_from_plreg <- function(txdat,
                                     ydat,
                                     tzdat,
                                     exdat,
                                     ezdat,
                                     bws,
                                     B,
                                     counts = NULL,
                                     ridge = 1.0e-12) {
  txdat <- toFrame(txdat)
  tzdat <- toFrame(tzdat)
  exdat <- toFrame(exdat)
  ezdat <- toFrame(ezdat)
  B <- as.integer(B)

  n <- nrow(txdat)
  neval <- nrow(exdat)
  p <- ncol(txdat)
  if (nrow(tzdat) != n)
    stop("plreg inid fast path requires aligned txdat/tzdat rows")
  if (nrow(ezdat) != neval)
    stop("plreg inid fast path requires aligned exdat/ezdat rows")
  if (n < 1L || neval < 1L || p < 1L || B < 1L)
    stop("invalid plreg inid fast path dimensions")

  y.num <- if (is.factor(ydat)) {
    ty <- adjustLevels(data.frame(ydat), bws$bw$yzbw$ydati)
    bws$bw$yzbw$ydati$all.dlev[[1L]][as.integer(ty[, 1L])]
  } else {
    as.double(ydat)
  }
  if (length(y.num) != n)
    stop("length of ydat must match training rows")

  x.num <- .np_plreg_numeric_x_matrix(txdat = txdat, exdat = exdat, bws = bws)
  x.train.num <- x.num$train
  x.eval.num <- x.num$eval

  counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)

  y.train <- .np_inid_boot_from_regression(
    xdat = tzdat,
    exdat = tzdat,
    bws = bws$bw$yzbw,
    ydat = y.num,
    B = B,
    counts = counts.mat,
    ridge = ridge
  )
  y.eval <- .np_inid_boot_from_regression(
    xdat = tzdat,
    exdat = ezdat,
    bws = bws$bw$yzbw,
    ydat = y.num,
    B = B,
    counts = counts.mat,
    ridge = ridge
  )

  x.train <- vector("list", p)
  x.eval <- vector("list", p)
  for (j in seq_len(p)) {
    x.train[[j]] <- .np_inid_boot_from_regression(
      xdat = tzdat,
      exdat = tzdat,
      bws = bws$bw[[j + 1L]],
      ydat = x.train.num[, j],
      B = B,
      counts = counts.mat,
      ridge = ridge
    )
    x.eval[[j]] <- .np_inid_boot_from_regression(
      xdat = tzdat,
      exdat = ezdat,
      bws = bws$bw[[j + 1L]],
      ydat = x.train.num[, j],
      B = B,
      counts = counts.mat,
      ridge = ridge
    )
  }

  xres.train0 <- matrix(0.0, nrow = n, ncol = p)
  xres.eval0 <- matrix(0.0, nrow = neval, ncol = p)
  for (j in seq_len(p)) {
    xres.train0[, j] <- x.train.num[, j] - as.double(x.train[[j]]$t0)
    xres.eval0[, j] <- x.eval.num[, j] - as.double(x.eval[[j]]$t0)
  }
  yres0 <- y.num - as.double(y.train$t0)
  beta0 <- .np_plreg_weighted_coef(
    X = xres.train0,
    y = yres0,
    w = rep.int(1.0, n),
    ridge = ridge
  )
  t0 <- as.double(y.eval$t0) + as.vector(xres.eval0 %*% beta0)

  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  xres.train.b <- matrix(0.0, nrow = n, ncol = p)
  xres.eval.b <- matrix(0.0, nrow = neval, ncol = p)

  for (b in seq_len(B)) {
    for (j in seq_len(p)) {
      xres.train.b[, j] <- x.train.num[, j] - x.train[[j]]$t[b, ]
      xres.eval.b[, j] <- x.eval.num[, j] - x.eval[[j]]$t[b, ]
    }
    yres.b <- y.num - y.train$t[b, ]
    beta.b <- .np_plreg_weighted_coef(
      X = xres.train.b,
      y = yres.b,
      w = counts.mat[, b],
      ridge = ridge
    )
    tmat[b, ] <- y.eval$t[b, ] + as.vector(xres.eval.b %*% beta.b)
  }

  if (any(!is.finite(t0)) || any(!is.finite(tmat)))
    stop("plreg inid fast path produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_boot_matrix_from_ksum <- function(ksum, B, nout, where = "ksum fast path") {
  if (is.null(dim(ksum))) {
    if (B == 1L && length(ksum) == nout)
      return(matrix(as.double(ksum), nrow = 1L))
    stop(sprintf("%s returned unexpected vector shape", where))
  }

  km <- as.matrix(ksum)
  if (nrow(km) == B && ncol(km) == nout)
    return(km)
  if (nrow(km) == nout && ncol(km) == B)
    return(t(km))
  if (B == 1L && length(km) == nout)
    return(matrix(as.double(km), nrow = 1L))

  stop(sprintf("%s returned unexpected matrix shape", where))
}

.np_inid_boot_from_ksum_unconditional <- function(xdat,
                                                  exdat,
                                                  bws,
                                                  B,
                                                  operator,
                                                  counts = NULL) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)

  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid unconditional inid bootstrap dimensions")

  npksum.fun <- .npRmpi_bootstrap_estimator("npksum")
  ones <- matrix(1.0, nrow = n, ncol = 1L)
  ksum0 <- npksum.fun(
    txdat = xdat,
    tydat = ones,
    exdat = exdat,
    bws = bws,
    weights = ones,
    operator = operator,
    bandwidth.divide = TRUE
  )$ksum
  t0 <- as.numeric(.np_boot_matrix_from_ksum(ksum0, B = 1L, nout = neval,
                                             where = "npksum unconditional baseline")[1L, ]) / n

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    ksum <- npksum.fun(
      txdat = xdat,
      tydat = ones,
      exdat = exdat,
      bws = bws,
      weights = counts.mat,
      operator = operator,
      bandwidth.divide = TRUE
    )$ksum
    tmat <- .np_boot_matrix_from_ksum(ksum, B = B, nout = neval,
                                      where = "npksum unconditional counts")
    return(list(t = tmat / n, t0 = t0))
  }

  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  prob <- rep.int(1 / n, n)
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
    ksum.chunk <- npksum.fun(
      txdat = xdat,
      tydat = ones,
      exdat = exdat,
      bws = bws,
      weights = counts.chunk,
      operator = operator,
      bandwidth.divide = TRUE
    )$ksum
    tmat[start:stopi, ] <- .np_boot_matrix_from_ksum(
      ksum.chunk, B = bsz, nout = neval, where = "npksum unconditional chunk"
    ) / n
    start <- stopi + 1L
  }

  list(t = tmat, t0 = t0)
}

.np_con_inid_ksum_eligible <- function(bws) {
  isTRUE(identical(bws$type, "fixed")) &&
    isTRUE(identical(bws$cxkertype, bws$cykertype)) &&
    isTRUE(identical(bws$cxkerorder, bws$cykerorder)) &&
    isTRUE(identical(bws$cxkerbound, bws$cykerbound)) &&
    isTRUE(identical(bws$uxkertype, bws$uykertype)) &&
    isTRUE(identical(bws$oxkertype, bws$oykertype))
}

.np_con_make_kbandwidth_x <- function(bws, xdat) {
  xdat <- toFrame(xdat)
  kbandwidth.numeric.fun <- .npRmpi_bootstrap_estimator("kbandwidth.numeric")
  untangle.fun <- .npRmpi_bootstrap_estimator("untangle")
  kbandwidth.numeric.fun(
    bw = bws$xbw,
    bwscaling = FALSE,
    bwtype = bws$type,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = if (!is.null(bws$cxkerlb)) bws$cxkerlb else NULL,
    ckerub = if (!is.null(bws$cxkerub)) bws$cxkerub else NULL,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype,
    nobs = nrow(xdat),
    xdati = untangle.fun(xdat),
    xnames = names(xdat)
  )
}

.np_con_make_kbandwidth_xy <- function(bws, xdat, ydat) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  xydat <- data.frame(xdat, ydat)
  kbandwidth.numeric.fun <- .npRmpi_bootstrap_estimator("kbandwidth.numeric")
  untangle.fun <- .npRmpi_bootstrap_estimator("untangle")
  ckerlb <- c(if (is.null(bws$cxkerlb)) numeric(0) else bws$cxkerlb,
              if (is.null(bws$cykerlb)) numeric(0) else bws$cykerlb)
  ckerub <- c(if (is.null(bws$cxkerub)) numeric(0) else bws$cxkerub,
              if (is.null(bws$cykerub)) numeric(0) else bws$cykerub)

  kbandwidth.numeric.fun(
    bw = c(bws$xbw, bws$ybw),
    bwscaling = FALSE,
    bwtype = bws$type,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = if (length(ckerlb)) ckerlb else NULL,
    ckerub = if (length(ckerub)) ckerub else NULL,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype,
    nobs = nrow(xydat),
    xdati = untangle.fun(xydat),
    xnames = names(xydat)
  )
}

.np_inid_boot_from_ksum_conditional <- function(xdat,
                                                ydat,
                                                exdat,
                                                eydat,
                                                bws,
                                                B,
                                                cdf,
                                                counts = NULL) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)

  if (nrow(ydat) != n || nrow(eydat) != neval)
    stop("conditional inid fast path requires aligned x/y training and evaluation rows")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid conditional inid bootstrap dimensions")
  if (!.np_con_inid_ksum_eligible(bws))
    return(NULL)
  npksum.fun <- .npRmpi_bootstrap_estimator("npksum")

  kbx <- tryCatch(.np_con_make_kbandwidth_x(bws = bws, xdat = xdat),
                  error = function(e) NULL)
  kbxy <- tryCatch(.np_con_make_kbandwidth_xy(bws = bws, xdat = xdat, ydat = ydat),
                   error = function(e) NULL)
  if (is.null(kbx) || is.null(kbxy))
    return(NULL)

  xop <- rep.int("normal", ncol(xdat))
  yop <- rep.int(if (cdf) "integral" else "normal", ncol(ydat))
  xyop <- c(xop, yop)

  xydat <- data.frame(xdat, ydat)
  exydat <- data.frame(exdat, eydat)
  ones <- matrix(1.0, nrow = n, ncol = 1L)

  den0 <- npksum.fun(
    txdat = xdat,
    tydat = ones,
    exdat = exdat,
    bws = kbx,
    weights = ones,
    operator = xop,
    bandwidth.divide = TRUE
  )$ksum
  num0 <- npksum.fun(
    txdat = xydat,
    tydat = ones,
    exdat = exydat,
    bws = kbxy,
    weights = ones,
    operator = xyop,
    bandwidth.divide = TRUE
  )$ksum
  den0 <- as.numeric(.np_boot_matrix_from_ksum(den0, B = 1L, nout = neval,
                                               where = "npksum conditional denominator baseline")[1L, ]) / n
  num0 <- as.numeric(.np_boot_matrix_from_ksum(num0, B = 1L, nout = neval,
                                               where = "npksum conditional numerator baseline")[1L, ]) / n
  t0 <- num0 / pmax(den0, .Machine$double.eps)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    den <- npksum.fun(
      txdat = xdat,
      tydat = ones,
      exdat = exdat,
      bws = kbx,
      weights = counts.mat,
      operator = xop,
      bandwidth.divide = TRUE
    )$ksum
    num <- npksum.fun(
      txdat = xydat,
      tydat = ones,
      exdat = exydat,
      bws = kbxy,
      weights = counts.mat,
      operator = xyop,
      bandwidth.divide = TRUE
    )$ksum
    den <- .np_boot_matrix_from_ksum(den, B = B, nout = neval,
                                     where = "npksum conditional denominator counts") / n
    num <- .np_boot_matrix_from_ksum(num, B = B, nout = neval,
                                     where = "npksum conditional numerator counts") / n
    return(list(t = num / pmax(den, .Machine$double.eps), t0 = t0))
  }

  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  prob <- rep.int(1 / n, n)
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)

    den <- npksum.fun(
      txdat = xdat,
      tydat = ones,
      exdat = exdat,
      bws = kbx,
      weights = counts.chunk,
      operator = xop,
      bandwidth.divide = TRUE
    )$ksum
    num <- npksum.fun(
      txdat = xydat,
      tydat = ones,
      exdat = exydat,
      bws = kbxy,
      weights = counts.chunk,
      operator = xyop,
      bandwidth.divide = TRUE
    )$ksum

    den <- .np_boot_matrix_from_ksum(den, B = bsz, nout = neval,
                                     where = "npksum conditional denominator chunk") / n
    num <- .np_boot_matrix_from_ksum(num, B = bsz, nout = neval,
                                     where = "npksum conditional numerator chunk") / n
    tmat[start:stopi, ] <- num / pmax(den, .Machine$double.eps)
    start <- stopi + 1L
  }

  list(t = tmat, t0 = t0)
}

gen.label = function(label, altlabel){
  paste(if (is.null(label)) altlabel else label)
}

gen.tflabel = function(condition, tlabel, flabel){
  paste(if (isTRUE(condition)) tlabel else flabel)
}

draw.error.bands = function(ex, ely, ehy, lty = 2, col = par("col")){
  lines(ex,ely,lty=lty,col=col)
  lines(ex,ehy,lty=lty,col=col)
}

draw.error.bars = function(ex, ely, ehy, hbar = TRUE, hbarscale = 0.3, lty = 2, col = par("col")){
  yy = double(3*length(ex))
  jj = seq_along(ex)*3

  yy[jj-2] = ely
  yy[jj-1] = ehy
  yy[jj] = NA
  
  xx = double(3*length(ex))
  xx[jj-2] = ex
  xx[jj-1] = ex
  xx[jj] = NA

  lines(xx,yy,lty=lty,col=col)

  if (hbar){
    ## hbars look silly if they are too wide in relation to their height
    ## this only matters in the limit of few points, since that is when
    ## hbardist may get relatively large

    golden = (1+sqrt(5))/2
    hbardist = abs(max(ex) - min(ex))/length(ex)*hbarscale

    yg = abs(yy[jj-2]-yy[jj-1])/golden
    htest = (hbardist >= yg)
    
    hdelta = pmin(yg, hbardist)/2
    xx[jj-2] = ex - hdelta
    xx[jj-1] = ex + hdelta
    
    ty = yy[jj-1]
    yy[jj-1] = yy[jj-2]

    lines(xx,yy,col=col)

    yy[jj-2] = ty
    yy[jj-1] = ty

    lines(xx,yy,col=col)
  }
}

draw.errors =
  function(ex, ely, ehy,
           plot.errors.style,
           plot.errors.bar,
           plot.errors.bar.num,
           lty,
           col = par("col")){
    if (plot.errors.style == "bar"){
      ei = seq(1,length(ex),length.out = min(length(ex),plot.errors.bar.num))
      draw.error.bars(ex = ex[ei],
                      ely = ely[ei],
                      ehy = ehy[ei],
                      hbar = (plot.errors.bar == "I"),
                      lty = lty,
                      col = col)
    } else if (plot.errors.style == "band") {
      draw.error.bands(ex = ex,
                       ely = ely,
                       ehy = ehy,
                       lty = lty,
                       col = col)
    }
  }

draw.all.error.types <- function(ex, center, all.err,
                                 plot.errors.style = "band",
                                 plot.errors.bar = "|",
                                 plot.errors.bar.num = min(length(ex), 25),
                                 lty = 2, add.legend = TRUE, legend.loc = "topleft",
                                 xi.factor = FALSE){
  if (is.null(all.err)) return(invisible(NULL))

  if (xi.factor) {
    plot.errors.style <- "bar"
    plot.errors.bar <- "I"
  }

  draw_one <- function(err, col) {
    if (is.null(err)) return(invisible(NULL))
    lower <- center - err[,1]
    upper <- center + err[,2]
    good <- complete.cases(ex, lower, upper)
    if (!any(good)) return(invisible(NULL))
    draw.errors(ex = ex[good], ely = lower[good], ehy = upper[good],
                plot.errors.style = plot.errors.style,
                plot.errors.bar = plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = lty, col = col)
  }

  draw_one(all.err$pointwise, "red")
  draw_one(all.err$simultaneous, "green3")
  draw_one(all.err$bonferroni, "blue")

  if (add.legend) {
    legend(legend.loc,
           legend = c("Pointwise","Simultaneous","Bonferroni"),
           lty = 2, col = c("red","green3","blue"), lwd = 2, bty = "n")
  }
}

plotFactor <- function(f, y, ...){
  dot.args <- list(...)
  dot.names <- names(dot.args)
  has.user.lty <- !is.null(dot.names) && any(dot.names == "lty")

  if (has.user.lty) {
    do.call(plot, c(list(x = f, y = y), dot.args))
  } else {
    plot(x = f, y = y, lty = "blank", ...)
  }

  l.f = rep(f, each=3)
  l.f[3*seq_along(f)] = NA

  l.y = unlist(lapply(y, function (p) { c(0,p,NA) }))

  lines(x = l.f, y = l.y, lty = 2)
  points(x = f, y = y)
}

.np_plot_panel_fun <- function(plot.bootstrap, plot.bxp) {
  if (plot.bootstrap && plot.bxp) bxp else plotFactor
}

.np_plot_resolve_xydat <- function(bws, xdat, ydat, miss.xy) {
  if (any(miss.xy) && !all(miss.xy))
    stop("one of, but not both, xdat and ydat was specified")

  if (all(miss.xy) && !is.null(bws$formula)) {
    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1, m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"), "term.labels"), drop = FALSE]
    return(list(xdat = xdat, ydat = ydat))
  }

  if (all(miss.xy) && !is.null(bws$call)) {
    xdat <- data.frame(.np_eval_bws_call_arg(bws, "xdat"))
    ydat <- .np_eval_bws_call_arg(bws, "ydat")
  }

  xdat <- toFrame(xdat)
  goodrows <- seq_len(nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (!is.null(rows.omit))
    goodrows[rows.omit] <- 0

  if (all(goodrows == 0))
    stop("Data has no rows without NAs")

  goodrows <- goodrows[goodrows != 0]
  list(xdat = xdat[goodrows, , drop = FALSE],
       ydat = ydat[goodrows])
}

.npRmpi_guard_bootstrap_plot_autodispatch <- function(plot.errors.method,
                                                      where = "plot()",
                                                      allow.direct.bootstrap = FALSE) {
  invisible(TRUE)
}

.npRmpi_with_local_bootstrap <- function(expr) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  force(expr)
}

.npRmpi_bootstrap_estimator <- function(name) {
  nprmpi.ns <- asNamespace("npRmpi")
  if (exists(name, envir = nprmpi.ns, mode = "function", inherits = FALSE))
    return(get(name, envir = nprmpi.ns, inherits = FALSE))

  get(name, mode = "function", envir = parent.frame(), inherits = TRUE)
}

.npRmpi_bootstrap_uses_local_namespace <- function(fun) {
  if (!is.function(fun))
    return(FALSE)
  env <- environment(fun)
  if (is.null(env))
    return(FALSE)
  environmentName(env) %in% c("np", "namespace:np", "npRmpi", "namespace:npRmpi")
}

.npRmpi_bootstrap_maybe_local <- function(expr, estimators = list()) {
  use.local <- length(estimators) > 0L &&
    all(vapply(estimators, .npRmpi_bootstrap_uses_local_namespace, logical(1)))
  if (use.local)
    return(.npRmpi_with_local_bootstrap(expr))
  force(expr)
}

.npRmpi_plot_behavior_for_rank <- function(plot.behavior) {
  if (!.npRmpi_autodispatch_called_from_bcast())
    return(plot.behavior)

  rank <- tryCatch(mpi.comm.rank(), error = function(e) NA_integer_)
  if (is.na(rank) || rank == 0L)
    return(plot.behavior)

  "data"
}

.npRmpi_profile_env <- local({
  env <- new.env(parent = emptyenv())
  env$last <- NULL
  env$history <- list()
  env$active_id <- NULL
  env$active <- list()
  env
})

.npRmpi_profile_enabled <- function() {
  !isFALSE(getOption("npRmpi.profile.enable", TRUE))
}

.npRmpi_profile_level <- function() {
  lvl <- as.character(getOption("npRmpi.profile.level", "basic"))[1L]
  if (is.na(lvl) || !(lvl %in% c("basic", "detailed")))
    lvl <- "basic"
  lvl
}

.npRmpi_profile_history_limit <- function() {
  lim <- suppressWarnings(as.integer(getOption("npRmpi.profile.history.max", 200L))[1L])
  if (is.na(lim) || lim < 1L)
    lim <- 200L
  lim
}

.np_nrows_safe <- function(x) {
  if (is.null(x))
    return(NA_integer_)
  if (is.atomic(x) && is.null(dim(x)))
    return(as.integer(length(x)))
  xf <- tryCatch(toFrame(x), error = function(e) NULL)
  if (is.null(xf))
    return(NA_integer_)
  as.integer(nrow(xf))
}

.npRmpi_profile_bootstrap_begin <- function(where, method = NA_character_,
                                            B = NA_integer_, ntrain = NA_integer_,
                                            neval = NA_integer_) {
  if (!isTRUE(.npRmpi_profile_enabled()))
    return(NULL)

  rank <- tryCatch(as.integer(mpi.comm.rank()), error = function(e) NA_integer_)
  size <- tryCatch(as.integer(mpi.comm.size()), error = function(e) NA_integer_)
  active <- list(
    where = where,
    level = .npRmpi_profile_level(),
    method = if (length(method)) as.character(method)[1L] else NA_character_,
    B = suppressWarnings(as.integer(B)[1L]),
    ntrain = suppressWarnings(as.integer(ntrain)[1L]),
    neval = suppressWarnings(as.integer(neval)[1L]),
    rank = rank,
    size = size,
    via_bcast = isTRUE(.npRmpi_autodispatch_called_from_bcast()),
    comm_elapsed_sec = 0.0,
    comm_calls = 0L,
    comm_notes = character(0),
    start_proc = proc.time(),
    start_wall = Sys.time()
  )
  active$id <- paste0(format(active$start_wall, "%Y%m%d%H%M%OS6"), "-", sample.int(1e8, 1L))
  .npRmpi_profile_env$active_id <- active$id
  .npRmpi_profile_env$active[[active$id]] <- active
  active
}

.npRmpi_profile_add_comm_elapsed <- function(elapsed_sec, where = NA_character_) {
  id <- .npRmpi_profile_env$active_id
  if (is.null(id))
    return(invisible(FALSE))
  cur <- .npRmpi_profile_env$active[[id]]
  if (is.null(cur))
    return(invisible(FALSE))

  elapsed <- suppressWarnings(as.double(elapsed_sec)[1L])
  if (is.na(elapsed) || elapsed < 0)
    elapsed <- 0.0

  cur$comm_elapsed_sec <- as.double(cur$comm_elapsed_sec) + elapsed
  cur$comm_calls <- as.integer(cur$comm_calls) + 1L
  if (!is.na(where) && nzchar(as.character(where)[1L]) &&
      identical(cur$level, "detailed")) {
    cur$comm_notes <- c(cur$comm_notes, as.character(where)[1L])
  }

  .npRmpi_profile_env$active[[id]] <- cur
  invisible(TRUE)
}

.npRmpi_profile_bootstrap_end <- function(value, ctx) {
  if (is.null(ctx))
    return(value)

  id <- ctx$id
  active <- if (!is.null(id)) .npRmpi_profile_env$active[[id]] else NULL
  if (is.list(active)) {
    ctx$comm_elapsed_sec <- as.double(active$comm_elapsed_sec)
    ctx$comm_calls <- as.integer(active$comm_calls)
    if (identical(active$level, "detailed"))
      ctx$comm_notes <- active$comm_notes
  }

  dt <- proc.time() - ctx$start_proc
  wall <- unname(as.double(dt[["elapsed"]]))
  comm <- suppressWarnings(as.double(ctx$comm_elapsed_sec)[1L])
  if (is.na(comm) || comm < 0)
    comm <- 0.0
  if (!is.finite(wall) || wall <= 0)
    wall <- 0.0
  compute <- max(0.0, wall - comm)
  denom <- comm + compute
  ratio <- if (denom > 0) min(1.0, max(0.0, comm / denom)) else NA_real_

  record <- list(
    where = ctx$where,
    level = ctx$level,
    method = ctx$method,
    B = ctx$B,
    ntrain = ctx$ntrain,
    neval = ctx$neval,
    rank = ctx$rank,
    size = ctx$size,
    via_bcast = ctx$via_bcast,
    wall_elapsed_sec = wall,
    comm_elapsed_sec = comm,
    compute_elapsed_sec = compute,
    comm_ratio = ratio,
    comm_calls = suppressWarnings(as.integer(ctx$comm_calls)[1L]),
    user_sec = unname(as.double(dt[["user.self"]])),
    system_sec = unname(as.double(dt[["sys.self"]])),
    timestamp_start = ctx$start_wall,
    timestamp_end = Sys.time()
  )
  if (identical(ctx$level, "detailed")) {
    record$comm_notes <- if (length(ctx$comm_notes)) ctx$comm_notes else character(0)
  }

  .npRmpi_profile_env$last <- record
  .npRmpi_profile_env$history <- c(.npRmpi_profile_env$history, list(record))
  if (!is.null(id)) {
    .npRmpi_profile_env$active[[id]] <- NULL
    if (identical(.npRmpi_profile_env$active_id, id))
      .npRmpi_profile_env$active_id <- NULL
  }
  keep <- .npRmpi_profile_history_limit()
  if (length(.npRmpi_profile_env$history) > keep) {
    .npRmpi_profile_env$history <- tail(.npRmpi_profile_env$history, keep)
  }

  if (is.list(value))
    value$timing.profile <- record
  value
}

.npRmpi_profile_finalize_bootstrap <- function(boot.err, bxp, boot.all.err, ctx) {
  out <- list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
  .npRmpi_profile_bootstrap_end(out, ctx)
}

.npRmpi_profile_last <- function() {
  .npRmpi_profile_env$last
}

.npRmpi_profile_history <- function(n = NULL) {
  h <- .npRmpi_profile_env$history
  if (is.null(n))
    return(h)
  n <- suppressWarnings(as.integer(n)[1L])
  if (is.na(n) || n < 1L)
    return(list())
  tail(h, n)
}

.npRmpi_profile_clear <- function() {
  .npRmpi_profile_env$last <- NULL
  .npRmpi_profile_env$history <- list()
  .npRmpi_profile_env$active_id <- NULL
  .npRmpi_profile_env$active <- list()
  invisible(NULL)
}

## Rank-based simultaneous confidence set helper, vendored from
## MCPAN::SCSrank (MCPAN 1.1-21, GPL-2; Schaarschmidt, Gerhard, Sill).
np.plot.SCSrank <- function(x, conf.level = 0.95, alternative = "two.sided", ...) {
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))

  DataMatrix <- x
  N <- nrow(DataMatrix)
  k <- round(conf.level * N, 0)
  RankDat <- apply(DataMatrix, 2, rank)

  switch(alternative,
    "two.sided" = {
      W1 <- apply(RankDat, 1, max)
      W2 <- N + 1 - apply(RankDat, 1, min)

      Wmat <- cbind(W1, W2)
      w <- apply(Wmat, 1, max)
      tstar <- round(sort(w)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(sortx[N + 1 - tstar], sortx[tstar])
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    },
    "less" = {
      W1 <- apply(RankDat, 1, max)
      tstar <- round(sort(W1)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(-Inf, sortx[tstar])
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    },
    "greater" = {
      W2 <- N + 1 - apply(RankDat, 1, min)
      tstar <- round(sort(W2)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(sortx[N + 1 - tstar], Inf)
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    }
  )

  colnames(SCS) <- c("lower", "upper")

  attr(SCS, which = "k") <- k
  attr(SCS, which = "N") <- N
  OUT <- list(conf.int = SCS, conf.level = conf.level, alternative = alternative)
  return(OUT)
}

compute.bootstrap.quantile.bounds <- function(boot.t, alpha, band.type) {
  neval <- ncol(boot.t)

  if (band.type == "pointwise") {
    probs <- c(alpha / 2.0, 1.0 - alpha / 2.0)
    return(t(apply(boot.t, 2, quantile, probs = probs)))
  }

  if (band.type == "bonferroni") {
    probs <- c(alpha / (2.0 * neval), 1.0 - alpha / (2.0 * neval))
    return(t(apply(boot.t, 2, quantile, probs = probs)))
  }

  if (band.type == "simultaneous") {
    return(np.plot.SCSrank(boot.t, conf.level = 1.0 - alpha)$conf.int)
  }

  if (band.type == "all") {
    return(list(
      pointwise = compute.bootstrap.quantile.bounds(boot.t, alpha, "pointwise"),
      bonferroni = compute.bootstrap.quantile.bounds(boot.t, alpha, "bonferroni"),
      simultaneous = compute.bootstrap.quantile.bounds(boot.t, alpha, "simultaneous")
    ))
  }

  stop("'band.type' must be one of pointwise, bonferroni, simultaneous, all")
}

compute.all.error.range <- function(center, all.err) {
  if (is.null(all.err)) {
    return(c(NA_real_, NA_real_))
  }
  lower <- c(center - all.err$pointwise[,1],
             center - all.err$simultaneous[,1],
             center - all.err$bonferroni[,1])
  upper <- c(center + all.err$pointwise[,2],
             center + all.err$simultaneous[,2],
             center + all.err$bonferroni[,2])
  c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
}

compute.default.error.range <- function(center, err) {
  lower <- c(center - err[,1], err[,3] - err[,1])
  upper <- c(center + err[,2], err[,3] + err[,2])
  c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
}

.np_plot_normalize_common_options <- function(plot.behavior,
                                             plot.errors.method,
                                             plot.errors.boot.method,
                                             plot.errors.boot.wild = c("rademacher", "mammen"),
                                             plot.errors.boot.blocklen,
                                             plot.errors.center,
                                             plot.errors.type,
                                             plot.errors.alpha,
                                             plot.errors.style,
                                             plot.errors.bar,
                                             xdat,
                                             common.scale,
                                             ylim,
                                             allow_asymptotic_quantile = TRUE) {
  scalar_choice <- function(value, default) {
    if (is.null(value) || length(value) < 1L || is.na(value[1L])) default else value[1L]
  }

  plot.behavior <- match.arg(
    scalar_choice(plot.behavior, "plot"),
    c("plot", "plot-data", "data")
  )
  plot.errors.method <- match.arg(
    scalar_choice(plot.errors.method, "none"),
    c("none", "bootstrap", "asymptotic")
  )
  plot.errors.boot.method <- match.arg(
    scalar_choice(plot.errors.boot.method, "wild"),
    c("wild", "inid", "fixed", "geom")
  )
  plot.errors.boot.wild <- match.arg(
    scalar_choice(plot.errors.boot.wild, "rademacher"),
    c("rademacher", "mammen")
  )
  plot.errors.center <- match.arg(
    scalar_choice(plot.errors.center, "estimate"),
    c("estimate", "bias-corrected")
  )
  plot.errors.type <- match.arg(
    scalar_choice(plot.errors.type, "simultaneous"),
    c("simultaneous", "pointwise", "bonferroni", "pmzsd", "all")
  )

  if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
      is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
    stop("the tail probability plot.errors.alpha must lie in (0,0.5)")

  plot.errors.style <- match.arg(
    scalar_choice(plot.errors.style, "band"),
    c("band", "bar")
  )
  plot.errors.bar <- match.arg(
    scalar_choice(plot.errors.bar, "|"),
    c("|", "I")
  )

  common.scale <- common.scale | (!is.null(ylim))

  if (plot.errors.method == "none" && plot.errors.type == "all") {
    warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
    plot.errors.method <- "bootstrap"
  }

  if (allow_asymptotic_quantile && plot.errors.method == "asymptotic") {
    if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      warning("bootstrap quantile bands cannot be calculated with asymptotics, calculating pmzsd errors")
      plot.errors.type <- "pmzsd"
    }

    if (plot.errors.center == "bias-corrected") {
      warning("no bias corrections can be calculated with asymptotics, centering on estimate")
      plot.errors.center <- "estimate"
    }
  }

  if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
      is.null(plot.errors.boot.blocklen))
    plot.errors.boot.blocklen <- b.star(xdat, round = TRUE)[1,1]

  list(
    plot.behavior = plot.behavior,
    plot.errors.method = plot.errors.method,
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.wild = plot.errors.boot.wild,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    plot.errors.style = plot.errors.style,
    plot.errors.bar = plot.errors.bar,
    common.scale = common.scale,
    plot.errors = (plot.errors.method != "none")
  )
}


compute.bootstrap.errors = function(...,bws){
  UseMethod("compute.bootstrap.errors",bws)
}

compute.bootstrap.errors.rbandwidth =
  function(xdat, ydat,
           exdat,
           gradients,
           gradient.order,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.rbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.rbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    npreg_fit <- .npRmpi_bootstrap_estimator("npreg.rbandwidth")
    npreghat_fit <- .npRmpi_bootstrap_estimator("npreghat.rbandwidth")
    boot.out <- NULL

    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL
    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"

    if (is.wild.hat && gradients) {
      cont.idx <- which(bws$xdati$icon)
      if (is.na(match(slice.index, cont.idx))) {
        warning("plot.errors.boot.method='wild' supports gradients only for continuous slices; using requested bootstrap method fallback")
        is.wild.hat <- FALSE
      }
    }

    inid.helper.ok <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      !isTRUE(gradients) &&
      identical(bws$type, "fixed")

    if (is.inid && !isTRUE(inid.helper.ok)) {
      warning("inid regression helper unavailable for this configuration; using explicit bootstrap fallback")
    }

    if (is.inid && isTRUE(inid.helper.ok)) {
      boot.out <- tryCatch(
        .np_inid_boot_from_regression(
          xdat = xdat,
          exdat = exdat,
          bws = bws,
          ydat = ydat,
          B = plot.errors.boot.num
        ),
        error = function(e) {
          stop(sprintf("inid regression helper failed in compute.bootstrap.errors.rbandwidth (%s)",
                       conditionMessage(e)),
               call. = FALSE)
        }
      )
    }

    if (is.null(boot.out) && is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.train <- suppressWarnings(npreg_fit(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        gradients = FALSE,
        warn.glp.gradient = FALSE
      ))

      s.vec <- NULL
      if (gradients) {
        cont.idx <- which(bws$xdati$icon)
        cpos <- match(slice.index, cont.idx)
        gorder <- if (length(gradient.order) == 1L) {
          rep.int(as.integer(gradient.order), length(cont.idx))
        } else {
          as.integer(gradient.order)
        }
        if (length(gorder) != length(cont.idx))
          gorder <- rep.int(1L, length(cont.idx))
        s.vec <- integer(length(cont.idx))
        s.vec[cpos] <- gorder[cpos]
      }

      H <- suppressWarnings(npreghat_fit(
        bws = bws,
        txdat = xdat,
        exdat = exdat,
        s = s.vec,
        output = "matrix"
      ))

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - fit.train$mean)
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- list(
        t = .np_wild_boot_t(
          H = H,
          fit.mean = fit.train$mean,
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
    }

    if (is.null(boot.out)) {
      if (is.inid) {
        boofun <- function(data, indices) {
          fit <- suppressWarnings(npreg_fit(
            txdat = data[indices, seq_len(ncol(data) - 1L), drop = FALSE],
            tydat = data[indices, ncol(data), drop = TRUE],
            exdat = exdat, bws = bws,
            gradients = gradients,
            gradient.order = gradient.order,
            warn.glp.gradient = FALSE
          ))
          if (gradients) fit$grad[, slice.index] else fit$mean
        }

        boot.out <- .npRmpi_bootstrap_maybe_local({
          boot(
            data = data.frame(xdat, ydat),
            statistic = boofun,
            R = plot.errors.boot.num
          )
        }, estimators = list(npreg_fit, npreghat_fit))
      } else {
        boofun <- function(tsb) {
          fit <- suppressWarnings(npreg_fit(
            txdat = tsb[, seq_len(ncol(tsb) - 1L), drop = FALSE],
            tydat = tsb[, ncol(tsb)],
            exdat = exdat, bws = bws,
            gradients = gradients,
            gradient.order = gradient.order,
            warn.glp.gradient = FALSE
          ))
          if (gradients) fit$grad[, slice.index] else fit$mean
        }

        boot.out <- .npRmpi_bootstrap_maybe_local({
          tsboot(
            tseries = data.frame(xdat, ydat),
            statistic = boofun,
            R = plot.errors.boot.num,
            l = plot.errors.boot.blocklen,
            sim = plot.errors.boot.method
          )
        }, estimators = list(npreg_fit, npreghat_fit))
      }
    }

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.scbandwidth =
  function(xdat, ydat, zdat,
           exdat, ezdat,
           gradients,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.scbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.scbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    miss.z <- missing(zdat)
    npscoef_fit <- .npRmpi_bootstrap_estimator("npscoef.scbandwidth")
    npscoefhat_fit <- .npRmpi_bootstrap_estimator("npscoefhat")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"
    boot.out <- NULL

    if (is.inid) {
      if (!isTRUE(.np_plot_inid_fastpath_enabled()))
        stop("inid bootstrap requires fastpath-enabled helper for smooth coefficient plots", call. = FALSE)
      if (isTRUE(gradients))
        stop("inid bootstrap for smooth coefficient gradients is not supported in helper mode", call. = FALSE)
      boot.out <- tryCatch(
        .np_inid_boot_from_scoef(
          txdat = xdat,
          ydat = ydat,
          tzdat = if (miss.z) NULL else zdat,
          exdat = exdat,
          ezdat = if (miss.z) NULL else ezdat,
          bws = bws,
          B = plot.errors.boot.num,
          leave.one.out = FALSE
        ),
        error = function(e) {
          stop(sprintf("inid smooth coefficient helper failed in compute.bootstrap.errors.scbandwidth (%s)",
                       conditionMessage(e)),
               call. = FALSE)
        }
      )
    }

    if (is.null(boot.out) && is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.args <- list(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        iterate = FALSE
      )
      hat.args <- list(
        bws = bws,
        txdat = xdat,
        exdat = exdat,
        output = "matrix",
        iterate = FALSE
      )
      if (!miss.z) {
        fit.args$tzdat <- zdat
        hat.args$tzdat <- zdat
        hat.args$ezdat <- ezdat
      }

      fit.train <- do.call(npscoef_fit, fit.args)
      H <- do.call(npscoefhat_fit, hat.args)

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - as.vector(fit.train$mean))
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- list(
        t = .np_wild_boot_t(
          H = H,
          fit.mean = as.vector(fit.train$mean),
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
    }

    if (is.null(boot.out)) {
      xcols <- seq_len(ncol(xdat))
      ycol <- ncol(xdat) + 1L
      zcols <- if (miss.z) integer(0) else (ycol + 1L):(ycol + ncol(zdat))

      boofun <- function(tsb) {
        npscoef_fit(
          txdat = tsb[, xcols, drop = FALSE],
          tydat = tsb[, ycol],
          tzdat = if (miss.z) NULL else tsb[, zcols, drop = FALSE],
          exdat = exdat,
          ezdat = if (miss.z) NULL else ezdat,
          bws = bws,
          iterate = FALSE
        )$mean
      }

      boot.data <- if (miss.z) data.frame(xdat, ydat) else data.frame(xdat, ydat, zdat)
      boot.out <- .npRmpi_bootstrap_maybe_local({
        tsboot(
          tseries = boot.data, statistic = boofun, R = plot.errors.boot.num,
          l = plot.errors.boot.blocklen, sim = plot.errors.boot.method
        )
      }, estimators = list(npscoef_fit, npscoefhat_fit))
    }

    all.bp <- list()

    if ((slice.index > 0) && (((slice.index <= ncol(xdat)) && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])) ||
                              ((slice.index > ncol(xdat)) && (bws$zdati$iord[slice.index-ncol(xdat)] || bws$zdati$iuno[slice.index-ncol(xdat)])))) {
      boot.frame <- as.data.frame(boot.out$t)

      if(slice.index <= ncol(xdat))
          u.lev <- bws$xdati$all.ulev[[slice.index]]
      else
          u.lev <- bws$zdati$all.ulev[[slice.index-ncol(xdat)]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))

      if(slice.index <= ncol(xdat))
          all.bp$names <- bws$xdati$all.lev[[slice.index]]
      else
          all.bp$names <- bws$zdati$all.lev[[slice.index-ncol(xdat)]]
      rm(boot.frame)
    }
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.plbandwidth =
  function(xdat, ydat, zdat,
           exdat, ezdat,
           gradients,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.plbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.plbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    npplreg_fit <- .npRmpi_bootstrap_estimator("npplreg.plbandwidth")
    npplreghat_fit <- .npRmpi_bootstrap_estimator("npplreghat")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"

    if (is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.train <- npplreg_fit(
        txdat = xdat,
        tydat = ydat,
        tzdat = zdat,
        bws = bws
      )
      H <- npplreghat_fit(
        bws = bws,
        txdat = xdat,
        tzdat = zdat,
        exdat = exdat,
        ezdat = ezdat,
        output = "matrix"
      )

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - as.vector(fit.train$mean))
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- list(
        t = .np_wild_boot_t(
          H = H,
          fit.mean = as.vector(fit.train$mean),
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
    } else {
      boot.out <- NULL
      if (is.inid) {
        if (!isTRUE(.np_plot_inid_fastpath_enabled()))
          stop("inid bootstrap requires fastpath-enabled helper for partially linear plots", call. = FALSE)
        boot.out <- tryCatch(
          .np_inid_boot_from_plreg(
            txdat = xdat,
            ydat = ydat,
            tzdat = zdat,
            exdat = exdat,
            ezdat = ezdat,
            bws = bws,
            B = plot.errors.boot.num
          ),
          error = function(e) {
            stop(sprintf("inid plreg helper failed in compute.bootstrap.errors.plbandwidth (%s)",
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      }

      if (is.null(boot.out)) {
      boofun <- function(tsb) {
        npplreg_fit(
          txdat = tsb[, seq_len(ncol(xdat)), drop = FALSE],
          tydat = tsb[, ncol(xdat) + 1L],
          tzdat = tsb[, (ncol(xdat) + 2L):ncol(tsb), drop = FALSE],
          exdat = exdat, ezdat = ezdat, bws = bws
        )$mean
      }

      boot.out <- .npRmpi_bootstrap_maybe_local({
        tsboot(tseries = data.frame(xdat, ydat, zdat), statistic = boofun,
               R = plot.errors.boot.num,
               l = plot.errors.boot.blocklen,
               sim = plot.errors.boot.method)
      }, estimators = list(npplreg_fit, npplreghat_fit))
      }
    }

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$zdati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.bandwidth =
  function(xdat, 
           exdat,
           cdf,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.bandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.bandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    npudens_fit <- .npRmpi_bootstrap_estimator("npudens.bandwidth")
    npudist_fit <- .npRmpi_bootstrap_estimator("npudist.dbandwidth")
    .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "unconditional density/distribution estimators")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"
    fast.inid <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.inid) &&
      isTRUE(identical(bws$type, "fixed"))

    boot.out <- NULL
    if (fast.inid) {
      op <- if (cdf) "integral" else "normal"
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_ksum_unconditional(
            xdat = xdat,
            exdat = exdat,
            bws = bws,
            B = plot.errors.boot.num,
            operator = op
          ),
          error = function(e) {
            warning(sprintf("inid ksum fast path failed in compute.bootstrap.errors.bandwidth (%s); using bootstrap fallback",
                            conditionMessage(e)))
            NULL
          }
        )
      })
    }

    if (is.null(boot.out)) {
      boofun <- if (is.inid) {
        function(data, indices) {
          fit <- if (cdf) {
            npudist_fit(tdat = xdat[indices, , drop = FALSE], edat = exdat, bws = bws)
          } else {
            npudens_fit(tdat = xdat[indices, , drop = FALSE], edat = exdat, bws = bws)
          }
          if (cdf) fit$dist else fit$dens
        }
      } else {
        function(tsb) {
          fit <- if (cdf) {
            npudist_fit(tdat = tsb, edat = exdat, bws = bws)
          } else {
            npudens_fit(tdat = tsb, edat = exdat, bws = bws)
          }
          if (cdf) fit$dist else fit$dens
        }
      }

      boot.out <- .npRmpi_with_local_bootstrap({
        if (is.inid) {
          boot(data = data.frame(xdat), statistic = boofun,
               R = plot.errors.boot.num)
        } else {
          tsboot(tseries = data.frame(xdat), statistic = boofun,
                 R = plot.errors.boot.num,
                 l = plot.errors.boot.blocklen,
                 sim = plot.errors.boot.method)
        }
      })
    }

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.dbandwidth =
  function(xdat, 
           exdat,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.dbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.dbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    npudist_fit <- .npRmpi_bootstrap_estimator("npudist.dbandwidth")
    .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "unconditional density/distribution estimators")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"
    fast.inid <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.inid) &&
      isTRUE(identical(bws$type, "fixed"))

    boot.out <- NULL
    if (fast.inid) {
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_ksum_unconditional(
            xdat = xdat,
            exdat = exdat,
            bws = bws,
            B = plot.errors.boot.num,
            operator = "integral"
          ),
          error = function(e) {
            warning(sprintf("inid ksum fast path failed in compute.bootstrap.errors.dbandwidth (%s); using bootstrap fallback",
                            conditionMessage(e)))
            NULL
          }
        )
      })
    }

    if (is.null(boot.out)) {
      boofun <- if (is.inid) {
        function(data, indices) {
          npudist_fit(tdat = xdat[indices, , drop = FALSE], edat = exdat, bws = bws)$dist
        }
      } else {
        function(tsb) {
          npudist_fit(tdat = tsb, edat = exdat, bws = bws)$dist
        }
      }

      boot.out <- .npRmpi_with_local_bootstrap({
        if (is.inid) {
          boot(data = data.frame(xdat), statistic = boofun,
               R = plot.errors.boot.num)
        } else {
          tsboot(tseries = data.frame(xdat), statistic = boofun,
                 R = plot.errors.boot.num,
                 l = plot.errors.boot.blocklen,
                 sim = plot.errors.boot.method)
        }
      })
    }

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.conbandwidth =
  function(xdat, ydat,
           exdat, eydat,
           cdf,
           quantreg,
           tau,
           gradients,
           gradient.index,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.conbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.conbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    npqreg_fit <- .npRmpi_bootstrap_estimator("npqreg.conbandwidth")
    npcdist_fit <- .npRmpi_bootstrap_estimator("npcdist")
    npcdens_fit <- .npRmpi_bootstrap_estimator("npcdens.conbandwidth")
    exdat = toFrame(exdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    tboo =
      if(quantreg) "quant"
      else if (cdf) "dist"
      else "dens"

    if (!identical(tboo, "quant")) {
      .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "conditional density/distribution estimators")
    }

    is.inid = plot.errors.boot.method=="inid"
    fast.inid <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.inid) &&
      isTRUE(!quantreg) &&
      isTRUE(!gradients) &&
      isTRUE(identical(bws$type, "fixed"))

    fit.cond <- function(tx, ty) {
      switch(
        tboo,
        quant = npqreg_fit(txdat = tx, tydat = ty, exdat = exdat, tau = tau, bws = bws, gradients = gradients),
        dist = npcdist_fit(txdat = tx, tydat = ty, exdat = exdat, eydat = eydat, bws = bws, gradients = gradients),
        dens = npcdens_fit(txdat = tx, tydat = ty, exdat = exdat, eydat = eydat, bws = bws, gradients = gradients)
      )
    }
    out.cond <- function(fit) {
      switch(
        tboo,
        quant = if (gradients) fit$yqgrad[, gradient.index] else fit$quantile,
        dist = if (gradients) fit$congrad[, gradient.index] else fit$condist,
        dens = if (gradients) fit$congrad[, gradient.index] else fit$condens
      )
    }
    boot.out <- NULL
    if (fast.inid) {
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_ksum_conditional(
            xdat = xdat,
            ydat = ydat,
            exdat = exdat,
            eydat = eydat,
            bws = bws,
            B = plot.errors.boot.num,
            cdf = cdf
          ),
          error = function(e) {
            warning(sprintf("inid ksum fast path failed in compute.bootstrap.errors.conbandwidth (%s); using bootstrap fallback",
                            conditionMessage(e)))
            NULL
          }
        )
      })
    }

    if (is.null(boot.out)) {
      boofun <- if (is.inid) {
        function(data, indices) out.cond(fit.cond(
          tx = xdat[indices, , drop = FALSE],
          ty = ydat[indices, , drop = FALSE]
        ))
      } else {
        function(tsb) out.cond(fit.cond(
          tx = tsb[, seq_len(ncol(xdat)), drop = FALSE],
          ty = tsb[, (ncol(xdat) + 1L):ncol(tsb), drop = FALSE]
        ))
      }
      boot.out <- .npRmpi_with_local_bootstrap({
        if (is.inid) {
          boot(data = data.frame(xdat, ydat), statistic = boofun,
               R = plot.errors.boot.num)
        } else {
          tsboot(tseries = data.frame(xdat, ydat), statistic = boofun,
                 R = plot.errors.boot.num,
                 l = plot.errors.boot.blocklen,
                 sim = plot.errors.boot.method)
        }
      })
    }

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$ydati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.condbandwidth =
  function(xdat, ydat,
           exdat, eydat,
           cdf,
           quantreg,
           tau,
           gradients,
           gradient.index,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.condbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.condbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    npqreg_fit <- .npRmpi_bootstrap_estimator("npqreg.condbandwidth")
    npcdist_fit <- .npRmpi_bootstrap_estimator("npcdist.condbandwidth")
    npcdens_fit <- .npRmpi_bootstrap_estimator("npcdens")
    exdat = toFrame(exdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    tboo =
      if(quantreg) "quant"
      else if (cdf) "dist"
      else "dens"

    if (!identical(tboo, "quant")) {
      .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "conditional density/distribution estimators")
    }

    is.inid = plot.errors.boot.method=="inid"
    fast.inid <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.inid) &&
      isTRUE(!quantreg) &&
      isTRUE(!gradients) &&
      isTRUE(identical(bws$type, "fixed"))

    fit.cond <- function(tx, ty) {
      switch(
        tboo,
        quant = npqreg_fit(txdat = tx, tydat = ty, exdat = exdat, tau = tau, bws = bws, gradients = gradients),
        dist = npcdist_fit(txdat = tx, tydat = ty, exdat = exdat, eydat = eydat, bws = bws, gradients = gradients),
        dens = npcdens_fit(txdat = tx, tydat = ty, exdat = exdat, eydat = eydat, bws = bws, gradients = gradients)
      )
    }
    out.cond <- function(fit) {
      switch(
        tboo,
        quant = if (gradients) fit$yqgrad[, gradient.index] else fit$quantile,
        dist = if (gradients) fit$congrad[, gradient.index] else fit$condist,
        dens = if (gradients) fit$congrad[, gradient.index] else fit$condens
      )
    }
    boot.out <- NULL
    if (fast.inid) {
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_ksum_conditional(
            xdat = xdat,
            ydat = ydat,
            exdat = exdat,
            eydat = eydat,
            bws = bws,
            B = plot.errors.boot.num,
            cdf = cdf
          ),
          error = function(e) {
            warning(sprintf("inid ksum fast path failed in compute.bootstrap.errors.condbandwidth (%s); using bootstrap fallback",
                            conditionMessage(e)))
            NULL
          }
        )
      })
    }

    if (is.null(boot.out)) {
      boofun <- if (is.inid) {
        function(data, indices) out.cond(fit.cond(
          tx = xdat[indices, , drop = FALSE],
          ty = ydat[indices, , drop = FALSE]
        ))
      } else {
        function(tsb) out.cond(fit.cond(
          tx = tsb[, seq_len(ncol(xdat)), drop = FALSE],
          ty = tsb[, (ncol(xdat) + 1L):ncol(tsb), drop = FALSE]
        ))
      }
      boot.out <- .npRmpi_with_local_bootstrap({
        if (is.inid) {
          boot(data = data.frame(xdat, ydat), statistic = boofun,
               R = plot.errors.boot.num)
        } else {
          tsboot(tseries = data.frame(xdat, ydat), statistic = boofun,
                 R = plot.errors.boot.num,
                 l = plot.errors.boot.blocklen,
                 sim = plot.errors.boot.method)
        }
      })
    }

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$ydati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.sibandwidth =
  function(xdat, ydat,
           gradients,
           plot.errors.boot.method,
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.sibandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.sibandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(xdat)
    )

    npindex_fit <- .npRmpi_bootstrap_estimator("npindex.sibandwidth")
    npindexhat_fit <- .npRmpi_bootstrap_estimator("npindexhat")
    boot.err = matrix(data = NA, nrow = nrow(xdat), ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method=="inid"

    if (is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.train <- npindex_fit(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        gradients = FALSE
      )
      H <- npindexhat_fit(
        bws = bws,
        txdat = xdat,
        exdat = xdat,
        output = "matrix",
        s = if (gradients) 1L else 0L
      )

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - as.vector(fit.train$mean))
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- list(
        t = .np_wild_boot_t(
          H = H,
          fit.mean = as.vector(fit.train$mean),
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
    } else if (is.inid) {
      inid.helper.ok <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
        !isTRUE(gradients) &&
        identical(bws$type, "fixed")
      if (!isTRUE(inid.helper.ok)) {
        warning("inid single-index helper unavailable for this configuration; using explicit bootstrap fallback")
        boot.out <- NULL
      } else {
        boot.out <- .npRmpi_bootstrap_maybe_local({
          tryCatch({
            tx.index <- data.frame(index = as.vector(toMatrix(xdat) %*% bws$beta))
            rbw <- .np_indexhat_rbw(bws = bws, idx.train = tx.index)
            .np_inid_boot_from_regression(
              xdat = tx.index,
              exdat = tx.index,
              bws = rbw,
              ydat = ydat,
              B = plot.errors.boot.num
            )
          }, error = function(e) {
            stop(sprintf("inid single-index helper failed in compute.bootstrap.errors.sibandwidth (%s)",
                         conditionMessage(e)),
                 call. = FALSE)
          })
        }, estimators = list(npindex_fit, npindexhat_fit))
      }
    }

    if (is.null(boot.out)) {
      ## beta[1] is always 1.0, so use first column of gradients matrix ...
      if (is.inid) {
        boofun.inid <- function(data, indices) {
          fit <- npindex_fit(
            txdat = data[indices, seq_len(ncol(data) - 1L), drop = FALSE],
            tydat = data[indices, ncol(data), drop = TRUE],
            exdat = xdat, bws = bws,
            gradients = gradients
          )
          if (gradients) fit$grad[,1] else fit$mean
        }
        boot.out <- .npRmpi_bootstrap_maybe_local({
          boot(
            data = data.frame(xdat, ydat),
            statistic = boofun.inid,
            R = plot.errors.boot.num
          )
        }, estimators = list(npindex_fit, npindexhat_fit))
      } else {
        boofun <- function(tsb) {
          fit <- npindex_fit(
            txdat = tsb[, seq_len(ncol(tsb) - 1L), drop = FALSE],
            tydat = tsb[, ncol(tsb)],
            exdat = xdat, bws = bws,
            gradients = gradients
          )
          if (gradients) fit$grad[,1] else fit$mean
        }
        boot.out <- .npRmpi_bootstrap_maybe_local({
          tsboot(tseries = data.frame(xdat, ydat), statistic = boofun,
                 R = plot.errors.boot.num,
                 l = plot.errors.boot.blocklen,
                 sim = plot.errors.boot.method)
        }, estimators = list(npindex_fit, npindexhat_fit))
      }
    }
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = list(),
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }


uocquantile <- function(x, prob) {
  if(anyNA(prob)) stop("'prob' contains missing values")
  if(any(prob < 0 | prob > 1, na.rm = TRUE)) stop("'prob' outside [0,1]")
  if(anyNA(x)) stop("missing values and NaN's not allowed")
  if (is.ordered(x)){
    x <- droplevels(x)
    tq = unclass(table(x))
    tq = tq / sum(tq)
    tq[length(tq)] <- 1.0
    bscape <- levels(x)
    tq <- cumsum(tq)
    j <- sapply(prob, function(p){ which(tq >= p)[1] })
    bscape[j]
  } else if (is.factor(x)) {
    ## just returns mode
    x <- droplevels(x)
    tq = unclass(table(x))
    j = which(tq == max(tq))[1]
    levels(x)[j]
  } else {
    quantile(x, probs = prob)
  }
}


trim.quantiles <- function(dat, trim){
  if (sign(trim) == sign(-1)){
    trim <- abs(trim)
    tq <- quantile(dat, probs = c(0.0, 0.0+trim, 1.0-trim,1.0))
    tq <- c(2.0*tq[1]-tq[2], 2.0*tq[4]-tq[3])
  }
  else {
    tq <- quantile(dat, probs = c(0.0+trim, 1.0-trim))
  }
  tq
}
