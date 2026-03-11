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

.np_plot_progress_enabled <- function() {
  isTRUE(getOption("np.messages", TRUE)) &&
    isTRUE(getOption("np.plot.progress", TRUE)) &&
    isTRUE(.np_progress_is_interactive())
}

.np_plot_progress_interval_sec <- function() {
  val <- suppressWarnings(as.numeric(getOption("np.plot.progress.interval.sec", 0.5))[1L])
  if (!is.finite(val) || is.na(val) || val < 0)
    val <- 0.5
  val
}

.np_plot_progress_begin <- function(total, label) {
  total <- as.integer(total)
  if (is.na(total) || total < 1L || !.np_plot_progress_enabled())
    return(NULL)

  label <- as.character(label)[1L]
  .np_progress_note(paste0(label, "..."))

  state <- .np_progress_begin(label = label, total = total, domain = "plot")
  state$enabled <- isTRUE(.np_plot_progress_enabled())
  state$throttle_sec <- .np_plot_progress_interval_sec()
  state$last_emit <- state$started - state$throttle_sec
  state
}

.np_plot_progress_tick <- function(state, done, force = FALSE) {
  if (is.null(state))
    return(state)

  done <- as.integer(done)
  if (is.na(done))
    done <- 0L
  done <- max(0L, min(state$total, done))

  if (isTRUE(force))
    state$last_emit <- -Inf

  .np_progress_step(state = state, done = done)
}

.np_plot_progress_end <- function(state) {
  if (is.null(state))
    return(invisible(NULL))

  if (isTRUE(state$known_total) && identical(state$last_done, state$total))
    return(invisible(NULL))

  state$last_emit <- -Inf
  .np_progress_end(state)
  invisible(NULL)
}

.np_plot_activity_begin <- function(label) {
  if (!.np_plot_progress_enabled())
    return(NULL)

  .np_progress_note(as.character(label)[1L])
  TRUE
}

.np_plot_activity_end <- function(console) {
  invisible(NULL)
}

.np_plot_capture_par <- function(names = character()) {
  names <- unique(as.character(names))
  if (!length(names))
    return(list())
  par(names)
}

.np_plot_restore_par <- function(oldpar, reset.new = TRUE) {
  if (isTRUE(reset.new))
    suppressWarnings(try(par(new = FALSE), silent = TRUE))
  if (!is.null(oldpar) && length(oldpar))
    suppressWarnings(try(par(oldpar), silent = TRUE))
  invisible(NULL)
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

  chunk.size <- .np_wild_chunk_size(n = n, B = B)
  out <- matrix(NA_real_, nrow = B, ncol = nrow(H))
  fit.mean <- as.double(fit.mean)
  residuals <- as.double(residuals)
  wild <- .np_plot_normalize_wild(wild)
  draw.fun <- if (identical(wild, "mammen")) .np_mammen_draws else .np_rademacher_draws
  progress <- .np_plot_progress_begin(total = B, label = "Plot bootstrap wild")
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    draws <- draw.fun(n = n, B = bsz)
    ystar <- residuals * draws
    ystar <- ystar + fit.mean
    out[start:stopi, ] <- t(H %*% ystar)
    progress <- .np_plot_progress_tick(state = progress, done = stopi)
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

.np_plot_normalize_wild <- function(wild = c("rademacher", "mammen")) {
  if (length(wild) > 1L)
    wild <- wild[1L]
  match.arg(wild, c("mammen", "rademacher"))
}

.np_plot_boot_from_hat_wild <- function(H, ydat, fit.mean, B, wild) {
  fit.mean <- as.vector(fit.mean)
  list(
    t = .np_wild_boot_t(
      H = H,
      fit.mean = fit.mean,
      residuals = as.double(ydat - fit.mean),
      B = as.integer(B),
      wild = wild
    ),
    t0 = as.vector(H %*% as.double(ydat))
  )
}

.np_plot_require_bws <- function(bws, where) {
  if (is.null(bws))
    stop(sprintf("required argument 'bws' is missing or NULL in %s", where))
  invisible(TRUE)
}

.np_plot_boot_factor_boxplots <- function(boot.t, tdati, ti, B) {
  all.bp <- list()
  ti <- as.integer(ti)[1L]
  if (is.na(ti) || ti < 1L)
    return(all.bp)
  if (ti > length(tdati$iord) || ti > length(tdati$iuno))
    return(all.bp)
  if (!(isTRUE(tdati$iord[ti]) || isTRUE(tdati$iuno[ti])))
    return(all.bp)

  boot.frame <- as.data.frame(boot.t)
  u.lev <- tdati$all.ulev[[ti]]
  stopifnot(length(u.lev) == ncol(boot.frame))

  all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
  all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

  for (i in seq_along(u.lev)) {
    t.bp <- boxplot.stats(boot.frame[, i])
    all.bp$stats[, i] <- t.bp$stats
    all.bp$conf[, i] <- t.bp$conf
    all.bp$out <- c(all.bp$out, t.bp$out)
    all.bp$group <- c(all.bp$group, rep.int(i, length(t.bp$out)))
  }

  all.bp$n <- rep.int(as.integer(B), length(u.lev))
  all.bp$names <- tdati$all.lev[[ti]]
  all.bp
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

.np_counts_to_indices <- function(counts.col) {
  counts.col <- as.integer(counts.col)
  rep.int(seq_along(counts.col), counts.col)
}

.np_block_counts_drawer <- function(n,
                                    B,
                                    blocklen,
                                    sim = c("fixed", "geom"),
                                    n.sim = n,
                                    endcorr = TRUE) {
  sim <- match.arg(sim)
  n <- as.integer(n)
  B <- as.integer(B)
  n.sim <- as.integer(n.sim)
  blocklen <- as.integer(blocklen)

  if (n < 1L || B < 1L || n.sim < 1L)
    stop("invalid block bootstrap dimensions")
  if (length(blocklen) != 1L || is.na(blocklen) || blocklen < 1L || blocklen > n)
    stop("invalid block length for block bootstrap")

  ts.array <- utils::getFromNamespace("ts.array", "boot")
  make.ends <- utils::getFromNamespace("make.ends", "boot")
  ts.draws <- ts.array(
    n = n,
    n.sim = n.sim,
    R = B,
    l = blocklen,
    sim = sim,
    endcorr = isTRUE(endcorr)
  )

  starts <- as.matrix(ts.draws$starts)
  lengths <- ts.draws$lengths

  function(start, stopi) {
    start <- as.integer(start)
    stopi <- as.integer(stopi)
    if (start < 1L || stopi < start || stopi > B)
      stop("invalid block bootstrap chunk bounds")

    idx <- seq.int(start, stopi)
    out <- matrix(0.0, nrow = n, ncol = length(idx))

    for (jj in seq_along(idx)) {
      rr <- idx[jj]
      ends <- if (identical(sim, "geom")) {
        cbind(starts[rr, ], lengths[rr, ])
      } else {
        cbind(starts[rr, ], lengths)
      }

      inds <- apply(ends, 1L, make.ends, n)
      inds <- if (is.list(inds)) {
        as.integer(unlist(inds)[seq_len(n.sim)])
      } else {
        as.integer(inds)[seq_len(n.sim)]
      }
      out[, jj] <- tabulate(inds, nbins = n)
    }

    out
  }
}

.np_inid_lc_boot_from_hat <- function(H, ydat, B, counts = NULL, counts.drawer = NULL) {
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

  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  prob <- rep.int(1 / n, n)
  tmat <- matrix(NA_real_, nrow = B, ncol = nrow(H))
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- if (!is.null(counts.drawer)) {
      .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
    } else {
      stats::rmultinom(n = bsz, size = n, prob = prob)
    }
    den <- crossprod(counts.chunk, W)
    num <- crossprod(counts.chunk, Wy)
    tmat[start:stopi, ] <- num / pmax(den, .Machine$double.eps)
    progress <- .np_plot_progress_tick(state = progress, done = stopi)
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

.np_inid_lp_predict_row <- function(A, z, rhs, ridge.grid) {
  ridge.grid <- as.double(ridge.grid)
  if (!length(ridge.grid) || anyNA(ridge.grid) || any(!is.finite(ridge.grid)) || any(ridge.grid < 0))
    stop("argument 'ridge.grid' must be a non-empty non-negative finite numeric vector")
  z <- as.double(z)
  if (!length(z))
    return(NA_real_)
  z1 <- z[1L]

  for (ridge in ridge.grid) {
    Ar <- A
    zr <- z
    if (ridge > 0)
      diag(Ar) <- diag(Ar) + ridge
    if (ridge > 0)
      zr[1L] <- z1 + ridge * z1 / NZD(Ar[1L, 1L])
    beta <- tryCatch(
      drop(solve(Ar, matrix(zr, ncol = 1L))),
      error = function(e) NULL
    )
    if (!is.null(beta) && all(is.finite(beta)))
      return(sum(rhs * beta))
  }

  NA_real_
}

.np_inid_lp_predict_chunk_general <- function(Mvals, Zvals, rhs, ridge.grid) {
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
      ridge.grid = ridge.grid
    )
  }

  out
}

.np_inid_lp_predict_chunk <- function(Mvals, Zvals, rhs, ridge.grid) {
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
        ridge.grid = ridge.grid
      )
    }
  }

  out
}

.np_inid_boot_from_regression_localpoly_fixed <- function(xdat,
                                                          exdat,
                                                          bws,
                                                          ydat,
                                                          B,
                                                          counts = NULL,
                                                          counts.drawer = NULL,
                                                          ridge = 1.0e-12,
                                                          gradients = FALSE,
                                                          gradient.order = 1L,
                                                          slice.index = 1L) {
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
  if (!identical(bws$type, "fixed"))
    stop("local-polynomial fixed regression helper requires bwtype='fixed'")

  regtype <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  if (identical(regtype, "lc"))
    stop("local-polynomial fixed regression helper requires regtype='ll' or 'lp'")

  ncon <- bws$ncon
  degree <- if (identical(regtype, "ll")) {
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

  gradient.vec <- NULL
  if (isTRUE(gradients)) {
    cont.idx <- which(bws$icon)
    cpos <- match(slice.index, cont.idx)
    if (is.na(cpos))
      stop("fixed regression gradient helper requires a continuous slice", call. = FALSE)

    gradient.vec <- integer(ncon)
    if (identical(regtype, "lp")) {
      gradient.order <- npValidateGlpGradientOrder(
        regtype = "lp",
        gradient.order = gradient.order,
        ncon = ncon
      )
      if (gradient.order[cpos] > degree[cpos]) {
        return(list(
          t = matrix(NA_real_, nrow = B, ncol = neval),
          t0 = rep(NA_real_, neval)
        ))
      }
      gradient.vec[cpos] <- gradient.order[cpos]
    } else {
      gradient.vec[cpos] <- 1L
    }
  }

  kw <- npksum(
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
    gradient.vec = gradient.vec,
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
  ridge.grid <- npRidgeSequenceFromBase(n.train = n, ridge.base = ridge, cap = 1.0)

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
        ridge.grid = ridge.grid
      )[1L]
    } else {
      .np_inid_lp_predict_chunk(
        Mvals = M0,
        Zvals = Z0,
        rhs = rhs[i, ],
        ridge.grid = ridge.grid
      )[1L]
    }
  }

  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  fill_chunk <- function(counts.chunk, start, stopi) {
    for (i in seq_len(neval)) {
      Mvals <- crossprod(counts.chunk, Mfeat[[i]])
      Zvals <- crossprod(counts.chunk, Zfeat[[i]])
      tmat[start:stopi, i] <<- if (p > 3L) {
        .np_inid_lp_predict_chunk_general(
          Mvals = Mvals,
          Zvals = Zvals,
          rhs = rhs[i, ],
          ridge.grid = ridge.grid
        )
      } else {
        .np_inid_lp_predict_chunk(
          Mvals = Mvals,
          Zvals = Zvals,
          rhs = rhs[i, ],
          ridge.grid = ridge.grid
        )
      }
    }
  }

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    fill_chunk(counts.chunk = counts.mat, start = 1L, stopi = B)
    progress <- .np_plot_progress_tick(state = progress, done = B, force = TRUE)
  } else {
    chunk.size <- .np_inid_chunk_size(n = n, B = B)
    start <- 1L
    while (start <= B) {
      stopi <- min(B, start + chunk.size - 1L)
      bsz <- stopi - start + 1L
      counts.chunk <- if (!is.null(counts.drawer)) {
        .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
      } else {
        stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
      }
      fill_chunk(counts.chunk = counts.chunk, start = start, stopi = stopi)
      progress <- .np_plot_progress_tick(state = progress, done = stopi)
      start <- stopi + 1L
    }
  }

  if (any(!is.finite(t0)) || any(!is.finite(tmat)))
    stop("inid regression helper path produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_regression <- function(xdat,
                                          exdat,
                                          bws,
                                          ydat,
                                          B,
                                          counts = NULL,
                                          counts.drawer = NULL,
                                          ridge = 1.0e-12,
                                          gradients = FALSE,
                                          gradient.order = 1L,
                                          slice.index = 1L) {
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
  if (!identical(bws$type, "fixed")) {
    return(.np_inid_boot_from_reghat_exact(
      xdat = xdat,
      exdat = exdat,
      bws = bws,
      ydat = ydat,
      B = B,
      counts = counts,
      counts.drawer = counts.drawer,
      gradients = gradients,
      gradient.order = gradient.order,
      slice.index = slice.index
    ))
  }

  if (isTRUE(gradients)) {
    if (!identical(regtype, "lc")) {
      return(.np_inid_boot_from_regression_localpoly_fixed(
        xdat = xdat,
        exdat = exdat,
        bws = bws,
        ydat = ydat,
        B = B,
        counts = counts,
        counts.drawer = counts.drawer,
        ridge = ridge,
        gradients = TRUE,
        gradient.order = gradient.order,
        slice.index = slice.index
      ))
    }

    return(.np_inid_boot_from_reghat_exact(
      xdat = xdat,
      exdat = exdat,
      bws = bws,
      ydat = ydat,
      B = B,
      counts = counts,
      counts.drawer = counts.drawer,
      gradients = TRUE,
      gradient.order = gradient.order,
      slice.index = slice.index
    ))
  }

  if (identical(regtype, "lc")) {
    H <- suppressWarnings(
      tryCatch(
        npreghat.rbandwidth(
          bws = bws,
          txdat = xdat,
          exdat = exdat,
          s = 0L,
          output = "matrix"
        ),
        error = function(e) NULL
      )
    )
    if (!is.null(H)) {
      if (!is.matrix(H))
        H <- matrix(as.double(H), nrow = neval, ncol = n)
      if (nrow(H) == neval && ncol(H) == n) {
        return(.np_inid_lc_boot_from_hat(
          H = H,
          ydat = ydat,
          B = B,
          counts = counts,
          counts.drawer = counts.drawer
        ))
      }
    }
  }

  .np_inid_boot_from_regression_localpoly_fixed(
    xdat = xdat,
    exdat = exdat,
    bws = bws,
    ydat = ydat,
    B = B,
    counts = counts,
    counts.drawer = counts.drawer,
    ridge = ridge,
    gradients = FALSE,
    gradient.order = gradient.order,
    slice.index = slice.index
  )
}

.np_inid_boot_from_reghat_exact <- function(xdat,
                                            exdat,
                                            bws,
                                            ydat,
                                            B,
                                            counts = NULL,
                                            counts.drawer = NULL,
                                            gradients = FALSE,
                                            gradient.order = 1L,
                                            slice.index = 1L) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  ydat <- as.double(ydat)
  B <- as.integer(B)

  n <- nrow(xdat)
  neval <- nrow(exdat)
  if (length(ydat) != n)
    stop("length of ydat must match training rows in exact regression bootstrap helper")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid exact regression bootstrap dimensions")

  fit_hat <- function(x.train, y.train) {
    fit <- .np_regression_direct(
      bws = bws,
      txdat = x.train,
      tydat = y.train,
      exdat = exdat,
      gradients = isTRUE(gradients),
      gradient.order = gradient.order
    )
    if (isTRUE(gradients))
      as.vector(fit$grad[, slice.index])
    else
      as.vector(fit$mean)
  }

  t0 <- fit_hat(x.train = xdat, y.train = ydat)
  tmat <- matrix(NA_real_, nrow = B, ncol = length(t0))
  counts.mat <- if (!is.null(counts)) .np_inid_counts_matrix(n = n, B = B, counts = counts) else NULL
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  start <- 1L
  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- if (!is.null(counts.mat)) {
      counts.mat[, start:stopi, drop = FALSE]
    } else if (!is.null(counts.drawer)) {
      .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
    } else {
      stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    }

    for (jj in seq_len(bsz)) {
      idx <- .np_counts_to_indices(counts.chunk[, jj])
      tmat[start + jj - 1L, ] <- fit_hat(
        x.train = xdat[idx, , drop = FALSE],
        y.train = ydat[idx]
      )
    }

    progress <- .np_plot_progress_tick(state = progress, done = stopi)
    start <- stopi + 1L
  }

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_index <- function(xdat,
                                     ydat,
                                     bws,
                                     B,
                                     counts = NULL,
                                     counts.drawer = NULL,
                                     gradients = FALSE) {
  if (!identical(bws$type, "fixed")) {
    return(.np_inid_boot_from_index_exact(
      xdat = xdat,
      ydat = ydat,
      bws = bws,
      B = B,
      counts = counts,
      counts.drawer = counts.drawer,
      gradients = gradients
    ))
  }

  if (isTRUE(gradients)) {
    return(.np_inid_boot_from_index_gradient_fixed(
      xdat = xdat,
      ydat = ydat,
      bws = bws,
      B = B,
      counts = counts,
      counts.drawer = counts.drawer
    ))
  }

  xdat <- toFrame(xdat)
  B <- as.integer(B)
  n <- nrow(xdat)

  if (length(ydat) != n)
    stop("length of ydat must match training rows in single-index inid helper")
  if (n < 1L || B < 1L)
    stop("invalid single-index inid helper dimensions")
  if (isTRUE(gradients))
    stop("single-index inid helper does not support gradients", call. = FALSE)

  idx.train <- data.frame(index = as.vector(toMatrix(xdat) %*% bws$beta))
  y.num <- if (is.factor(ydat)) {
    yadj <- adjustLevels(data.frame(ydat), bws$ydati)
    bws$ydati$all.dlev[[1L]][as.integer(yadj[, 1L])]
  } else {
    as.double(ydat)
  }

  spec <- .npindex_resolve_spec(bws, where = "single-index inid helper")
  regtype <- spec$regtype.engine
  kbw <- .np_indexhat_kbw(bws = bws, idx.train = idx.train)
  kw <- .np_kernel_weights_direct(
    bws = kbw,
    txdat = idx.train,
    exdat = idx.train,
    bandwidth.divide = TRUE,
    kernel.pow = 1.0
  )

  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = n)
  if (nrow(kw) != n || ncol(kw) != n)
    stop("single-index inid helper kernel-weight matrix shape mismatch")

  if (identical(regtype, "lc")) {
    H <- sweep(t(kw), 1L, pmax(colSums(kw), .Machine$double.eps), "/",
               check.margin = FALSE)
    return(.np_inid_lc_boot_from_hat(
      H = H,
      ydat = y.num,
      B = B,
      counts = counts,
      counts.drawer = counts.drawer
    ))
  }

  degree <- if (identical(regtype, "ll")) {
    1L
  } else {
    spec$degree.engine
  }

  W <- W.lp(
    xdat = idx.train,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  W.eval <- W.lp(
    xdat = idx.train,
    exdat = idx.train,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  W <- as.matrix(W)
  W.eval <- as.matrix(W.eval)

  p <- ncol(W)
  mcols <- p * (p + 1L) / 2L
  ridge.grid <- npRidgeSequenceFromBase(n.train = n, ridge.base = 1.0e-12, cap = 1.0)
  rhs <- W.eval
  ones <- matrix(1.0, nrow = n, ncol = 1L)

  Mfeat <- vector("list", n)
  Zfeat <- vector("list", n)
  t0 <- numeric(n)

  for (i in seq_len(n)) {
    k <- as.double(kw[, i])
    WK <- W * k
    Zfeat[[i]] <- WK * y.num

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
        ridge.grid = ridge.grid
      )[1L]
    } else {
      .np_inid_lp_predict_chunk(
        Mvals = M0,
        Zvals = Z0,
        rhs = rhs[i, ],
        ridge.grid = ridge.grid
      )[1L]
    }
  }

  tmat <- matrix(NA_real_, nrow = B, ncol = n)
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  fill_chunk <- function(counts.chunk, start, stopi) {
    for (i in seq_len(n)) {
      Mvals <- crossprod(counts.chunk, Mfeat[[i]])
      Zvals <- crossprod(counts.chunk, Zfeat[[i]])
      tmat[start:stopi, i] <<- if (p > 3L) {
        .np_inid_lp_predict_chunk_general(
          Mvals = Mvals,
          Zvals = Zvals,
          rhs = rhs[i, ],
          ridge.grid = ridge.grid
        )
      } else {
        .np_inid_lp_predict_chunk(
          Mvals = Mvals,
          Zvals = Zvals,
          rhs = rhs[i, ],
          ridge.grid = ridge.grid
        )
      }
    }
  }

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    fill_chunk(counts.chunk = counts.mat, start = 1L, stopi = B)
    progress <- .np_plot_progress_tick(state = progress, done = B, force = TRUE)
  } else {
    chunk.size <- .np_inid_chunk_size(n = n, B = B)
    start <- 1L
    while (start <= B) {
      stopi <- min(B, start + chunk.size - 1L)
      bsz <- stopi - start + 1L
      counts.chunk <- if (!is.null(counts.drawer)) {
        .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
      } else {
        stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
      }
      fill_chunk(counts.chunk = counts.chunk, start = start, stopi = stopi)
      progress <- .np_plot_progress_tick(state = progress, done = stopi)
      start <- stopi + 1L
    }
  }

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_index_gradient_fixed <- function(xdat,
                                                    ydat,
                                                    bws,
                                                    B,
                                                    counts = NULL,
                                                    counts.drawer = NULL) {
  xdat <- toFrame(xdat)
  B <- as.integer(B)
  n <- nrow(xdat)

  if (length(ydat) != n)
    stop("length of ydat must match training rows in single-index fixed gradient helper")
  if (n < 1L || B < 1L)
    stop("invalid single-index fixed gradient helper dimensions")

  idx.train <- data.frame(index = as.vector(toMatrix(xdat) %*% bws$beta))
  rbw <- .np_indexhat_rbw(bws = bws, idx.train = idx.train)

  .np_inid_boot_from_regression(
    xdat = idx.train,
    exdat = idx.train,
    bws = rbw,
    ydat = ydat,
    B = B,
    counts = counts,
    counts.drawer = counts.drawer,
    gradients = TRUE,
    gradient.order = 1L,
    slice.index = 1L
  )
}

.np_inid_boot_from_index_exact <- function(xdat,
                                           ydat,
                                           bws,
                                           B,
                                           counts = NULL,
                                           counts.drawer = NULL,
                                           gradients = FALSE) {
  xdat <- toFrame(xdat)
  B <- as.integer(B)
  n <- nrow(xdat)

  if (length(ydat) != n)
    stop("length of ydat must match training rows in exact single-index bootstrap helper")
  if (n < 1L || B < 1L)
    stop("invalid exact single-index bootstrap dimensions")
  if (isTRUE(gradients))
    stop("exact single-index bootstrap helper does not support gradients", call. = FALSE)

  fit_hat <- function(x.train, y.train) {
    as.vector(npindexhat(
      bws = bws,
      txdat = x.train,
      exdat = xdat,
      y = y.train,
      output = "apply",
      s = 0L
    ))
  }

  t0 <- fit_hat(x.train = xdat, y.train = ydat)
  tmat <- matrix(NA_real_, nrow = B, ncol = length(t0))
  counts.mat <- if (!is.null(counts)) .np_inid_counts_matrix(n = n, B = B, counts = counts) else NULL
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  start <- 1L
  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- if (!is.null(counts.mat)) {
      counts.mat[, start:stopi, drop = FALSE]
    } else if (!is.null(counts.drawer)) {
      .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
    } else {
      stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    }

    for (jj in seq_len(bsz)) {
      idx <- .np_counts_to_indices(counts.chunk[, jj])
      tmat[start + jj - 1L, ] <- fit_hat(
        x.train = xdat[idx, , drop = FALSE],
        y.train = ydat[idx]
      )
    }

    progress <- .np_plot_progress_tick(state = progress, done = stopi)
    start <- stopi + 1L
  }

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

.np_inid_scoef_predict_row <- function(mrow, zrow, rhs, ridge.grid) {
  A <- .np_inid_lp_unpack_sym_row(mrow = mrow, p = length(rhs))
  tyw <- as.double(zrow)
  nc <- ncol(A)
  ridge.grid <- as.double(ridge.grid)
  if (!length(ridge.grid) || anyNA(ridge.grid) || any(!is.finite(ridge.grid)) || any(ridge.grid < 0))
    stop("argument 'ridge.grid' must be a non-empty non-negative finite numeric vector")

  maxPenalty <- sqrt(.Machine$double.xmax)
  coef.ii <- rep(maxPenalty, nc)
  for (ridge in ridge.grid) {
    ridge.val <- ridge * tyw[1L] / NZD(A[1L, 1L])
    coef.try <- tryCatch(
      solve(
        A + diag(rep(ridge, nc)),
        tyw + c(ridge.val, rep(0, nc - 1L))
      ),
      error = function(e) e
    )
    if (!(inherits(coef.try, "error") || any(!is.finite(coef.try))))
      return(sum(as.double(rhs) * as.double(coef.try)))
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
                                     counts.drawer = NULL,
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

  kw <- .npscoef_weight_matrix(
    bws = bws,
    tzdat = tzdat,
    ezdat = ezdat,
    leave.one.out = leave.one.out
  )
  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(txdat))

  n <- nrow(W.train)
  neval <- nrow(W.eval)
  if (nrow(kw) != n || ncol(kw) != neval)
    stop("smooth coefficient inid helper kernel-weight matrix shape mismatch")

  p <- ncol(W.train)
  mcols <- p * (p + 1L) / 2L
  ones <- matrix(1.0, nrow = n, ncol = 1L)
  ridge.grid <- npRidgeSequenceAdditive(n.train = n, cap = 1.0)

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
        ridge.grid = ridge.grid
      )
    }
    t0[i] <- t0i
  }

  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

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
            ridge.grid = ridge.grid
          )
        }
      }
      tmat[start:stopi, i] <<- out
    }
  }

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    fill_chunk(counts.chunk = counts.mat, start = 1L, stopi = B)
    progress <- .np_plot_progress_tick(state = progress, done = B, force = TRUE)
  } else {
    chunk.size <- .np_inid_chunk_size(n = n, B = B)
    start <- 1L
    while (start <= B) {
      stopi <- min(B, start + chunk.size - 1L)
      bsz <- stopi - start + 1L
      counts.chunk <- if (!is.null(counts.drawer)) {
        .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
      } else {
        stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
      }
      fill_chunk(counts.chunk = counts.chunk, start = start, stopi = stopi)
      progress <- .np_plot_progress_tick(state = progress, done = stopi)
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
    stop("plreg inid helper path requires at least one linear regressor")
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

  XtWX <- crossprod(X, X * w)
  XtWy <- drop(crossprod(X, y * w))
  XtWy0 <- XtWy[1L]
  ridge.grid <- npRidgeSequenceFromBase(
    n.train = nrow(X),
    ridge.base = max(0.0, as.double(ridge)),
    cap = 1.0
  )

  for (ridge.try in ridge.grid) {
    A <- XtWX
    ridge <- ridge.try
    z <- XtWy
    if (ridge > 0)
      diag(A) <- diag(A) + ridge
    if (ridge > 0)
      z[1L] <- XtWy0 + ridge * XtWy0 / NZD(A[1L, 1L])
    beta <- tryCatch(
      drop(solve(A, matrix(z, ncol = 1L))),
      error = function(e) NULL
    )
    if (!is.null(beta) && all(is.finite(beta)))
      return(as.double(beta))
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
                                     counts.drawer = NULL,
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
    stop("plreg inid helper path requires aligned txdat/tzdat rows")
  if (nrow(ezdat) != neval)
    stop("plreg inid helper path requires aligned exdat/ezdat rows")
  if (n < 1L || neval < 1L || p < 1L || B < 1L)
    stop("invalid plreg inid helper path dimensions")

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

  counts.mat <- if (!is.null(counts)) {
    .np_inid_counts_matrix(n = n, B = B, counts = counts)
  } else if (!is.null(counts.drawer)) {
    .np_inid_counts_matrix(n = n, B = B, counts = counts.drawer(1L, B))
  } else {
    .np_inid_counts_matrix(n = n, B = B)
  }

  y.train <- .np_inid_boot_from_regression(
    xdat = tzdat,
    exdat = tzdat,
    bws = bws$bw$yzbw,
    ydat = y.num,
    B = B,
    counts = counts.mat,
    counts.drawer = counts.drawer,
    ridge = ridge
  )
  y.eval <- .np_inid_boot_from_regression(
    xdat = tzdat,
    exdat = ezdat,
    bws = bws$bw$yzbw,
    ydat = y.num,
    B = B,
    counts = counts.mat,
    counts.drawer = counts.drawer,
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
      counts.drawer = counts.drawer,
      ridge = ridge
    )
    x.eval[[j]] <- .np_inid_boot_from_regression(
      xdat = tzdat,
      exdat = ezdat,
      bws = bws$bw[[j + 1L]],
      ydat = x.train.num[, j],
      B = B,
      counts = counts.mat,
      counts.drawer = counts.drawer,
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
    stop("plreg inid helper path produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_boot_matrix_from_ksum <- function(ksum, B, nout, where = "ksum helper path") {
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

.np_ksum_unconditional_eval_exact <- function(xdat, exdat, bws, operator) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  ones <- matrix(1.0, nrow = nrow(xdat), ncol = 1L)
  as.numeric(npksum(
    txdat = xdat,
    tydat = ones,
    exdat = exdat,
    bws = bws,
    weights = ones,
    operator = operator,
    bandwidth.divide = TRUE
  )$ksum) / nrow(xdat)
}

.np_operator_matrix_from_ksum <- function(ksum, nrow.out, ncol.out, where) {
  km <- as.matrix(ksum)

  if (nrow(km) == nrow.out && ncol(km) == ncol.out)
    return(km)
  if (nrow(km) == ncol.out && ncol(km) == nrow.out)
    return(t(km))

  stop(sprintf("%s returned unexpected operator shape", where))
}

.np_ksum_unconditional_operator_fixed <- function(xdat, exdat, bws, operator) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  n <- nrow(xdat)

  K <- npksum(
    txdat = xdat,
    exdat = exdat,
    bws = bws,
    tydat = diag(n),
    operator = operator,
    bandwidth.divide = TRUE
  )$ksum

  .np_operator_matrix_from_ksum(
    ksum = K,
    nrow.out = nrow(exdat),
    ncol.out = n,
    where = "npksum unconditional operator"
  ) / n
}

.np_apply_operator_counts <- function(K, counts) {
  t(K %*% counts)
}

.np_inid_boot_from_ksum_unconditional_exact <- function(xdat,
                                                        exdat,
                                                        bws,
                                                        B,
                                                        operator,
                                                        counts = NULL,
                                                        counts.drawer = NULL) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)

  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid unconditional exact bootstrap dimensions")

  fit_one <- function(x.train) {
    .np_ksum_unconditional_eval_exact(
      xdat = x.train,
      exdat = exdat,
      bws = bws,
      operator = operator
    )
  }

  t0 <- fit_one(x.train = xdat)
  tmat <- matrix(NA_real_, nrow = B, ncol = length(t0))
  counts.mat <- if (!is.null(counts)) .np_inid_counts_matrix(n = n, B = B, counts = counts) else NULL
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  start <- 1L
  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- if (!is.null(counts.mat)) {
      counts.mat[, start:stopi, drop = FALSE]
    } else if (!is.null(counts.drawer)) {
      .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
    } else {
      stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    }

    for (jj in seq_len(bsz)) {
      idx <- .np_counts_to_indices(counts.chunk[, jj])
      tmat[start + jj - 1L, ] <- fit_one(x.train = xdat[idx, , drop = FALSE])
    }

    progress <- .np_plot_progress_tick(state = progress, done = stopi)
    start <- stopi + 1L
  }

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_ksum_unconditional <- function(xdat,
                                                  exdat,
                                                  bws,
                                                  B,
                                                  operator,
                                                  counts = NULL,
                                                  counts.drawer = NULL) {
  if (!identical(bws$type, "fixed")) {
    return(.np_inid_boot_from_ksum_unconditional_exact(
      xdat = xdat,
      exdat = exdat,
      bws = bws,
      B = B,
      operator = operator,
      counts = counts,
      counts.drawer = counts.drawer
    ))
  }

  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)

  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid unconditional inid bootstrap dimensions")

  K <- .np_ksum_unconditional_operator_fixed(
    xdat = xdat,
    exdat = exdat,
    bws = bws,
    operator = operator
  )
  t0 <- rowSums(K)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    return(list(t = .np_apply_operator_counts(K = K, counts = counts.mat), t0 = t0))
  }

  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- if (!is.null(counts.drawer)) {
      .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
    } else {
      stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    }
    tmat[start:stopi, ] <- .np_apply_operator_counts(K = K, counts = counts.chunk)
    progress <- .np_plot_progress_tick(state = progress, done = stopi)
    start <- stopi + 1L
  }

  list(t = tmat, t0 = t0)
}

.np_con_inid_ksum_eligible <- function(bws) {
  isTRUE(identical(bws$cxkertype, bws$cykertype)) &&
    isTRUE(identical(bws$cxkerorder, bws$cykerorder)) &&
    isTRUE(identical(bws$cxkerbound, bws$cykerbound)) &&
    isTRUE(identical(bws$uxkertype, bws$uykertype)) &&
    isTRUE(identical(bws$oxkertype, bws$oykertype))
}

.np_con_xregtype <- function(bws) {
  regtype <- if (is.null(bws$regtype.engine)) bws$regtype else bws$regtype.engine
  if (is.null(regtype)) "lc" else as.character(regtype)
}

.np_con_make_kbandwidth_x <- function(bws, xdat) {
  xdat <- toFrame(xdat)
  kbandwidth.numeric(
    bw = bws$xbw,
    bwscaling = FALSE,
    # npksum helper constructors require raw bandwidths; bwscaling flags are
    # non-fit-defining here and are intentionally normalized to FALSE.
    bwtype = bws$type,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = if (!is.null(bws$cxkerlb)) bws$cxkerlb else NULL,
    ckerub = if (!is.null(bws$cxkerub)) bws$cxkerub else NULL,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype,
    nobs = nrow(xdat),
    xdati = untangle(xdat),
    xnames = names(xdat)
  )
}

.np_con_make_kbandwidth_xy <- function(bws, xdat, ydat) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  xydat <- data.frame(xdat, ydat)
  ckerlb <- c(if (is.null(bws$cxkerlb)) numeric(0) else bws$cxkerlb,
              if (is.null(bws$cykerlb)) numeric(0) else bws$cykerlb)
  ckerub <- c(if (is.null(bws$cxkerub)) numeric(0) else bws$cxkerub,
              if (is.null(bws$cykerub)) numeric(0) else bws$cykerub)

  kbandwidth.numeric(
    bw = c(bws$xbw, bws$ybw),
    bwscaling = FALSE,
    # npksum helper constructors require raw bandwidths; bwscaling flags are
    # non-fit-defining here and are intentionally normalized to FALSE.
    bwtype = bws$type,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = if (length(ckerlb)) ckerlb else NULL,
    ckerub = if (length(ckerub)) ckerub else NULL,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype,
    nobs = nrow(xydat),
    xdati = untangle(xydat),
    xnames = names(xydat)
  )
}

.np_ksum_conditional_eval_exact <- function(xdat,
                                            ydat,
                                            exdat,
                                            eydat,
                                            kbx,
                                            kbxy,
                                            cdf) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  xop <- rep.int("normal", ncol(xdat))
  yop <- rep.int(if (cdf) "integral" else "normal", ncol(ydat))
  xyop <- c(xop, yop)
  ones <- matrix(1.0, nrow = nrow(xdat), ncol = 1L)
  xydat <- data.frame(xdat, ydat)
  exydat <- data.frame(exdat, eydat)

  den <- as.numeric(npksum(
    txdat = xdat,
    tydat = ones,
    exdat = exdat,
    bws = kbx,
    weights = ones,
    operator = xop,
    bandwidth.divide = TRUE
  )$ksum) / nrow(xdat)
  num <- as.numeric(npksum(
    txdat = xydat,
    tydat = ones,
    exdat = exydat,
    bws = kbxy,
    weights = ones,
    operator = xyop,
    bandwidth.divide = TRUE
  )$ksum) / nrow(xydat)

  num / pmax(den, .Machine$double.eps)
}

.np_ksum_conditional_operator_fixed <- function(xdat,
                                                ydat,
                                                exdat,
                                                eydat,
                                                bws,
                                                cdf) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  n <- nrow(xdat)

  kbx <- .np_con_make_kbandwidth_x(bws = bws, xdat = xdat)
  kbxy <- .np_con_make_kbandwidth_xy(bws = bws, xdat = xdat, ydat = ydat)
  xop <- rep.int("normal", ncol(xdat))
  yop <- rep.int(if (cdf) "integral" else "normal", ncol(ydat))
  xyop <- c(xop, yop)

  Kden <- npksum(
    txdat = xdat,
    exdat = exdat,
    bws = kbx,
    tydat = diag(n),
    operator = xop,
    bandwidth.divide = TRUE
  )$ksum
  Knum <- npksum(
    txdat = data.frame(xdat, ydat),
    exdat = data.frame(exdat, eydat),
    bws = kbxy,
    tydat = diag(n),
    operator = xyop,
    bandwidth.divide = TRUE
  )$ksum

  list(
    den = .np_operator_matrix_from_ksum(
      ksum = Kden,
      nrow.out = nrow(exdat),
      ncol.out = n,
      where = "npksum conditional denominator operator"
    ) / n,
    num = .np_operator_matrix_from_ksum(
      ksum = Knum,
      nrow.out = nrow(exdat),
      ncol.out = n,
      where = "npksum conditional numerator operator"
    ) / n
  )
}

.np_inid_boot_from_conditional_localpoly_fixed <- function(xdat,
                                                           ydat,
                                                           exdat,
                                                           eydat,
                                                           bws,
                                                           B,
                                                           cdf,
                                                           counts = NULL,
                                                           counts.drawer = NULL) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  B <- as.integer(B)

  n <- nrow(xdat)
  neval <- nrow(exdat)
  if (nrow(ydat) != n || nrow(eydat) != neval)
    stop("conditional localpoly bootstrap helper requires aligned x/y training and evaluation rows")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid conditional localpoly bootstrap dimensions")
  if (!identical(bws$type, "fixed"))
    stop("conditional localpoly bootstrap helper supports only fixed bandwidths")

  xbw <- .npcdhat_make_xbw(bws = bws, txdat = xdat)
  ybw <- .npcdhat_make_ybw(bws = bws, tydat = ydat)
  Gy <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = ydat,
    exdat = eydat,
    operator = rep.int(if (cdf) "integral" else "normal", ncol(ydat))
  )
  Gy <- .np_operator_matrix_from_ksum(
    ksum = Gy,
    nrow.out = neval,
    ncol.out = n,
    where = "conditional localpoly y-operator"
  )

  counts.mat <- if (!is.null(counts)) {
    .np_inid_counts_matrix(n = n, B = B, counts = counts)
  } else if (!is.null(counts.drawer)) {
    .np_inid_counts_matrix(n = n, B = B, counts = counts.drawer(1L, B))
  } else {
    .np_inid_counts_matrix(n = n, B = B)
  }

  t0 <- numeric(neval)
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = neval, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  old.progress <- getOption("np.plot.progress", TRUE)
  options(np.plot.progress = FALSE)
  on.exit({
    options(np.plot.progress = old.progress)
  }, add = TRUE)

  # Exact duplicate-sample refits for fixed conditional ll/lp reduce to
  # reusing the fixed x-side local-polynomial bootstrap helper on y-kernel
  # pseudo-responses, one evaluation row at a time.
  for (i in seq_len(neval)) {
    out <- .np_inid_boot_from_regression(
      xdat = xdat,
      exdat = exdat[i, , drop = FALSE],
      bws = xbw,
      ydat = as.numeric(Gy[i, ]),
      B = B,
      counts = counts.mat
    )
    t0[i] <- out$t0[1L]
    tmat[, i] <- out$t[, 1L]
    progress <- .np_plot_progress_tick(state = progress, done = i)
  }

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_ksum_conditional_exact <- function(xdat,
                                                      ydat,
                                                      exdat,
                                                      eydat,
                                                      bws,
                                                      B,
                                                      cdf,
                                                      counts = NULL,
                                                      counts.drawer = NULL) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  B <- as.integer(B)
  n <- nrow(xdat)

  if (nrow(ydat) != n || nrow(exdat) != nrow(eydat))
    stop("conditional exact bootstrap helper requires aligned x/y training and evaluation rows")
  if (n < 1L || nrow(exdat) < 1L || B < 1L)
    stop("invalid conditional exact bootstrap dimensions")
  if (!.np_con_inid_ksum_eligible(bws))
    return(NULL)

  kbx <- tryCatch(.np_con_make_kbandwidth_x(bws = bws, xdat = xdat),
                  error = function(e) NULL)
  kbxy <- tryCatch(.np_con_make_kbandwidth_xy(bws = bws, xdat = xdat, ydat = ydat),
                   error = function(e) NULL)
  if (is.null(kbx) || is.null(kbxy))
    return(NULL)

  fit_one <- function(x.train, y.train) {
    .np_ksum_conditional_eval_exact(
      xdat = x.train,
      ydat = y.train,
      exdat = exdat,
      eydat = eydat,
      kbx = kbx,
      kbxy = kbxy,
      cdf = cdf
    )
  }

  t0 <- fit_one(x.train = xdat, y.train = ydat)
  tmat <- matrix(NA_real_, nrow = B, ncol = length(t0))
  counts.mat <- if (!is.null(counts)) .np_inid_counts_matrix(n = n, B = B, counts = counts) else NULL
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  start <- 1L
  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- if (!is.null(counts.mat)) {
      counts.mat[, start:stopi, drop = FALSE]
    } else if (!is.null(counts.drawer)) {
      .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
    } else {
      stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    }

    for (jj in seq_len(bsz)) {
      idx <- .np_counts_to_indices(counts.chunk[, jj])
      tmat[start + jj - 1L, ] <- fit_one(
        x.train = xdat[idx, , drop = FALSE],
        y.train = ydat[idx, , drop = FALSE]
      )
    }

    progress <- .np_plot_progress_tick(state = progress, done = stopi)
    start <- stopi + 1L
  }

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_ksum_conditional <- function(xdat,
                                                ydat,
                                                exdat,
                                                eydat,
                                                bws,
                                                B,
                                                cdf,
                                                counts = NULL,
                                                counts.drawer = NULL) {
  regtype <- .np_con_xregtype(bws)

  if (!identical(bws$type, "fixed")) {
    return(.np_inid_boot_from_ksum_conditional_exact(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      bws = bws,
      B = B,
      cdf = cdf,
      counts = counts,
      counts.drawer = counts.drawer
    ))
  }

  if (!identical(regtype, "lc")) {
    return(.np_inid_boot_from_conditional_localpoly_fixed(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      bws = bws,
      B = B,
      cdf = cdf,
      counts = counts,
      counts.drawer = counts.drawer
    ))
  }

  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)

  if (nrow(ydat) != n || nrow(eydat) != neval)
    stop("conditional inid helper path requires aligned x/y training and evaluation rows")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid conditional inid bootstrap dimensions")
  if (!.np_con_inid_ksum_eligible(bws))
    return(NULL)

  ops <- tryCatch(
    .np_ksum_conditional_operator_fixed(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      bws = bws,
      cdf = cdf
    ),
    error = function(e) NULL
  )
  if (is.null(ops))
    return(NULL)
  t0 <- rowSums(ops$num) / pmax(rowSums(ops$den), .Machine$double.eps)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    den <- t(ops$den %*% counts.mat)
    num <- t(ops$num %*% counts.mat)
    return(list(t = num / pmax(den, .Machine$double.eps), t0 = t0))
  }

  chunk.size <- .np_inid_chunk_size(n = n, B = B)
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  progress.label <- if (!is.null(counts.drawer)) "Plot bootstrap block" else "Plot bootstrap inid"
  progress <- .np_plot_progress_begin(total = B, label = progress.label)
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    counts.chunk <- if (!is.null(counts.drawer)) {
      .np_inid_counts_matrix(n = n, B = bsz, counts = counts.drawer(start, stopi))
    } else {
      stats::rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    }
    den <- t(ops$den %*% counts.chunk)
    num <- t(ops$num %*% counts.chunk)
    tmat[start:stopi, ] <- num / pmax(den, .Machine$double.eps)
    progress <- .np_plot_progress_tick(state = progress, done = stopi)
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

.np_plot_regression_eval <- function(bws,
                                     xdat,
                                     ydat,
                                     exdat,
                                     gradients = FALSE,
                                     gradient.order = 1L,
                                     need.asymptotic = FALSE) {
  if (isTRUE(need.asymptotic)) {
    return(npreg(
      txdat = xdat,
      tydat = ydat,
      exdat = exdat,
      bws = bws,
      gradients = gradients,
      gradient.order = gradient.order,
      warn.glp.gradient = FALSE
    ))
  }

  fit <- .np_regression_direct(
    bws = bws,
    txdat = xdat,
    tydat = ydat,
    exdat = exdat,
    gradients = gradients,
    gradient.order = gradient.order,
    local.mode = identical(bws$type, "generalized_nn")
  )

  neval <- length(fit$mean)
  fit$merr <- rep(NA_real_, neval)
  if (isTRUE(gradients))
    fit$gerr <- matrix(NA_real_, nrow = neval, ncol = NCOL(fit$grad))

  fit
}

.np_plot_singleindex_asymptotic_eval <- function(bws,
                                                 txdat,
                                                 tydat,
                                                 exdat = NULL,
                                                 gradients = FALSE) {
  no.ex <- is.null(exdat)
  gradients <- npValidateScalarLogical(gradients, "gradients")

  txdat <- toFrame(txdat)
  if (!no.ex) {
    exdat <- toFrame(exdat)
    if (!(txdat %~% exdat))
      stop("'txdat' and 'exdat' are not similar data frames!")
  }

  if (is.factor(tydat)) {
    tydat <- adjustLevels(data.frame(tydat), bws$ydati)[, 1L]
    tydat <- (bws$ydati$all.dlev[[1L]])[as.integer(tydat)]
  } else {
    tydat <- as.double(tydat)
  }

  txdat <- adjustLevels(txdat, bws$xdati)
  if (!no.ex)
    exdat <- adjustLevels(exdat, bws$xdati)

  txmat <- toMatrix(txdat)
  exmat <- if (no.ex) txmat else toMatrix(exdat)
  index <- as.vector(txmat %*% bws$beta)
  index.eval <- as.vector(exmat %*% bws$beta)

  spec <- .npindex_resolve_spec(bws, where = "plot.singleindex")
  regtype <- spec$regtype.engine
  npreg.args <- list(
    txdat = data.frame(index = index),
    tydat = tydat,
    exdat = data.frame(index = index.eval),
    bws = bws$bw,
    bwtype = bws$type,
    ckertype = bws$ckertype,
    ckerorder = bws$ckerorder,
    regtype = regtype,
    gradients = gradients,
    warn.glp.gradient = FALSE
  )
  if (identical(regtype, "lp")) {
    npreg.args$basis <- spec$basis.engine
    npreg.args$degree <- spec$degree.engine
    npreg.args$bernstein.basis <- spec$bernstein.basis.engine
  }

  fit <- do.call(npreg, npreg.args)
  out <- list(
    index = index.eval,
    mean = as.vector(fit$mean),
    merr = as.double(fit$merr)
  )

  if (gradients) {
    grad.index <- as.vector(fit$grad[, 1L])
    gerr.index <- as.vector(fit$gerr[, 1L])
    beta <- as.vector(bws$beta)
    out$grad.index <- grad.index
    out$gerr.index <- gerr.index
    out$grad <- grad.index %o% beta
    out$gerr <- gerr.index %o% abs(beta)
  }

  out
}

.np_plot_plreg_asymptotic_fit <-
  function(bws,
           xdat,
           ydat,
           zdat,
           exdat,
           ezdat) {

    xdat <- toFrame(xdat)
    zdat <- toFrame(zdat)

    keep.rows <- rep_len(TRUE, nrow(xdat))
    rows.omit <- attr(na.omit(data.frame(xdat, ydat, zdat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    xdat <- xdat[keep.rows, , drop = FALSE]
    ydat <- ydat[keep.rows]
    zdat <- zdat[keep.rows, , drop = FALSE]

    exdat <- toFrame(exdat)
    ezdat <- toFrame(ezdat)

    keep.eval <- rep_len(TRUE, nrow(exdat))
    rows.omit <- attr(na.omit(data.frame(exdat, ezdat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.eval[as.integer(rows.omit)] <- FALSE

    if (!any(keep.eval))
      stop("Evaluation data has no rows without NAs")

    exdat <- exdat[keep.eval, , drop = FALSE]
    ezdat <- ezdat[keep.eval, , drop = FALSE]

    fit <- npplreg(
      bws = bws,
      txdat = xdat,
      tydat = ydat,
      tzdat = zdat,
      exdat = exdat,
      ezdat = ezdat
    )

    yfit <- npreg(
      txdat = zdat,
      tydat = ydat,
      exdat = ezdat,
      bws = bws$bw$yzbw
    )

    neval <- nrow(exdat)
    resx.eval <- matrix(0.0, nrow = neval, ncol = ncol(xdat))
    x.err.sq <- numeric(neval)

    for (j in seq_len(ncol(xdat))) {
      xfit <- npreg(
        txdat = zdat,
        tydat = xdat[, j],
        exdat = ezdat,
        bws = bws$bw[[j + 1L]]
      )

      if (is.factor(xdat[1L, j])) {
        tmp.eval <- adjustLevels(exdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati, allowNewCells = TRUE)
        x.num.eval <- (bws$bw[[j + 1L]]$ydati$all.dlev[[1L]])[as.integer(tmp.eval[, 1L])]
      } else {
        x.num.eval <- as.double(exdat[, j])
      }

      resx.eval[, j] <- x.num.eval - xfit$mean
      x.err.sq <- x.err.sq + (fit$xcoef[j]^2) * (as.double(xfit$merr)^2)
    }

    beta.vcov.term <- rowSums((resx.eval %*% fit$xcoefvcov) * resx.eval)
    fit$merr <- sqrt(pmax(as.double(yfit$merr)^2 + x.err.sq + beta.vcov.term, 0.0))
    fit
  }

.np_plot_unconditional_eval <- function(xdat,
                                        exdat,
                                        bws,
                                        cdf = FALSE,
                                        need.asymptotic = FALSE) {
  if (isTRUE(need.asymptotic)) {
    return(if (isTRUE(cdf)) {
      npudist(tdat = xdat, edat = exdat, bws = bws)
    } else {
      npudens(tdat = xdat, edat = exdat, bws = bws)
    })
  }

  est <- .np_ksum_unconditional_eval_exact(
    xdat = xdat,
    exdat = exdat,
    bws = bws,
    operator = if (isTRUE(cdf)) "integral" else "normal"
  )

  if (isTRUE(cdf)) {
    list(dist = est, derr = rep(NA_real_, length(est)))
  } else {
    list(dens = est, derr = rep(NA_real_, length(est)))
  }
}

.np_plot_quantile_eval <- function(bws,
                                   txdat,
                                   tydat,
                                   exdat,
                                   tau = 0.5,
                                   gradients = FALSE,
                                   ftol = 1.490116e-07,
                                   tol = 1.490116e-04,
                                   small = 1.490116e-05,
                                   itmax = 10000,
                                   lbc.dir = 0.5,
                                   dfc.dir = 3,
                                   cfac.dir = 2.5 * (3.0 - sqrt(5)),
                                   initc.dir = 1.0,
                                   lbd.dir = 0.1,
                                   hbd.dir = 1.0,
                                   dfac.dir = 0.25 * (3.0 - sqrt(5)),
                                   initd.dir = 1.0,
                                   ...) {
  fit.start <- proc.time()[3]
  gradients <- npValidateScalarLogical(gradients, "gradients")

  if (!is.numeric(itmax) || length(itmax) != 1L || is.na(itmax) ||
      !is.finite(itmax) || itmax < 1 || itmax != floor(itmax))
    stop("'itmax' must be a positive integer")
  if (!is.numeric(ftol) || length(ftol) != 1L || is.na(ftol) ||
      !is.finite(ftol) || ftol <= 0)
    stop("'ftol' must be a positive finite numeric scalar")
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0)
    stop("'tol' must be a positive finite numeric scalar")
  if (!is.numeric(small) || length(small) != 1L || is.na(small) ||
      !is.finite(small) || small <= 0)
    stop("'small' must be a positive finite numeric scalar")

  itmax <- as.integer(itmax)
  ftol <- as.double(ftol)
  tol <- as.double(tol)
  small <- as.double(small)

  no.ex <- missing(exdat)

  xdat <- toFrame(txdat)
  ydat <- toFrame(tydat)

  if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau <= 0 || tau >= 1)
    stop("'tau' must be a single numeric value in (0,1)")
  if (ncol(ydat) != 1L)
    stop("'tydat' has more than one column")

  if (!no.ex) {
    exdat <- toFrame(exdat)
    if (!xdat %~% exdat)
      stop("'txdat' and 'exdat' are not similar data frames!")
  }

  if (gradients)
    stop("gradients not currently supported for this object")

  if (length(bws$xbw) != length(xdat))
    stop("length of bandwidth vector does not match number of columns of 'txdat'")
  if (length(bws$ybw) != 1L)
    stop("length of bandwidth vector does not match number of columns of 'tydat'")

  if (any(bws$iyord) || any(bws$iyuno) || coarseclass(ydat[, 1L]) != "numeric")
    stop("'tydat' is not continuous")

  if ((any(bws$ixcon) &&
       !all(vapply(xdat[, bws$ixcon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
      (any(bws$ixord) &&
       !all(vapply(xdat[, bws$ixord, drop = FALSE], inherits, logical(1), "ordered"))) ||
      (any(bws$ixuno) &&
       !all(vapply(xdat[, bws$ixuno, drop = FALSE], inherits, logical(1), "factor"))))
    stop("supplied bandwidths do not match 'txdat' in type")

  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE
  if (!any(keep.rows))
    stop("Training data has no rows without NAs")

  xdat <- xdat[keep.rows, , drop = FALSE]
  ydat <- ydat[keep.rows, , drop = FALSE]

  if (!no.ex) {
    keep.eval <- rep_len(TRUE, nrow(exdat))
    rows.omit <- attr(na.omit(exdat), "na.action")
    if (length(rows.omit) > 0L)
      keep.eval[as.integer(rows.omit)] <- FALSE
    exdat <- exdat[keep.eval, , drop = FALSE]
  }

  tnrow <- nrow(xdat)
  enrow <- if (no.ex) tnrow else nrow(exdat)

  xdat <- adjustLevels(xdat, bws$xdati)
  ydat <- adjustLevels(ydat, bws$ydati)
  if (!no.ex)
    exdat <- adjustLevels(exdat, bws$xdati)

  txeval <- if (no.ex) xdat else exdat
  xdat.df <- xdat
  ydat.df <- ydat
  if (!no.ex)
    exdat.df <- exdat

  ydat <- toMatrix(ydat)
  xdat <- toMatrix(xdat)

  xuno <- xdat[, bws$ixuno, drop = FALSE]
  xcon <- xdat[, bws$ixcon, drop = FALSE]
  xord <- xdat[, bws$ixord, drop = FALSE]

  if (!no.ex) {
    exdat <- toMatrix(exdat)
    exuno <- exdat[, bws$ixuno, drop = FALSE]
    excon <- exdat[, bws$ixcon, drop = FALSE]
    exord <- exdat[, bws$ixord, drop = FALSE]
  } else {
    exuno <- data.frame()
    excon <- data.frame()
    exord <- data.frame()
  }

  myopti <- list(
    num_obs_train = tnrow,
    num_obs_eval = enrow,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_den_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN
    ),
    int_MINIMIZE_IO = if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
    xkerneval = switch(bws$cxkertype,
      gaussian = CKER_GAUSS + bws$cxkerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$cxkerorder / 2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS
    ),
    ykerneval = switch(bws$cykertype,
      gaussian = CKER_GAUSS + bws$cykerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$cykerorder / 2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS
    ),
    uxkerneval = switch(bws$uxkertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR
    ),
    uykerneval = switch(bws$uykertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR
    ),
    oxkerneval = switch(bws$oxkertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR,
      racineliyan = OKER_RLY
    ),
    oykerneval = switch(bws$oykertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY
    ),
    num_yuno = bws$ynuno,
    num_yord = bws$ynord,
    num_ycon = bws$yncon,
    num_xuno = bws$xnuno,
    num_xord = bws$xnord,
    num_xcon = bws$xncon,
    no.ex = no.ex,
    gradients = gradients,
    itmax = itmax,
    xmcv.numRow = attr(bws$xmcv, "num.row"),
    nmulti = itmax,
    dfc.dir = dfc.dir
  )

  myoptd <- list(
    ftol = ftol,
    tol = tol,
    small = small,
    lbc.dir = lbc.dir,
    cfac.dir = cfac.dir,
    initc.dir = initc.dir,
    lbd.dir = lbd.dir,
    hbd.dir = hbd.dir,
    dfac.dir = dfac.dir,
    initd.dir = initd.dir
  )

  myout <- .Call(
    "C_np_quantile_conditional",
    as.double(ydat),
    as.double(xuno), as.double(xord), as.double(xcon),
    as.double(exuno), as.double(exord), as.double(excon),
    as.double(tau),
    as.double(c(
      bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
      bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
      bws$xbw[bws$ixuno], bws$xbw[bws$ixord]
    )),
    as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
    as.double(bws$nconfac), as.double(bws$ncatfac), as.double(bws$sdev),
    as.integer(myopti),
    as.double(myoptd),
    as.integer(enrow),
    as.integer(bws$xndim),
    as.logical(gradients),
    PACKAGE = "np"
  )[c("yq", "yqerr", "yqgrad")]

  if (all(!is.finite(myout$yqerr) | myout$yqerr <= 0.0)) {
    dens.obj <- tryCatch(
      .np_plot_conditional_eval(
        bws = bws,
        xdat = xdat.df,
        ydat = ydat.df,
        exdat = if (no.ex) xdat.df else exdat.df,
        eydat = setNames(data.frame(as.double(myout$yq)), names(ydat.df)),
        cdf = FALSE,
        gradients = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(dens.obj)) {
      dens.q <- as.double(dens.obj$condens)
      qvar <- tau * (1.0 - tau) / (tnrow * NZD(dens.q)^2)
      myout$yqerr <- sqrt(pmax(qvar, 0.0))
      myout$yqerr[!is.finite(myout$yqerr)] <- NA_real_
    }
  }

  if (gradients) {
    myout$yqgrad <- matrix(data = myout$yqgrad, nrow = enrow, ncol = bws$xndim, byrow = FALSE)
    rorder <- numeric(bws$xndim)
    xidx <- seq_len(bws$xndim)
    rorder[c(xidx[bws$ixcon], xidx[bws$ixuno], xidx[bws$ixord])] <- xidx
    myout$yqgrad <- myout$yqgrad[, rorder, drop = FALSE]
  } else {
    myout$yqgrad <- NA
  }

  fit.elapsed <- proc.time()[3] - fit.start
  optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
  total.time <- fit.elapsed + if (is.na(optim.time)) 0.0 else optim.time

  qregression(
    bws = bws,
    xeval = txeval,
    tau = tau,
    quantile = myout$yq,
    quanterr = myout$yqerr,
    quantgrad = myout$yqgrad,
    ntrain = tnrow,
    trainiseval = no.ex,
    gradients = gradients,
    timing = bws$timing,
    total.time = total.time,
    optim.time = optim.time,
    fit.time = fit.elapsed
  )
}

.np_plot_conditional_eval <- function(bws,
                                      xdat,
                                      ydat,
                                      exdat,
                                      eydat,
                                      cdf = FALSE,
                                      gradients = FALSE,
                                      proper = FALSE,
                                      proper.method = NULL,
                                      proper.control = list()) {
  fit.start <- proc.time()[3]
  proper <- npValidateScalarLogical(proper, "proper")

  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)

  if (nrow(xdat) != nrow(ydat))
    stop("conditional plot helper requires aligned training rows")
  if (nrow(exdat) != nrow(eydat))
    stop("conditional plot helper requires aligned evaluation rows")
  if (!xdat %~% exdat)
    stop("'xdat' and 'exdat' are not similar data frames!")
  if (!ydat %~% eydat)
    stop("'ydat' and 'eydat' are not similar data frames!")

  xdat <- adjustLevels(xdat, bws$xdati)
  ydat <- adjustLevels(ydat, bws$ydati)
  exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
  eydat <- adjustLevels(eydat, bws$ydati, allowNewCells = TRUE)
  npKernelBoundsCheckEval(exdat, bws$ixcon, bws$cxkerlb, bws$cxkerub, argprefix = "cxker")
  npKernelBoundsCheckEval(eydat, bws$iycon, bws$cykerlb, bws$cykerub, argprefix = "cyker")

  txeval <- exdat
  tyeval <- eydat
  tnrow <- nrow(xdat)
  enrow <- nrow(exdat)

  ydat <- toMatrix(ydat)
  tyuno <- ydat[, bws$iyuno, drop = FALSE]
  tycon <- ydat[, bws$iycon, drop = FALSE]
  tyord <- ydat[, bws$iyord, drop = FALSE]

  xdat <- toMatrix(xdat)
  txuno <- xdat[, bws$ixuno, drop = FALSE]
  txcon <- xdat[, bws$ixcon, drop = FALSE]
  txord <- xdat[, bws$ixord, drop = FALSE]

  eydat <- toMatrix(eydat)
  eyuno <- eydat[, bws$iyuno, drop = FALSE]
  eycon <- eydat[, bws$iycon, drop = FALSE]
  eyord <- eydat[, bws$iyord, drop = FALSE]

  exdat <- toMatrix(exdat)
  exuno <- exdat[, bws$ixuno, drop = FALSE]
  excon <- exdat[, bws$ixcon, drop = FALSE]
  exord <- exdat[, bws$ixord, drop = FALSE]

  reg.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  basis.engine <- if (is.null(bws$basis.engine)) {
    if (is.null(bws$basis)) "glp" else bws$basis
  } else {
    bws$basis.engine
  }
  degree.engine <- if (is.null(bws$degree.engine)) {
    if (bws$xncon > 0L) {
      if (identical(reg.engine, "lc")) rep.int(0L, bws$xncon) else npValidateGlpDegree(
        regtype = "lp",
        degree = bws$degree,
        ncon = bws$xncon
      )
    } else {
      integer(0)
    }
  } else {
    as.integer(bws$degree.engine)
  }
  bernstein.engine <- if (is.null(bws$bernstein.basis.engine)) {
    isTRUE(bws$bernstein.basis)
  } else {
    isTRUE(bws$bernstein.basis.engine)
  }

  reg.c <- npRegtypeToC(
    regtype = if (identical(reg.engine, "lp")) "lp" else "lc",
    degree = degree.engine,
    ncon = bws$xncon,
    context = if (isTRUE(cdf)) "npcdist" else "npcdens"
  )
  degree.c <- if (bws$xncon > 0L) {
    as.integer(if (is.null(reg.c$degree)) rep.int(0L, bws$xncon) else reg.c$degree)
  } else {
    integer(0)
  }
  basis.code <- as.integer(npLpBasisCode(basis.engine))

  myopti <- list(
    num_obs_train = tnrow,
    num_obs_eval = enrow,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_den_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN
    ),
    int_MINIMIZE_IO = if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
    xkerneval = switch(bws$cxkertype,
      gaussian = CKER_GAUSS + bws$cxkerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$cxkerorder / 2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS
    ),
    ykerneval = switch(bws$cykertype,
      gaussian = CKER_GAUSS + bws$cykerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$cykerorder / 2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS
    ),
    uxkerneval = switch(bws$uxkertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR
    ),
    uykerneval = switch(bws$uykertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR
    ),
    oxkerneval = switch(bws$oxkertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY
    ),
    oykerneval = switch(bws$oykertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY
    ),
    num_yuno = bws$ynuno,
    num_yord = bws$ynord,
    num_ycon = bws$yncon,
    num_xuno = bws$xnuno,
    num_xord = bws$xnord,
    num_xcon = bws$xncon,
    no.exy = FALSE,
    gradients = gradients,
    ymcv.numRow = attr(bws$ymcv, "num.row"),
    xmcv.numRow = attr(bws$xmcv, "num.row"),
    densOrDist = if (isTRUE(cdf)) NP_DO_DIST else NP_DO_DENS,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO
  )

  cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
  cyker.bounds.c <- npKernelBoundsMarshal(bws$cykerlb[bws$iycon], bws$cykerub[bws$iycon])

  myout <- .Call(
    "C_np_density_conditional",
    as.double(tyuno), as.double(tyord), as.double(tycon),
    as.double(txuno), as.double(txord), as.double(txcon),
    as.double(eyuno), as.double(eyord), as.double(eycon),
    as.double(exuno), as.double(exord), as.double(excon),
    as.double(c(
      bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
      bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
      bws$xbw[bws$ixuno], bws$xbw[bws$ixord]
    )),
    as.double(bws$ymcv), as.double(attr(bws$ymcv, "pad.num")),
    as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
    as.double(bws$nconfac), as.double(bws$ncatfac), as.double(bws$sdev),
    as.integer(myopti),
    as.integer(enrow),
    as.integer(bws$xndim),
    as.double(cxker.bounds.c$lb),
    as.double(cxker.bounds.c$ub),
    as.double(cyker.bounds.c$lb),
    as.double(cyker.bounds.c$ub),
    as.integer(reg.c$code),
    as.integer(degree.c),
    as.integer(bernstein.engine),
    basis.code,
    PACKAGE = "np"
  )

  if (isTRUE(cdf))
    names(myout)[1L] <- "condist"

  if (isTRUE(gradients)) {
    xidx <- seq_len(bws$xndim)
    rorder <- numeric(bws$xndim)
    rorder[c(xidx[bws$ixcon], xidx[bws$ixuno], xidx[bws$ixord])] <- xidx
    myout$congrad <- matrix(data = myout$congrad, nrow = enrow, ncol = bws$xndim, byrow = FALSE)
    myout$congrad <- myout$congrad[, rorder, drop = FALSE]
    myout$congerr <- matrix(data = myout$congerr, nrow = enrow, ncol = bws$xndim, byrow = FALSE)
    myout$congerr <- myout$congerr[, rorder, drop = FALSE]
  } else {
    myout$congrad <- NA
    myout$congerr <- NA
  }

  fit.elapsed <- proc.time()[3] - fit.start
  optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
  total.time <- fit.elapsed + if (is.na(optim.time)) 0.0 else optim.time

  if (isTRUE(cdf)) {
    out <- condistribution(
      bws = bws,
      xeval = txeval,
      yeval = tyeval,
      condist = myout$condist,
      conderr = myout$conderr,
      congrad = myout$congrad,
      congerr = myout$congerr,
      ntrain = tnrow,
      trainiseval = FALSE,
      gradients = gradients,
      rows.omit = integer(0),
      timing = bws$timing,
      total.time = total.time,
      optim.time = optim.time,
      fit.time = fit.elapsed
    )

    if (isTRUE(proper)) {
      .np_condist_finalize_proper_object(
        object = out,
        proper = TRUE,
        proper.method = proper.method,
        proper.control = proper.control,
        where = "plot()"
      )
    } else {
      out
    }
  } else {
    out <- condensity(
      bws = bws,
      xeval = txeval,
      yeval = tyeval,
      condens = myout$condens,
      conderr = myout$conderr,
      congrad = myout$congrad,
      congerr = myout$congerr,
      ll = myout$log_likelihood,
      ntrain = tnrow,
      trainiseval = FALSE,
      gradients = gradients,
      rows.omit = integer(0),
      timing = bws$timing,
      total.time = total.time,
      optim.time = optim.time,
      fit.time = fit.elapsed
    )

    if (isTRUE(proper)) {
      .np_condens_finalize_proper_object(
        object = out,
        proper = TRUE,
        proper.method = proper.method,
        proper.control = proper.control,
        where = "plot()"
      )
    } else {
      out
    }
  }
}

## Rank-based simultaneous confidence set helper, vendored from
## MCPAN::SCSrank (MCPAN 1.1-21, GPL-2; Schaarschmidt, Gerhard, Sill).
np.plot.SCSrank <- function(x, conf.level = 0.95, alternative = "two.sided", ...) {
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))

  DataMatrix <- x
  N <- nrow(DataMatrix)
  K <- ncol(DataMatrix)
  k <- round(conf.level * N, 0)
  row.max <- rep.int(-Inf, N)
  row.min <- rep.int(Inf, N)

  for (j in seq_len(K)) {
    ranks.j <- rank(DataMatrix[, j])
    row.max <- pmax(row.max, ranks.j)
    row.min <- pmin(row.min, ranks.j)
  }

  SCS <- matrix(NA_real_, nrow = K, ncol = 2L)

  switch(alternative,
    "two.sided" = {
      tstar <- round(sort.int(pmax(row.max, N + 1 - row.min))[k], 0)
      lower.idx <- N + 1 - tstar
      upper.idx <- tstar
      for (j in seq_len(K)) {
        sortx <- sort.int(DataMatrix[, j])
        SCS[j, ] <- c(sortx[lower.idx], sortx[upper.idx])
      }
    },
    "less" = {
      tstar <- round(sort.int(row.max)[k], 0)
      for (j in seq_len(K)) {
        sortx <- sort.int(DataMatrix[, j])
        SCS[j, ] <- c(-Inf, sortx[tstar])
      }
    },
    "greater" = {
      tstar <- round(sort.int(N + 1 - row.min)[k], 0)
      lower.idx <- N + 1 - tstar
      for (j in seq_len(K)) {
        sortx <- sort.int(DataMatrix[, j])
        SCS[j, ] <- c(sortx[lower.idx], Inf)
      }
    }
  )

  colnames(SCS) <- c("lower", "upper")

  attr(SCS, which = "k") <- k
  attr(SCS, which = "N") <- N
  OUT <- list(conf.int = SCS, conf.level = conf.level, alternative = alternative)
  return(OUT)
}

.np_plot_quantile_type7_sorted <- function(sorted.x, probs) {
  n <- length(sorted.x)
  probs <- pmin(pmax(as.double(probs), 0), 1)
  if (n < 1L)
    return(rep.int(NA_real_, length(probs)))

  h <- 1 + (n - 1) * probs
  lo <- floor(h)
  hi <- ceiling(h)
  q <- sorted.x[lo] + (h - lo) * (sorted.x[hi] - sorted.x[lo])
  q[probs <= 0] <- sorted.x[1L]
  q[probs >= 1] <- sorted.x[n]
  q
}

.np_plot_quantile_bounds_multi <- function(boot.t, probs.list) {
  probs.list <- unclass(probs.list)
  if (!length(probs.list))
    return(list())

  neval <- ncol(boot.t)
  row.names <- colnames(boot.t)
  out <- lapply(probs.list, function(probs) {
    prob.names <- names(stats::quantile(c(0, 1), probs = probs))
    matrix(NA_real_, nrow = neval, ncol = length(probs),
           dimnames = list(row.names, prob.names))
  })
  names(out) <- names(probs.list)

  for (j in seq_len(neval)) {
    sorted.j <- sort.int(boot.t[, j])
    for (nm in names(probs.list)) {
      out[[nm]][j, ] <- .np_plot_quantile_type7_sorted(sorted.j, probs.list[[nm]])
    }
  }

  out
}

.np_plot_quantile_bounds_single <- function(boot.t, probs) {
  t(apply(boot.t, 2L, quantile, probs = probs))
}

compute.bootstrap.quantile.bounds <- function(boot.t, alpha, band.type, warn.coverage = TRUE) {
  B <- nrow(boot.t)
  neval <- ncol(boot.t)

  .np_plot_bootstrap_tail_warning <- function(B, alpha, band.type, neval) {
    m <- if (identical(band.type, "pointwise")) 1L else max(1L, as.integer(neval))
    min.B <- ceiling((2.0 * m) / alpha - 1.0)
    if (B >= min.B)
      return(invisible(NULL))

    m.desc <- if (m == 1L) {
      "m=1 (pointwise tails)"
    } else {
      sprintf("m=n.eval=%d (Bonferroni-conservative tails)", neval)
    }
    .np_warning(sprintf(
      paste0("plot.errors.boot.num=%d is too small for plot.errors.type='%s' ",
             "(alpha=%g). Minimum recommended is %d using ",
             "B >= ceiling(2*m/alpha - 1), with %s. ",
             "For 2D perspective plots on a full neval x neval grid, m=neval^2."),
      B, band.type, alpha, min.B, m.desc
    ), call. = FALSE)
  }

  if (isTRUE(warn.coverage) && band.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
    warn.type <- if (identical(band.type, "all")) "bonferroni/simultaneous/all" else band.type
    .np_plot_bootstrap_tail_warning(B = B, alpha = alpha, band.type = warn.type, neval = neval)
  }

  if (band.type == "pointwise") {
    return(.np_plot_quantile_bounds_single(
      boot.t = boot.t,
      probs = c(alpha / 2.0, 1.0 - alpha / 2.0)
    ))
  }

  if (band.type == "bonferroni") {
    return(.np_plot_quantile_bounds_single(
      boot.t = boot.t,
      probs = c(alpha / (2.0 * neval), 1.0 - alpha / (2.0 * neval))
    ))
  }

  if (band.type == "simultaneous") {
    return(np.plot.SCSrank(boot.t, conf.level = 1.0 - alpha)$conf.int)
  }

  if (band.type == "all") {
    quantile.bounds <- .np_plot_quantile_bounds_multi(
      boot.t = boot.t,
      probs.list = list(
        pointwise = c(alpha / 2.0, 1.0 - alpha / 2.0),
        bonferroni = c(alpha / (2.0 * neval), 1.0 - alpha / (2.0 * neval))
      )
    )
    return(list(
      pointwise = quantile.bounds$pointwise,
      bonferroni = quantile.bounds$bonferroni,
      simultaneous = compute.bootstrap.quantile.bounds(boot.t, alpha, "simultaneous", warn.coverage = FALSE)
    ))
  }

  stop("'band.type' must be one of pointwise, bonferroni, simultaneous, all")
}

.np_plot_bootstrap_interval_summary <- function(boot.t, t0, alpha, band.type) {
  if (!(band.type %in% c("pointwise", "bonferroni", "simultaneous", "all")))
    stop("'band.type' must be one of pointwise, bonferroni, simultaneous, all")

  activity <- .np_plot_activity_begin("Constructing bootstrap bands...")
  on.exit(.np_plot_activity_end(activity), add = TRUE)

  if (identical(band.type, "all")) {
    all.bounds <- compute.bootstrap.quantile.bounds(
      boot.t = boot.t,
      alpha = alpha,
      band.type = "all"
    )
    bounds <- all.bounds$pointwise
    all.err <- lapply(all.bounds, function(bb) {
      cbind(t0 - bb[, 1L], bb[, 2L] - t0)
    })
  } else {
    bounds <- compute.bootstrap.quantile.bounds(
      boot.t = boot.t,
      alpha = alpha,
      band.type = band.type
    )
    all.err <- NULL
  }

  list(
    err = cbind(t0 - bounds[, 1L], bounds[, 2L] - t0),
    all.err = all.err
  )
}

.np_plot_asymptotic_error_from_se <- function(se, alpha, band.type, m = length(se)) {
  se <- as.numeric(se)
  n <- length(se)
  m <- max(1L, as.integer(m[1L]))

  if (!length(se))
    return(list(err = matrix(numeric(0), nrow = 0L, ncol = 2L), all.err = NULL))

  make_err <- function(mult) {
    cbind(mult * se, mult * se)
  }

  if (band.type == "all") {
    err.pointwise <- make_err(qnorm(alpha / 2.0, lower.tail = FALSE))
    err.bonf <- make_err(qnorm(alpha / (2.0 * m), lower.tail = FALSE))
    err.sim <- matrix(NA_real_, nrow = n, ncol = 2L)
    return(list(
      err = err.pointwise,
      all.err = list(
        pointwise = err.pointwise,
        simultaneous = err.sim,
        bonferroni = err.bonf
      )
    ))
  }

  if (band.type == "simultaneous")
    return(list(err = matrix(NA_real_, nrow = n, ncol = 2L), all.err = NULL))

  mult <- if (band.type == "bonferroni") {
    qnorm(alpha / (2.0 * m), lower.tail = FALSE)
  } else if (band.type %in% c("pmzsd", "pointwise")) {
    qnorm(alpha / 2.0, lower.tail = FALSE)
  } else {
    stop("unsupported asymptotic interval type")
  }

  list(err = make_err(mult), all.err = NULL)
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
  rng <- c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
  if (all(is.finite(rng)))
    return(rng)

  center <- center[is.finite(center)]
  if (!length(center))
    return(c(NA_real_, NA_real_))
  range(center, finite = TRUE)
}

compute.default.error.range <- function(center, err) {
  lower <- c(center - err[,1], err[,3] - err[,1])
  upper <- c(center + err[,2], err[,3] + err[,2])
  rng <- c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
  if (all(is.finite(rng)))
    return(rng)

  center <- center[is.finite(center)]
  if (!length(center))
    return(c(NA_real_, NA_real_))
  range(center, finite = TRUE)
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
    .np_warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
    plot.errors.method <- "bootstrap"
  }

  if (allow_asymptotic_quantile && plot.errors.method == "asymptotic") {
    if (plot.errors.type == "simultaneous")
      .np_warning("asymptotic simultaneous confidence bands are unavailable here; returning NA interval limits")
    if (plot.errors.type == "all")
      .np_warning("asymptotic simultaneous confidence bands are unavailable here; 'all' returns pointwise/bonferroni and NA simultaneous")

    if (plot.errors.center == "bias-corrected") {
      .np_warning("no bias corrections can be calculated with asymptotics, centering on estimate")
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
    boot.out <- NULL

    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL
    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))

    if (is.wild.hat && gradients) {
      cont.idx <- which(bws$xdati$icon)
      if (is.na(match(slice.index, cont.idx))) {
        stop("plot.errors.boot.method='wild' supports gradients only for continuous slices in compute.bootstrap.errors.rbandwidth", call. = FALSE)
      }
    }

    inid.helper.ok <- TRUE
    block.helper.ok <- TRUE

    if ((is.inid && isTRUE(inid.helper.ok)) || (is.block && isTRUE(block.helper.ok))) {
      counts.drawer <- if (is.block) {
        .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      } else {
        NULL
      }
      boot.out <- tryCatch(
        if (identical(bws$type, "fixed")) {
          .np_inid_boot_from_regression(
            xdat = xdat,
            exdat = exdat,
            bws = bws,
            ydat = ydat,
            B = plot.errors.boot.num,
            counts.drawer = counts.drawer,
            gradients = gradients,
            gradient.order = gradient.order,
            slice.index = slice.index
          )
        } else {
          .np_inid_boot_from_reghat_exact(
            xdat = xdat,
            exdat = exdat,
            bws = bws,
            ydat = ydat,
            B = plot.errors.boot.num,
            counts.drawer = counts.drawer,
            gradients = gradients,
            gradient.order = gradient.order,
            slice.index = slice.index
          )
        },
        error = function(e) {
          stop(sprintf("%s regression helper failed in compute.bootstrap.errors.rbandwidth (%s)",
                       if (is.block) plot.errors.boot.method else "inid",
                       conditionMessage(e)),
               call. = FALSE)
        }
      )
    }

    if (is.null(boot.out) && is.wild.hat) {
      plot.errors.boot.wild <- .np_plot_normalize_wild(plot.errors.boot.wild)

      fit.mean <- as.vector(suppressWarnings(npreghat(
        bws = bws,
        txdat = xdat,
        exdat = xdat,
        y = ydat,
        output = "apply"
      )))

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

      H <- suppressWarnings(npreghat(
        bws = bws,
        txdat = xdat,
        exdat = exdat,
        s = s.vec,
        output = "matrix"
      ))

      boot.out <- .np_plot_boot_from_hat_wild(
        H = H,
        ydat = ydat,
        fit.mean = fit.mean,
        B = plot.errors.boot.num,
        wild = plot.errors.boot.wild
      )
    }

    if (is.null(boot.out))
      stop(sprintf("unresolved bootstrap execution path for method '%s' in compute.bootstrap.errors.rbandwidth", plot.errors.boot.method), call. = FALSE)

    all.bp <- .np_plot_boot_factor_boxplots(
      boot.t = boot.out$t,
      tdati = bws$xdati,
      ti = slice.index,
      B = plot.errors.boot.num
    )
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      interval.summary <- .np_plot_bootstrap_interval_summary(
        boot.t = boot.out$t,
        t0 = boot.out$t0,
        alpha = plot.errors.alpha,
        band.type = plot.errors.type
      )
      boot.err[,1:2] <- interval.summary$err
      boot.all.err <- interval.summary$all.err
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
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
    miss.z <- missing(zdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
    regtype <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
    boot.out <- NULL

    if (is.inid) {
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
    if (is.null(boot.out) && is.block) {
      if (isTRUE(gradients))
        stop("fixed/geom bootstrap for smooth coefficient gradients is not supported in helper mode", call. = FALSE)
      counts.drawer <- .np_block_counts_drawer(
        n = nrow(xdat),
        B = plot.errors.boot.num,
        blocklen = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method
      )
      boot.out <- tryCatch(
        .np_inid_boot_from_scoef(
          txdat = xdat,
          ydat = ydat,
          tzdat = if (miss.z) NULL else zdat,
          exdat = exdat,
          ezdat = if (miss.z) NULL else ezdat,
          bws = bws,
          B = plot.errors.boot.num,
          counts.drawer = counts.drawer,
          leave.one.out = FALSE
        ),
        error = function(e) {
          stop(sprintf("%s smooth coefficient helper failed in compute.bootstrap.errors.scbandwidth (%s)",
                       plot.errors.boot.method,
                       conditionMessage(e)),
               call. = FALSE)
        }
      )
    }

    if (is.null(boot.out) && is.wild.hat) {
      plot.errors.boot.wild <- .np_plot_normalize_wild(plot.errors.boot.wild)

      hat.eval.args <- list(
        bws = bws,
        txdat = xdat,
        exdat = exdat,
        output = "matrix",
        iterate = FALSE
      )
      hat.train.args <- list(
        bws = bws,
        txdat = xdat,
        exdat = xdat,
        y = ydat,
        output = "apply",
        iterate = FALSE
      )
      if (!miss.z) {
        hat.eval.args$tzdat <- zdat
        hat.eval.args$ezdat <- ezdat
        hat.train.args$tzdat <- zdat
        hat.train.args$ezdat <- zdat
      }

      fit.mean <- as.vector(do.call(npscoefhat, hat.train.args))
      H <- do.call(npscoefhat, hat.eval.args)

      boot.out <- .np_plot_boot_from_hat_wild(
        H = H,
        ydat = ydat,
        fit.mean = fit.mean,
        B = plot.errors.boot.num,
        wild = plot.errors.boot.wild
      )
    }

    if (is.null(boot.out))
      stop(sprintf("unresolved bootstrap execution path for method '%s' in compute.bootstrap.errors.scbandwidth", plot.errors.boot.method), call. = FALSE)

    tdati <- if (slice.index <= ncol(xdat)) bws$xdati else bws$zdati
    ti <- if (slice.index <= ncol(xdat)) slice.index else slice.index - ncol(xdat)
    all.bp <- .np_plot_boot_factor_boxplots(
      boot.t = boot.out$t,
      tdati = tdati,
      ti = ti,
      B = plot.errors.boot.num
    )
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      interval.summary <- .np_plot_bootstrap_interval_summary(
        boot.t = boot.out$t,
        t0 = boot.out$t0,
        alpha = plot.errors.alpha,
        band.type = plot.errors.type
      )
      boot.err[,1:2] <- interval.summary$err
      boot.all.err <- interval.summary$all.err
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
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
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))

    if (is.wild.hat) {
      plot.errors.boot.wild <- .np_plot_normalize_wild(plot.errors.boot.wild)

      fit.mean <- as.vector(npplreghat(
        bws = bws,
        txdat = xdat,
        tzdat = zdat,
        exdat = xdat,
        ezdat = zdat,
        y = ydat,
        output = "apply"
      ))
      H <- npplreghat(
        bws = bws,
        txdat = xdat,
        tzdat = zdat,
        exdat = exdat,
        ezdat = ezdat,
        output = "matrix"
      )

      boot.out <- .np_plot_boot_from_hat_wild(
        H = H,
        ydat = ydat,
        fit.mean = fit.mean,
        B = plot.errors.boot.num,
        wild = plot.errors.boot.wild
      )
    } else {
      boot.out <- NULL
      if (is.inid) {
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
      } else if (is.block) {
        counts.drawer <- .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
        boot.out <- tryCatch(
          .np_inid_boot_from_plreg(
            txdat = xdat,
            ydat = ydat,
            tzdat = zdat,
            exdat = exdat,
            ezdat = ezdat,
            bws = bws,
            B = plot.errors.boot.num,
            counts.drawer = counts.drawer
          ),
          error = function(e) {
            stop(sprintf("%s plreg helper failed in compute.bootstrap.errors.plbandwidth (%s)",
                         plot.errors.boot.method,
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      }

      if (is.null(boot.out))
        stop(sprintf("unresolved bootstrap execution path for method '%s' in compute.bootstrap.errors.plbandwidth", plot.errors.boot.method), call. = FALSE)
    }

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$zdati
      ti <- slice.index - bws$xndim
    }
    all.bp <- .np_plot_boot_factor_boxplots(
      boot.t = boot.out$t,
      tdati = tdati,
      ti = ti,
      B = plot.errors.boot.num
    )

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      interval.summary <- .np_plot_bootstrap_interval_summary(
        boot.t = boot.out$t,
        t0 = boot.out$t0,
        alpha = plot.errors.alpha,
        band.type = plot.errors.type
      )
      boot.err[,1:2] <- interval.summary$err
      boot.all.err <- interval.summary$all.err
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
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
    .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "unconditional density/distribution estimators")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
    fast.inid <- isTRUE(is.inid)
    fast.block <- isTRUE(is.block)

    if (is.inid && !isTRUE(fast.inid))
      stop("inid unconditional helper unavailable for this configuration; no alternate fallback is permitted", call. = FALSE)
    if (is.block && !isTRUE(fast.block))
      stop(sprintf("%s unconditional helper unavailable for this configuration; no alternate fallback is permitted", plot.errors.boot.method), call. = FALSE)

    boot.out <- NULL
    if (fast.inid || fast.block) {
      op <- if (cdf) "integral" else "normal"
      counts.drawer <- if (fast.block) {
        .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      } else {
        NULL
      }
      boot.out <- tryCatch(
        .np_inid_boot_from_ksum_unconditional(
          xdat = xdat,
          exdat = exdat,
          bws = bws,
          B = plot.errors.boot.num,
          operator = op,
          counts.drawer = counts.drawer
        ),
        error = function(e) {
          stop(sprintf("%s ksum helper failed in compute.bootstrap.errors.bandwidth (%s)",
                       if (fast.block) plot.errors.boot.method else "inid",
                       conditionMessage(e)),
               call. = FALSE)
        }
      )
    }

    if (is.null(boot.out))
      stop("no canonical helper path available for this unconditional bootstrap configuration", call. = FALSE)

    all.bp <- .np_plot_boot_factor_boxplots(
      boot.t = boot.out$t,
      tdati = bws$xdati,
      ti = slice.index,
      B = plot.errors.boot.num
    )

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      interval.summary <- .np_plot_bootstrap_interval_summary(
        boot.t = boot.out$t,
        t0 = boot.out$t0,
        alpha = plot.errors.alpha,
        band.type = plot.errors.type
      )
      boot.err[,1:2] <- interval.summary$err
      boot.all.err <- interval.summary$all.err
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
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
    compute.bootstrap.errors.bandwidth(
      xdat = xdat,
      exdat = exdat,
      cdf = TRUE,
      slice.index = slice.index,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      ...,
      bws = bws
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
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
    fast.inid <- isTRUE(is.inid) &&
      isTRUE(!quantreg) &&
      isTRUE(!gradients) &&
      isTRUE(.np_con_inid_ksum_eligible(bws))
    fast.block <- isTRUE(is.block) &&
      isTRUE(!quantreg) &&
      isTRUE(!gradients) &&
      isTRUE(.np_con_inid_ksum_eligible(bws))

    if (is.inid && !isTRUE(fast.inid))
      stop("inid conditional helper unavailable for this configuration; no alternate fallback is permitted", call. = FALSE)
    if (is.block && !isTRUE(fast.block))
      stop(sprintf("%s conditional helper unavailable for this configuration; no alternate fallback is permitted", plot.errors.boot.method), call. = FALSE)

    boot.out <- NULL
    if (fast.inid || fast.block) {
      counts.drawer <- if (fast.block) {
        .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      } else {
        NULL
      }
      boot.out <- tryCatch(
        .np_inid_boot_from_ksum_conditional(
          xdat = xdat,
          ydat = ydat,
          exdat = exdat,
          eydat = eydat,
          bws = bws,
          B = plot.errors.boot.num,
          cdf = cdf,
          counts.drawer = counts.drawer
        ),
        error = function(e) {
          stop(sprintf("%s ksum helper failed in compute.bootstrap.errors.conbandwidth (%s)",
                       if (fast.block) plot.errors.boot.method else "inid",
                       conditionMessage(e)),
               call. = FALSE)
        }
      )
    }

    if (is.null(boot.out))
      stop("no canonical helper path available for this conditional bootstrap configuration", call. = FALSE)

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$ydati
      ti <- slice.index - bws$xndim
    }
    all.bp <- .np_plot_boot_factor_boxplots(
      boot.t = boot.out$t,
      tdati = tdati,
      ti = ti,
      B = plot.errors.boot.num
    )

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      interval.summary <- .np_plot_bootstrap_interval_summary(
        boot.t = boot.out$t,
        t0 = boot.out$t0,
        alpha = plot.errors.alpha,
        band.type = plot.errors.type
      )
      boot.err[,1:2] <- interval.summary$err
      boot.all.err <- interval.summary$all.err
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
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
    compute.bootstrap.errors.conbandwidth(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      cdf = cdf,
      quantreg = quantreg,
      tau = tau,
      gradients = gradients,
      gradient.index = gradient.index,
      slice.index = slice.index,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      ...,
      bws = bws
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

    boot.err = matrix(data = NA, nrow = nrow(xdat), ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method=="inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))

    if (is.wild.hat) {
      plot.errors.boot.wild <- .np_plot_normalize_wild(plot.errors.boot.wild)

      fit.mean <- as.vector(npindexhat(
        bws = bws,
        txdat = xdat,
        exdat = xdat,
        y = ydat,
        output = "apply",
        s = 0L
      ))
      H <- npindexhat(
        bws = bws,
        txdat = xdat,
        exdat = xdat,
        output = "matrix",
        s = if (gradients) 1L else 0L
      )

      boot.out <- .np_plot_boot_from_hat_wild(
        H = H,
        ydat = ydat,
        fit.mean = fit.mean,
        B = plot.errors.boot.num,
        wild = plot.errors.boot.wild
      )
    } else if (is.inid) {
      inid.helper.ok <- (!isTRUE(gradients)) || identical(bws$type, "fixed")
      if (!isTRUE(inid.helper.ok)) {
        stop("inid bootstrap requires helper mode with gradients=FALSE in compute.bootstrap.errors.sibandwidth", call. = FALSE)
      } else {
        boot.out <- tryCatch({
          .np_inid_boot_from_index(
            xdat = xdat,
            ydat = ydat,
            bws = bws,
            B = plot.errors.boot.num,
            gradients = gradients
          )
        }, error = function(e) {
          stop(sprintf("inid single-index helper failed in compute.bootstrap.errors.sibandwidth (%s)",
                       conditionMessage(e)),
               call. = FALSE)
        })
      }
    } else if (is.block) {
      block.helper.ok <- (!isTRUE(gradients)) || identical(bws$type, "fixed")
      if (!isTRUE(block.helper.ok)) {
        stop(sprintf("%s bootstrap requires helper mode with gradients=FALSE in compute.bootstrap.errors.sibandwidth", plot.errors.boot.method), call. = FALSE)
      } else {
        boot.out <- tryCatch({
          counts.drawer <- .np_block_counts_drawer(
            n = nrow(xdat),
            B = plot.errors.boot.num,
            blocklen = plot.errors.boot.blocklen,
            sim = plot.errors.boot.method
          )
          .np_inid_boot_from_index(
            xdat = xdat,
            ydat = ydat,
            bws = bws,
            B = plot.errors.boot.num,
            counts.drawer = counts.drawer,
            gradients = gradients
          )
        }, error = function(e) {
          stop(sprintf("%s single-index helper failed in compute.bootstrap.errors.sibandwidth (%s)",
                       plot.errors.boot.method,
                       conditionMessage(e)),
               call. = FALSE)
        })
      }
    }

    if (is.null(boot.out))
      stop(sprintf("unresolved bootstrap execution path for method '%s' in compute.bootstrap.errors.sibandwidth", plot.errors.boot.method), call. = FALSE)
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      interval.summary <- .np_plot_bootstrap_interval_summary(
        boot.t = boot.out$t,
        t0 = boot.out$t0,
        alpha = plot.errors.alpha,
        band.type = plot.errors.type
      )
      boot.err[,1:2] <- interval.summary$err
      boot.all.err <- interval.summary$all.err
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = list(), boot.all.err = boot.all.err)
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
