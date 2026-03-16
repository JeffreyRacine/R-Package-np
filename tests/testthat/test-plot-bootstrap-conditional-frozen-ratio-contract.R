test_that("npRmpi conditional frozen helper preserves numerator/denominator recombination", {
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  options(npRmpi.autodispatch = FALSE, np.messages = FALSE)

  ns <- asNamespace("npRmpi")
  frozen.fun <- get(".np_inid_boot_from_hat_conditional_frozen", envir = ns)
  make_xkbw <- get(".npcdhat_make_xkbw", envir = ns)
  make_ybw <- get(".npcdhat_make_ybw", envir = ns)
  make_kmat <- get(".npcdhat_make_kernel_matrix", envir = ns)
  make_kbx <- get(".np_con_make_kbandwidth_x", envir = ns)
  make_kbxy <- get(".np_con_make_kbandwidth_xy", envir = ns)
  exact_state_build <- get(".np_ksum_exact_state_build", envir = ns)
  exact_state_eval <- get(".np_ksum_eval_exact_state", envir = ns)
  bind_fast <- get(".np_bind_data_frames_fast", envir = ns)

  xdat <- data.frame(x = c(0.05, 0.15, 0.30, 0.45, 0.70, 0.90))
  ydat <- data.frame(y = c(0.10, 0.20, 0.35, 0.50, 0.72, 0.88))
  counts <- cbind(
    c(2, 0, 1, 1, 1, 1),
    c(0, 2, 1, 1, 1, 1),
    c(1, 1, 2, 0, 1, 1)
  )
  exdat <- data.frame(x = c(0.10, 0.35, 0.80))
  eydat <- data.frame(y = c(0.12, 0.40, 0.82))

  manual_boot <- function(bw, cdf) {
    xkbw <- make_xkbw(bws = bw, txdat = xdat)
    ykbw <- make_ybw(bws = bw, tydat = ydat)
    Kx <- make_kmat(
      kbw = xkbw,
      txdat = xdat,
      exdat = exdat,
      operator = rep.int("normal", ncol(xdat))
    )
    Ky <- make_kmat(
      kbw = ykbw,
      txdat = ydat,
      exdat = eydat,
      operator = rep.int(if (isTRUE(cdf)) "integral" else "normal", ncol(ydat))
    )

    n <- nrow(xdat)
    den.op <- Kx / n
    num.op <- (Kx * Ky) / n
    den <- t(den.op %*% counts)
    num <- t(num.op %*% counts)

    list(
      t = num / pmax(den, .Machine$double.eps),
      t0 = rowSums(num.op) / pmax(rowSums(den.op), .Machine$double.eps)
    )
  }

  weighted_state_boot <- function(bw, cdf) {
    kbx <- make_kbx(bws = bw, xdat = xdat)
    kbxy <- make_kbxy(bws = bw, xdat = xdat, ydat = ydat)
    den.state <- exact_state_build(
      bws = kbx,
      exdat = exdat,
      operator = rep.int("normal", ncol(xdat))
    )
    num.state <- exact_state_build(
      bws = kbxy,
      exdat = bind_fast(exdat, eydat),
      operator = c(
        rep.int("normal", ncol(xdat)),
        rep.int(if (isTRUE(cdf)) "integral" else "normal", ncol(ydat))
      )
    )

    n.total <- colSums(counts)
    tmat <- t(vapply(seq_len(ncol(counts)), function(j) {
      den <- as.numeric(exact_state_eval(
        state = den.state,
        txdat = xdat,
        weights = matrix(as.double(counts[, j]), ncol = 1L)
      )) / n.total[j]
      num <- as.numeric(exact_state_eval(
        state = num.state,
        txdat = bind_fast(xdat, ydat),
        weights = matrix(as.double(counts[, j]), ncol = 1L)
      )) / n.total[j]
      num / pmax(den, .Machine$double.eps)
    }, numeric(nrow(exdat))))

    den0 <- as.numeric(exact_state_eval(
      state = den.state,
      txdat = xdat,
      weights = matrix(1.0, nrow = nrow(xdat), ncol = 1L)
    )) / nrow(xdat)
    num0 <- as.numeric(exact_state_eval(
      state = num.state,
      txdat = bind_fast(xdat, ydat),
      weights = matrix(1.0, nrow = nrow(xdat), ncol = 1L)
    )) / nrow(xdat)

    list(
      t = tmat,
      t0 = num0 / pmax(den0, .Machine$double.eps)
    )
  }

  cases <- expand.grid(
    bwtype = c("generalized_nn", "adaptive_nn"),
    cdf = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for (ii in seq_len(nrow(cases))) {
    bwtype.ii <- cases$bwtype[ii]
    cdf.ii <- isTRUE(cases$cdf[ii])

    bw <- if (cdf.ii) {
      suppressWarnings(npcdistbw(
        xdat = xdat,
        ydat = ydat,
        bws = c(3L, 3L),
        bwtype = bwtype.ii,
        bandwidth.compute = FALSE
      ))
    } else {
      suppressWarnings(npcdensbw(
        xdat = xdat,
        ydat = ydat,
        bws = c(3L, 3L),
        bwtype = bwtype.ii,
        bandwidth.compute = FALSE
      ))
    }

    helper.out <- frozen.fun(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      bws = bw,
      B = ncol(counts),
      cdf = cdf.ii,
      counts = counts
    )
    manual.out <- manual_boot(bw = bw, cdf = cdf.ii)
    weighted.out <- weighted_state_boot(bw = bw, cdf = cdf.ii)

    expect_equal(
      helper.out$t0,
      manual.out$t0,
      tolerance = 1e-12,
      info = paste(bwtype.ii, if (cdf.ii) "dist" else "dens", "t0")
    )
    expect_equal(
      helper.out$t,
      manual.out$t,
      tolerance = 1e-12,
      info = paste(bwtype.ii, if (cdf.ii) "dist" else "dens", "t")
    )
    expect_equal(
      helper.out$t0,
      weighted.out$t0,
      tolerance = 1e-12,
      info = paste(bwtype.ii, if (cdf.ii) "dist" else "dens", "weighted-state t0")
    )
    expect_equal(
      helper.out$t,
      weighted.out$t,
      tolerance = 1e-12,
      info = paste(bwtype.ii, if (cdf.ii) "dist" else "dens", "weighted-state t")
    )
  }
})
