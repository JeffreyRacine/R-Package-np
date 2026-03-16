library(np)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

manual_conditional_frozen_boot <- function(bw, xdat, ydat, exdat, eydat, counts, cdf) {
  ns <- asNamespace("np")
  make_xkbw <- get(".npcdhat_make_xkbw", envir = ns)
  make_ybw <- get(".npcdhat_make_ybw", envir = ns)
  make_kmat <- get(".npcdhat_make_kernel_matrix", envir = ns)

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

test_that("conditional frozen helper preserves numerator/denominator recombination for nonfixed dens/dist", {
  ns <- asNamespace("np")
  frozen.fun <- get(".np_inid_boot_from_hat_conditional_frozen", envir = ns)

  xdat <- data.frame(x = c(0.05, 0.15, 0.30, 0.45, 0.70, 0.90))
  ydat <- data.frame(y = c(0.10, 0.20, 0.35, 0.50, 0.72, 0.88))
  counts <- cbind(
    c(2, 0, 1, 1, 1, 1),
    c(0, 2, 1, 1, 1, 1),
    c(1, 1, 2, 0, 1, 1)
  )
  exdat <- data.frame(x = c(0.10, 0.35, 0.80))
  eydat <- data.frame(y = c(0.12, 0.40, 0.82))

  cases <- expand.grid(
    bwtype = c("generalized_nn", "adaptive_nn"),
    cdf = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for (ii in seq_len(nrow(cases))) {
    bw <- quiet_eval(
      if (isTRUE(cases$cdf[ii])) {
        npcdistbw(
          xdat = xdat,
          ydat = ydat,
          bws = c(3L, 3L),
          bwtype = cases$bwtype[ii],
          bandwidth.compute = FALSE
        )
      } else {
        npcdensbw(
          xdat = xdat,
          ydat = ydat,
          bws = c(3L, 3L),
          bwtype = cases$bwtype[ii],
          bandwidth.compute = FALSE
        )
      }
    )

    helper.out <- frozen.fun(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      bws = bw,
      B = ncol(counts),
      cdf = cases$cdf[ii],
      counts = counts
    )
    manual.out <- manual_conditional_frozen_boot(
      bw = bw,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      counts = counts,
      cdf = cases$cdf[ii]
    )

    expect_equal(
      helper.out$t0,
      manual.out$t0,
      tolerance = 1e-12,
      info = paste(cases$bwtype[ii], if (isTRUE(cases$cdf[ii])) "dist" else "dens", "t0")
    )
    expect_equal(
      helper.out$t,
      manual.out$t,
      tolerance = 1e-12,
      info = paste(cases$bwtype[ii], if (isTRUE(cases$cdf[ii])) "dist" else "dens", "t")
    )
  }
})
