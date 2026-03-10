test_that("fixed conditional ll/lp density-distribution helper matches duplicate-sample refits", {
  skip_if_not_installed("np")

  helper <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")

  run_case <- function(family = c("dens", "dist"), regtype = c("ll", "lp")) {
    family <- match.arg(family)
    regtype <- match.arg(regtype)

    set.seed(switch(
      paste(family, regtype, sep = "_"),
      dens_ll = 603301L,
      dist_ll = 603302L,
      dens_lp = 603303L,
      dist_lp = 603304L
    ))

    n <- 55L
    B <- 7L
    tx <- data.frame(x = sort(runif(n)))
    ty <- data.frame(y = 0.5 * tx$x + rnorm(n, sd = 0.2))
    ex <- data.frame(x = seq(min(tx$x), max(tx$x), length.out = 13L))
    ey <- data.frame(y = seq(as.numeric(quantile(ty$y, 0.15)), as.numeric(quantile(ty$y, 0.85)), length.out = 13L))
    counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

    bw.args <- list(
      xdat = tx,
      ydat = ty,
      bws = c(0.28, 0.28),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = regtype
    )
    if (identical(regtype, "lp")) {
      bw.args$basis <- "glp"
      bw.args$degree <- 2L
    }

    bws <- if (identical(family, "dens")) {
      do.call(npcdensbw, bw.args)
    } else {
      do.call(npcdistbw, bw.args)
    }

    fast <- helper(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = bws,
      B = B,
      cdf = identical(family, "dist"),
      counts = counts
    )
    fast.drawer <- helper(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = bws,
      B = B,
      cdf = identical(family, "dist"),
      counts.drawer = function(start, stop) counts[, start:stop, drop = FALSE]
    )

    explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      explicit[b, ] <- if (identical(family, "dens")) {
        npcdens(
          txdat = tx[idx, , drop = FALSE],
          tydat = ty[idx, , drop = FALSE],
          exdat = ex,
          eydat = ey,
          bws = bws
        )$condens
      } else {
        npcdist(
          txdat = tx[idx, , drop = FALSE],
          tydat = ty[idx, , drop = FALSE],
          exdat = ex,
          eydat = ey,
          bws = bws
        )$condist
      }
    }

    t0 <- if (identical(family, "dens")) {
      npcdens(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = bws)$condens
    } else {
      npcdist(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = bws)$condist
    }

    expect_equal(fast$t, explicit, tolerance = 1e-10, info = paste(family, regtype))
    expect_equal(fast$t0, t0, tolerance = 1e-12, info = paste(family, regtype))
    expect_equal(fast.drawer$t, explicit, tolerance = 1e-10, info = paste(family, regtype, "drawer"))
    expect_equal(fast.drawer$t0, t0, tolerance = 1e-12, info = paste(family, regtype, "drawer"))
  }

  run_case("dens", "ll")
  run_case("dist", "ll")
  run_case("dens", "lp")
  run_case("dist", "lp")
})
