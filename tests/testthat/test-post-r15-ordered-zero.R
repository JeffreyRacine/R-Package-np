test_that("ordered zero initializes lookup before cached-category reuse", {
  withr::local_options(np.messages = FALSE)
  set.seed(19312)
  x <- data.frame(x = runif(36, -1, 1), o = ordered(rep(0:1, 18)))
  y <- x$x + .6 * as.integer(x$o) + rnorm(36, sd = .3)
  for (family in c("npreg", "npcdens", "npcdist")) {
    conditional <- family != "npreg"
    for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
      h <- if (type == "fixed") .7 else 24
      b <- do.call(get(paste0(family, "bw")), list(xdat = x,
        ydat = if (conditional) data.frame(y = y) else y,
        bws = c(if (conditional) h, h, .3), bwtype = type, bandwidth.compute = FALSE))
      ex <- x[c(1, 3, 2, 4), , drop = FALSE]; ex$x <- ex$x + .037
      run <- function(xx, gradients = FALSE) {
        a <- list(bws = b, txdat = x, tydat = if (conditional) data.frame(y = y) else y,
                  exdat = xx, gradients = gradients, se = gradients)
        if (conditional) a$eydat <- data.frame(y = rep(.5, nrow(xx)))
        do.call(get(family), a)
      }
      lo <- hi <- ex
      lo$o <- ordered(rep(0, 4), levels = 0:1)
      hi$o <- ordered(rep(1, 4), levels = 0:1)
      oracle <- fitted(run(hi)) - fitted(run(lo))
      batch <- run(ex, TRUE)
      expect_equal(as.numeric(gradients(batch)[, 2]), as.numeric(oracle), tolerance = 1e-11,
        info = paste(family, type))
      expect_equal(as.numeric(gradients(run(ex[1, , drop = FALSE], TRUE))[, 2]),
                   as.numeric(oracle[1]), tolerance = 1e-11)
      perm <- c(3, 1, 4, 2)
      expect_equal(gradients(run(ex[perm, , drop = FALSE], TRUE))[order(perm), 2],
                   gradients(batch)[, 2], tolerance = 1e-11)
    }
  }
})
