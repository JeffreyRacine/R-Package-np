test_that("missing responses are identity columns in weighted kernel sums", {
  x <- data.frame(x = c(.03, .11, .27, .39, .56, .72, .91))
  for (neval in c(1L, 3L)) {
    ex <- data.frame(x = c(.13, .45, .83)[seq_len(neval)])
    for (nw in c(1L, 2L, 3L)) {
      w <- outer(seq_len(nrow(x)), seq_len(nw), function(i, j) i / 9 + j)
      for (power in c(1, 2)) {
        a <- list(txdat = x, exdat = ex, bws = .24, weights = w,
                  kernel.pow = power, permutation.operator = "derivative")
        implicit <- do.call(npksum, a)
        explicit <- do.call(npksum, c(a, list(tydat = rep(1, nrow(x)))))
        expect_identical(implicit$ksum, explicit$ksum)
        expect_identical(implicit$p.ksum, explicit$p.ksum)
        k <- outer(x$x, ex$x, function(t, e) dnorm((e - t) / .24))^power
        expected <- crossprod(w, k)
        expect_equal(as.vector(implicit$ksum), as.vector(expected), tolerance = 1e-13)
        if (nw > 1L) expect_identical(dim(implicit$ksum), c(nw, neval))
      }
    }
    a <- list(txdat = x, exdat = ex, bws = .24)
    expect_identical(do.call(npksum, a)$ksum,
                     do.call(npksum, c(a, list(tydat = rep(1, nrow(x)))) )$ksum)
    # Public axes are weight, response, evaluation, including singleton axes.
    for (ny in c(1L, 2L)) {
      y <- outer(seq_len(nrow(x)), seq_len(ny), function(i, j) i - j / 2)
      w <- cbind(seq_len(nrow(x)), -seq_len(nrow(x)))
      fit <- npksum(txdat = x, exdat = ex, tydat = y, weights = w, bws = .24)$ksum
      k <- outer(x$x, ex$x, function(t, e) dnorm((e - t) / .24))
      expected <- array(0, c(2L, ny, neval))
      for (i in seq_len(neval)) expected[, , i] <- crossprod(w, y * k[, i])
      expect_equal(as.vector(fit), as.vector(expected), tolerance = 1e-13)
      expect_identical(dim(fit), if (ny == 1L) c(2L, neval) else c(2L, ny, neval))
    }
  }
})
