test_that("ordered contrasts validate evaluation indices before native work", {
  withr::local_options(np.messages = FALSE)
  x <- data.frame(o = ordered(rep(0:1, 10)))
  y <- seq_len(nrow(x))
  for (kernel in c("wangvanryzin", "liracine", "racine")) {
    b <- npregbw(xdat = x, ydat = y, bws = .3,
                 okertype = kernel, bandwidth.compute = FALSE)
    for (values in list(c(2, 0), c(0, 2), c(2, 2))) {
      ex <- data.frame(o = ordered(values, levels = 0:2))
      expect_error(npksum(txdat = x, exdat = ex, bws = b, compute.ocg = TRUE),
                   "ordered evaluation values in the bandwidth category table")
      expect_true(all(is.finite(npksum(txdat = x, exdat = ex, bws = b)$ksum)))
    }
    ex <- x[c(2, 1, 2), , drop = FALSE]
    a <- npksum(txdat = x, exdat = ex, bws = b, compute.ocg = TRUE)
    z <- npksum(txdat = x, exdat = ex[3:1, , drop = FALSE], bws = b,
               compute.ocg = TRUE)
    expect_equal(as.numeric(a$ksum), rev(as.numeric(z$ksum)))
    expect_equal(as.numeric(a$p.ksum), rev(as.numeric(z$p.ksum)))
  }
})

test_that("ordered contrast indexing preserves declared unused categories", {
  withr::local_options(np.messages = FALSE)
  x <- data.frame(o = ordered(rep(c(0, 2), 10), levels = 0:2))
  b <- npregbw(xdat = x, ydat = seq_len(20), bws = .3,
               bandwidth.compute = FALSE)
  ex <- data.frame(o = ordered(c(1, 0), levels = 0:2))
  expect_no_error(npksum(txdat = x, exdat = ex, bws = b, compute.ocg = TRUE))
})
