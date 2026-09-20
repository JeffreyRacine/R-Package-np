test_that("density equality cross normalization does not overflow integer row counts", {
  owner <- getFromNamespace("npdeneqtest", "np")
  local <- new.env(parent = environment(owner))
  local$.npksum_power12 <- function(txdat, exdat = NULL, ...) {
    index <- if (!is.null(exdat)) 3L else if (txdat[[1L]][1L] == 1) 1L else 2L
    list(ksum = c(3.25, 7.5, 1.125)[index],
         ksum.power2 = c(4.5, 8.75, 1.75)[index])
  }
  statistic <- NULL
  for (expr in as.list(body(owner))[-1L])
    if (is.call(expr) && identical(expr[[1L]], quote(`<-`)) &&
        identical(expr[[2L]], quote(teststat)))
      statistic <- eval(expr[[3L]], local)
  expect_true(is.function(statistic))
  for (sizes in list(c(46340L, 46340L), c(46341L, 46341L),
                    c(43000L, 50000L), c(50000L, 43000L))) {
    x <- data.frame(x = rep.int(1, sizes[1L]))
    y <- data.frame(x = rep.int(2, sizes[2L]))
    n <- as.double(sizes[1L]); m <- as.double(sizes[2L])
    In <- 3.25/(n*(n-1)) + 7.5/(m*(m-1)) - 2*1.125/(n*m)
    variance <- 2*(4.5/(n^2*(n-1)^2) + 8.75/(m^2*(m-1)^2) +
                    2*1.75/(n^2*m^2))
    expect_identical(statistic(x, y, .4, .4), list(Tn = In/sqrt(variance), In = In))
  }
})
