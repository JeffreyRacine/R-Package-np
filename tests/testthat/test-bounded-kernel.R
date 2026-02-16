library(np)

test_that("bounded gaussian kernel matches direct normalization at fixed bandwidth", {
  x <- sort(c(5, 11, 21, 31, 46, 75, 98, 122, 145, 165, 195, 224, 245, 293, 321, 330, 350, 420))
  h <- 100
  a <- min(x)
  b <- max(x)

  bw <- npudensbw(dat = data.frame(x = x),
                  bws = h,
                  bandwidth.compute = FALSE,
                  ckertype = "gaussian",
                  ckerbound = "range")
  fit <- npudens(tdat = data.frame(x = x), bws = bw)

  fdirect <- sapply(x, function(xx) {
    mean(dnorm((x - xx)/h)/(h * (pnorm((b - xx)/h) - pnorm((a - xx)/h))))
  })

  expect_equal(as.numeric(fit$dens), as.numeric(fdirect), tolerance = 1e-8)
})

test_that("cv bandwidth with range bounds is larger than unbounded for stress example", {
  x <- sort(c(5, 11, 21, 31, 46, 75, 98, 122, 145, 165, 195, 224, 245, 293, 321, 330, 350, 420))

  bw_bound <- npudensbw(dat = data.frame(x = x),
                        bwmethod = "cv.ml",
                        ckertype = "gaussian",
                        ckerbound = "range",
                        nmulti = 1)

  bw_inf <- npudensbw(dat = data.frame(x = x),
                      bwmethod = "cv.ml",
                      ckertype = "gaussian",
                      ckerbound = "fixed",
                      ckerlb = -Inf,
                      ckerub = Inf,
                      nmulti = 1)

  expect_gt(as.numeric(bw_bound$bw[1]), as.numeric(bw_inf$bw[1]))
})
