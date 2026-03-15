suppressPackageStartupMessages(library(np))

exact_train_kmax <- function(x) {
  u <- sort(unique(as.numeric(x)))
  if (length(u) <= 1L) {
    return(0L)
  }

  caps <- vapply(seq_along(u), function(i) {
    length(unique(abs(u[-i] - u[i])))
  }, integer(1))

  min(caps)
}

test_that("exact current-core kmax matches the toy unique-radius oracle", {
  x <- data.frame(x = c(0, 1, 2, 3))
  y <- c(0, 1, 2, 3)
  kmax <- exact_train_kmax(x$x)

  expect_identical(kmax, 2L)

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    bw.ok <- npregbw(
      xdat = x,
      ydat = y,
      regtype = "lc",
      bwtype = bwtype,
      bws = kmax,
      bandwidth.compute = FALSE
    )

    bw.bad <- npregbw(
      xdat = x,
      ydat = y,
      regtype = "lc",
      bwtype = bwtype,
      bws = kmax + 1L,
      bandwidth.compute = FALSE
    )

    expect_s3_class(npreg(bws = bw.ok, txdat = x, tydat = y), "npregression")
    expect_error(
      npreg(bws = bw.bad, txdat = x, tydat = y),
      "invalid bandwidth",
      info = bwtype
    )
  }
})

test_that("nonfixed regression search respects exact current-core kmax on repeated data", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  data(cps71, package = "np")

  x <- data.frame(age = cps71$age)
  y <- cps71$logwage
  kmax <- exact_train_kmax(x$age)

  expect_identical(kmax, 22L)

  bw.adap <- npregbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bwtype = "adaptive_nn",
    nmulti = 1
  )
  bw.gen <- npregbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bwtype = "generalized_nn",
    nmulti = 1
  )

  expect_gte(as.integer(round(unname(bw.adap$bw))), 1L)
  expect_lte(as.integer(round(unname(bw.adap$bw))), kmax)
  expect_gte(as.integer(round(unname(bw.gen$bw))), 1L)
  expect_lte(as.integer(round(unname(bw.gen$bw))), kmax)
})
