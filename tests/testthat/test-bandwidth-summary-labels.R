library(np)

make_regression_summary_label_fixture <- function(type) {
  untangle <- getFromNamespace("untangle", "np")
  rbandwidth <- getFromNamespace("rbandwidth", "np")
  npregression <- getFromNamespace("npregression", "np")

  xdat <- data.frame(x = c(0.1, 0.4, 0.9))
  ydat <- data.frame(y = c(1, 2, 3))
  xdati <- untangle(xdat)
  ydati <- untangle(ydat)

  bws <- rbandwidth(
    bw = 2,
    regtype = "lc",
    bwmethod = "cv.ls",
    bwscaling = FALSE,
    bwtype = type,
    ckertype = "gaussian",
    ckerorder = 2,
    ukertype = "aitchisonaitken",
    okertype = "liracine",
    nobs = nrow(xdat),
    xdati = xdati,
    ydati = ydati,
    xnames = names(xdat),
    ynames = names(ydat),
    sfactor = 1,
    bandwidth = if (identical(type, "fixed")) 0.25 else 2,
    nconfac = 1L,
    ncatfac = 0L,
    sdev = stats::sd(xdat$x),
    bandwidth.compute = FALSE
  )

  npregression(
    bws = bws,
    eval = as.matrix(xdat),
    mean = ydat$y,
    ntrain = nrow(xdat),
    trainiseval = TRUE
  )
}

make_conditional_polynomial_label_fixture <- function(kind) {
  untangle <- getFromNamespace("untangle", "np")
  ctor <- getFromNamespace(
    switch(kind,
           density = "conbandwidth",
           distribution = "condbandwidth"),
    "np"
  )

  xdat <- data.frame(x = c(0.1, 0.4, 0.9, 1.2, 1.7))
  ydat <- data.frame(y = c(0.2, 0.3, 0.8, 1.1, 1.5))
  xdati <- untangle(xdat)
  ydati <- untangle(ydat)

  ctor(
    xbw = 0.5,
    ybw = 0.4,
    bwmethod = "manual",
    bwtype = "fixed",
    nobs = nrow(xdat),
    xdati = xdati,
    ydati = ydati,
    xnames = names(xdat),
    ynames = names(ydat),
    sfactor = list(x = 2.8, y = 2.6),
    bandwidth = list(x = 0.5, y = 0.4),
    nconfac = 0.75,
    ncatfac = 0.5,
    sdev = c(xcon = stats::sd(xdat$x), ycon = stats::sd(ydat$y)),
    bandwidth.compute = FALSE,
    regtype = "lp",
    pregtype = "Local-Polynomial",
    basis = "glp",
    degree = 2L
  )
}

make_quantile_polynomial_label_fixture <- function() {
  qregression <- getFromNamespace("qregression", "np")
  qregression(
    bws = make_conditional_polynomial_label_fixture("distribution"),
    xeval = matrix(seq_len(5), ncol = 1L),
    tau = 0.5,
    quantile = rep(0.5, 5),
    ntrain = 5L,
    trainiseval = TRUE
  )
}

test_that("fixed mixed-data bandwidth summaries label continuous scale factors correctly", {
  set.seed(123)
  n <- 40L
  xdat <- data.frame(
    x = runif(n),
    z1 = runif(n),
    z2 = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  ydat <- xdat$x + xdat$z1 + as.numeric(xdat$z2 == "b") + rnorm(n)

  bw <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "lc",
    bwmethod = "cv.aic",
    nmulti = 1
  )

  s <- paste(capture.output(summary(bw)), collapse = "\n")

  expect_match(s, "Exp\\. Var\\. Name: x\\s+Bandwidth: .*Scale Factor:")
  expect_match(s, "Exp\\. Var\\. Name: z1\\s+Bandwidth: .*Scale Factor:")
  expect_match(s, "Exp\\. Var\\. Name: z2\\s+Lambda: .*Lambda Max:")
})

test_that("regression summaries distinguish nearest-neighbor selectors from fixed bandwidths", {
  fixed <- paste(capture.output(summary(make_regression_summary_label_fixture("fixed"))), collapse = "\n")
  adaptive <- paste(capture.output(summary(make_regression_summary_label_fixture("adaptive_nn"))), collapse = "\n")
  generalized <- paste(capture.output(summary(make_regression_summary_label_fixture("generalized_nn"))), collapse = "\n")

  expect_match(fixed, "Kernel Regression Estimator:\\s+Local-Constant")
  expect_match(fixed, "Bandwidth\\(s\\):\\s+2")
  expect_match(adaptive, "Type\\s+NN Index")
  expect_match(adaptive, "Value\\s+2")
  expect_match(generalized, "Type\\s+NN Index")
  expect_match(generalized, "Value\\s+2")
})

test_that("density and distribution bandwidth summaries label polynomial type", {
  regression <- paste(capture.output(summary(make_regression_summary_label_fixture("fixed")$bws)), collapse = "\n")
  condensity <- paste(capture.output(summary(make_conditional_polynomial_label_fixture("density"))), collapse = "\n")
  condist <- paste(capture.output(summary(make_conditional_polynomial_label_fixture("distribution"))), collapse = "\n")
  quantile <- paste(capture.output(summary(make_quantile_polynomial_label_fixture())), collapse = "\n")

  expect_true(grepl("Regression Type: Local-Constant", regression, fixed = TRUE))
  expect_true(grepl("Polynomial Type: Local-Polynomial", condensity, fixed = TRUE))
  expect_true(grepl("Polynomial Type: Local-Polynomial", condist, fixed = TRUE))
  expect_true(grepl("Polynomial Type: Local-Polynomial", quantile, fixed = TRUE))
  expect_false(grepl("Regression Type:", condensity, fixed = TRUE))
  expect_false(grepl("Regression Type:", condist, fixed = TRUE))
  expect_false(grepl("Kernel Regression Estimator:", quantile, fixed = TRUE))
})
