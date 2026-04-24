library(npRmpi)

make_regression_summary_label_fixture <- function(type) {
  untangle <- getFromNamespace("untangle", "npRmpi")
  rbandwidth <- getFromNamespace("rbandwidth", "npRmpi")
  npregression <- getFromNamespace("npregression", "npRmpi")

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

test_that("fixed mixed-data bandwidth summaries label continuous scale factors correctly", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

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
  expect_match(s, "Exp\\. Var\\. Name: z2\\s+Bandwidth: .*Lambda Max:")
})

test_that("regression summaries distinguish nearest-neighbor selectors from fixed bandwidths", {
  fixed <- paste(capture.output(summary(make_regression_summary_label_fixture("fixed"))), collapse = "\n")
  adaptive <- paste(capture.output(summary(make_regression_summary_label_fixture("adaptive_nn"))), collapse = "\n")
  generalized <- paste(capture.output(summary(make_regression_summary_label_fixture("generalized_nn"))), collapse = "\n")

  expect_match(fixed, "Bandwidth\\(s\\):\\s+2")
  expect_match(adaptive, "Bandwidth Nearest Neighbor\\(s\\):\\s+2")
  expect_match(generalized, "Bandwidth Nearest Neighbor\\(s\\):\\s+2")
})
