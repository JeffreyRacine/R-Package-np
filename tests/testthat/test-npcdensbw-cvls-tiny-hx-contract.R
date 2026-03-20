library(np)

test_that("public npcdensbw cv.ls bounded LP does not reward numerically singular tiny hx", {
  set.seed(42)
  n <- 350L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rbeta(n, 1, 1))

  bw.bad <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.0610566549, 0.0009802178),
    bwtype = "fixed",
    bwmethod = "cv.ls",
    regtype = "lp",
    degree = 3,
    basis = "glp",
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )
  bw.good <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(1.20411811, 0.05536272),
    bwtype = "fixed",
    bwmethod = "cv.ls",
    regtype = "lp",
    degree = 3,
    basis = "glp",
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )

  obj.bad <- np:::.npcdensbw_eval_only(x, y, bw.bad)$objective
  obj.good <- np:::.npcdensbw_eval_only(x, y, bw.good)$objective

  expect_true(is.finite(obj.bad))
  expect_true(is.finite(obj.good))
  expect_lte(obj.bad, obj.good + 1e-10)
})
