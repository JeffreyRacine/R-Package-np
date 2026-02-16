test_that("bounded conditional kernels change npcdens on bounded-support DGP", {
  library(np)
  set.seed(20260216)

  n <- 120
  x <- runif(n)
  y <- rbeta(n, shape1 = 1 + 4 * (x - 0.5)^2, shape2 = 2 - 4 * (x - 0.5)^2)
  dat <- data.frame(y = y, x = x)

  b_emp <- npcdensbw(
    y ~ x,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cykerbound = "range"
  )
  b_inf <- npcdensbw(
    y ~ x,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cykerbound = "fixed",
    cykerlb = -Inf,
    cykerub = Inf
  )

  f_emp <- npcdens(bws = b_emp)
  f_inf <- npcdens(bws = b_inf)

  grid <- seq(min(y), max(y), length.out = 120)
  nd <- data.frame(y = grid, x = rep(0.5, length(grid)))
  p_emp <- predict(f_emp, newdata = nd)
  p_inf <- predict(f_inf, newdata = nd)

  expect_true(max(abs(p_emp - p_inf)) > 1e-6)
})
