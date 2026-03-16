library(np)

test_that("lc first-derivative owner matches npreg and apply on common cells", {
  run_case <- function(bwtype, bwval, tol = 1e-10) {
    train <- data.frame(x = x)
    frame <- data.frame(train, y = y)
    bw <- npregbw(y ~ x,
                  data = frame,
                  regtype = "lc",
                  bwtype = bwtype,
                  bandwidth.compute = FALSE,
                  bws = bwval)

    H.train <- npreghat(bws = bw, txdat = train, s = 1L)
    H.eval <- npreghat(bws = bw, txdat = train, exdat = eval, s = 1L)
    a.train <- npreghat(bws = bw, txdat = train, y = y, s = 1L, output = "apply")
    a.eval <- npreghat(bws = bw, txdat = train, exdat = eval, y = y, s = 1L, output = "apply")
    g.train <- npreg(bws = bw, gradients = TRUE)
    g.eval <- npreg(bws = bw, exdat = eval, gradients = TRUE)

    expect_equal(drop(H.train %*% y), a.train, tolerance = tol)
    expect_equal(drop(H.eval %*% y), a.eval, tolerance = tol)
    expect_equal(drop(H.train %*% y), g.train$grad[, 1L], tolerance = tol)
    expect_equal(drop(H.eval %*% y), g.eval$grad[, 1L], tolerance = tol)
  }

  set.seed(410)
  n <- 60L
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.2)
  eval <- data.frame(x = seq(min(x), max(x), length.out = 31L))

  run_case("fixed", 0.45)
  run_case("generalized_nn", 9)
  run_case("adaptive_nn", 9)
})
