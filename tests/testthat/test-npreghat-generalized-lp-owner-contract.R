library(np)

test_that("generalized higher-degree lp owner matches npreg on evaluation data", {
  run_case <- function(train_df, y, eval_df, degree, basis = "glp", bern = FALSE, tol = 1e-10) {
    frame <- train_df
    frame$y <- y
    fml <- as.formula(paste("y~", paste(names(train_df), collapse = "+")))
    bw <- npregbw(fml,
                  data = frame,
                  regtype = "lp",
                  degree = degree,
                  basis = basis,
                  bernstein.basis = bern,
                  bwtype = "generalized_nn")

    H <- npreghat(bws = bw, txdat = train_df, exdat = eval_df)
    a <- npreghat(bws = bw, txdat = train_df, exdat = eval_df, y = y, output = "apply")
    g <- npreg(bws = bw, exdat = eval_df)

    expect_equal(drop(H %*% y), a, tolerance = tol)
    expect_equal(drop(H %*% y), g$mean, tolerance = tol)
  }

  set.seed(101)
  n <- 80
  x <- runif(n, -1, 1)
  y <- x^2 + rnorm(n, sd = 0.15)
  xe <- data.frame(x = seq(-1, 1, length.out = 41))
  run_case(data.frame(x = x), y, xe, degree = 2L)
  run_case(data.frame(x = x), y, xe, degree = 3L)

  set.seed(102)
  n <- 70
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  y <- x1^2 - x2 + x1 * x2 + rnorm(n, sd = 0.1)
  xe <- expand.grid(x1 = seq(-1, 1, length.out = 6),
                    x2 = seq(-1, 1, length.out = 6))
  for (b in c("glp", "additive", "tensor")) {
    run_case(data.frame(x1 = x1, x2 = x2), y, xe, degree = c(2L, 2L), basis = b)
  }

  set.seed(103)
  n <- 70
  x1 <- c(0, 1, runif(n - 2L))
  x2 <- c(1, 0, runif(n - 2L))
  y <- sin(pi * x1) + x2^2 + rnorm(n, sd = 0.1)
  xe <- expand.grid(x1 = seq(0.1, 0.9, length.out = 6),
                    x2 = seq(0.1, 0.9, length.out = 6))
  run_case(data.frame(x1 = x1, x2 = x2),
           y,
           xe,
           degree = c(2L, 2L),
           basis = "tensor",
           bern = TRUE,
           tol = 1e-9)

  set.seed(104)
  n <- 70
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  u <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  o <- ordered(sample(1:3, n, replace = TRUE))
  y <- x1^2 + x2 - 0.5 * (u == "b") + 0.3 * (o == 3) + rnorm(n, sd = 0.1)
  xe <- data.frame(
    x1 = seq(-1, 1, length.out = 12),
    x2 = seq(-1, 1, length.out = 12),
    u = factor(rep(c("a", "b", "c"), length.out = 12), levels = levels(u)),
    o = ordered(rep(1:3, length.out = 12), levels = levels(o))
  )
  run_case(data.frame(x1 = x1, x2 = x2, u = u, o = o),
           y,
           xe,
           degree = c(2L, 2L),
           basis = "tensor")
})

test_that("tree-enabled generalized higher-degree lp owner stays exact via core fallback", {
  old_tree <- getOption("np.tree")
  options(np.tree = TRUE)
  on.exit(options(np.tree = old_tree), add = TRUE)

  set.seed(203)
  n <- 80
  x <- runif(n, -1, 1)
  y <- x^2 + rnorm(n, sd = 0.2)
  d <- data.frame(x = x)
  bw <- npregbw(y ~ x,
                data = data.frame(d, y = y),
                regtype = "lp",
                degree = 2L,
                basis = "glp",
                bwtype = "generalized_nn")

  H <- npreghat(bws = bw, txdat = d)
  a <- npreghat(bws = bw, txdat = d, y = y, output = "apply")
  g <- npreg(bws = bw)

  expect_equal(drop(H %*% y), a, tolerance = 1e-10)
  expect_equal(drop(H %*% y), g$mean, tolerance = 1e-10)
})

test_that("generalized higher-degree lp derivative owner matches npreg on evaluation data", {
  run_case <- function(train_df, y, eval_df, degree, s, basis = "glp", bern = FALSE, tol = 1e-9) {
    frame <- train_df
    frame$y <- y
    fml <- as.formula(paste("y~", paste(names(train_df), collapse = "+")))
    bw <- npregbw(fml,
                  data = frame,
                  regtype = "lp",
                  degree = degree,
                  basis = basis,
                  bernstein.basis = bern,
                  bwtype = "generalized_nn")

    H <- npreghat(bws = bw, txdat = train_df, exdat = eval_df, s = s)
    a <- npreghat(bws = bw, txdat = train_df, exdat = eval_df, y = y, output = "apply", s = s)
    g <- npreg(bws = bw, exdat = eval_df, gradients = TRUE)
    target.col <- which(as.integer(s) == 1L)[1L]

    expect_equal(drop(H %*% y), a, tolerance = tol)
    expect_equal(drop(H %*% y), g$grad[, target.col], tolerance = tol)
  }

  set.seed(301)
  n <- 80
  x <- runif(n, -1, 1)
  y <- x^2 + rnorm(n, sd = 0.15)
  xe <- data.frame(x = seq(-1, 1, length.out = 41))
  run_case(data.frame(x = x), y, xe, degree = 2L, s = 1L)
  run_case(data.frame(x = x), y, xe, degree = 3L, s = 1L)

  set.seed(302)
  n <- 70
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  y <- x1^2 - x2 + x1 * x2 + rnorm(n, sd = 0.1)
  xe <- expand.grid(x1 = seq(-1, 1, length.out = 6),
                    x2 = seq(-1, 1, length.out = 6))
  run_case(data.frame(x1 = x1, x2 = x2),
           y,
           xe,
           degree = c(2L, 2L),
           s = c(1L, 0L),
           basis = "tensor")
})
