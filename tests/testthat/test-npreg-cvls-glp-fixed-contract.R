fixed_cvls_objective <- function(xdat,
                                 ydat,
                                 regtype,
                                 degree,
                                 bws,
                                 basis = "glp",
                                 bernstein.basis = FALSE) {
  ns <- asNamespace("npRmpi")
  xdat <- as.data.frame(xdat)
  bw0 <- suppressMessages(suppressPackageStartupMessages(
    npRmpi::npregbw(
      xdat = xdat,
      ydat = ydat,
      bws = bws,
      regtype = regtype,
      degree = degree,
      basis = basis,
      bernstein.basis = bernstein.basis,
      bwmethod = "cv.ls",
      bandwidth.compute = FALSE
    )
  ))
  out <- get(".npregbw_eval_only", ns)(xdat, ydat, bw0, localize = TRUE)
  unname(out$objective[1L])
}

test_that("fixed cv.ls is exact on in-class local polynomial fixtures", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  tol <- 1e-20

  set.seed(20260308)
  x1 <- sort(runif(30))
  y1 <- 1 + 2 * x1
  obj.ll <- fixed_cvls_objective(
    xdat = data.frame(x = x1),
    ydat = y1,
    regtype = "ll",
    degree = 1L,
    bws = 0.25
  )
  obj.lp1 <- fixed_cvls_objective(
    xdat = data.frame(x = x1),
    ydat = y1,
    regtype = "lp",
    degree = 1L,
    bws = 0.25
  )
  expect_lt(abs(obj.ll), tol)
  expect_lt(abs(obj.lp1), tol)
  expect_equal(obj.ll, obj.lp1, tolerance = 1e-24)

  set.seed(20260309)
  x2 <- runif(36)
  z2 <- runif(36)
  y2 <- 0.5 + 1.2 * x2 - 0.7 * z2 + 0.8 * x2^2 - 0.3 * x2 * z2 + 0.6 * z2^2
  obj.lp2 <- fixed_cvls_objective(
    xdat = data.frame(x = x2, z = z2),
    ydat = y2,
    regtype = "lp",
    degree = c(2L, 2L),
    bws = c(0.3, 0.3)
  )
  expect_lt(abs(obj.lp2), tol)

  set.seed(20260310)
  x3 <- runif(40)
  z3 <- runif(40)
  w3 <- runif(40)
  y3 <- 0.2 + 0.9 * x3 + 0.6 * z3 - 0.4 * w3 + 0.7 * x3^2
  obj.lp211 <- fixed_cvls_objective(
    xdat = data.frame(x = x3, z = z3, w = w3),
    ydat = y3,
    regtype = "lp",
    degree = c(2L, 1L, 1L),
    bws = c(0.35, 0.35, 0.35)
  )
  expect_lt(abs(obj.lp211), tol)
})

test_that("fixed cv.ls keeps ll and canonical lp degree-1 aligned off-model", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260311)
  x <- sort(runif(48))
  y <- sin(2 * pi * x) + 0.15 * x

  obj.ll <- fixed_cvls_objective(
    xdat = data.frame(x = x),
    ydat = y,
    regtype = "ll",
    degree = 1L,
    bws = 0.22
  )
  obj.lp <- fixed_cvls_objective(
    xdat = data.frame(x = x),
    ydat = y,
    regtype = "lp",
    degree = 1L,
    bws = 0.22
  )

  expect_equal(obj.ll, obj.lp, tolerance = 1e-12)

  set.seed(20260312)
  x1 <- runif(52)
  x2 <- runif(52)
  u <- factor(sample(c("a", "b", "c"), 52, replace = TRUE))
  o <- ordered(sample(1:3, 52, replace = TRUE))
  y2 <- sin(2 * pi * x1) + 0.2 * x2 + as.numeric(u) / 9 + as.numeric(o) / 11

  obj2.ll <- fixed_cvls_objective(
    xdat = data.frame(x1 = x1, x2 = x2, u = u, o = o),
    ydat = y2,
    regtype = "ll",
    degree = c(1L, 1L),
    bws = c(0.28, 0.31, 0.45, 0.55)
  )
  obj2.lp <- fixed_cvls_objective(
    xdat = data.frame(x1 = x1, x2 = x2, u = u, o = o),
    ydat = y2,
    regtype = "lp",
    degree = c(1L, 1L),
    bws = c(0.28, 0.31, 0.45, 0.55)
  )

  expect_equal(obj2.ll, obj2.lp, tolerance = 1e-12)
})

test_that("fixed width-three LP tree uses canonical support and objective", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260809)
  n <- 128L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(2 * pi * x$x1) * cos(2 * pi * x$x2) + rnorm(n, sd = 0.1)
  bw <- npRmpi::npregbw(
    xdat = x,
    ydat = y,
    bws = c(0.24, 0.21),
    regtype = "ll",
    degree = c(1L, 1L),
    ckertype = "epanechnikov",
    bwmethod = "cv.ls",
    bandwidth.compute = FALSE
  )

  old_opts <- options(np.tree = TRUE, np.categorical.compress = FALSE)
  old_oracle <- Sys.getenv("NP_LP_TREE_ORACLE", unset = NA_character_)
  on.exit({
    options(old_opts)
    if(is.na(old_oracle)) Sys.unsetenv("NP_LP_TREE_ORACLE")
    else Sys.setenv(NP_LP_TREE_ORACLE = old_oracle)
  }, add = TRUE)

  Sys.setenv(NP_LP_TREE_ORACLE = "1")
  tree_message <- capture.output(
    tree_objective <- npRmpi:::.npregbw_eval_only(
      x, y, bw, localize = TRUE
    )$objective,
    type = "message"
  )
  expect_true(any(grepl(
    "NP_LP_TREE_ORACLE .*missing=0 mismatch=0",
    tree_message
  )))

  Sys.unsetenv("NP_LP_TREE_ORACLE")
  options(np.tree = "auto")
  auto_objective <- npRmpi:::.npregbw_eval_only(
    x, y, bw, localize = TRUE
  )$objective
  expect_identical(auto_objective, tree_objective)

  options(np.tree = FALSE)
  dense_objective <- npRmpi:::.npregbw_eval_only(
    x, y, bw, localize = TRUE
  )$objective
  expect_equal(tree_objective, dense_objective, tolerance = 1e-12)
})

test_that("fixed sparse width-six LP agrees with the independent hat oracle", {
  skip_if_not(spawn_mpi_slaves(1), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260812)
  n <- 128L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(4 * pi * x$x1) * cos(2 * pi * x$x2) + rnorm(n, sd = 0.1)
  evaluate <- getFromNamespace(".npregbw_eval_only", "npRmpi")
  old_opts <- options(np.tree = TRUE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  for (bernstein in c(FALSE, TRUE)) {
    bw <- npRmpi::npregbw(
      xdat = x,
      ydat = y,
      bws = c(0.22, 0.19),
      regtype = "lp",
      degree = c(2L, 3L),
      basis = "additive",
      bernstein.basis = bernstein,
      ckertype = "epanechnikov",
      bwmethod = "cv.ls",
      bandwidth.compute = FALSE
    )

    native <- evaluate(x, y, bw, localize = TRUE)$objective
    hat <- suppressWarnings(npRmpi::npreghat(
      bws = bw, txdat = x, output = "matrix"
    ))
    fitted <- drop(hat %*% y)
    reference <- mean(((y - fitted) / (1.0 - diag(hat)))^2)

    expect_equal(native, reference, tolerance = 2e-10,
                 info = paste("Bernstein", bernstein))
  }
})
