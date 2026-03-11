fixed_cvls_objective <- function(xdat,
                                 ydat,
                                 regtype,
                                 degree,
                                 bws,
                                 basis = "glp",
                                 bernstein.basis = FALSE) {
  ns <- asNamespace("npRmpi")
  xdat <- as.data.frame(xdat)
  bw0 <- npRmpi::npregbw(
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

  xmat <- get("toMatrix", ns)(xdat)
  runo <- xmat[, bw0$iuno, drop = FALSE]
  rcon <- xmat[, bw0$icon, drop = FALSE]
  rord <- xmat[, bw0$iord, drop = FALSE]
  mysd <- get("EssDee", ns)(rcon)
  n <- nrow(xdat)
  nconfac <- n^(-1.0 / (2.0 * bw0$ckerorder + bw0$ncon))
  ncatfac <- n^(-2.0 / (2.0 * bw0$ckerorder + bw0$ncon))
  reg.c <- get("npRegtypeToC", ns)(
    regtype = bw0$regtype,
    degree = bw0$degree,
    ncon = bw0$ncon,
    context = "npregbw"
  )

  myopti <- as.integer(list(
    n,
    get("IMULTI_FALSE", ns),
    1L,
    get("USE_START_YES", ns),
    if (bw0$scaling) get("SF_NORMAL", ns) else get("SF_ARB", ns),
    get("BW_FIXED", ns),
    0L,
    get("RE_MIN_FALSE", ns),
    get("IO_MIN_TRUE", ns),
    get("BWM_CVLS", ns),
    get("CKER_GAUSS", ns) + bw0$ckerorder / 2 - 1,
    get("UKER_AIT", ns),
    get("OKER_WANG", ns),
    bw0$nuno,
    bw0$nord,
    bw0$ncon,
    reg.c$code,
    get("DO_TREE_NO", ns),
    FALSE,
    3L,
    FALSE
  ))
  myoptd <- as.double(list(
    1.490116e-07,
    1.490116e-04,
    1.490116e-05,
    0.5,
    2.5 * (3.0 - sqrt(5)),
    1.0,
    0.1,
    1,
    0.25 * (3.0 - sqrt(5)),
    1.0,
    0.1,
    2.0,
    0.5,
    0.1,
    0.9,
    0.375,
    nconfac,
    ncatfac
  ))
  ck <- get("npKernelBoundsMarshal", ns)(bw0$ckerlb[bw0$icon], bw0$ckerub[bw0$icon])

  out <- .Call(
    "C_np_regression_bw_eval",
    as.double(runo),
    as.double(rord),
    as.double(rcon),
    as.double(ydat),
    as.double(mysd),
    myopti,
    myoptd,
    as.double(c(bw0$bw[bw0$icon], bw0$bw[bw0$iuno], bw0$bw[bw0$iord])),
    as.integer(1L),
    as.integer(1L),
    as.double(10),
    as.integer(if (bw0$ncon > 0) bw0$degree else integer(1)),
    as.integer(isTRUE(bw0$bernstein.basis)),
    as.integer(get("npLpBasisCode", ns)(bw0$basis)),
    as.double(ck$lb),
    as.double(ck$ub),
    PACKAGE = "npRmpi"
  )

  unname(out$fval[1])
}

test_that("fixed cv.ls is exact on in-class local polynomial fixtures", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
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
})

test_that("fixed cv.ls keeps ll and canonical lp degree-1 aligned off-model", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
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
