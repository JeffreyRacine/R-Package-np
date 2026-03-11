test_that("npRmpi helper constructors forward kernel options and normalize bwscaling", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(7101)
  n <- 50
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(y = runif(n))

  np.ns <- asNamespace("npRmpi")
  for (builder in c("npcdensbw", "npcdistbw")) {
    bw <- do.call(builder, list(
      xdat = xdat,
      ydat = ydat,
      bws = c(0.22, 0.22),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      bwscaling = TRUE,
      cxkertype = "epanechnikov",
      cykertype = "epanechnikov",
      cxkerbound = "fixed",
      cykerbound = "fixed",
      cxkerlb = 0.0,
      cykerlb = 0.0,
      cxkerub = 1.0,
      cykerub = 1.0
    ))

    cap <- new.env(parent = emptyenv())
    cap$calls <- list()

    trace(
      what = "kbandwidth.numeric",
      where = np.ns,
      tracer = bquote({
        assign(
          "calls",
          c(
            get("calls", envir = .(cap)),
            list(list(
              bwscaling = bwscaling,
              ckertype = ckertype,
              ckerorder = ckerorder,
              ckerbound = ckerbound
            ))
          ),
          envir = .(cap)
        )
      }),
      print = FALSE
    )
    on.exit(try(untrace("kbandwidth.numeric", where = np.ns), silent = TRUE), add = TRUE)

    make.kx <- getFromNamespace(".np_con_make_kbandwidth_x", "npRmpi")
    make.kxy <- getFromNamespace(".np_con_make_kbandwidth_xy", "npRmpi")

    obj1 <- make.kx(bws = bw, xdat = xdat)
    obj2 <- make.kxy(bws = bw, xdat = xdat, ydat = ydat)

    expect_false(is.null(obj1))
    expect_false(is.null(obj2))
    expect_true(length(cap$calls) >= 2L)

    for (call in cap$calls) {
      expect_identical(isTRUE(call$bwscaling), FALSE)
      expect_identical(as.character(call$ckertype), as.character(bw$cxkertype))
      expect_identical(as.character(call$ckerorder), as.character(bw$cxkerorder))
      expect_identical(as.character(call$ckerbound), as.character(bw$cxkerbound))
    }
  }
})

test_that("npRmpi semihat regbw args forward index LP/kernel options with bound collapse", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(7102)
  n <- 64
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- rnorm(n)

  bw <- npindexbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.2, 0.2, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = c(0, 0),
    ckerub = c(1, 1)
  )
  idx.train <- data.frame(index = as.vector(as.matrix(xdat) %*% bw$beta))

  make.args <- getFromNamespace(".np_semihat_make_regbw_args", "npRmpi")
  args <- make.args(
    source = bw,
    xdat = idx.train,
    ydat = rep.int(0.0, nrow(idx.train)),
    bw = bw$bw
  )

  expect_identical(args$regtype, as.character(bw$regtype))
  expect_identical(args$basis, bw$basis)
  expect_equal(args$degree, bw$degree)
  expect_identical(isTRUE(args$bernstein.basis), isTRUE(bw$bernstein.basis))
  expect_identical(args$bwtype, bw$type)
  expect_identical(args$ckertype, bw$ckertype)
  expect_identical(args$ckerorder, bw$ckerorder)
  expect_identical(args$ckerbound, bw$ckerbound)
  expect_equal(as.double(args$ckerlb), 0)
  expect_equal(as.double(args$ckerub), 1)
  expect_identical(args$bandwidth.compute, FALSE)
})

test_that("npRmpi semihat regbw args forward smooth-coef LP/kernel/scaling options", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(7103)
  n <- 58
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  ydat <- rnorm(n)

  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = 0.25,
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwscaling = TRUE,
    ckertype = "epanechnikov",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  make.args <- getFromNamespace(".np_semihat_make_regbw_args", "npRmpi")
  args <- make.args(
    source = bw,
    xdat = zdat,
    ydat = rep.int(0.0, nrow(zdat)),
    bw = bw$bw
  )

  expect_identical(args$regtype, as.character(bw$regtype))
  expect_identical(args$basis, bw$basis)
  expect_equal(args$degree, bw$degree)
  expect_identical(isTRUE(args$bernstein.basis), isTRUE(bw$bernstein.basis))
  expect_identical(args$bwscaling, isTRUE(bw$scaling))
  expect_identical(args$bwtype, bw$type)
  expect_identical(args$ckertype, bw$ckertype)
  expect_identical(args$ckerorder, bw$ckerorder)
  expect_identical(args$ckerbound, bw$ckerbound)
  expect_equal(as.double(args$ckerlb), 0)
  expect_equal(as.double(args$ckerub), 1)
  expect_identical(args$bandwidth.compute, FALSE)

  if (!is.null(bw$method) && bw$method %in% c("cv.ls", "cv.aic")) {
    expect_identical(args$bwmethod, bw$method)
  }
})
