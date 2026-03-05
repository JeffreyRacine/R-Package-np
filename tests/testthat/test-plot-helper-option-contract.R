test_that("npRmpi helper constructors forward kernel options and normalize bwscaling", {
  skip_if_not_installed("np")

  set.seed(7101)
  n <- 50
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(y = runif(n))

  bw <- np::npcdensbw(
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
  )

  np.ns <- asNamespace("npRmpi")
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
  on.exit(untrace("kbandwidth.numeric", where = np.ns), add = TRUE)

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
})
