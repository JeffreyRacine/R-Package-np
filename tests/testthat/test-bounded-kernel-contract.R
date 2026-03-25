test_that("fixed +/-Inf bounds are parity-equivalent to none at fixed bandwidth", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(80)
  dat <- data.frame(x = x)

  bw.none <- npudensbw(
    dat = dat,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "none"
  )

  bw.inf <- npudensbw(
    dat = dat,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "fixed",
    ckerlb = -Inf,
    ckerub = Inf
  )

  fit.none <- npudens(bws = bw.none, tdat = dat)
  fit.inf <- npudens(bws = bw.inf, tdat = dat)

  expect_equal(as.numeric(fit.none$dens), as.numeric(fit.inf$dens), tolerance = 1e-12)
})

test_that("scalar fixed bounds recycle over multiple continuous variables", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  dat <- data.frame(x1 = runif(64), x2 = runif(64))

  bw.scalar <- npudensbw(
    dat = dat,
    bws = c(0.2, 0.25),
    bandwidth.compute = FALSE,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  bw.vector <- npudensbw(
    dat = dat,
    bws = c(0.2, 0.25),
    bandwidth.compute = FALSE,
    ckerbound = "fixed",
    ckerlb = c(0, 0),
    ckerub = c(1, 1)
  )

  fit.scalar <- npudens(bws = bw.scalar, tdat = dat)
  fit.vector <- npudens(bws = bw.vector, tdat = dat)

  expect_equal(as.numeric(fit.scalar$dens), as.numeric(fit.vector$dens), tolerance = 1e-12)
})

test_that("invalid fixed bounds are rejected with clear diagnostics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(50)
  dat <- data.frame(x = x)

  expect_error(
    npudensbw(
      dat = dat,
      bws = 0.2,
      bandwidth.compute = FALSE,
      ckerbound = "fixed",
      ckerlb = 1,
      ckerub = 1
    ),
    "lower < upper"
  )

  expect_error(
    npudensbw(
      dat = dat,
      bws = 0.2,
      bandwidth.compute = FALSE,
      ckerbound = "fixed",
      ckerlb = 0.2,
      ckerub = 0.8
    ),
    "Violations: x"
  )
})

test_that("bounded generalized_nn is available for certified core public routes", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(36)
  y <- cos(2 * pi * x) + rnorm(36, sd = 0.05)
  xy <- data.frame(x = x)
  yy <- data.frame(y = runif(36))

  bw.ud.cvml <- npudensbw(
    dat = xy,
    bwmethod = "cv.ml",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  bw.ud.cvls <- npudensbw(
    dat = xy,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.ud <- npudens(bws = bw.ud.cvls, tdat = xy)

  bw.dist <- npudistbw(
    dat = xy,
    bwmethod = "cv.cdf",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.dist <- npudist(bws = bw.dist, tdat = xy)

  bw.reg <- npregbw(
    xdat = xy,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.reg <- npreg(bws = bw.reg, txdat = xy, tydat = y)

  bw.cd.cvml <- npcdensbw(
    xdat = xy,
    ydat = yy,
    bwmethod = "cv.ml",
    bwtype = "generalized_nn",
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  bw.cd.cvls <- npcdensbw(
    xdat = xy,
    ydat = yy,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  fit.cd <- npcdens(bws = bw.cd.cvls, txdat = xy, tydat = yy)

  bw.cdist <- npcdistbw(
    xdat = xy,
    ydat = yy,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  fit.cdist <- npcdist(bws = bw.cdist, txdat = xy, tydat = yy)

  expect_true(all(is.finite(as.numeric(bw.ud.cvml$bw))))
  expect_true(is.finite(bw.ud.cvml$fval))
  expect_true(all(is.finite(as.numeric(bw.ud.cvls$bw))))
  expect_true(is.finite(bw.ud.cvls$fval))
  expect_true(all(is.finite(as.numeric(fit.ud$dens))))

  expect_true(all(is.finite(as.numeric(bw.dist$bw))))
  expect_true(is.finite(bw.dist$fval))
  expect_true(all(is.finite(as.numeric(fit.dist$dist))))

  expect_true(all(is.finite(as.numeric(bw.reg$bw))))
  expect_true(is.finite(bw.reg$fval))
  expect_true(all(is.finite(as.numeric(fit.reg$mean))))

  expect_true(all(is.finite(as.numeric(bw.cd.cvml$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cd.cvml$ybw))))
  expect_true(is.finite(bw.cd.cvml$fval))
  expect_true(all(is.finite(as.numeric(bw.cd.cvls$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cd.cvls$ybw))))
  expect_true(is.finite(bw.cd.cvls$fval))
  expect_true(all(is.finite(as.numeric(fit.cd$condens))))

  expect_true(all(is.finite(as.numeric(bw.cdist$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cdist$ybw))))
  expect_true(is.finite(bw.cdist$fval))
  expect_true(all(is.finite(as.numeric(fit.cdist$condist))))
})

test_that("bounded adaptive_nn is available for certified non-distribution routes", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260325)
  x <- runif(36)
  y <- sin(2 * pi * x) + rnorm(36, sd = 0.05)
  xy <- data.frame(x = x)
  yy <- data.frame(y = runif(36))

  bw.ud <- npudensbw(
    dat = xy,
    bwmethod = "cv.ml",
    bwtype = "adaptive_nn",
    ckerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )
  fit.ud <- npudens(bws = bw.ud, tdat = xy)

  bw.reg <- npregbw(
    xdat = xy,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    bwtype = "adaptive_nn",
    ckerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )
  fit.reg <- npreg(bws = bw.reg, txdat = xy, tydat = y, gradients = TRUE)

  bw.cd <- npcdensbw(
    xdat = xy,
    ydat = yy,
    bwmethod = "cv.ml",
    bwtype = "adaptive_nn",
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )
  fit.cd <- npcdens(bws = bw.cd, txdat = xy, tydat = yy)

  expect_true(all(is.finite(as.numeric(bw.ud$bw))))
  expect_true(is.finite(bw.ud$fval))
  expect_true(all(is.finite(as.numeric(fit.ud$dens))))

  expect_true(all(is.finite(as.numeric(bw.reg$bw))))
  expect_true(is.finite(bw.reg$fval))
  expect_true(all(is.finite(as.numeric(fit.reg$mean))))
  expect_true(all(is.finite(as.numeric(fit.reg$gradients))))

  expect_true(all(is.finite(as.numeric(bw.cd$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cd$ybw))))
  expect_true(is.finite(bw.cd$fval))
  expect_true(all(is.finite(as.numeric(fit.cd$condens))))
  expect_true(all(as.numeric(fit.cd$condens) >= 0))
})

test_that("npplreg, npindex, and npscoef bounded generalized_nn are available", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(32)
  y <- cos(2 * pi * x)
  y_sc <- (1 + sin(2 * pi * x)) * x + rnorm(32, sd = 0.05)
  xy <- data.frame(x = x)

  bw.pl <- npplregbw(
    xdat = xy,
    ydat = y,
    zdat = xy,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.pl <- npplreg(bws = bw.pl, txdat = xy, tydat = y, tzdat = xy)

  bw.index <- npindexbw(
    xdat = data.frame(x1 = x, x2 = x^2),
    ydat = y,
    method = "ichimura",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.index <- npindex(bws = bw.index, txdat = data.frame(x1 = x, x2 = x^2), tydat = y)

  bw.sc <- npscoefbw(
    xdat = xy,
    ydat = y_sc,
    zdat = xy,
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.sc <- npscoef(bws = bw.sc, txdat = xy, tydat = y_sc, tzdat = xy, iterate = FALSE)

  expect_true(all(is.finite(as.numeric(unlist(lapply(bw.pl$bw, function(obj) obj$bw))))))
  expect_true(all(is.finite(as.numeric(fit.pl$mean))))
  expect_true(all(is.finite(as.numeric(bw.index$beta))))
  expect_true(is.finite(as.numeric(bw.index$bw)))
  expect_true(all(is.finite(as.numeric(fit.index$mean))))
  expect_true(all(is.finite(as.numeric(bw.sc$bw))))
  expect_true(all(is.finite(as.numeric(fit.sc$mean))))
})

.bounded_adaptive_udist_fitonly_fixture <- function() {
  dat <- structure(list(x = c(0.0293716583400965, 0.0337911003734916,
    0.0688598649576306, 0.146281852386892, 0.15000450075604, 0.168401539325714,
    0.173830970656127, 0.201029423624277, 0.210880494909361, 0.216087364125997,
    0.261414775624871, 0.332859222311527, 0.340369544224814, 0.398810388986021,
    0.412044096272439, 0.460934564704075, 0.522144361166283, 0.529954593162984,
    0.559345822315663, 0.573270461056381, 0.644053619122133, 0.668163585709408,
    0.686013063415885, 0.689600103069097, 0.764754617586732, 0.80543372547254,
    0.814246771391481, 0.830958560109138, 0.849165873136371, 0.849887294694781,
    0.936625305796042, 0.97392991813831)), class = "data.frame", row.names = c(NA,
    -32L))
  bw <- structure(list(bw = 1, method = "cv.cdf", pmethod = "Least Squares Cross-Validation",
    fval = 0.183985330074149, ifval = 1, num.feval = 75, fval.history = 0.183985330074149,
    eval.history = 75, invalid.history = 7, scaling = FALSE,
    pscaling = "Bandwidth Nearest Neighbor(s)", type = "adaptive_nn",
    ptype = "Adaptive Nearest Neighbour", ckertype = "gaussian",
    ckerorder = 2, ckerbound = "range", ckerlb = 0.0293716583400965,
    ckerub = 0.97392991813831, pckertype = "Second-Order Gaussian (bounded/range)",
    ukertype = "aitchisonaitken", pukertype = "Aitchison and Aitken",
    okertype = "liracine", pokertype = "Li and Racine (normalized)",
    nobs = 32L, ndim = 1L, ncon = 1L, nuno = 0L, nord = 0L, icon = c(x = TRUE),
    iuno = c(x = FALSE), iord = c(x = FALSE), xnames = "x", xdati = list(
      iord = c(x = FALSE), iuno = c(x = FALSE), icon = c(x = TRUE),
      inumord = c(x = FALSE), all.lev = list(x = NULL), all.ulev = list(
        x = NULL), all.dlev = list(x = NULL), all.nlev = list(
        x = NULL), all.min = list(x = 0.0293716583400965),
      all.max = list(x = 0.97392991813831)), sfactor = list(
      x = 10.9247445704317), bandwidth = list(x = 1), nconfac = 0.314980262473718,
    ncatfac = 0.0992125657480125, sdev = 0.290606529376362, sumNum = list(
      x = 10.9247445704317), xmcv = structure(numeric(0), dim = c(0L,
    0L), num.row = 0L, pad.num = 0.565213999172708), dati = list(
      x = list(iord = c(x = FALSE), iuno = c(x = FALSE), icon = c(x = TRUE),
        inumord = c(x = FALSE), all.lev = list(x = NULL),
        all.ulev = list(x = NULL), all.dlev = list(x = NULL),
        all.nlev = list(x = NULL), all.min = list(x = 0.0293716583400965),
        all.max = list(x = 0.97392991813831))), varnames = list(
      x = "x"), vartitle = list(x = ""), vartitleabb = list(
      x = ""), rows.omit = NA, nobs.omit = 0, timing = 5.4e-05,
    total.time = c(elapsed = 0.00700000000000001), klist = list(
      x = list(ckertype = "gaussian", pckertype = "Second-Order Gaussian (bounded/range)",
        ukertype = "aitchisonaitken", pukertype = "Aitchison and Aitken",
        okertype = "liracine", pokertype = "Li and Racine (normalized)")),
    call = npudistbw.NULL(dat = dat, ... = pairlist(bwtype = "adaptive_nn", bwmethod = "cv.cdf", ckerbound = "range", nmulti = 1))), class = "dbandwidth")
  list(dat = dat, bw = bw)
}

.bounded_adaptive_cdist_fitonly_fixture <- function() {
  x <- structure(list(x = c(0.0293716583400965, 0.0337911003734916,
    0.146281852386892, 0.15000450075604, 0.168401539325714, 0.173830970656127,
    0.201029423624277, 0.210880494909361, 0.261414775624871, 0.332859222311527,
    0.340369544224814, 0.398810388986021, 0.412044096272439, 0.460934564704075,
    0.522144361166283, 0.529954593162984, 0.573270461056381, 0.644053619122133,
    0.668163585709408, 0.686013063415885, 0.689600103069097, 0.764754617586732,
    0.80543372547254, 0.814246771391481, 0.830958560109138, 0.849887294694781,
    0.936625305796042, 0.97392991813831)), class = "data.frame", row.names = c(NA,
    -28L))
  y <- structure(list(y = c(0.0647218793164939, 0.0688598649576306,
    0.0753922816365957, 0.116926148533821, 0.141229883302003, 0.146168121369556,
    0.185237535042688, 0.216087364125997, 0.237273563398048, 0.241650115931407,
    0.280670147156343, 0.282998539274558, 0.41492543485947, 0.444256805581972,
    0.452224941225722, 0.501736348960549, 0.506873786216602, 0.554796958575025,
    0.559345822315663, 0.603990532224998, 0.627316205063835, 0.627725793980062,
    0.714035893790424, 0.715121431509033, 0.849165873136371, 0.861479766201228,
    0.869792656507343, 0.889108955627307)), class = "data.frame", row.names = c(NA,
    -28L))
  bw <- structure(list(xbw = 7, ybw = 1, method = "cv.ls",
    pmethod = "Least Squares Cross-Validation", fval = 0.0556322703376515,
    ifval = 1, num.feval = 179, num.feval.fast = 0, fval.history = 0.0556322703376515,
    eval.history = 179, invalid.history = 10, scaling = FALSE,
    pscaling = "Bandwidth Nearest Neighbor(s)", type = "adaptive_nn",
    ptype = "Adaptive Nearest Neighbour", cxkertype = "gaussian",
    cykertype = "gaussian", cxkerorder = 2, cykerorder = 2, cxkerbound = "range",
    cxkerlb = 0.0293716583400965, cxkerub = 0.97392991813831,
    cykerbound = "range", cykerlb = 0.0647218793164939, cykerub = 0.889108955627307,
    pcxkertype = "Second-Order Gaussian (bounded/range)", pcykertype = "Second-Order Gaussian (bounded/range)",
    uxkertype = "aitchisonaitken", uykertype = "aitchisonaitken",
    puxkertype = "Aitchison and Aitken", puykertype = "Aitchison and Aitken",
    oxkertype = "liracine", oykertype = "liracine", poxkertype = "Li and Racine",
    poykertype = "Li and Racine (normalized)", nobs = 28L, xndim = 1L,
    yndim = 1L, ndim = 2L, xncon = 1L, xnuno = 0L, xnord = 0L,
    yncon = 1L, ynuno = 0L, ynord = 0L, ncon = 2L, ixcon = c(x = TRUE),
    ixuno = c(x = FALSE), ixord = c(x = FALSE), iycon = c(y = TRUE),
    iyuno = c(y = FALSE), iyord = c(y = FALSE), xnames = "x",
    ynames = "y", xdati = list(iord = c(x = FALSE), iuno = c(x = FALSE),
      icon = c(x = TRUE), inumord = c(x = FALSE), all.lev = list(
        x = NULL), all.ulev = list(x = NULL), all.dlev = list(
        x = NULL), all.nlev = list(x = NULL), all.min = list(
        x = 0.0293716583400965), all.max = list(x = 0.97392991813831)),
    ydati = list(iord = c(y = FALSE), iuno = c(y = FALSE), icon = c(y = TRUE),
      inumord = c(y = FALSE), all.lev = list(y = NULL), all.ulev = list(
        y = NULL), all.dlev = list(y = NULL), all.nlev = list(
        y = NULL), all.min = list(y = 0.0647218793164939),
      all.max = list(y = 0.889108955627307)), xmcv = structure(numeric(0), dim = c(0L,
    0L), num.row = 0L, pad.num = -0.577228845077323), ymcv = structure(numeric(0), dim = c(0L,
    0L), num.row = 0L, pad.num = 0.302938417806345), sfactor = list(
      x = 42.3889978414751, y = 6.47952560003134), bandwidth = list(
      x = 7, y = 1), nconfac = 0.573861375250308, ncatfac = 0.329316878004175,
    sdev = c(xcon = 0.287764950469254, ycon = 0.268936528850516),
    sumNum = list(x = 42.3889978414751, y = 6.47952560003134),
    dati = list(x = list(iord = c(x = FALSE), iuno = c(x = FALSE),
      icon = c(x = TRUE), inumord = c(x = FALSE), all.lev = list(
        x = NULL), all.ulev = list(x = NULL), all.dlev = list(
        x = NULL), all.nlev = list(x = NULL), all.min = list(
        x = 0.0293716583400965), all.max = list(x = 0.97392991813831)),
      y = list(iord = c(y = FALSE), iuno = c(y = FALSE), icon = c(y = TRUE),
        inumord = c(y = FALSE), all.lev = list(y = NULL),
        all.ulev = list(y = NULL), all.dlev = list(y = NULL),
        all.nlev = list(y = NULL), all.min = list(y = 0.0647218793164939),
        all.max = list(y = 0.889108955627307))), varnames = list(
      x = "x", y = "y"), vartitle = list(x = "Explanatory",
      y = "Dependent"), vartitleabb = list(x = "Exp.", y = "Dep."),
    rows.omit = NA, nobs.omit = 0, timing = 8.2e-05, total.time = c(elapsed = 0.017),
    regtype = "lc", pregtype = "Local-Constant", basis = "glp",
    degree = 0L, bernstein.basis = FALSE, regtype.engine = "lc",
    basis.engine = "glp", degree.engine = 0L, bernstein.basis.engine = FALSE,
    klist = list(x = list(ckertype = "gaussian", ckerbound = "range",
      ckerlb = 0.0293716583400965, ckerub = 0.97392991813831,
      pckertype = "Second-Order Gaussian (bounded/range)",
      ukertype = "aitchisonaitken", pukertype = "Aitchison and Aitken",
      okertype = "liracine", pokertype = "Li and Racine"),
      y = list(ckertype = "gaussian", ckerbound = "range",
        ckerlb = 0.0647218793164939, ckerub = 0.889108955627307,
        pckertype = "Second-Order Gaussian (bounded/range)",
        ukertype = "aitchisonaitken", pukertype = "Aitchison and Aitken",
        okertype = "liracine", pokertype = "Li and Racine (normalized)")),
    call = npcdistbw.NULL(xdat = x, ydat = y, ... = pairlist(bwtype = "adaptive_nn", bwmethod = "cv.ls", cxkerbound = "range", cykerbound = "range", nmulti = 1))), class = "condbandwidth")
  list(x = x, y = y, bw = bw)
}

test_that("bounded adaptive_nn fit-only remains available for distribution routes", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  ud.fixture <- .bounded_adaptive_udist_fitonly_fixture()
  fit.ud <- npudist(bws = ud.fixture$bw, tdat = ud.fixture$dat)

  cd.fixture <- .bounded_adaptive_cdist_fitonly_fixture()
  fit.cdist <- npcdist(bws = cd.fixture$bw, txdat = cd.fixture$x, tydat = cd.fixture$y)

  expect_true(all(is.finite(as.numeric(fit.ud$dist))))
  expect_true(all(is.finite(as.numeric(fit.cdist$condist))))
})

test_that("deferred bounded public families remain blocked", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(32)
  y <- cos(2 * pi * x)
  xy <- data.frame(x = x)
  yy <- data.frame(y = y)

  expect_error(
    npudistbw(
      dat = xy,
      bwmethod = "cv.cdf",
      bwtype = "adaptive_nn",
      ckerbound = "range",
      nmulti = 1
    ),
    "bounded adaptive_nn remains unsupported for npudistbw\\(\\) in npRmpi"
  )

  expect_error(
    npcdistbw(
      xdat = xy,
      ydat = yy,
      bwmethod = "cv.ls",
      bwtype = "adaptive_nn",
      cxkerbound = "range",
      cykerbound = "range",
      nmulti = 1
    ),
    "bounded adaptive_nn remains unsupported for npcdistbw\\(\\) in npRmpi"
  )

  expect_error(
    npplregbw(
      xdat = xy,
      ydat = y,
      zdat = xy,
      bwmethod = "cv.ls",
      bwtype = "generalized_nn",
      ckerbound = "range",
      nmulti = 1
    ),
    "finite continuous kernel bounds require bwtype = \"fixed\""
  )
})

test_that("evaluation support violations are caught before native execution", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(70)
  dat <- data.frame(x = x)

  bw <- npudensbw(
    dat = dat,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "range"
  )

  expect_error(
    npudens(
      bws = bw,
      tdat = dat,
      edat = data.frame(x = max(x) + 0.05)
    ),
    "x >"
  )
})

test_that("predict paths enforce bounded eval checks with variable diagnostics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(80)
  y <- runif(80)
  dat.x <- data.frame(x = x)
  dat.y <- data.frame(y = y)

  bw.den <- npudensbw(
    dat = dat.x,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "range"
  )
  fit.den <- npudens(bws = bw.den, tdat = dat.x)
  expect_error(
    predict(fit.den, edat = data.frame(x = max(x) + 0.01)),
    "Evaluation data violate 'ckerbound' bounds: x >"
  )

  bw.dist <- npudistbw(
    dat = dat.x,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "range"
  )
  fit.dist <- npudist(bws = bw.dist, tdat = dat.x)
  expect_error(
    predict(fit.dist, edat = data.frame(x = max(x) + 0.01)),
    "Evaluation data violate 'ckerbound' bounds: x >"
  )

  bw.reg <- npregbw(
    xdat = dat.x,
    ydat = y,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "range"
  )
  fit.reg <- npreg(bws = bw.reg, txdat = dat.x, tydat = y)
  expect_error(
    predict(fit.reg, exdat = data.frame(x = max(x) + 0.01)),
    "Evaluation data violate 'ckerbound' bounds: x >"
  )

  bw.cd <- npcdensbw(
    xdat = dat.x,
    ydat = dat.y,
    bws = c(0.2, 0.2),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )
  fit.cd <- npcdens(bws = bw.cd, txdat = dat.x, tydat = dat.y)
  expect_error(
    predict(fit.cd, exdat = data.frame(x = max(x) + 0.01), eydat = data.frame(y = mean(y))),
    "Evaluation data violate 'cxkerbound' bounds: x >"
  )
  expect_error(
    predict(fit.cd, exdat = data.frame(x = mean(x)), eydat = data.frame(y = max(y) + 0.01)),
    "Evaluation data violate 'cykerbound' bounds: y >"
  )

  bw.cdist <- npcdistbw(
    xdat = dat.x,
    ydat = dat.y,
    bws = c(0.2, 0.2),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )
  fit.cdist <- npcdist(bws = bw.cdist, txdat = dat.x, tydat = dat.y)
  expect_error(
    predict(fit.cdist, exdat = data.frame(x = max(x) + 0.01), eydat = data.frame(y = mean(y))),
    "Evaluation data violate 'cxkerbound' bounds: x >"
  )
  expect_error(
    predict(fit.cdist, exdat = data.frame(x = mean(x)), eydat = data.frame(y = max(y) + 0.01)),
    "Evaluation data violate 'cykerbound' bounds: y >"
  )
})
