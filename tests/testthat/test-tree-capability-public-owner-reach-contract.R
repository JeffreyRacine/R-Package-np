test_that("fixed density and distribution objectives preserve tree parity", {
  old_options <- options(np.messages = FALSE, np.largeh = FALSE,
                         np.macMseries.accelerate = FALSE)
  on.exit(options(old_options), add = TRUE)

  set.seed(20260818L)
  n <- 43L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y = runif(n))
  x_range <- vapply(x, function(value) diff(range(value)), numeric(1L))
  y_range <- diff(range(y$y))
  x_bw <- 1.25 * x_range / sqrt(5)
  y_bw <- 1.25 * y_range / sqrt(5)

  udens_ml <- npudensbw(
    dat = x, bws = x_bw, bandwidth.compute = FALSE,
    bwmethod = "cv.ml", ckertype = "epanechnikov"
  )
  udens_ls <- npudensbw(
    dat = x, bws = x_bw, bandwidth.compute = FALSE,
    bwmethod = "cv.ls", ckertype = "epanechnikov"
  )
  udist <- npudistbw(
    dat = x, bws = x_bw, bandwidth.compute = FALSE,
    bwmethod = "cv.cdf", ckertype = "epanechnikov"
  )
  cdens_ml <- npcdensbw(
    xdat = x, ydat = y, bws = c(y_bw, x_bw),
    bandwidth.compute = FALSE, bwmethod = "cv.ml", regtype = "lc",
    cxkertype = "epanechnikov", cykertype = "epanechnikov"
  )

  evaluate <- function(tree) {
    options(np.tree = tree)
    list(
      udens_ml = npudensbw(
        dat = x, bws = udens_ml, bandwidth.compute = TRUE,
        eval.only = TRUE
      )$fval,
      udens_ls = npudensbw(
        dat = x, bws = udens_ls, bandwidth.compute = TRUE,
        eval.only = TRUE
      )$fval,
      udist = npudistbw(
        dat = x, bws = udist, bandwidth.compute = TRUE,
        eval.only = TRUE, do.full.integral = TRUE
      )$fval,
      cdens_ml = np:::.npcdensbw_eval_only(x, y, cdens_ml)$objective
    )
  }

  dense <- evaluate(FALSE)
  tree <- evaluate(TRUE)
  for (field in names(dense))
    expect_equal(tree[[field]], dense[[field]], tolerance = 2e-11,
                 info = field)
})

test_that("fixed estimator descendants preserve prepared-order scatter", {
  old_options <- options(np.messages = FALSE, np.largeh = FALSE,
                         np.macMseries.accelerate = FALSE)
  on.exit(options(old_options), add = TRUE)

  set.seed(20260819L)
  n <- 41L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y = 0.4 + x$x1 - 0.6 * x$x2 + rnorm(n, sd = 0.2))
  x_range <- vapply(x, function(value) diff(range(value)), numeric(1L))
  y_range <- diff(range(y$y))
  x_bw <- 1.25 * x_range / sqrt(5)
  y_bw <- 1.25 * y_range / sqrt(5)

  udens_bw <- npudensbw(
    dat = x, bws = x_bw, bandwidth.compute = FALSE,
    ckertype = "epanechnikov"
  )
  udist_bw <- npudistbw(
    dat = x, bws = x_bw, bandwidth.compute = FALSE,
    ckertype = "epanechnikov"
  )
  cdens_bw <- npcdensbw(
    xdat = x, ydat = y, bws = c(y_bw, x_bw),
    bandwidth.compute = FALSE, regtype = "lc",
    cxkertype = "epanechnikov", cykertype = "epanechnikov"
  )
  cdist_bw <- npcdistbw(
    xdat = x, ydat = y, bws = c(y_bw, x_bw),
    bandwidth.compute = FALSE, regtype = "lc",
    cxkertype = "epanechnikov", cykertype = "epanechnikov"
  )
  reg_lc_bw <- npregbw(
    xdat = x, ydat = y$y, bws = x_bw, bandwidth.compute = FALSE,
    regtype = "lc", ckertype = "epanechnikov"
  )
  reg_lp_bw <- npregbw(
    xdat = x, ydat = y$y, bws = x_bw, bandwidth.compute = FALSE,
    regtype = "lp", degree = c(1L, 1L), ckertype = "epanechnikov"
  )

  evaluate <- function(tree) {
    options(np.tree = tree)
    list(
      udens = npudens(bws = udens_bw, tdat = x),
      udist = npudist(bws = udist_bw, tdat = x),
      cdens = npcdens(
        bws = cdens_bw, txdat = x, tydat = y,
        gradients = TRUE
      ),
      cdist = npcdist(
        bws = cdist_bw, txdat = x, tydat = y,
        gradients = TRUE
      ),
      reg_lc = npreg(bws = reg_lc_bw, txdat = x, tydat = y$y),
      reg_lp = npreg(
        bws = reg_lp_bw, txdat = x, tydat = y$y,
        se = TRUE, gradients = TRUE
      ),
      reghat = unclass(npreghat(bws = reg_lc_bw, txdat = x))
    )
  }

  dense <- evaluate(FALSE)
  tree <- evaluate(TRUE)
  fields <- list(
    udens = c("dens", "derr"),
    udist = c("dist", "derr"),
    cdens = c("condens", "conderr", "congrad", "congerr"),
    cdist = c("condist", "conderr", "congrad", "congerr"),
    reg_lc = "mean",
    reg_lp = c("mean", "merr", "grad", "gerr")
  )

  for (owner in names(fields))
    for (field in fields[[owner]])
      expect_equal(
        tree[[owner]][[field]], dense[[owner]][[field]],
        tolerance = 2e-11, info = paste(owner, field)
      )
  expect_equal(tree$reghat, dense$reghat, tolerance = 2e-11)
})
