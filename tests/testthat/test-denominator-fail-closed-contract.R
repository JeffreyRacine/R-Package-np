denominator_fail_closed_run_local <- function(package, expression) {
  if (identical(package, "npRmpi"))
    getFromNamespace(".npRmpi_with_local_regression", package)(expression)
  else
    expression
}

test_that("regression fits and hats fail closed for exact zero mass", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  denominator_fail_closed_run_local(package, {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(-0.2, -0.1, 0, 0.1, 0.2))
  y <- as.double(seq_len(nrow(x)))
  ex <- data.frame(x = c(0, 10))
  bw <- npregbw(
    xdat = x, ydat = y, bws = 0.05,
    bandwidth.compute = FALSE, bwmethod = "cv.ls",
    regtype = "lc", ckertype = "epanechnikov"
  )

  fit <- npreg(bws = bw, exdat = ex, se = TRUE, gradients = TRUE)
  hat <- npreghat(bws = bw, exdat = ex)

  expect_equal(fit$mean[[1L]], 3, tolerance = 1e-15)
  expect_equal(rowSums(hat)[[1L]], 1, tolerance = 1e-15)
  expect_true(is.na(fit$mean[[2L]]))
  expect_true(is.na(fit$merr[[2L]]))
  expect_true(all(is.na(fit$grad[2L, ])))
  expect_true(all(is.na(fit$gerr[2L, ])))
  expect_true(all(is.na(hat[2L, ])))
  })
})

test_that("conditional fits and hats inherit row-local zero-mass missingness", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  denominator_fail_closed_run_local(package, {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(-0.2, -0.1, 0, 0.1, 0.2))
  y <- data.frame(y = c(-0.2, -0.1, 0, 0.1, 0.2))
  ex <- data.frame(x = c(0, 10))
  ey <- data.frame(y = c(0, 0))

  dbw <- npcdensbw(
    xdat = x, ydat = y, bws = c(0.1, 0.05),
    bandwidth.compute = FALSE, bwmethod = "cv.ml",
    cxkertype = "epanechnikov", cykertype = "epanechnikov"
  )
  dens <- npcdens(
    bws = dbw, txdat = x, tydat = y, exdat = ex, eydat = ey,
    gradients = TRUE
  )
  dhat <- npcdenshat(
    bws = dbw, txdat = x, tydat = y, exdat = ex, eydat = ey
  )

  expect_true(is.finite(dens$condens[[1L]]))
  expect_true(is.na(dens$condens[[2L]]))
  expect_true(is.na(dens$conderr[[2L]]))
  expect_true(all(is.na(dens$congrad[2L, ])))
  expect_true(all(is.na(dens$congerr[2L, ])))
  expect_true(all(is.na(dhat[2L, ])))

  cbw <- npcdistbw(
    xdat = x, ydat = y, bws = c(0.1, 0.05),
    bandwidth.compute = FALSE, bwmethod = "cv.ls",
    cxkertype = "epanechnikov", cykertype = "epanechnikov"
  )
  dist <- npcdist(
    bws = cbw, txdat = x, tydat = y, exdat = ex, eydat = ey,
    gradients = FALSE
  )
  chat <- npcdisthat(
    bws = cbw, txdat = x, tydat = y, exdat = ex, eydat = ey
  )

  expect_true(is.finite(dist$condist[[1L]]))
  expect_true(is.na(dist$condist[[2L]]))
  expect_true(is.na(dist$conderr[[2L]]))
  expect_true(all(is.na(chat[2L, ])))
  })
})

test_that("positive-small required mass remains valid and unperturbed", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  denominator_fail_closed_run_local(package, {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = factor(rep("a", 5L), levels = c("a", "b")))
  ex <- data.frame(x = factor("b", levels = c("a", "b")))
  response <- as.double(seq_len(5L))
  tiny <- .Machine$double.xmin

  rbw <- npregbw(
    xdat = x, ydat = response, bws = tiny,
    bandwidth.compute = FALSE, bwtype = "fixed", bwscaling = FALSE,
    ukertype = "aitchisonaitken", bwmethod = "cv.ls", regtype = "lc"
  )
  expect_identical(
    npreg(bws = rbw, exdat = ex, se = FALSE, gradients = FALSE)$mean,
    3
  )

  y <- data.frame(y = c(-0.2, -0.1, 0, 0.1, 0.2))
  cbw <- npcdensbw(
    xdat = x, ydat = y, bws = c(0.1, tiny),
    bandwidth.compute = FALSE, bwtype = "fixed", bwscaling = FALSE,
    uxkertype = "aitchisonaitken", cykertype = "epanechnikov",
    bwmethod = "cv.ml", regtype = "lc"
  )
  fit <- npcdens(
    bws = cbw, txdat = x, tydat = y, exdat = ex,
    eydat = data.frame(y = 0), gradients = FALSE
  )
  hat <- npcdenshat(
    bws = cbw, txdat = x, tydat = y, exdat = ex,
    eydat = data.frame(y = 0)
  )

  expect_identical(fit$condens, as.vector(rowSums(hat)))
  })
})

test_that("zero numerators remain legitimate when required mass is positive", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  denominator_fail_closed_run_local(package, {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(-0.2, -0.1, 0, 0.1, 0.2))
  y <- data.frame(y = c(-0.2, -0.1, 0, 0.1, 0.2))
  bw <- npcdensbw(
    xdat = x, ydat = y, bws = c(0.05, 0.05),
    bandwidth.compute = FALSE, bwmethod = "cv.ml",
    cxkertype = "epanechnikov", cykertype = "epanechnikov"
  )
  fit <- npcdens(
    bws = bw, txdat = x, tydat = y,
    exdat = data.frame(x = 0), eydat = data.frame(y = 10),
    gradients = FALSE
  )

  expect_identical(fit$condens, 0)
  })
})
