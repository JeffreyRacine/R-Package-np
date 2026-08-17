adaptive_conditional_cvml_gaussian2 <- function(u) stats::dnorm(u)

adaptive_conditional_cvml_epanechnikov2 <- function(u) {
  edge <- sqrt(5)
  value <- numeric(length(u))
  inside <- abs(u) < edge
  value[inside] <- 3 / (4 * edge) * (1 - u[inside]^2 / 5)
  value
}

adaptive_conditional_cvml_fold_radius <- function(value, held_out, donor, k) {
  distance <- sort(abs(value[-c(held_out, donor)] - value[[donor]]),
                   method = "quick")
  lookup <- min(k, length(distance))
  distance[[lookup]] * k / lookup
}

adaptive_conditional_cvml_design <- function(delta, degree) {
  if (all(degree == 0L))
    return(matrix(1, nrow = nrow(delta), ncol = 1L))
  if (ncol(delta) == 1L)
    return(vapply(0:degree[[1L]], function(power) delta[, 1L]^power,
                  numeric(nrow(delta))))
  stopifnot(all(degree == 1L))
  cbind(1, delta)
}

adaptive_conditional_cvml_literal <- function(x, y, kx, ky, degree, kernel) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  fit <- numeric(n)

  for (held_out in seq_len(n)) {
    donor <- setdiff(seq_len(n), held_out)
    xweight <- rep.int(1, length(donor))
    for (coordinate in seq_len(ncol(x))) {
      radius <- vapply(donor, function(index) {
        adaptive_conditional_cvml_fold_radius(
          x[, coordinate], held_out, index, kx[[coordinate]])
      }, numeric(1L))
      xweight <- xweight * kernel(
        (x[held_out, coordinate] - x[donor, coordinate]) / radius) /
        radius
    }
    yradius <- vapply(donor, function(index) {
      adaptive_conditional_cvml_fold_radius(y, held_out, index, ky)
    }, numeric(1L))
    ykernel <- kernel((y[[held_out]] - y[donor]) / yradius) / yradius
    design <- adaptive_conditional_cvml_design(
      sweep(x[donor, , drop = FALSE], 2L, x[held_out, ], "-"), degree)
    coefficient <- solve(
      crossprod(design, xweight * design),
      c(1, rep.int(0, ncol(design) - 1L)))
    influence <- drop(xweight * (design %*% coefficient))
    fit[[held_out]] <- sum(influence * ykernel)
  }

  stopifnot(all(is.finite(fit)), all(fit > 0))
  sum(log(fit))
}

adaptive_conditional_cvml_bw <- function(x, y, kx, ky, degree, kernel) {
  x <- as.data.frame(x)
  names(x) <- paste0("x", seq_len(ncol(x)))
  call <- list(
    xdat = x, ydat = data.frame(y = y), bws = as.double(c(ky, kx)),
    bandwidth.compute = FALSE, bwmethod = "cv.ml",
    bwtype = "adaptive_nn", bwscaling = FALSE,
    regtype = if (all(degree == 0L)) "lc" else "lp",
    cxkertype = kernel, cxkerorder = 2L,
    cykertype = kernel, cykerorder = 2L)
  if (!all(degree == 0L)) {
    call$degree <- as.integer(degree)
    call$degree.select <- "manual"
    call$basis <- "glp"
  }
  do.call(npcdensbw, call)
}

adaptive_conditional_cvml_objective <- function(x, y, bw,
                                                invalid.penalty = "baseline") {
  x <- as.data.frame(x)
  names(x) <- paste0("x", seq_len(ncol(x)))
  as.numeric(getFromNamespace(".npRmpi_with_local_regression", "npRmpi")(
    npRmpi:::.npcdensbw_eval_only(
      x, data.frame(y = y), bw,
      invalid.penalty = invalid.penalty)
  )$objective[[1L]])
}

test_that("adaptive conditional CVML uses one exact fold geometry on X and Y", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old_options <- options(np.messages = FALSE, np.largeh = FALSE)
  on.exit(options(old_options), add = TRUE)

  index <- seq_len(17L)
  fixtures <- list(
    lc_p1 = list(
      x = matrix(seq(-1.9, 2.1, length.out = 13L), ncol = 1L),
      kx = 4L, ky = 5L, degree = 0L),
    lp_p2 = list(
      x = cbind(seq(-1.7, 2.3, length.out = 17L),
                cos(index * 0.47) - index / 31),
      kx = c(7L, 8L), ky = 7L, degree = c(1L, 1L)))

  for (fixture in fixtures) {
    x <- fixture$x
    y <- 0.35 + 0.62 * x[, 1L] +
      if (ncol(x) > 1L) -0.27 * x[, 2L] else 0 +
      sin(seq_len(nrow(x)) * 1.17) / 7
    for (kernel_name in c("gaussian", "epanechnikov")) {
      kernel <- if (kernel_name == "gaussian") {
        adaptive_conditional_cvml_gaussian2
      } else {
        adaptive_conditional_cvml_epanechnikov2
      }
      expected <- adaptive_conditional_cvml_literal(
        x, y, fixture$kx, fixture$ky, fixture$degree, kernel)
      observed <- unlist(lapply(c(FALSE, TRUE), function(tree) {
        vapply(c(FALSE, TRUE), function(accelerate) {
          options(np.tree = tree, np.macMseries.accelerate = accelerate)
          bw <- adaptive_conditional_cvml_bw(
            x, y, fixture$kx, fixture$ky, fixture$degree, kernel_name)
          adaptive_conditional_cvml_objective(x, y, bw)
        }, numeric(1L))
      }))
      expect_equal(observed, rep(expected, length(observed)),
                   tolerance = 5e-9)
    }
  }
})

test_that("adaptive conditional CVML fold geometry composes categorical sides", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old_options <- options(np.messages = FALSE, np.largeh = FALSE)
  on.exit(options(old_options), add = TRUE)

  n <- 17L
  z <- seq(-2, 2, length.out = n)
  xu <- factor(rep(c("a", "b", "c"), length.out = n))
  y <- 0.4 + 0.6 * z + sin(seq_len(n)) / 9
  yo <- ordered(rep(1:3, length.out = n))
  kx <- 7L
  ky <- 6L
  lambda_x <- 0.2
  lambda_y <- 0.25
  fit <- vapply(seq_len(n), function(held_out) {
    donor <- setdiff(seq_len(n), held_out)
    xradius <- vapply(donor, function(index) {
      adaptive_conditional_cvml_fold_radius(z, held_out, index, kx)
    }, numeric(1L))
    yradius <- vapply(donor, function(index) {
      adaptive_conditional_cvml_fold_radius(y, held_out, index, ky)
    }, numeric(1L))
    xcategorical <- ifelse(
      xu[donor] == xu[[held_out]], 1 - lambda_x,
      lambda_x / (nlevels(xu) - 1L))
    ydistance <- abs(as.integer(yo[donor]) - as.integer(yo[[held_out]]))
    ycategorical <- lambda_y^ydistance *
      (1 - lambda_y) / (1 + lambda_y)
    xweight <- stats::dnorm((z[[held_out]] - z[donor]) / xradius) /
      xradius * xcategorical
    ykernel <- stats::dnorm((y[[held_out]] - y[donor]) / yradius) /
      yradius * ycategorical
    sum(xweight / sum(xweight) * ykernel)
  }, numeric(1L))
  expected <- sum(log(fit))
  xdat <- data.frame(z = z, xu = xu)
  ydat <- data.frame(y = y, yo = yo)

  observed <- unlist(lapply(c(FALSE, TRUE), function(tree) {
    vapply(c(FALSE, TRUE), function(accelerate) {
      options(np.tree = tree, np.macMseries.accelerate = accelerate)
      bw <- npcdensbw(
        xdat = xdat, ydat = ydat,
        bws = c(ky, lambda_y, kx, lambda_x),
        bandwidth.compute = FALSE, bwmethod = "cv.ml",
        bwtype = "adaptive_nn", bwscaling = FALSE, regtype = "lc",
        cxkertype = "gaussian", cykertype = "gaussian",
        uxkertype = "aitchisonaitken", oykertype = "liracine")
      as.numeric(getFromNamespace(".npRmpi_with_local_regression", "npRmpi")(
        npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw)
      )$objective[[1L]])
    }, numeric(1L))
  }))
  expect_equal(observed, rep(expected, length(observed)), tolerance = 2e-10)
})

test_that("adaptive conditional CVML honors extension, replay, and failures", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old_options <- options(np.messages = FALSE, np.tree = FALSE,
                         np.macMseries.accelerate = FALSE,
                         np.extendednn = TRUE)
  on.exit(options(old_options), add = TRUE)

  x <- data.frame(x = c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1))
  y <- c(-2.2, -1.1, -0.4, 0.15, 0.8, 1.5, 2.6, 3.9)
  n <- nrow(x)
  admitted <- adaptive_conditional_cvml_bw(
    x, y, n - 2L, n - 2L, 0L, "gaussian")
  replay <- unserialize(serialize(admitted, NULL, version = 3L))
  expect_identical(
    adaptive_conditional_cvml_objective(x, y, replay, "dbmax"),
    adaptive_conditional_cvml_objective(x, y, admitted, "dbmax"))

  saturated <- admitted
  saturated$xbw[] <- n - 1L
  saturated$ybw[] <- n - 1L
  expect_equal(
    adaptive_conditional_cvml_objective(x, y, saturated, "dbmax"),
    adaptive_conditional_cvml_literal(
      x, y, n - 1L, n - 1L, 0L,
      adaptive_conditional_cvml_gaussian2),
    tolerance = 5e-9)

  beyond <- admitted
  beyond$xbw[] <- n + 2L
  beyond$ybw[] <- n + 2L
  expect_equal(
    adaptive_conditional_cvml_objective(x, y, beyond, "dbmax"),
    adaptive_conditional_cvml_literal(
      x, y, n + 2L, n + 2L, 0L,
      adaptive_conditional_cvml_gaussian2),
    tolerance = 5e-9)

  zero_x <- data.frame(x = c(0, 0, 0, 0, 1, 2))
  zero_y <- c(0, 0, 0, 0, 1, 2)
  zero_bw <- adaptive_conditional_cvml_bw(
    zero_x, zero_y, 2L, 2L, 0L, "gaussian")
  expect_identical(
    adaptive_conditional_cvml_objective(zero_x, zero_y, zero_bw, "dbmax"),
    -.Machine$double.xmax)
})
