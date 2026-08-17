adaptive_conditional_cdist_gaussian_density <- function(u) stats::dnorm(u)
adaptive_conditional_cdist_gaussian_cdf <- function(u) stats::pnorm(u)

adaptive_conditional_cdist_epan_density <- function(u) {
  edge <- sqrt(5)
  value <- numeric(length(u))
  inside <- abs(u) < edge
  value[inside] <- 3 / (4 * edge) * (1 - u[inside]^2 / 5)
  value
}

adaptive_conditional_cdist_epan_cdf <- function(u) {
  edge <- sqrt(5)
  value <- numeric(length(u))
  value[u >= edge] <- 1
  inside <- abs(u) < edge
  value[inside] <- 0.5 + 3 * u[inside] / (4 * edge) -
    u[inside]^3 / (20 * edge)
  value
}

adaptive_conditional_cdist_fold_radius <- function(value, held_out, donor, k) {
  sort(abs(value[-c(held_out, donor)] - value[[donor]]),
       method = "quick")[[k]]
}

adaptive_conditional_cdist_design <- function(train, evaluation, degree) {
  delta <- sweep(as.matrix(train), 2L, evaluation, "-")
  if (all(degree == 0L))
    return(matrix(1, nrow = nrow(delta), ncol = 1L))
  if (ncol(delta) == 1L)
    return(vapply(0:degree[[1L]], function(power) delta[, 1L]^power,
                  numeric(nrow(delta))))
  stopifnot(all(degree == 1L))
  cbind(1, delta)
}

adaptive_conditional_cdist_unordered <- function(target, donor, lambda,
                                                  nlevels) {
  ifelse(target == donor, 1 - lambda, lambda / (nlevels - 1L))
}

adaptive_conditional_cdist_ordered_cdf <- function(target, donor, lambda) {
  distance <- abs(as.integer(target) - as.integer(donor))
  geometric <- lambda^distance / (1 + lambda)
  ifelse(as.integer(target) < as.integer(donor),
         geometric, 1 - lambda * geometric)
}

adaptive_conditional_cdist_literal <- function(
  x, y, kx, ky, degree, density_kernel, cdf_kernel,
  xu = NULL, lambda_x = NULL, yo = NULL, lambda_y = NULL
) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  total <- 0

  for (held_out in seq_len(n)) {
    donor <- setdiff(seq_len(n), held_out)
    xweight <- rep.int(1, length(donor))
    for (coordinate in seq_len(ncol(x))) {
      radius <- vapply(donor, function(index) {
        adaptive_conditional_cdist_fold_radius(
          x[, coordinate], held_out, index, kx[[coordinate]])
      }, numeric(1L))
      xweight <- xweight * density_kernel(
        (x[held_out, coordinate] - x[donor, coordinate]) / radius) /
        radius
    }
    if (!is.null(xu)) {
      xweight <- xweight * adaptive_conditional_cdist_unordered(
        xu[[held_out]], xu[donor], lambda_x, nlevels(xu))
    }
    design <- adaptive_conditional_cdist_design(
      x[donor, , drop = FALSE], x[held_out, ], degree)
    coefficient <- solve(
      crossprod(design, xweight * design),
      c(1, rep.int(0, ncol(design) - 1L)))
    influence <- drop(xweight * (design %*% coefficient))

    yradius <- lapply(seq_len(ncol(y)), function(coordinate) {
      vapply(donor, function(index) {
        adaptive_conditional_cdist_fold_radius(
          y[, coordinate], held_out, index, ky[[coordinate]])
      }, numeric(1L))
    })
    for (query in seq_len(n)) {
      if (query == held_out)
        next
      yrow <- rep.int(1, length(donor))
      for (coordinate in seq_len(ncol(y))) {
        yrow <- yrow * cdf_kernel(
          (y[query, coordinate] - y[donor, coordinate]) /
            yradius[[coordinate]])
      }
      if (!is.null(yo)) {
        yrow <- yrow * adaptive_conditional_cdist_ordered_cdf(
          yo[[query]], yo[donor], lambda_y)
      }
      fit <- sum(influence * yrow)
      indicator <- as.integer(all(y[held_out, ] <= y[query, ]))
      if (!is.null(yo)) {
        indicator <- indicator *
          as.integer(as.integer(yo[[held_out]]) <= as.integer(yo[[query]]))
      }
      total <- total + (indicator - fit)^2
    }
  }
  total / (n * (n - 1L))
}

adaptive_conditional_cdist_state <- function(
  x, y, kx, ky, degree, kernel_name,
  xu = NULL, lambda_x = NULL, yo = NULL, lambda_y = NULL
) {
  xdat <- as.data.frame(x)
  names(xdat) <- paste0("x", seq_len(ncol(xdat)))
  ydat <- as.data.frame(y)
  names(ydat) <- paste0("y", seq_len(ncol(ydat)))
  bandwidths <- as.double(c(ky, kx))
  if (!is.null(xu)) {
    xdat$xu <- xu
    bandwidths <- c(ky, lambda_y, kx, lambda_x)
  }
  if (!is.null(yo))
    ydat$yo <- yo
  call <- list(
    xdat = xdat, ydat = ydat, bws = bandwidths,
    bandwidth.compute = FALSE, bwmethod = "cv.ls",
    bwtype = "adaptive_nn", bwscaling = FALSE,
    regtype = if (all(degree == 0L)) "lc" else "lp",
    cxkertype = kernel_name, cxkerorder = 2L,
    cykertype = kernel_name, cykerorder = 2L)
  if (!all(degree == 0L)) {
    call$degree <- as.integer(degree)
    call$degree.select <- "manual"
    call$basis <- "glp"
  }
  if (!is.null(xu)) {
    call$uxkertype <- "aitchisonaitken"
    call$oykertype <- "liracine"
  }
  list(xdat = xdat, ydat = ydat, bw = do.call(npcdistbw, call))
}

adaptive_conditional_cdist_objective <- function(state,
                                                  invalid.penalty = "baseline") {
  as.numeric(getFromNamespace(".npRmpi_with_local_cdist_eval", "npRmpi")(
    npRmpi:::.npcdistbw_eval_only(
      state$xdat, state$ydat, bws = state$bw,
      do.full.integral = TRUE,
      invalid.penalty = invalid.penalty)
  )$objective[[1L]])
}

test_that("adaptive empirical conditional-CDF CV uses exact donor folds", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old <- options(np.messages = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  index <- seq_len(13L)
  fixtures <- list(
    lc_p1_y1 = list(
      x = matrix(seq(-1.8, 2.2, length.out = 13L), ncol = 1L),
      y2 = FALSE, kx = 4L, ky = 5L, degree = 0L),
    lp_p2_y2 = list(
      x = cbind(seq(-1.7, 2.3, length.out = 13L),
                cos(index * 0.93) - index / 31),
      y2 = TRUE, kx = c(6L, 7L), ky = c(5L, 6L),
      degree = c(1L, 1L)),
    lp_p3_y1 = list(
      x = cbind(seq(-1.7, 2.3, length.out = 13L),
                cos(index * 0.93) - index / 31,
                sin(index * 1.37) + index / 37),
      y2 = FALSE, kx = c(9L, 9L, 9L), ky = 5L,
      degree = c(1L, 1L, 1L)))

  for (fixture in fixtures) {
    y1 <- 0.35 + 0.58 * fixture$x[, 1L] +
      (if (ncol(fixture$x) > 1L) -0.24 * fixture$x[, 2L] else 0) +
      sin(index) / 9
    y <- if (fixture$y2) {
      cbind(y1, -0.31 * fixture$x[, 1L] + cos(index * 0.61) / 7)
    } else {
      cbind(y1)
    }
    for (kernel_name in c("gaussian", "epanechnikov")) {
      density_kernel <- if (kernel_name == "gaussian") {
        adaptive_conditional_cdist_gaussian_density
      } else {
        adaptive_conditional_cdist_epan_density
      }
      cdf_kernel <- if (kernel_name == "gaussian") {
        adaptive_conditional_cdist_gaussian_cdf
      } else {
        adaptive_conditional_cdist_epan_cdf
      }
      expected <- adaptive_conditional_cdist_literal(
        fixture$x, y, fixture$kx, fixture$ky, fixture$degree,
        density_kernel, cdf_kernel)
      observed <- unlist(lapply(c(FALSE, TRUE), function(tree) {
        vapply(c(FALSE, TRUE), function(accelerate) {
          options(np.tree = tree, np.macMseries.accelerate = accelerate)
          state <- adaptive_conditional_cdist_state(
            fixture$x, y, fixture$kx, fixture$ky, fixture$degree,
            kernel_name)
          adaptive_conditional_cdist_objective(state)
        }, numeric(1L))
      }))
      expect_equal(observed, rep(expected, length(observed)), tolerance = 5e-9)
    }
  }
})

test_that("adaptive empirical conditional-CDF folds compose categorical sides", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old <- options(np.messages = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  index <- seq_len(13L)
  x <- matrix(seq(-1.8, 2.2, length.out = 13L), ncol = 1L)
  y <- cbind(0.4 + 0.55 * x[, 1L] + sin(index * 1.13) / 8)
  xu <- factor(rep(c("a", "b", "c"), length.out = nrow(x)))
  yo <- ordered(rep(1:3, length.out = nrow(x)))
  expected <- adaptive_conditional_cdist_literal(
    x, y, 5L, 5L, 0L,
    adaptive_conditional_cdist_gaussian_density,
    adaptive_conditional_cdist_gaussian_cdf,
    xu = xu, lambda_x = 0.2, yo = yo, lambda_y = 0.25)
  observed <- unlist(lapply(c(FALSE, TRUE), function(tree) {
    vapply(c(FALSE, TRUE), function(accelerate) {
      options(np.tree = tree, np.macMseries.accelerate = accelerate)
      state <- adaptive_conditional_cdist_state(
        x, y, 5L, 5L, 0L, "gaussian",
        xu = xu, lambda_x = 0.2, yo = yo, lambda_y = 0.25)
      adaptive_conditional_cdist_objective(state)
    }, numeric(1L))
  }))
  expect_equal(observed, rep(expected, length(observed)), tolerance = 5e-9)
})

test_that("adaptive empirical conditional-CDF honors saturation and replay", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)

  x <- matrix(c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1), ncol = 1L)
  y <- cbind(c(-2.2, -1.1, -0.4, 0.15, 0.8, 1.5, 2.6, 3.9))
  n <- nrow(x)
  state <- adaptive_conditional_cdist_state(
    x, y, n - 2L, n - 2L, 0L, "gaussian")
  replay <- state
  replay$bw <- unserialize(serialize(state$bw, NULL, version = 3L))
  expected <- adaptive_conditional_cdist_literal(
    x, y, n - 2L, n - 2L, 0L,
    adaptive_conditional_cdist_gaussian_density,
    adaptive_conditional_cdist_gaussian_cdf)
  expect_equal(adaptive_conditional_cdist_objective(replay, "dbmax"),
               expected, tolerance = 5e-9)
  expect_identical(adaptive_conditional_cdist_objective(replay, "dbmax"),
                   adaptive_conditional_cdist_objective(state, "dbmax"))

  saturated <- state
  saturated$bw$xbw[] <- n - 1L
  saturated$bw$ybw[] <- n - 1L
  expect_identical(adaptive_conditional_cdist_objective(saturated, "dbmax"),
                   .Machine$double.xmax)

  zero <- adaptive_conditional_cdist_state(
    matrix(c(0, 0, 0, 0, 1, 2), ncol = 1L),
    cbind(c(0, 0, 0, 0, 1, 2)), 2L, 2L, 0L, "gaussian")
  expect_identical(adaptive_conditional_cdist_objective(zero, "dbmax"),
                   .Machine$double.xmax)
})
