adaptive_conditional_cvls_gaussian2 <- function(u) stats::dnorm(u)

adaptive_conditional_cvls_epanechnikov2 <- function(u) {
  edge <- sqrt(5)
  value <- numeric(length(u))
  inside <- abs(u) < edge
  value[inside] <- 3 / (4 * edge) * (1 - u[inside]^2 / 5)
  value
}

adaptive_conditional_cvls_fold_radius <- function(value, held_out, donor, k) {
  sort(abs(value[-c(held_out, donor)] - value[[donor]]),
       method = "quick")[[k]]
}

adaptive_conditional_cvls_design <- function(train, evaluation, degree) {
  delta <- sweep(as.matrix(train), 2L, evaluation, "-")
  if (all(degree == 0L))
    return(matrix(1, nrow = nrow(delta), ncol = 1L))
  if (ncol(delta) == 1L)
    return(vapply(0:degree[[1L]], function(power) delta[, 1L]^power,
                  numeric(nrow(delta))))
  stopifnot(all(degree == 1L))
  cbind(1, delta)
}

adaptive_conditional_cvls_continuous_overlap <- function(
  first, second, hfirst, hsecond, kernel_name, kernel
) {
  if (identical(kernel_name, "gaussian")) {
    scale <- sqrt(hfirst^2 + hsecond^2)
    return(stats::dnorm((first - second) / scale) / scale)
  }
  edge <- sqrt(5)
  lower <- max(first - edge * hfirst, second - edge * hsecond)
  upper <- min(first + edge * hfirst, second + edge * hsecond)
  if (lower >= upper)
    return(0)
  stats::integrate(
    function(value) {
      kernel((value - first) / hfirst) / hfirst *
        kernel((value - second) / hsecond) / hsecond
    }, lower = lower, upper = upper,
    subdivisions = 200L, rel.tol = 1e-12, abs.tol = 1e-14,
    stop.on.error = TRUE)$value
}

adaptive_conditional_cvls_unordered <- function(target, donor, lambda,
                                                 nlevels) {
  ifelse(target == donor, 1 - lambda, lambda / (nlevels - 1L))
}

adaptive_conditional_cvls_ordered <- function(target, donor, lambda) {
  lambda^abs(as.integer(target) - as.integer(donor)) *
    (1 - lambda) / (1 + lambda)
}

adaptive_conditional_cvls_ordered_overlap <- function(first, second,
                                                       lambda) {
  distance <- abs(as.integer(first) - as.integer(second))
  scale <- ((1 - lambda) / (1 + lambda))^2
  scale * lambda^distance *
    ((1 + lambda^2) / (1 - lambda^2) + distance)
}

adaptive_conditional_cvls_literal <- function(
  x, y, kx, ky, degree, kernel_name, kernel,
  xu = NULL, lambda_x = NULL, yo = NULL, lambda_y = NULL
) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  contribution <- numeric(n)

  for (held_out in seq_len(n)) {
    donor <- setdiff(seq_len(n), held_out)
    xweight <- rep.int(1, length(donor))
    for (coordinate in seq_len(ncol(x))) {
      radius <- vapply(donor, function(index) {
        adaptive_conditional_cvls_fold_radius(
          x[, coordinate], held_out, index, kx[[coordinate]])
      }, numeric(1L))
      xweight <- xweight * kernel(
        (x[held_out, coordinate] - x[donor, coordinate]) / radius) /
        radius
    }
    if (!is.null(xu)) {
      xweight <- xweight * adaptive_conditional_cvls_unordered(
        xu[[held_out]], xu[donor], lambda_x, nlevels(xu))
    }
    design <- adaptive_conditional_cvls_design(
      x[donor, , drop = FALSE], x[held_out, ], degree)
    coefficient <- solve(
      crossprod(design, xweight * design),
      c(1, rep.int(0, ncol(design) - 1L)))
    influence <- drop(xweight * (design %*% coefficient))

    yradius <- vapply(donor, function(index) {
      adaptive_conditional_cvls_fold_radius(y, held_out, index, ky)
    }, numeric(1L))
    ynormal <- kernel((y[[held_out]] - y[donor]) / yradius) / yradius
    if (!is.null(yo)) {
      ynormal <- ynormal * adaptive_conditional_cvls_ordered(
        yo[[held_out]], yo[donor], lambda_y)
    }

    overlap <- outer(seq_along(donor), seq_along(donor),
                     Vectorize(function(first, second) {
      value <- adaptive_conditional_cvls_continuous_overlap(
        y[donor[[first]]], y[donor[[second]]],
        yradius[[first]], yradius[[second]], kernel_name, kernel)
      if (!is.null(yo)) {
        value <- value * adaptive_conditional_cvls_ordered_overlap(
          yo[donor[[first]]], yo[donor[[second]]], lambda_y)
      }
      value
    }))
    linear <- sum(influence * ynormal)
    quadratic <- sum((influence %o% influence) * overlap)
    contribution[[held_out]] <- quadratic - 2 * linear
  }

  stopifnot(all(is.finite(contribution)))
  -mean(contribution)
}

adaptive_conditional_cvls_bw <- function(
  x, y, kx, ky, degree, kernel_name,
  xu = NULL, lambda_x = NULL, yo = NULL, lambda_y = NULL
) {
  xdat <- as.data.frame(x)
  names(xdat) <- paste0("x", seq_len(ncol(xdat)))
  ydat <- data.frame(y = y)
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
  list(xdat = xdat, ydat = ydat, bw = do.call(npcdensbw, call))
}

adaptive_conditional_cvls_objective <- function(state,
                                                 invalid.penalty = "baseline") {
  as.numeric(getFromNamespace(".npRmpi_with_local_regression", "npRmpi")(
    npRmpi:::.npcdensbw_eval_only(
      state$xdat, state$ydat, state$bw,
      invalid.penalty = invalid.penalty)
  )$objective[[1L]])
}

test_that("adaptive conditional CVLS uses one exact fold for I1 and I2", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old <- options(np.messages = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  index <- seq_len(13L)
  fixtures <- list(
    lc_p1 = list(
      x = matrix(seq(-1.8, 2.2, length.out = 13L), ncol = 1L),
      kx = 4L, ky = 5L, degree = 0L),
    lp_p2 = list(
      x = cbind(seq(-1.7, 2.3, length.out = 13L),
                cos(index * 0.47) - index / 31),
      kx = c(5L, 6L), ky = 5L, degree = c(1L, 1L)),
    lp_p3 = list(
      x = cbind(seq(-1.7, 2.3, length.out = 13L),
                cos(index * 0.93) - index / 31,
                sin(index * 1.37) + index / 37),
      kx = c(9L, 9L, 9L), ky = 5L, degree = c(1L, 1L, 1L)))

  for (fixture in fixtures) {
    x <- fixture$x
    y <- 0.35 + 0.58 * x[, 1L] +
      (if (ncol(x) > 1L) -0.24 * x[, 2L] else 0) +
      (if (ncol(x) > 2L) 0.18 * x[, 3L] else 0) + sin(index) / 9
    for (kernel_name in c("gaussian", "epanechnikov")) {
      kernel <- if (kernel_name == "gaussian") {
        adaptive_conditional_cvls_gaussian2
      } else {
        adaptive_conditional_cvls_epanechnikov2
      }
      expected <- adaptive_conditional_cvls_literal(
        x, y, fixture$kx, fixture$ky, fixture$degree,
        kernel_name, kernel)
      observed <- unlist(lapply(c(FALSE, TRUE), function(tree) {
        vapply(c(FALSE, TRUE), function(accelerate) {
          options(np.tree = tree, np.macMseries.accelerate = accelerate)
          state <- adaptive_conditional_cvls_bw(
            x, y, fixture$kx, fixture$ky, fixture$degree, kernel_name)
          adaptive_conditional_cvls_objective(state)
        }, numeric(1L))
      }))
      expect_equal(observed, rep(expected, length(observed)), tolerance = 5e-9)
    }
  }
})

test_that("adaptive conditional CVLS exact folds compose categorical sides", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old <- options(np.messages = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 13L
  x <- matrix(seq(-1.8, 2.2, length.out = n), ncol = 1L)
  y <- 0.4 + 0.55 * x[, 1L] + sin(seq_len(n) * 1.13) / 8
  xu <- factor(rep(c("a", "b", "c"), length.out = n))
  yo <- ordered(rep(1:3, length.out = n))
  expected <- adaptive_conditional_cvls_literal(
    x, y, 5L, 5L, 0L, "gaussian", adaptive_conditional_cvls_gaussian2,
    xu = xu, lambda_x = 0.2, yo = yo, lambda_y = 0.25)
  observed <- unlist(lapply(c(FALSE, TRUE), function(tree) {
    vapply(c(FALSE, TRUE), function(accelerate) {
      options(np.tree = tree, np.macMseries.accelerate = accelerate)
      state <- adaptive_conditional_cvls_bw(
        x, y, 5L, 5L, 0L, "gaussian",
        xu = xu, lambda_x = 0.2, yo = yo, lambda_y = 0.25)
      adaptive_conditional_cvls_objective(state)
    }, numeric(1L))
  }))
  expect_equal(observed, rep(expected, length(observed)), tolerance = 5e-9)
})

test_that("adaptive conditional CVLS honors saturation, replay, and failures", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)

  x <- matrix(c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1), ncol = 1L)
  y <- c(-2.2, -1.1, -0.4, 0.15, 0.8, 1.5, 2.6, 3.9)
  n <- nrow(x)
  state <- adaptive_conditional_cvls_bw(
    x, y, n - 2L, n - 2L, 0L, "gaussian")
  replay <- state
  replay$bw <- unserialize(serialize(state$bw, NULL, version = 3L))
  expect_identical(adaptive_conditional_cvls_objective(replay, "dbmax"),
                   adaptive_conditional_cvls_objective(state, "dbmax"))

  saturated <- state
  saturated$bw$xbw[] <- n - 1L
  saturated$bw$ybw[] <- n - 1L
  expect_identical(adaptive_conditional_cvls_objective(saturated, "dbmax"),
                   -.Machine$double.xmax)

  zero <- adaptive_conditional_cvls_bw(
    matrix(c(0, 0, 0, 0, 1, 2), ncol = 1L),
    c(0, 0, 0, 0, 1, 2), 2L, 2L, 0L, "gaussian")
  expect_identical(adaptive_conditional_cvls_objective(zero, "dbmax"),
                   -.Machine$double.xmax)
})
