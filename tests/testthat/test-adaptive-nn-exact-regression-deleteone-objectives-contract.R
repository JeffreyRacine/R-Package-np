adaptive_deleteone_gaussian2 <- function(u) stats::dnorm(u)

adaptive_deleteone_epanechnikov2 <- function(u) {
  edge <- sqrt(5)
  value <- numeric(length(u))
  inside <- abs(u) < edge
  value[inside] <- 3 / (4 * edge) * (1 - u[inside]^2 / 5)
  value
}

adaptive_deleteone_literal_fit <- function(x, response, k, kernel) {
  x <- as.matrix(x)
  n <- nrow(x)
  fitted <- numeric(n)
  for (held_out in seq_len(n)) {
    numerator <- denominator <- 0
    for (donor in setdiff(seq_len(n), held_out)) {
      weight <- 1
      for (coordinate in seq_len(ncol(x))) {
        radius <- sort(
          abs(x[-c(held_out, donor), coordinate] - x[donor, coordinate]),
          method = "quick")[k[coordinate]]
        weight <- weight * kernel(
          (x[held_out, coordinate] - x[donor, coordinate]) / radius) /
          radius
      }
      numerator <- numerator + weight * response[donor]
      denominator <- denominator + weight
    }
    fitted[held_out] <- numerator / denominator
  }
  fitted
}

adaptive_deleteone_literal_lp1_fit <- function(x, response, k, kernel) {
  x <- as.matrix(x)
  n <- nrow(x)
  fitted <- numeric(n)
  for (held_out in seq_len(n)) {
    donors <- setdiff(seq_len(n), held_out)
    weight <- rep(1, length(donors))
    for (coordinate in seq_len(ncol(x))) {
      radius <- vapply(donors, function(donor) {
        sort(abs(x[-c(held_out, donor), coordinate] -
                   x[donor, coordinate]), method = "quick")[k[coordinate]]
      }, numeric(1L))
      weight <- weight * kernel(
        (x[held_out, coordinate] - x[donors, coordinate]) / radius) /
        radius
    }
    design <- cbind(
      1, sweep(x[donors,, drop = FALSE], 2L, x[held_out, ], "-"))
    coefficient <- solve(
      crossprod(design, weight * design),
      crossprod(design, weight * response[donors]))
    fitted[held_out] <- coefficient[1L]
  }
  fitted
}

adaptive_deleteone_package_ks <- function(x, y, k, kernel) {
  x <- as.data.frame(x)
  names(x) <- paste0("x", seq_len(ncol(x)))
  bw <- npregbw(
    xdat = x, ydat = y, bws = as.double(k), regtype = "lc",
    bwmethod = "cv.ls", bwtype = "adaptive_nn", bwscaling = FALSE,
    ckertype = kernel, ckerorder = 2L, bandwidth.compute = FALSE)
  as.numeric(np:::.npregbw_eval_only(x, y, bw, objective = "ks")$objective)
}

adaptive_deleteone_package_check <- function(x, y, scale, k, tau, delta,
                                             kernel, regtype = "lc",
                                             degree = NULL) {
  x <- as.data.frame(x)
  names(x) <- paste0("x", seq_len(ncol(x)))
  arguments <- list(
    xdat = x, ydat = y, bws = as.double(k), scale = scale,
    tau = tau, delta = delta, delta.bounds = c(0.1, 0.9),
    regtype = regtype, bwtype = "adaptive_nn", bwscaling = FALSE,
    ckertype = kernel, ckerorder = 2L, bandwidth.compute = FALSE)
  if (identical(regtype, "lp")) {
    arguments$degree <- as.integer(degree)
    arguments$degree.select <- "manual"
    arguments$basis <- "glp"
  }
  bw <- do.call(nplsqregbw, arguments)
  as.numeric(bw$objective)
}

adaptive_deleteone_package_lp <- function(x, y, k, kernel, objective) {
  x <- as.data.frame(x)
  names(x) <- paste0("x", seq_len(ncol(x)))
  bw <- npregbw(
    xdat = x, ydat = y, bws = as.double(k), regtype = "lp",
    degree = rep.int(1L, ncol(x)), degree.select = "manual",
    basis = "glp", bwmethod = "cv.ls", bwtype = "adaptive_nn",
    bwscaling = FALSE, ckertype = kernel, ckerorder = 2L,
    bandwidth.compute = FALSE)
  as.numeric(np:::.npregbw_eval_only(
    x, y, bw, objective = objective)$objective)
}

test_that("adaptive scalar KS and check loss use exact delete-one radii", {
  old_options <- options(np.tree = FALSE, np.macMseries.accelerate = FALSE,
                         np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  x <- cbind(seq(-1.7, 2.3, length.out = 13L),
             sin(seq_len(13L) * 0.73) + seq_len(13L) / 17)
  k <- c(4L, 5L)
  latent <- 0.35 + 0.8 * x[, 1L] - 0.45 * x[, 2L]
  y_binary <- as.double(
    latent + sin(seq_len(nrow(x)) * 1.7) / 3 > median(latent))
  y_check <- latent + cos(seq_len(nrow(x)) * 0.8) / 5
  scale <- 0.7 + seq_len(nrow(x)) / (3 * nrow(x))
  tau <- 0.37
  delta <- 0.58

  for (kernel_name in c("gaussian", "epanechnikov")) {
    kernel <- if (kernel_name == "gaussian") {
      adaptive_deleteone_gaussian2
    } else {
      adaptive_deleteone_epanechnikov2
    }

    fit_ks <- adaptive_deleteone_literal_fit(x, y_binary, k, kernel)
    floor <- sqrt(.Machine$double.eps)
    probability <- pmin(1 - floor, pmax(floor, fit_ks))
    expected_ks <- mean(-(
      y_binary * log(probability) +
        (1 - y_binary) * log1p(-probability)))

    quantile_response <- y_check + scale * stats::qnorm(delta)
    fit_check <- adaptive_deleteone_literal_fit(
      x, quantile_response, k, kernel)
    residual <- y_check - fit_check
    expected_check <- mean(
      residual * (tau - as.double(residual < 0)))

    expect_equal(adaptive_deleteone_package_ks(
      x, y_binary, k, kernel_name), expected_ks,
      tolerance = 4096 * .Machine$double.eps)
    expect_equal(adaptive_deleteone_package_check(
      x, y_check, scale, k, tau, delta, kernel_name), expected_check,
      tolerance = 4096 * .Machine$double.eps)
  }
})

test_that("adaptive LP objectives share exact delete-one radii", {
  old_options <- options(np.tree = TRUE, np.macMseries.accelerate = TRUE,
                         np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  x <- cbind(seq(-1.8, 2.2, length.out = 15L),
             sin(seq_len(15L) * 0.61) + seq_len(15L) / 27)
  k <- c(6L, 7L)
  latent <- 0.25 + 0.65 * x[, 1L] - 0.3 * x[, 2L]
  y <- latent + sin(seq_len(nrow(x)) * 1.13) / 8
  y_binary <- as.double(
    latent + cos(seq_len(nrow(x)) * 1.41) / 3 > median(latent))
  scale <- 0.8 + seq_len(nrow(x)) / (4 * nrow(x))
  tau <- 0.41
  delta <- 0.57

  for (kernel_name in c("gaussian", "epanechnikov")) {
    kernel <- if (kernel_name == "gaussian") {
      adaptive_deleteone_gaussian2
    } else {
      adaptive_deleteone_epanechnikov2
    }

    fit_ls <- adaptive_deleteone_literal_lp1_fit(x, y, k, kernel)
    expected_ls <- mean((y - fit_ls)^2)
    fit_ks <- adaptive_deleteone_literal_lp1_fit(x, y_binary, k, kernel)
    floor <- sqrt(.Machine$double.eps)
    probability <- pmin(1 - floor, pmax(floor, fit_ks))
    expected_ks <- mean(-(
      y_binary * log(probability) +
        (1 - y_binary) * log1p(-probability)))
    fit_check <- adaptive_deleteone_literal_lp1_fit(
      x, y + scale * stats::qnorm(delta), k, kernel)
    residual <- y - fit_check
    expected_check <- mean(
      residual * (tau - as.double(residual < 0)))

    expect_equal(adaptive_deleteone_package_lp(
      x, y, k, kernel_name, "ls"), expected_ls,
      tolerance = 2e-9)
    expect_equal(adaptive_deleteone_package_lp(
      x, y_binary, k, kernel_name, "ks"), expected_ks,
      tolerance = 2e-9)
    expect_equal(adaptive_deleteone_package_check(
      x, y, scale, k, tau, delta, kernel_name,
      regtype = "lp", degree = c(1L, 1L)), expected_check,
      tolerance = 2e-9)
  }
})
