adaptive_density_gaussian2 <- function(u) stats::dnorm(u)

adaptive_density_epanechnikov2 <- function(u) {
  edge <- sqrt(5)
  value <- numeric(length(u))
  inside <- abs(u) < edge
  value[inside] <- 3 / (4 * edge) * (1 - u[inside]^2 / 5)
  value
}

adaptive_density_fold_fit <- function(x, k, kernel) {
  x <- as.matrix(x)
  n <- nrow(x)
  fit <- numeric(n)
  for (held_out in seq_len(n)) {
    donors <- setdiff(seq_len(n), held_out)
    weight <- rep(1, length(donors))
    for (coordinate in seq_len(ncol(x))) {
      radius <- vapply(donors, function(donor) {
        distance <- sort(abs(x[-c(held_out, donor), coordinate] -
                             x[donor, coordinate]), method = "quick")
        lookup <- min(k[coordinate], length(distance))
        distance[[lookup]] * k[coordinate] / lookup
      }, numeric(1L))
      weight <- weight * kernel(
        (x[held_out, coordinate] - x[donors, coordinate]) / radius) /
        radius
    }
    fit[held_out] <- sum(weight) / (n - 1)
  }
  fit
}

adaptive_density_full_radius <- function(x, donor, k) {
  distance <- sort(abs(x[-donor] - x[donor]), method = "quick")
  lookup <- min(k, length(distance))
  distance[[lookup]] * k / lookup
}

adaptive_density_gaussian_i1 <- function(x, k) {
  x <- as.matrix(x)
  n <- nrow(x)
  bandwidth <- matrix(NA_real_, nrow = n, ncol = ncol(x))
  for (donor in seq_len(n))
    for (coordinate in seq_len(ncol(x)))
      bandwidth[donor, coordinate] <- adaptive_density_full_radius(
        x[, coordinate], donor, k[coordinate])

  total <- 0
  for (left in seq_len(n)) {
    for (right in seq_len(n)) {
      pair <- 1
      for (coordinate in seq_len(ncol(x))) {
        combined <- sqrt(bandwidth[left, coordinate]^2 +
                         bandwidth[right, coordinate]^2)
        pair <- pair * stats::dnorm(
          (x[left, coordinate] - x[right, coordinate]) / combined) /
          combined
      }
      total <- total + pair
    }
  }
  total / n^2
}

adaptive_density_package_objective <- function(x, k, method, kernel) {
  dat <- as.data.frame(x)
  names(dat) <- paste0("x", seq_len(ncol(dat)))
  bw <- npudensbw(
    dat = dat, bws = as.double(k), bwmethod = method,
    bwtype = "adaptive_nn", bwscaling = FALSE,
    ckertype = kernel, ckerorder = 2L, bandwidth.compute = FALSE)
  evaluated <- np:::npudensbw.bandwidth(
    dat = dat, bws = bw, bandwidth.compute = TRUE,
    nmulti = 1L, powell.remin = FALSE, bwsolver = "powell",
    eval.only = TRUE)
  -as.numeric(evaluated$fval)
}

adaptive_density_package_objective_from_bw <- function(dat, bw) {
  evaluated <- np:::npudensbw.bandwidth(
    dat = dat, bws = bw, bandwidth.compute = TRUE,
    nmulti = 1L, powell.remin = FALSE, bwsolver = "powell",
    eval.only = TRUE)
  -as.numeric(evaluated$fval)
}

test_that("adaptive density CVML uses exact delete-one donor radii", {
  old_options <- options(np.tree = TRUE, np.macMseries.accelerate = TRUE,
                         np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  x <- cbind(seq(-1.7, 2.3, length.out = 13L),
             sin(seq_len(13L) * 0.73) + seq_len(13L) / 17)
  k <- c(4L, 5L)

  for (kernel_name in c("gaussian", "epanechnikov")) {
    kernel <- if (kernel_name == "gaussian") {
      adaptive_density_gaussian2
    } else {
      adaptive_density_epanechnikov2
    }
    fit <- adaptive_density_fold_fit(x, k, kernel)
    expect_true(all(is.finite(fit) & fit > 0))
    expect_equal(
      adaptive_density_package_objective(x, k, "cv.ml", kernel_name),
      sum(-log(fit)), tolerance = 1e-10)
  }
})

test_that("adaptive density CVLS changes only its delete-one cross term", {
  old_options <- options(np.tree = FALSE, np.macMseries.accelerate = FALSE,
                         np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  x <- matrix(c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1),
              ncol = 1L)
  k <- 3L
  fit <- adaptive_density_fold_fit(x, k, adaptive_density_gaussian2)
  expected <- adaptive_density_gaussian_i1(x, k) - 2 * mean(fit)

  expect_equal(
    adaptive_density_package_objective(x, k, "cv.ls", "gaussian"),
    expected, tolerance = 2e-9)
})

test_that("adaptive density fold geometry composes with categorical kernels", {
  old_options <- options(np.tree = TRUE, np.macMseries.accelerate = FALSE,
                         np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  z <- c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1, 5, 6.2)
  group <- factor(c("a", "b", "a", "c", "b", "c", "a", "b", "c", "a"))
  dat <- data.frame(z = z, group = group)
  k <- 3L
  lambda <- 0.2
  fit <- vapply(seq_along(z), function(held_out) {
    donors <- setdiff(seq_along(z), held_out)
    sum(vapply(donors, function(donor) {
      radius <- sort(abs(z[-c(held_out, donor)] - z[donor]),
                     method = "quick")[k]
      categorical <- if (group[held_out] == group[donor]) {
        1 - lambda
      } else {
        lambda / (nlevels(group) - 1)
      }
      stats::dnorm((z[held_out] - z[donor]) / radius) /
        radius * categorical
    }, numeric(1L))) / (length(z) - 1)
  }, numeric(1L))
  bw <- npudensbw(
    dat = dat, bws = c(k, lambda), bwmethod = "cv.ml",
    bwtype = "adaptive_nn", bwscaling = FALSE,
    ckertype = "gaussian", ukertype = "aitchisonaitken",
    bandwidth.compute = FALSE)

  expect_equal(adaptive_density_package_objective_from_bw(dat, bw),
               sum(-log(fit)), tolerance = 1e-10)
})

test_that("adaptive density CV honors extended folds and serialized replay", {
  old_options <- options(np.tree = FALSE, np.macMseries.accelerate = FALSE,
                         np.messages = FALSE, np.extendednn = TRUE)
  on.exit(options(old_options), add = TRUE)

  dat <- data.frame(
    x = c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1))
  for (method in c("cv.ml", "cv.ls")) {
    admitted <- npudensbw(
      dat = dat, bws = nrow(dat) - 2L, bwmethod = method,
      bwtype = "adaptive_nn", bwscaling = FALSE,
      ckertype = "gaussian", bandwidth.compute = FALSE)
    replay <- unserialize(serialize(admitted, NULL, version = 3L))
    expect_identical(replay$bw, admitted$bw)
    expect_identical(replay$type, admitted$type)
    expect_identical(replay$method, admitted$method)
    expect_equal(adaptive_density_package_objective_from_bw(dat, replay),
                 adaptive_density_package_objective_from_bw(dat, admitted),
                 tolerance = 0)

    saturated <- admitted
    saturated$bw[] <- nrow(dat) - 1L
    expected <- adaptive_density_gaussian_i1(
      dat$x, saturated$bw) -
      2 * mean(adaptive_density_fold_fit(
        dat$x, saturated$bw, adaptive_density_gaussian2))
    if (identical(method, "cv.ml"))
      expected <- sum(-log(adaptive_density_fold_fit(
        dat$x, saturated$bw, adaptive_density_gaussian2)))
    expect_equal(
      adaptive_density_package_objective_from_bw(dat, saturated),
      expected,
      tolerance = 2e-9
    )

    beyond <- admitted
    beyond$bw[] <- nrow(dat) + 2L
    expected_beyond <- adaptive_density_gaussian_i1(
      dat$x, beyond$bw) -
      2 * mean(adaptive_density_fold_fit(
        dat$x, beyond$bw, adaptive_density_gaussian2))
    if (identical(method, "cv.ml"))
      expected_beyond <- sum(-log(adaptive_density_fold_fit(
        dat$x, beyond$bw, adaptive_density_gaussian2)))
    expect_equal(
      adaptive_density_package_objective_from_bw(dat, beyond),
      expected_beyond,
      tolerance = 2e-9
    )

    options(np.extendednn = FALSE)
    expect_error(
      adaptive_density_package_objective_from_bw(dat, saturated),
      "invalid bandwidth")
    options(np.extendednn = TRUE)
  }
})
