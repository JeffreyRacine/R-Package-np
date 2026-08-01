expect_beta_dual_square_equivalent <- function(actual, expected) {
  expect_equal(actual, expected, tolerance = 16 * .Machine$double.eps)
}

test_that("beta dual-power rows equal independently evaluated powers", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  dual_sum <- getFromNamespace(".npksum_power12", "np")

  training <- data.frame(
    x = c(.017, .08, .21, .47, .73, .94, .989),
    z = c(.9, .18, .67, .31, .79, .42, .06)
  )
  evaluation <- data.frame(
    x = c(.11, .37, .82),
    z = c(.2, .7, .44)
  )

  for (mode in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(mode, "fixed")) c(.14, .18) else c(4, 4)
    for (order in c(2L, 4L, 6L, 8L)) {
      arguments <- list(
        bws = bandwidth, txdat = training, exdat = evaluation,
        bwtype = mode, bandwidth.divide = TRUE,
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
      )
      dual <- do.call(dual_sum, arguments)
      power_one <- do.call(npksum, arguments)
      power_two <- do.call(npksum, c(arguments, list(kernel.pow = 2L)))

      expect_identical(dual$ksum, power_one$ksum)
      expect_beta_dual_square_equivalent(
        dual$ksum.power2, power_two$ksum
      )
    }
  }
})

test_that("beta dual-power rows preserve deletion and weighted semantics", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)
  dual_sum <- getFromNamespace(".npksum_power12", "np")
  weighted_sum <- getFromNamespace(".npksum_power12_weighted", "np")

  training <- data.frame(
    x = c(.03, .14, .29, .46, .63, .78, .92),
    z = c(.86, .2, .72, .35, .66, .43, .09)
  )
  counts <- c(1, 3, 2, 4, 1, 2, 3)
  arguments <- list(
    bws = c(.16, .2), txdat = training, leave.one.out = TRUE,
    bandwidth.divide = FALSE,
    ckertype = "beta", ckerorder = 8L,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )

  dual <- do.call(dual_sum, arguments)
  power_one <- do.call(npksum, arguments)
  power_two <- do.call(npksum, c(arguments, list(kernel.pow = 2L)))
  expect_identical(dual$ksum, power_one$ksum)
  expect_beta_dual_square_equivalent(dual$ksum.power2, power_two$ksum)

  weighted_arguments <- arguments
  weighted_arguments$leave.one.out <- NULL
  weighted <- do.call(
    weighted_sum,
    c(weighted_arguments, list(counts = counts))
  )
  weighted_one <- do.call(
    npksum,
    c(weighted_arguments, list(tydat = counts))
  )
  weighted_two <- do.call(
    npksum,
    c(weighted_arguments, list(tydat = counts, kernel.pow = 2L))
  )
  expect_identical(weighted$ksum, weighted_one$ksum)
  expect_beta_dual_square_equivalent(
    weighted$ksum.power2, weighted_two$ksum
  )
})

test_that("second-order beta dual powers match a direct scalar oracle", {
  dual_sum <- getFromNamespace(".npksum_power12", "np")
  training <- data.frame(x = c(.02, .11, .27, .51, .74, .93))
  evaluation <- data.frame(x = c(.08, .36, .81))
  bandwidth <- .17
  tau <- (1 / bandwidth)^2
  arguments <- list(
    bws = bandwidth, txdat = training, exdat = evaluation,
    ckertype = "beta", ckerorder = 2L,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )

  dual <- do.call(dual_sum, arguments)
  raw <- vapply(evaluation$x, function(target) {
    dbeta(training$x,
          shape1 = 1 + target * tau,
          shape2 = 1 + (1 - target) * tau)
  }, numeric(nrow(training)))

  expect_equal(as.double(dual$ksum), colSums(raw), tolerance = 2e-12)
  expect_equal(as.double(dual$ksum.power2), colSums(raw^2),
               tolerance = 2e-11)
})

test_that("beta dual powers honor empirical range resolution", {
  dual_sum <- getFromNamespace(".npksum_power12", "np")
  training <- data.frame(x = c(-2.1, -1.4, -.2, .8, 2.7, 4.3))
  evaluation <- data.frame(x = c(-1.7, .4, 3.6))
  arguments <- list(
    bws = .5, txdat = training, exdat = evaluation,
    ckertype = "beta", ckerorder = 6L, ckerbound = "range"
  )

  dual <- do.call(dual_sum, arguments)
  power_one <- do.call(npksum, arguments)
  power_two <- do.call(npksum, c(arguments, list(kernel.pow = 2L)))
  expect_identical(dual$ksum, power_one$ksum)
  expect_beta_dual_square_equivalent(dual$ksum.power2, power_two$ksum)
})

test_that("beta dual powers share the canonical mixed categorical row", {
  old <- options(np.messages = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  dual_sum <- getFromNamespace(".npksum_power12", "np")
  training <- data.frame(
    x = seq(.01, .99, length.out = 24L),
    u = factor(rep(letters[1:3], each = 8L)),
    o = ordered(rep(1:4, length.out = 24L), levels = 1:4)
  )
  evaluation <- training[c(2L, 7L, 13L, 20L), , drop = FALSE]

  for (mode in c("fixed", "generalized_nn", "adaptive_nn")) {
    continuous_bandwidth <- if (identical(mode, "fixed")) .14 else 12
    for (order in c(2L, 4L, 6L, 8L)) {
      bandwidth <- np:::kbandwidth(
        bw = c(continuous_bandwidth, .22, .31),
        xdati = np:::untangle(training),
        xnames = names(training),
        bwtype = mode,
        ckertype = "gaussian"
      )
      bandwidth[["ckertype"]] <- "beta"
      bandwidth[["ckerorder"]] <- order
      bandwidth[["ckerbound"]] <- "fixed"
      bandwidth[["ckerlb"]][bandwidth[["icon"]]] <- 0
      bandwidth[["ckerub"]][bandwidth[["icon"]]] <- 1
      arguments <- list(
        bws = bandwidth,
        txdat = training,
        exdat = evaluation
      )

      dense <- do.call(dual_sum, arguments)
      power_one <- do.call(npksum, arguments)
      power_two <- do.call(
        npksum, c(arguments, list(kernel.pow = 2L))
      )
      options(np.categorical.compress = TRUE)
      compressed <- do.call(dual_sum, arguments)
      options(np.categorical.compress = FALSE)

      expect_identical(dense$ksum, power_one$ksum)
      expect_beta_dual_square_equivalent(
        dense$ksum.power2, power_two$ksum
      )
      expect_identical(compressed$ksum, dense$ksum)
      expect_identical(compressed$ksum.power2, dense$ksum.power2)
    }
  }
})
