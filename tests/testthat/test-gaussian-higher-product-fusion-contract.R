library(np)

test_that("higher-order Gaussian product fusion is narrow and bounded", {
  src <- readLines(np_test_source_path("src", "jksum.c"), warn = FALSE)
  text <- paste(src, collapse = "\n")

  expect_match(
    text, "static int NP_NOINLINE np_accel_gauss_higher_product_try",
    fixed = TRUE
  )
  expect_match(
    text,
    "(kernel[d] < 1) || (kernel[d] > 3)",
    fixed = TRUE
  )
  expect_match(
    text,
    "np_accel_gauss_scratch_ensure(num_xt)",
    fixed = TRUE
  )
  expect_match(
    text,
    "static int NP_NOINLINE np_accel_gauss_product_kind",
    fixed = TRUE
  )
  expect_match(
    text,
    "else if((fused_gaussian_product_eligible == 2) && (!any_cont_largeh))",
    fixed = TRUE
  )
  expect_match(
    text,
    paste0(
      "np_accel_gauss_higher_product_try(\n",
      "            KERNEL_reg_np, xtc, xc, m, num_reg_continuous"
    ),
    fixed = TRUE
  )
  expect_match(
    text,
    "np_ckernelv(KERNEL_reg_np[i]",
    fixed = TRUE
  )
})

test_that("higher-order Gaussian fusion matches independent formulas", {
  old_options <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.largeh = FALSE,
    np.macMseries.accelerate = TRUE
  )
  on.exit(options(old_options), add = TRUE)

  gaussian_higher <- function(z, order) {
    z2 <- z * z
    polynomial <- switch(
      as.character(order),
      `4` = 1.5 - 0.5 * z2,
      `6` = 1.875 + z2 * (z2 * 0.125 - 1.25),
      `8` = 2.1875 + z2 *
        (-2.1875 + z2 * (0.4375 - z2 * 0.02083333333))
    )
    (1 / sqrt(2 * pi)) * polynomial * exp(-0.5 * z2)
  }
  adaptive_radius <- function(train, k) {
    vapply(seq_along(train), function(i) {
      distance <- abs(train - train[[i]])
      positive <- sort(distance[distance > 0])
      duplicate_count <- sum(distance == 0) - 1L
      if (duplicate_count >= k) positive[[1L]] else
        positive[[max(1L, k - duplicate_count)]]
    }, numeric(1L))
  }
  generalized_radius <- function(train, eval, k) {
    vapply(eval, function(value) sort(abs(train - value))[[k]], numeric(1L))
  }

  set.seed(2026073119L)
  ntrain <- 260L
  neval <- 256L
  xtrain <- data.frame(
    x1 = runif(ntrain, -0.9, 0.9),
    x2 = runif(ntrain, -0.8, 0.85)
  )
  xeval <- data.frame(
    x1 = runif(neval, -0.87, 0.88),
    x2 = runif(neval, -0.78, 0.82)
  )
  fixed_h <- c(0.31, 0.37)
  k <- c(79L, 83L)

  for (order in c(4L, 6L, 8L)) {
    for (topology in c("fixed", "generalized_nn", "adaptive_nn")) {
      bws <- if (topology == "fixed") fixed_h else k
      bandwidth <- switch(
        topology,
        fixed = lapply(fixed_h, rep, times = neval),
        generalized_nn = Map(generalized_radius, xtrain, xeval, k),
        adaptive_nn = Map(adaptive_radius, xtrain, k)
      )
      oracle <- matrix(1, nrow = ntrain, ncol = neval)
      denominator <- matrix(1, nrow = ntrain, ncol = neval)
      for (d in 1:2) {
        hmat <- if (topology == "adaptive_nn") {
          matrix(bandwidth[[d]], nrow = ntrain, ncol = neval)
        } else {
          matrix(bandwidth[[d]], nrow = ntrain, ncol = neval, byrow = TRUE)
        }
        z <- outer(
          xtrain[[d]], xeval[[d]], function(train, eval) eval - train
        ) / hmat
        oracle <- oracle * gaussian_higher(z, order)
        denominator <- denominator * hmat
      }

      raw <- npksum(
        bws = bws,
        txdat = xtrain,
        exdat = xeval,
        bwtype = topology,
        ckertype = "gaussian",
        ckerorder = order,
        return.kernel.weights = TRUE
      )
      divided <- npksum(
        bws = bws,
        txdat = xtrain,
        exdat = xeval,
        bwtype = topology,
        ckertype = "gaussian",
        ckerorder = order,
        bandwidth.divide = TRUE,
        return.kernel.weights = TRUE
      )

      expect_equal(as.matrix(raw$kw), oracle, tolerance = 2e-12)
      expect_equal(
        as.numeric(raw$ksum), colSums(oracle), tolerance = 2e-11
      )
      expect_equal(
        as.numeric(divided$ksum),
        colSums(oracle / denominator),
        tolerance = 2e-10
      )
      if (topology == "adaptive_nn") {
        expect_equal(
          as.matrix(divided$kw), oracle / denominator, tolerance = 2e-10
        )
      }
    }
  }
})
