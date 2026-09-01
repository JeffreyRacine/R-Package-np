.t4_ratio_package <- function() {
  if ("package:npRmpi" %in% search()) "npRmpi" else "np"
}

.t4_ratio_owner <- function(Kx, Ky) {
  package <- .t4_ratio_package()
  owner <- getFromNamespace(".npcdhat_ratio_matrix", package)
  testthat::with_mocked_bindings(
    owner(
      bws = list(),
      txdat = data.frame(x = seq_len(ncol(Kx))),
      tydat = data.frame(y = seq_len(ncol(Kx))),
      exdat = data.frame(x = seq_len(nrow(Kx))),
      eydat = data.frame(y = seq_len(nrow(Kx))),
      operator = "normal"
    ),
    .npcdhat_make_xkbw = function(...) list(owner = "x"),
    .npcdhat_make_ybw = function(...) list(owner = "y"),
    .npcdhat_make_kernel_matrix = function(kbw, ...) {
      if (identical(kbw$owner, "x")) Kx else Ky
    },
    .package = package
  )
}

test_that("conditional-hat ratios use every finite nonzero denominator", {
  small <- 2^-60
  subnormal <- 2^-1072
  Kx <- rbind(
    ordinary = c(2, 0),
    positive.small = c(small, small),
    negative.small = c(-small, -small),
    positive.subnormal = c(subnormal, subnormal),
    exact.zero = c(0, 0),
    missing = c(NA_real_, 1),
    positive.infinity = c(Inf, 1),
    negative.infinity = c(-Inf, -1)
  )
  Ky <- matrix(
    c(4, 8, 8, 16, 8, 16, 8, 16,
      4, 8, 4, 8, 4, 8, 4, 8),
    nrow = nrow(Kx), byrow = FALSE
  )
  denom <- rowSums(Kx) / ncol(Kx)
  numerator <- (Kx * Ky) / ncol(Kx)
  finite.nonzero <- is.finite(denom) & denom != 0.0

  mixed <- .t4_ratio_owner(Kx, Ky)
  direct <- sweep(
    numerator[finite.nonzero, , drop = FALSE],
    1L, denom[finite.nonzero], "/"
  )
  incumbent.exceptional <- sweep(
    numerator[!finite.nonzero, , drop = FALSE],
    1L,
    pmax(denom[!finite.nonzero], .Machine$double.eps),
    "/"
  )

  expect_true(is.matrix(mixed))
  expect_identical(
    as.vector(mixed[finite.nonzero, , drop = FALSE]),
    as.vector(direct)
  )
  expect_identical(
    as.vector(mixed[!finite.nonzero, , drop = FALSE]),
    as.vector(incumbent.exceptional)
  )

  fast <- .t4_ratio_owner(
    Kx[finite.nonzero, , drop = FALSE],
    Ky[finite.nonzero, , drop = FALSE]
  )
  expect_identical(as.vector(fast), as.vector(direct))
})

test_that("adaptive conditional hats preserve reachable positive-small support", {
  if (identical(.t4_ratio_package(), "npRmpi")) {
    if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }

  n <- 20L
  x <- data.frame(z = seq(0, 1, length.out = n))
  for (j in 1:10)
    x[[paste0("f", j)]] <- factor(rep("A", n), levels = c("A", "B"))
  y <- data.frame(y = seq(-1, 1, length.out = n))
  ex <- x[1L, , drop = FALSE]
  ex$z <- 0.5
  for (j in 1:10)
    ex[[paste0("f", j)]] <- factor("B", levels = c("A", "B"))

  bw <- npcdistbw(
    xdat = x, ydat = y,
    bws = c(3, 3, rep.int(0.01, 10L)),
    bandwidth.compute = FALSE,
    bwtype = "adaptive_nn",
    cxkertype = "gaussian",
    regtype = "lc"
  )
  package <- .t4_ratio_package()
  make.xbw <- getFromNamespace(".npcdhat_make_xkbw", package)
  make.kernel <- getFromNamespace(".npcdhat_make_kernel_matrix", package)
  xbw <- make.xbw(bws = bw, txdat = x)
  Kx <- make.kernel(
    kbw = xbw, txdat = x, exdat = ex,
    operator = rep.int("normal", ncol(x))
  )
  denominator <- sum(Kx) / n

  H <- npcdisthat(
    bws = bw, txdat = x, tydat = y,
    exdat = ex, eydat = data.frame(y = 10)
  )
  expect_gt(denominator, 0)
  expect_lt(denominator, .Machine$double.eps)
  expect_true(all(is.finite(H)))
  expect_equal(sum(H), 1, tolerance = 64 * .Machine$double.eps)
})
