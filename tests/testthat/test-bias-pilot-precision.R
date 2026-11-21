test_that("P1 higher derivatives reproduce moments in concentrated weights", {
  u <- seq(-1,1,length.out=60L)
  Q <- outer(u,0:5,`^`)
  at <- .9; h <- .15
  weight <- exp(-(u-at)^2/(2*h^2))
  ratios <- np:::.np_bias_gaussian_derivatives((u-at)/h,h,3L)
  values <- t(vapply(0:3,function(j) vapply(0:5,function(k)
    if(k<j)0 else factorial(k)/factorial(k-j)*at^(k-j),numeric(1)),numeric(6)))
  row <- np:::.np_bias_pilot_higher(Q,weight,values,ratios)
  expect_equal(drop(row %*% Q),values[4L,],tolerance=1e-9)
  expect_error(np:::.np_bias_pilot_higher(cbind(Q,Q[,2L]),weight,
    cbind(values,values[,2L]),ratios),"numerically unidentified",fixed=TRUE)
  # A fixed penalty enters the zeroth Gram matrix only. Its derivative is zero.
  P <- Q[,1:3]; V <- values[,1:3]
  ordinary <- np:::.np_bias_gaussian_derivatives((u-at)/.8,.8,3L)
  a <- np:::.np_bias_pilot_higher(P,rep(1,length(u)),V,ordinary,
    penalty=rep(.1,3))
  b <- np:::.np_bias_weighted_derivative(P,rep(1,length(u)),V,ordinary,
    penalty=rep(.1,3))
  expect_equal(a,b,tolerance=1e-10)
})

test_that("univariate P1 training and first derivative rows agree at both edges", {
  x <- data.frame(x = seq(0, 1, length.out = 96L))
  n <- nrow(x)
  b <- 0.06
  bw <- npregbw(xdat = x, ydat = x$x, bws = b, bandwidth.compute = FALSE)
  pilot <- list(x = x, xcon = x, numeric.cols = 1L, d = 1L, q = 5L,
                factor.cols = integer(), groups = list(all = seq_len(n)),
                b = b, bernstein = TRUE,
                pilot.bws = .np_bias_regbw(bw, x, b, 5L, TRUE))
  indices <- c(1L, n)
  means <- .np_bias_pilot_owner(pilot)(indices)
  for (j in seq_along(indices)) {
    at <- x[indices[j], , drop = FALSE]
    u <- (x$x - at$x) / b
    Q <- outer(u, 0:5, "^")
    root <- exp(-u^2 / 4)
    s <- svd(Q * root)
    B <- sweep(s$v %*% (t(s$u) / s$d), 2L, root, "*")
    w <- B[1L, ] * u / b
    expected <- B[2L, ] / b + w - drop((w %*% Q) %*% B)
    target <- .np_bias_pilot_rows(pilot, at, derivatives = TRUE)
    expect_lt(max(abs(means[j, ] - B[1L, ])), 1e-8)
    expect_lt(max(abs(target$derivatives[[1L]] - expected)), 1e-8)
  }
  # The incumbent solver decides identifiability before the QR calculation.
  Q <- outer(x$x, 0:5, "^")
  expect_error(.np_bias_pilot_qr(Q, as.double(seq_len(n) == 1L), Q[1L, ], x$x),
               "requires regularization", fixed = TRUE)
})

test_that("univariate P1 preserves the canonical solve on tied numerical support", {
  x <- data.frame(x = c(rep(c(.11, .12), each = 10L), 1:12))
  b <- .005
  bw <- npregbw(xdat = x, ydat = x$x, bws = b, bandwidth.compute = FALSE)
  kb <- .np_bias_regbw(bw, x, b, 2L, TRUE)
  p <- list(x = x, xcon = x, groups = list(all = seq_len(nrow(x))), numeric.cols = 1L,
            factor.cols = integer(), d = 1L, q = 2L, b = b, bernstein = TRUE, pilot.bws = kb)
  at <- x[1L, , drop = FALSE]
  u <- data.frame(x = (x$x - at$x) / b)
  k <- drop(.np_kernel_weights_direct(kb, x, at, bandwidth.divide = TRUE))
  Q <- W.lp(u, degree = 2L, basis = "glp", bernstein.basis = TRUE)
  v <- drop(W.lp(u, exdat = data.frame(x = 0), degree = 2L, basis = "glp", bernstein.basis = TRUE))
  canonical <- tryCatch(.np_bias_solve(Q, k, v, "P1 auxiliary pilot"), error = identity)
  if (inherits(canonical, "error")) {
    expect_error(.np_bias_pilot_rows(p, at, FALSE), conditionMessage(canonical), fixed = TRUE)
  } else {
    expect_identical(.np_bias_pilot_rows(p, at, FALSE)$mean, canonical)
  }
})
