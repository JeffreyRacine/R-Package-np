test_that("density-family ordered mass survives operator conversion", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  x <- data.frame(o = ordered(rep(0:2, 4), levels = 0:2))
  e <- x[c(1, 2), , drop = FALSE]
  y <- seq_len(nrow(x))^2
  delta <- outer(as.numeric(e$o), as.numeric(x$o), "-")
  mass <- .3^abs(delta) * .7 / 1.3 / nrow(x)
  cdf <- ifelse(delta < 0, .3^abs(delta) / 1.3,
                1 - .3^(abs(delta) + 1) / 1.3) / nrow(x)
  for (target in c("density", "distribution")) {
    bw <- if (target == "density")
      npudensbw(dat = x, bws = .3, bandwidth.compute = FALSE) else
      npudistbw(dat = x, bws = .3, bandwidth.compute = FALSE)
    hat <- if (target == "density") npudenshat else npudisthat
    expected <- if (target == "density") mass else cdf
    expect_equal(unname(as.matrix(hat(bw, x, e)))[, ], expected,
                 ignore_attr = TRUE, tolerance = 1e-12)
    expect_equal(hat(bw, x, e, y = y, output = "apply"),
                 drop(expected %*% y), tolerance = 1e-12)
    expect_identical(kbandwidth(bw)$okertype, "nliracine")
  }
  raw <- kbandwidth.numeric(.3, xdati = untangle(x), xnames = "o")
  expect_identical(.np_kbandwidth_okertype(raw), "liracine")
  # Copula marginal adapters deliberately receive plain density-role lists.
  expect_identical(.np_make_kbandwidth_unconditional(
    list(bw = .3, type = "fixed", ckertype = "gaussian", ckerorder = 2,
         ckerbound = "none", ukertype = "aitchisonaitken", okertype = "liracine"),
    x)$okertype, "nliracine")
})
