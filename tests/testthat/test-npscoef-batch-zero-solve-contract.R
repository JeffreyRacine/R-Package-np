test_that("npscoef native zero-ridge batch solve matches base LAPACK solves", {
  p <- 4L
  neval <- 5L
  tww <- array(0.0, dim = c(p, p, neval))
  tyw <- matrix(0.0, nrow = p, ncol = neval)

  for (ii in seq_len(neval)) {
    design <- matrix(
      sin(seq_len((p + 2L) * p) + ii / 7),
      nrow = p + 2L,
      ncol = p
    )
    tww[, , ii] <- crossprod(design) + diag(0.25 + ii / 20, p)
    tyw[, ii] <- cos(seq_len(p) + ii / 11)
  }

  actual <- .Call(
    "C_np_npscoef_batch_zero_solve",
    tww,
    tyw,
    PACKAGE = "npRmpi"
  )
  expected <- vapply(
    seq_len(neval),
    function(ii) solve(tww[, , ii], tyw[, ii]),
    numeric(p)
  )

  expect_equal(actual, expected, tolerance = 1e-13)
})

test_that("npscoef native zero-ridge batch solve returns difficult rows to R", {
  p <- 4L
  neval <- 3L
  tww <- array(rep(diag(p), neval), dim = c(p, p, neval))
  tyw <- matrix(as.double(seq_len(p * neval)), nrow = p, ncol = neval)

  singular <- tww
  singular[, , 2L] <- 0
  expect_null(.Call(
    "C_np_npscoef_batch_zero_solve",
    singular,
    tyw,
    PACKAGE = "npRmpi"
  ))

  nonfinite <- tww
  nonfinite[1L, 1L, 2L] <- NaN
  expect_null(.Call(
    "C_np_npscoef_batch_zero_solve",
    nonfinite,
    tyw,
    PACKAGE = "npRmpi"
  ))

  ill_conditioned <- tww
  ill_conditioned[, , 2L] <- diag(c(1, rep(.Machine$double.eps^2, p - 1L)))
  expect_null(.Call(
    "C_np_npscoef_batch_zero_solve",
    ill_conditioned,
    tyw,
    PACKAGE = "npRmpi"
  ))
})

test_that("npscoef width-one and width-two systems cannot enter batch solve", {
  for (p in 1:2) {
    tww <- array(diag(p), dim = c(p, p, 1L))
    tyw <- matrix(as.double(seq_len(p)), nrow = p, ncol = 1L)
    expect_null(npRmpi:::.npscoef_batch_zero_solve(tyw = tyw, tww = tww))
    expect_error(
      .Call(
        "C_np_npscoef_batch_zero_solve",
        tww,
        tyw,
        PACKAGE = "npRmpi"
      ),
      "incompatible dimensions"
    )
  }
})
