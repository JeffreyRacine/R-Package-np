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

test_that("npscoef batch projection matches incumbent row crossproducts numerically", {
  set.seed(163L)
  for (nbasis in c(1L, 3L, 6L, 9L, 15L, 26L)) {
    neval <- 37L
    pcoef <- 3L
    theta <- matrix(
      rnorm(nbasis * pcoef * neval),
      nrow = nbasis * pcoef,
      ncol = neval
    )
    Wz.eval <- matrix(
      runif(neval * nbasis, -1, 1),
      nrow = neval,
      ncol = nbasis
    )
    expected <- vapply(seq_len(neval), function(ii) {
      as.vector(crossprod(
        Wz.eval[ii, ],
        matrix(theta[, ii], nrow = nbasis, ncol = pcoef)
      ))
    }, numeric(pcoef))
    actual <- npRmpi:::.npscoef_batch_project(theta, Wz.eval)
    expect_equal(actual, expected, tolerance = 1e-14)
    expect_lte(max(abs(actual - expected)), 1e-14)
  }
})

test_that("npscoef batch projection rejects incompatible internal shapes", {
  expect_error(
    .Call(
      "C_np_npscoef_batch_project",
      matrix(1.0, nrow = 5L, ncol = 3L),
      matrix(1.0, nrow = 3L, ncol = 2L),
      PACKAGE = "npRmpi"
    ),
    "incompatible dimensions"
  )
})

test_that("npscoef native dimension attributes remain protected", {
  candidates <- unique(c(
    test_path("..", "..", "src", "npscoef_batch_solve.c"),
    test_path("..", "..", "..", "src", "npscoef_batch_solve.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src",
              "npscoef_batch_solve.c"),
    file.path(getwd(), "src", "npscoef_batch_solve.c"),
    file.path(getwd(), "..", "src", "npscoef_batch_solve.c")
  ))
  source_file <- candidates[nzchar(candidates) & file.exists(candidates)][1L]
  skip_if(is.na(source_file), "source file unavailable")
  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "tww_dim = PROTECT(getAttrib(tww_r, R_DimSymbol));",
    fixed = TRUE
  )
  expect_match(
    source,
    "tyw_dim = PROTECT(getAttrib(tyw_r, R_DimSymbol));",
    fixed = TRUE
  )
  expect_match(
    source,
    "theta_dim = PROTECT(getAttrib(theta_r, R_DimSymbol));",
    fixed = TRUE
  )
  expect_match(
    source,
    "wz_dim = PROTECT(getAttrib(wz_r, R_DimSymbol));",
    fixed = TRUE
  )
})
