test_that("complete LP term counts are basis-neutral", {
  fixtures <- list(
    list(basis = "glp", degree = c(2L, 2L), nterms = 6),
    list(basis = "additive", degree = c(2L, 3L), nterms = 6),
    list(basis = "tensor", degree = c(1L, 2L), nterms = 6),
    list(basis = "glp", degree = c(4L, 4L, 4L), nterms = 35),
    list(basis = "additive", degree = rep.int(10L, 10L), nterms = 101),
    list(basis = "tensor", degree = rep.int(4L, 5L), nterms = 3125)
  )

  for (fixture in fixtures) {
    actual <- np:::npLpCompleteTermStatus(
      fixture$basis,
      fixture$degree,
      nobs = fixture$nterms + 1
    )
    expect_identical(actual$status, "valid")
    expect_identical(as.double(actual$nterms), as.double(fixture$nterms))
    expect_true(actual$capacity.ok)
    expect_true(actual$rank.ok)
  }
})

test_that("LP statistical rank and native capacity are distinct", {
  degree <- list(
    glp = c(2L, 2L),
    additive = c(2L, 3L),
    tensor = c(1L, 2L)
  )
  for (basis in names(degree)) {
    expect_identical(
      np:::npLpCompleteTermStatus(basis, degree[[basis]], 5L)$status,
      "rank_inadmissible"
    )
    expect_identical(
      np:::npLpCompleteTermStatus(basis, degree[[basis]], 6L)$status,
      "rank_inadmissible"
    )
    expect_identical(
      np:::npLpCompleteTermStatus(basis, degree[[basis]], 7L)$status,
      "valid"
    )
  }

  expect_identical(
    np:::npLpCompleteTermStatus(
      "additive", rep.int(25L, 4001L), 200000L
    )$status,
    "capacity_exceeded"
  )
  expect_identical(
    np:::npLpCompleteTermStatus(
      "tensor", rep.int(10L, 5L), 200000L
    )$status,
    "capacity_exceeded"
  )
})

test_that("public LP bandwidth families reject nterms equal to n", {
  x <- data.frame(x = c(0, 1))
  y <- c(0, 1)
  admission <- "LP basis dimension \\(2\\) exceeds nobs - 1 \\(1\\)"

  expect_error(
    npregbw(
      xdat = x, ydat = y, regtype = "lp", degree = 1L,
      degree.select = "manual", bandwidth.compute = FALSE, bws = 0.5
    ),
    admission
  )
  expect_error(
    npcdensbw(
      xdat = x, ydat = y, regtype = "lp", degree = 1L,
      degree.select = "manual", bandwidth.compute = FALSE,
      bws = c(0.5, 0.5)
    ),
    admission
  )
  expect_error(
    npcdistbw(
      xdat = x, ydat = y, regtype = "lp", degree = 1L,
      degree.select = "manual", bandwidth.compute = FALSE,
      bws = c(0.5, 0.5)
    ),
    admission
  )
  expect_error(
    nplsqregbw(
      xdat = x, ydat = y, regtype = "lp", degree = 1L,
      degree.select = "manual", bandwidth.compute = FALSE,
      bws = 0.5, tau = 0.5
    ),
    admission
  )
})
