test_that("npValidateScalarLogical enforces scalar logical flags", {
  expect_identical(npValidateScalarLogical(TRUE, "x"), TRUE)
  expect_identical(npValidateScalarLogical(1, "x"), TRUE)
  expect_identical(npValidateScalarLogical(0, "x"), FALSE)
  expect_error(npValidateScalarLogical(c(TRUE, FALSE), "x"), "'x' must be TRUE or FALSE")
  expect_error(npValidateScalarLogical("foo", "x"), "'x' must be TRUE or FALSE")
})

test_that("core estimator methods reject non-scalar logical control flags", {
  tx <- data.frame(x = c(0.1, 0.2))
  ty <- c(0.3, 0.4)
  tydf <- data.frame(y = ty)

  expect_error(
    npcdens.conbandwidth(structure(list(), class = "conbandwidth"),
                         txdat = tx, tydat = tydf, gradients = c(TRUE, FALSE)),
    "'gradients' must be TRUE or FALSE"
  )
  expect_error(
    npcdist.condbandwidth(structure(list(), class = "condbandwidth"),
                          txdat = tx, tydat = tydf, gradients = c(TRUE, FALSE)),
    "'gradients' must be TRUE or FALSE"
  )
  expect_error(
    npqreg.condbandwidth(structure(list(), class = "condbandwidth"),
                         txdat = tx, tydat = ty, gradients = c(TRUE, FALSE)),
    "'gradients' must be TRUE or FALSE"
  )
  expect_error(
    npqreg.condbandwidth(structure(list(), class = "condbandwidth"),
                         txdat = tx, tydat = ty, itmax = 0),
    "'itmax' must be a positive integer"
  )
  expect_error(
    npqreg.condbandwidth(structure(list(), class = "condbandwidth"),
                         txdat = tx, tydat = ty, ftol = 0),
    "'ftol' is no longer accepted by npqreg"
  )
  expect_error(
    npindex.sibandwidth(structure(list(), class = "sibandwidth"),
                        txdat = tx, tydat = ty, gradients = c(TRUE, FALSE)),
    "'gradients' must be TRUE or FALSE"
  )
  expect_error(
    npindex.sibandwidth(structure(list(), class = "sibandwidth"),
                        txdat = tx, tydat = ty, boot.num = 0),
    "'boot.num' must be a positive integer"
  )
  expect_error(
    npscoef.scbandwidth(structure(list(), class = "scbandwidth"),
                        txdat = tx, tydat = ty, iterate = c(TRUE, FALSE)),
    "'iterate' must be TRUE or FALSE"
  )
  expect_error(
    npscoef.scbandwidth(structure(list(), class = "scbandwidth"),
                        txdat = tx, tydat = ty, maxiter = 0),
    "'maxiter' must be a positive integer"
  )
  expect_error(
    npscoef.scbandwidth(structure(list(), class = "scbandwidth"),
                        txdat = tx, tydat = ty, tol = -1),
    "'tol' must be a finite numeric scalar >= 0"
  )
  expect_error(
    npplreg.plbandwidth(structure(list(), class = "plbandwidth"),
                        txdat = tx, tydat = ty, tzdat = tx, residuals = c(TRUE, FALSE)),
    "'residuals' must be TRUE or FALSE"
  )
})
