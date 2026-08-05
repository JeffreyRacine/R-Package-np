capture_np_error <- function(expr) {
  tryCatch(force(expr), error = identity)
}

test_that("conditional beta errors name the failing public kernel side", {
  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  y <- data.frame(y = seq(0.08, 0.92, length.out = 24L))
  common <- list(
    xdat = x, ydat = y, bws = c(0.2, 0.2),
    bandwidth.compute = FALSE, regtype = "lp", degree = 1L
  )

  x.invalid <- c(common, list(
    cxkertype = "beta", cykertype = "beta", cykerbound = "range"
  ))
  y.invalid <- c(common, list(
    cxkertype = "gaussian", cykertype = "beta", cxkerbound = "range"
  ))

  expected.x <- paste(
    "conditional density explanatory beta kernel requires",
    "cxkerbound = \"fixed\" or \"range\" with finite cxkerlb and cxkerub"
  )
  expected.y <- paste(
    "conditional distribution dependent beta kernel requires",
    "cykerbound = \"fixed\" or \"range\" with finite cykerlb and cykerub"
  )

  x.condition <- capture_np_error(do.call(npcdensbw, x.invalid))
  y.condition <- capture_np_error(do.call(npcdistbw, y.invalid))

  expect_s3_class(x.condition, "simpleError")
  expect_s3_class(y.condition, "simpleError")
  expect_identical(conditionMessage(x.condition), expected.x)
  expect_identical(conditionMessage(y.condition), expected.y)
  expect_null(conditionCall(x.condition))
  expect_null(conditionCall(y.condition))
})

test_that("conditional beta validation retains X-before-Y failure precedence", {
  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  y <- data.frame(y = seq(0.08, 0.92, length.out = 24L))

  condition <- capture_np_error(npcdensbw(
    xdat = x, ydat = y, bws = c(0.2, 0.2),
    bandwidth.compute = FALSE, regtype = "lp", degree = 1L,
    cxkertype = "beta", cykertype = "beta"
  ))

  expect_identical(
    conditionMessage(condition),
    paste(
      "conditional density explanatory beta kernel requires",
      "cxkerbound = \"fixed\" or \"range\" with finite cxkerlb and cxkerub"
    )
  )
})
