test_that("plot contract: bootstrap args require explicit bootstrap mode across engine families", {
  skip_if_not_installed("np")

  msg <- "bootstrap controls require errors = \"bootstrap\""

  set.seed(811)
  n <- 60
  x <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + 0.3 * z + rnorm(n, sd = 0.08)

  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)
  idxdat <- data.frame(y = y, x = x, x2 = x2)

  rbw <- npregbw(xdat = xdat, ydat = y, bws = 0.25, bandwidth.compute = FALSE)
  ubw <- npudensbw(dat = xdat, bws = 0.25, bandwidth.compute = FALSE)
  dbw <- npudistbw(dat = xdat, bws = 0.25, bandwidth.compute = FALSE)
  cbw <- npcdensbw(xdat = xdat, ydat = ydat, bws = c(0.30, 0.30), bandwidth.compute = FALSE)
  cdbw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.30, 0.30), bandwidth.compute = FALSE)
  pbw <- npplregbw(
    xdat = data.frame(x = z),
    zdat = data.frame(z = x),
    ydat = y,
    bws = matrix(c(0.30, 0.30), nrow = 2),
    bandwidth.compute = FALSE
  )
  sbw <- npindexbw(
    y ~ x + x2,
    data = idxdat,
    bws = c(1, 0.30, 0.30),
    bandwidth.compute = FALSE,
    method = "ichimura"
  )
  scbw <- npscoefbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    bws = 0.20,
    bandwidth.compute = FALSE
  )

  expect_error(
    suppressWarnings(plot(rbw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(ubw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(dbw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(cbw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(cdbw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      plot(
        pbw,
        xdat = data.frame(x = z),
        ydat = y,
        zdat = data.frame(z = x),
        output = "data",
        perspective = FALSE,
        B = 5
      )
    ),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      plot(
        sbw,
        xdat = data.frame(x = x, x2 = x2),
        ydat = y,
        output = "data",
        perspective = FALSE,
        B = 5
      )
    ),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      plot(
        scbw,
        xdat = data.frame(x = x),
        ydat = y,
        zdat = data.frame(z = z),
        output = "data",
        perspective = FALSE,
        B = 5
      )
    ),
    msg,
    fixed = TRUE
  )
})


test_that("estimate centering is not a bootstrap control for asymptotics", {
  normalize <- getFromNamespace(".np_plot_normalize_public_dots", "np")
  for (key in c("center", "plot.errors.center")) {
    args <- list(errors = "asymptotic")
    args[[key]] <- "estimate"
    expect_identical(normalize(args)$plot.errors.center, "estimate")
    args$B <- 19L
    expect_error(normalize(args), 'bootstrap controls require errors', fixed = TRUE)
  }
  expect_error(normalize(list(errors = "asymptotic", center = "bias-corrected")),
               'bootstrap controls require errors', fixed = TRUE)
  expect_identical(normalize(list(errors = "asymptotic", center = "bias-corrected"),
                             allow_asymptotic_bias = TRUE)$plot.errors.center,
                   "bias-corrected")
})
