with_plot_state <- function(expr) {
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  before <- list(
    usr = par("usr"),
    mfrow = par("mfrow"),
    mar = par("mar"),
    new = par("new")
  )

  value <- suppressWarnings(force(expr))

  after <- list(
    usr = par("usr"),
    mfrow = par("mfrow"),
    mar = par("mar"),
    new = par("new")
  )

  list(value = value, before = before, after = after)
}

expect_plot_device_clean <- function(state, expect_mar_restore = FALSE, info = NULL) {
  expect_false(isTRUE(state$after$new), info = info)
  expect_equal(state$after$mfrow, state$before$mfrow, info = info)
  if (isTRUE(expect_mar_restore))
    expect_equal(state$after$mar, state$before$mar, info = info)
  expect_false(isTRUE(all.equal(state$after$usr, state$before$usr)), info = info)
}

test_that("plot methods keep device coordinates active while restoring layout state", {
  set.seed(20260310)
  n <- 60
  x <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + 0.3 * z + rnorm(n, sd = 0.08)

  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)
  idxdat <- data.frame(y = y, x = x, x2 = x2)

  rfit <- npreg(
    txdat = xdat,
    tydat = y,
    bws = npregbw(xdat = xdat, ydat = y, bws = 0.25, bandwidth.compute = FALSE)
  )
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

  cases <- list(
    list(
      label = "npregression fit",
      expr = quote(plot(rfit))
    ),
    list(
      label = "npudensbw",
      expr = quote(plot(ubw, xdat = xdat))
    ),
    list(
      label = "npudistbw",
      expr = quote(plot(dbw, xdat = xdat))
    ),
    list(
      label = "npcdensbw",
      expr = quote(plot(cbw, xdat = xdat, ydat = ydat, perspective = FALSE))
    ),
    list(
      label = "npcdistbw",
      expr = quote(plot(cdbw, xdat = xdat, ydat = ydat, perspective = FALSE))
    ),
    list(
      label = "npplregbw",
      expr = quote(plot(pbw, xdat = data.frame(x = z), ydat = y, zdat = data.frame(z = x), perspective = FALSE))
    ),
    list(
      label = "npindexbw",
      expr = quote(plot(sbw, xdat = data.frame(x = x, x2 = x2), ydat = y, perspective = FALSE))
    ),
    list(
      label = "npscoefbw",
      expr = quote(plot(scbw, xdat = data.frame(x = x), ydat = y, zdat = data.frame(z = z), perspective = FALSE))
    )
  )

  for (case in cases) {
    state <- with_plot_state(eval(case$expr, envir = environment()))
    expect_plot_device_clean(state, info = case$label)
  }
})

test_that("np.pairs restores layout parameters without resetting plot coordinates", {
  set.seed(20260310 + 1L)
  dat <- data.frame(a = rnorm(40), b = runif(40), c = rnorm(40))

  pair_list <- np.pairs(y_vars = c("a", "b", "c"), y_dat = dat)
  state <- with_plot_state(np.pairs.plot(pair_list))

  expect_plot_device_clean(state, expect_mar_restore = TRUE, info = "np.pairs")
})
