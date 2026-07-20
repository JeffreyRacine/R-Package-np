test_that("plot geometry normalizes non-finite coordinates without mutating inputs", {
  input <- matrix(
    c(-Inf, 1, NA_real_, NaN, Inf, 2),
    nrow = 2L,
    dimnames = list(c("a", "b"), c("x", "y", "z"))
  )
  original <- input

  geometry <- .np_plot_geometry_values(input)

  expect_identical(input, original)
  expect_identical(dim(geometry), dim(input))
  expect_identical(dimnames(geometry), dimnames(input))
  expect_true(all(is.na(geometry[!is.finite(input)])))
  expect_identical(geometry[is.finite(input)], input[is.finite(input)])
  expect_identical(
    .np_plot_finite_range(c(-Inf, NA_real_, 2, 5, Inf)),
    c(2, 5)
  )
  expect_identical(
    .np_plot_finite_range(c(-Inf, NA_real_, NaN, Inf)),
    c(-1, 1)
  )
})

test_that("beta endpoint infinities remain in returns but not plot geometry", {
  training <- data.frame(x = c(0, .12, .3, .55, .82, 1))
  response <- c(4, 1, 2, 3, 5, -1)
  bw <- npregbw(
    xdat = training, ydat = response,
    bws = .16, bandwidth.compute = FALSE, regtype = "lc",
    ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )

  direct <- NULL
  expect_warning(
    direct <- npreg(
      bws = bw, txdat = training, tydat = response,
      exdat = data.frame(x = c(0, .5, 1)), gradients = TRUE
    ),
    "infinite endpoint"
  )
  expect_true(all(is.infinite(gradients(direct)[c(1L, 3L), 1L])))
  expect_true(is.finite(gradients(direct)[2L, 1L]))

  plot.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot.file)
  on.exit({
    grDevices::dev.off()
    unlink(plot.file)
  }, add = TRUE)

  for (common in c(FALSE, TRUE)) {
    plotted <- NULL
    expect_warning(
      plotted <- plot(
        bw, xdat = training, ydat = response,
        gradients = TRUE, errors = "none", output = "plot-data",
        common_scale = common, neval = 11L, xtrim = 0
      ),
      "infinite endpoint"
    )
    returned <- as.double(gradients(plotted[[1L]]))
    expect_true(is.infinite(returned[1L]))
    expect_true(is.infinite(returned[length(returned)]))
    expect_true(all(is.finite(returned[-c(1L, length(returned))])))
  }

  for (common in c(FALSE, TRUE)) {
    plotted <- NULL
    expect_warning(
      plotted <- plot(
        bw, xdat = training, ydat = response,
        gradients = TRUE, errors = "asymptotic", output = "plot-data",
        common_scale = common, neval = 2L, xtrim = 0
      ),
      "infinite endpoint"
    )
    expect_true(all(is.infinite(as.double(gradients(plotted[[1L]])))))
  }
})

test_that("non-finite interval geometry is omitted rather than sent to graphics", {
  plot.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot.file)
  on.exit({
    grDevices::dev.off()
    unlink(plot.file)
  }, add = TRUE)
  graphics::plot(1, 1, type = "n", xlim = c(0, 4), ylim = c(-2, 2))

  expect_silent(draw.errors(
    ex = c(1, 2, 3, 4),
    ely = c(-Inf, -1, NA_real_, 0),
    ehy = c(Inf, 1, 2, NaN),
    plot.errors.style = "band",
    plot.errors.bar = "|",
    plot.errors.bar.num = 4L,
    lty = 2L
  ))
})

test_that("conditional beta endpoint gradients plot and retain infinities", {
  training.x <- data.frame(x = c(0, .12, .3, .55, .82, 1))
  training.y <- data.frame(y = c(.05, .21, .37, .62, .79, .96))
  plot.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot.file)
  on.exit({
    grDevices::dev.off()
    unlink(plot.file)
  }, add = TRUE)

  for (distribution in c(FALSE, TRUE)) {
    bwfun <- if (distribution) npcdistbw else npcdensbw
    bw <- bwfun(
      xdat = training.x, ydat = training.y,
      bws = c(.16, .18), bandwidth.compute = FALSE,
      cxkertype = "beta", cxkerorder = 8,
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykertype = "beta", cykerorder = 8,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1
    )

    plotted <- NULL
    expect_warning(
      plotted <- plot(
        bw, xdat = training.x, ydat = training.y,
        gradients = TRUE, errors = "none", output = "plot-data",
        common_scale = TRUE, neval = 11L, xtrim = 0, ytrim = 0
      ),
      "infinite endpoint"
    )
    returned <- as.double(gradients(plotted[[1L]]))
    expect_true(is.infinite(returned[1L]))
    expect_true(is.infinite(returned[length(returned)]))
    expect_true(all(is.finite(returned[-c(1L, length(returned))])))
  }
})
