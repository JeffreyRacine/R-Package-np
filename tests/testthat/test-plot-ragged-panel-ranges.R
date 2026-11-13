test_that("ragged panel ranges accept only declared trailing NA padding", {
  all.range <- getFromNamespace("compute.all.error.range", "npRmpi")
  default.range <- getFromNamespace("compute.default.error.range", "npRmpi")
  for (n in 1:5) {
    centre <- seq_len(n) / 7
    e <- cbind(seq_len(n) / 10, seq_len(n) / 8, centre + .03)
    families <- list(pointwise = e, simultaneous = 1.2 * e, bonferroni = 1.4 * e)
    expected <- range(unlist(lapply(families, function(x)
      c(centre - x[, 1L], centre + x[, 2L]))))
    expect_identical(all.range(centre, families), expected)
    expect_identical(all.range(c(centre, NA_real_, NA_real_), families), expected)
    expect_identical(default.range(c(centre, NA_real_), e), default.range(centre, e))
    expect_error(all.range(c(centre, 7), families), "rows do not match")
    bad <- families; bad$simultaneous <- rbind(e, 1)
    expect_error(all.range(centre, bad), "incompatible row counts")
  }
})

test_that("partially linear all-band ragged panels draw without recycling", {
  withr::local_options(np.messages = FALSE)
  set.seed(19345)
  d <- data.frame(x = runif(36), z = rnorm(36), u = factor(rep(letters[1:2], 18)),
                  o = ordered(rep(1:3, 12)))
  d$y <- d$x + .3 * d$z + rnorm(36)
  b <- npplregbw(y ~ z | x + u + o, data = d,
    bws = matrix(c(.6, .6, .2, .2, .3, .3), 2, 3), bandwidth.compute = FALSE)
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  for (common in c(FALSE, TRUE)) {
    set.seed(19346)
    expect_warning(plot(b, errors = "bootstrap", band = "all", B = 199L,
      common.scale = common, data.overlay = FALSE, neval = 4L, perspective = FALSE), NA)
  }
})

test_that("all-band rendering uses the same unpadded interval rows as range calculation", {
  draw <- getFromNamespace("draw.all.error.types", "npRmpi")
  rows <- list()
  testthat::local_mocked_bindings(draw.errors = function(ex, ely, ehy, ...) {
    rows[[length(rows)+1L]] <<- list(ex = ex, lower = ely, upper = ehy)
  }, .package = "npRmpi")
  e <- cbind(c(.1,.2), c(.2,.3), c(.3,.5))
  bands <- list(pointwise=e, simultaneous=e*1.2, bonferroni=e*1.4)
  draw(1:2, c(.3,.5), bands, xi.factor=TRUE, add.legend=FALSE,
       plot.errors.method="bootstrap")
  expected <- rows; rows <- list()
  expect_warning(draw(1:2,c(.3,.5,NA,NA,NA),bands,xi.factor=TRUE,
    add.legend=FALSE,plot.errors.method="bootstrap"),NA)
  expect_identical(rows,expected)
  expect_error(draw(1:2,c(.3,.5,9),bands,xi.factor=TRUE,
    add.legend=FALSE,plot.errors.method="bootstrap"),"rows do not match")
})
