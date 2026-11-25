test_that("compact CDF rows, signed range and process admission agree with direct integration", {
  x <- c(-.9, -.4, -.15, .22, .57, .82)
  g <- c(-2, -.38, .29, 2)
  h <- .25
  state <- list(type = "fixed", x = data.frame(x = x), n = length(x),
    q = 1L, s = 0L, numeric = TRUE, h = h, category = list(), max.bytes = 256*1024^2,
    metadata = list(kernel = "epanechnikov2"))
  grid <- data.frame(x = g)
  rows <- .np_udist_bc_rows(state, grid)
  t <- outer(g, x, "-") / (h*sqrt(5))
  kernel <- function(u) (45-150*u^2+105*u^4)/32
  primitive <- function(u) {
    if (u <= -1) return(0)
    if (u >= 1) return(1)
    integrate(kernel, -1, u)$value
  }
  reference <- matrix(vapply(t, primitive, 0.0), nrow(grid))
  expect_equal(rows$corrected, reference, tolerance = 1e-13)
  expect_equal(rows$center, rowMeans(reference), tolerance = 1e-13)
  centered <- sweep(reference, 1L, rowMeans(reference), "-") / length(x)
  expect_equal(tcrossprod(rows$factor), tcrossprod(centered), tolerance = 1e-13)
  expect_equal(rows$se[c(1L,4L)], c(0,0))
  expect_equal(unname(.np_udist_marginal_range(1L, "epanechnikov2")),
    c(primitive(-sqrt(3/7)), primitive(sqrt(3/7))), tolerance = 1e-13)
  for (law in c("inid","fixed","geom"))
    expect_error(.np_udist_joint_inference(state, grid, law = law), "Epanechnikov")
  expect_error(.np_udist_joint_inference(state, grid, calibration = "score-envelope"), "Epanechnikov")
})

test_that("compact CDF admission follows the native large-bandwidth boundary", {
  old <- options(np.largeh = TRUE, np.largeh.rel.tol = .001)
  on.exit(options(old))
  state <- list(type = "fixed", x = data.frame(x = c(0, 1)), n = 2L,
    q = 1L, s = 0L, numeric = TRUE, h = 1/.038, category = list(),
    metadata = list(kernel = "epanechnikov2"))
  grid <- data.frame(x = .5)
  expect_silent(.np_udist_bc_grid(state, grid))
  state$h <- 1/.030
  expect_error(.np_udist_bc_grid(state, grid), "large-h/global-limit", fixed = TRUE)
})
