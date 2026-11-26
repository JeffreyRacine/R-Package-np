test_that("the local-constant design envelope has a separate explicit control", {
  expect_identical(np_conditional_bias_control("pairs-score")$version, 1L)
  expect_identical(np_conditional_bias_control("pairs-score", simultaneous = "rms-envelope")$version, 3L)
  ctrl <- np_conditional_bias_control("pairs-score", simultaneous = "design-envelope")
  expect_identical(ctrl$version, 4L)
  expect_error(np_conditional_bias_control("analytic", simultaneous = "design-envelope"), "requires pairs-score")
  expect_error(np_conditional_bias_control("pairs-score", "local-wilson", "design-envelope"),
               "sharp-score")
})

test_that("design-envelope Bonferroni diagnostics identify the actual finite ranks", {
  check <- function(B, m) .np_conditional_design_resolution(B, m, .05)
  dense <- check(9999L, 2500L)
  expect_true(dense$limited)
  expect_equal(dense$expected.tail.draws, .09999)
  expect_identical(c(dense$lower.index, dense$upper.index), c(1L, 9999L))
  for (params in list(c(9999L,250L), c(1999L,50L), c(399L,9L))) {
    v <- check(params[1L], params[2L])
    expect_false(v$limited)
    ranks <- quantile(seq_len(v$B), c(v$nominal.tail, 1-v$nominal.tail), type=1L, names=FALSE)
    expect_equal(c(v$lower.index,v$upper.index), ranks)
    expect_equal(c(v$lower.rank.tail,v$upper.rank.tail), c(ranks[1L],v$B+1-ranks[2L])/(v$B+1))
  }
})

test_that("the local-constant design envelope preserves its sharp-score marginal and count roots", {
  x <- seq(-1, 1, length.out = 67L)
  d <- data.frame(x = x, y = sin(3*x) + .2*cos(seq_along(x)))
  bws <- npcdistbw(y ~ x, data = d, bws = c(.31, .28),
                   regtype = "lc", bandwidth.compute = FALSE)
  model <- npcdist(bws)
  args <- list(x = model, center = "bias-corrected", errors = "bootstrap",
               bootstrap = "inid", interval = "percentile", band = "all",
               neval = 5L, B = 131L, xtrim = .1, output = "data", random.seed = 17L)
  retrieve <- function(ctrl) suppressWarnings(do.call(plot, c(args, list(bias.control = ctrl))))
  base <- retrieve(np_conditional_bias_control("pairs-score"))
  alternative <- retrieve(np_conditional_bias_control("pairs-score", simultaneous = "design-envelope"))
  matrix.for <- function(obj, scope) as.data.frame(obj, panel = 1L, target = 1L, band = scope)
  point <- matrix.for(base, "pointwise")
  expect_identical(point, matrix.for(alternative, "pointwise"))
  expect_match(alternative$inference[[1L]]$endpoint_construction, "predictor-scaled")
  joint <- matrix.for(alternative, "simultaneous")
  bonf <- matrix.for(alternative, "bonferroni")
  expect_true(all(joint$lower <= point$lower & joint$upper >= point$upper))
  expect_true(all(bonf$lower <= joint$lower & bonf$upper >= joint$upper))
  expect_false(identical(joint$lower, matrix.for(base, "simultaneous")$lower))
})

test_that("the local-constant design option refuses polynomial conditional fits", {
  x <- seq(-1, 1, length.out = 45L)
  d <- data.frame(x = x, y = x^2 + cos(seq_along(x)))
  ctrl <- np_conditional_bias_control("pairs-score", simultaneous = "design-envelope")
  for (degree in 1:3) {
    b <- npcdistbw(y ~ x, data = d, bws = c(.3, .4), regtype = "lp",
                   degree = degree, bernstein.basis = TRUE, bandwidth.compute = FALSE)
    expect_error(.np_conditional_control(ctrl, b, TRUE, "bootstrap", "bias-corrected",
                                        "default", "inid", FALSE, FALSE, FALSE),
                 "local-constant fits")
  }
})
