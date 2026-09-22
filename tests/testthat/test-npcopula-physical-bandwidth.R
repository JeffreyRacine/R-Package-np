test_that("copula marginal owners receive physical bandwidths", {
  helper <- getFromNamespace(".npcopula_marginal_bw", "np")
  args.helper <- getFromNamespace(".npcopula_marginal_bw_args", "np")
  set.seed(170)
  d <- data.frame(x = rnorm(32, sd = 2), z = runif(32),
                  o = ordered(rep(1:4, 8)))
  for (target in c("density", "distribution"))
    for (type in c("fixed", "generalized_nn", "adaptive_nn"))
      for (scaled in c(FALSE, TRUE)) {
        constructor <- if(target == "density") npudensbw else npudistbw
        b <- constructor(dat = d,
          bws = c(if(type == "fixed") 1.1 else 18,
                  if(type == "fixed") .8 else 20, .2),
          bwtype = type, bwscaling = scaled, bandwidth.compute = FALSE)
        saved <- serialize(b, NULL)
        for (j in seq_along(d)) {
          expected <- unname(b$bandwidth$x[j])
          args <- args.helper(b, d, j, target)
          expect_identical(unname(args$bws), expected)
          kb <- helper(b, d, j, target, kbandwidth = TRUE)
          mb <- helper(b, d, j, target, kbandwidth = FALSE)
          expect_identical(unname(kb$bw), expected)
          expect_identical(unname(mb$bandwidth$x), expected)
        }
        expect_identical(serialize(b, NULL), saved)
      }
})

test_that("scaled copula samples, grids, uncertainty and bootstrap match physical twins", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(1701)
  d <- data.frame(x = rnorm(32, sd = 2), z = rnorm(32, sd = .4))
  u <- data.frame(x = c(.25, .6), z = c(.3, .7))
  for (target in c("density", "distribution")) {
    constructor <- if(target == "density") npudensbw else npudistbw
    scaled <- constructor(dat = d, bws = c(.8, 1.2), bwscaling = TRUE,
                            bandwidth.compute = FALSE)
    physical <- constructor(dat = d, bws = scaled$bandwidth$x,
                              bandwidth.compute = FALSE)
    for (grid in c(FALSE, TRUE)) {
      args <- if(grid) list(u = u, n.quasi.inv = 40L) else list()
      a <- do.call(npcopula, c(list(bws = scaled, data = d, se = TRUE), args))
      b <- do.call(npcopula, c(list(bws = physical, data = d, se = TRUE), args))
      expect_identical(as.data.frame(a), as.data.frame(b))
      expect_identical(fitted(a), fitted(b))
      expect_identical(se(a), se(b))
      expect_identical(predict(a, u = u, n.quasi.inv = 40L),
                       predict(b, u = u, n.quasi.inv = 40L))
      if(grid) {
        if(target == "distribution")
          expect_identical(plot(a, view = "all", output = "data"),
                           plot(b, view = "all", output = "data"))
        set.seed(1702)
        expect_warning(pa <- plot(a, errors = "bootstrap", B = 9,
          band = "all", output = "data"), "B=9 is too small", fixed = TRUE)
        seed <- .Random.seed
        set.seed(1702)
        expect_warning(pb <- plot(b, errors = "bootstrap", B = 9,
          band = "all", output = "data"), "B=9 is too small", fixed = TRUE)
        expect_identical(pa, pb)
        expect_identical(.Random.seed, seed)
      }
    }
  }
})
