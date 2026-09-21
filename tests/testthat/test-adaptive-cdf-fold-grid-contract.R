.ann_cdf_grid_objective <- function(dat, state, grid = NULL) {
  getFromNamespace("npudistbw.dbandwidth", "np")(
    dat = dat, bws = state, bandwidth.compute = TRUE, eval.only = TRUE,
    do.full.integral = is.null(grid), gdat = grid, nmulti = 1L,
    invalid.penalty = "dbmax")$fval
}

.ann_cdf_grid_refit <- function(dat, state, grid = NULL) {
  n <- nrow(dat)
  empirical <- is.null(grid)
  if (empirical) grid <- dat
  loss <- vapply(seq_len(n), function(held) {
    prediction <- as.numeric(fitted(npudist(
      bws = state, tdat = dat[-held, , drop = FALSE], edat = grid)))
    indicator <- rep(TRUE, nrow(grid))
    for (j in seq_along(dat))
      indicator <- indicator & (dat[[j]][held] <= grid[[j]])
    error <- (indicator - prediction)^2
    if (empirical) error <- error[-held]
    mean(error)
  }, numeric(1L))
  mean(loss)
}

test_that("ANN CDF CV deletes donors on empirical and external grids", {
  old <- options(np.messages = FALSE, np.largeh = FALSE,
                 np.largelambda = FALSE, np.extendednn = TRUE)
  on.exit(options(old))
  set.seed(471)
  n <- 17L
  for (p in c(1L, 2L, 3L, 5L)) for (kernel in c("gaussian", "beta")) {
    dat <- as.data.frame(matrix(runif(n * p, .03, .97), ncol = p))
    dat[2, ] <- dat[1, ] # delete occurrences, not all equal values
    grid <- as.data.frame(matrix(runif(7L * p), ncol = p))
    names(grid) <- names(dat)
    for (mixed in c(FALSE, TRUE)) {
      if (mixed) {
        dat$o <- ordered(rep(c(1, 3, 7), length.out = n), levels = c(1, 3, 7))
        grid$o <- ordered(rep(c(1, 3, 7), length.out = 7L), levels = c(1, 3, 7))
      }
      for (order in c(2L, 4L)) for (k in c(4L, n - 1L)) {
        args <- list(dat = dat, bws = c(rep(k, p), if (mixed) .3),
          bwtype = "adaptive_nn", ckertype = kernel, ckerorder = order,
          okertype = "racineliyan", bandwidth.compute = FALSE)
        if (kernel == "beta") args <- c(args, list(ckerbound = "fixed",
          ckerlb = rep(0, p), ckerub = rep(1, p)))
        state <- do.call(npudistbw, args)
        for (external in c(FALSE, TRUE)) {
          query <- if (external) grid else NULL
          expected <- .ann_cdf_grid_refit(dat, state, query)
          expect_equal(as.numeric(.ann_cdf_grid_objective(dat, state, query)),
            expected, tolerance = 3e-12,
            info = paste(p, kernel, mixed, order, k, external))
        }
      }
    }
  }
})

test_that("ANN CDF CV agrees with independent beta and Gaussian ranks", {
  old <- options(np.messages = FALSE, np.largeh = FALSE,
                 np.largelambda = FALSE, np.extendednn = FALSE)
  on.exit(options(old))
  x <- c(.03, .08, .15, .29, .41, .53, .69, .81, .92)
  grid <- c(0, .12, .46, .78, 1)
  n <- length(x); k <- 3L
  for (kernel in c("gaussian", "beta")) {
    args <- list(dat = data.frame(x), bws = k, bwtype = "adaptive_nn",
      bandwidth.compute = FALSE, ckertype = kernel)
    if (kernel == "beta") args <- c(args, list(ckerbound = "fixed",
      ckerlb = 0, ckerub = 1))
    state <- do.call(npudistbw, args)
    for (empirical in c(TRUE, FALSE)) {
      query <- if (empirical) x else grid
      errors <- vapply(seq_len(n), function(i) {
        donors <- setdiff(seq_len(n), i)
        h <- vapply(donors, function(donor)
          sort(abs(x[setdiff(donors, donor)] - x[donor]))[k], 0)
        fit <- vapply(query, function(q) {
          if (kernel == "gaussian") mean(pnorm((q-x[donors])/h)) else
            mean(pbeta(q, 1+x[donors]/h^2, 1+(1-x[donors])/h^2))
        }, numeric(1L))
        loss <- ((x[i] <= query) - fit)^2
        if (empirical) loss <- loss[-i]
        mean(loss)
      }, numeric(1L))
      expect_equal(as.numeric(.ann_cdf_grid_objective(data.frame(x), state,
        if (empirical) NULL else data.frame(x = grid))), mean(errors),
        # The legacy Gaussian CDF approximation has an established 2e-10
        # oracle contract; beta uses the tighter incomplete-beta comparison.
        tolerance = if (kernel == "gaussian") 2e-10 else 2e-13)
    }
  }
})

test_that("ANN interval ranks retain scaled inclusive ties", {
  old <- options(np.messages = FALSE, np.largeh = FALSE,
                 np.largelambda = FALSE, np.extendednn = TRUE)
  on.exit(options(old))
  dat <- data.frame(x = c(0, 0, .125, .25, .5, .75, .875, 1),
                    z = c(.5, .25, .75, .125, 1, 0, .875, .5))
  for (scale in c(1, 7)) {
    values <- 13 + scale * dat
    b <- npudistbw(dat = values, bws = c(7, 7),
      bwtype = "adaptive_nn", bandwidth.compute = FALSE)
    expect_equal(as.numeric(.ann_cdf_grid_objective(values, b)),
      .ann_cdf_grid_refit(values, b), tolerance = 3e-12)
  }
})
