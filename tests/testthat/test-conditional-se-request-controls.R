test_that("conditional SE requests do not change point or gradient results", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(24L)
  x <- data.frame(x = .05 + .9*i/25)
  y <- data.frame(y = .04 + .92*((i*7L) %% 29L)/29)
  ex <- x[c(3L, 11L, 19L), , drop = FALSE]
  ey <- y[c(3L, 11L, 19L), , drop = FALSE]
  for (cdf in c(FALSE, TRUE)) {
    bw.fun <- if (cdf) npcdistbw else npcdensbw
    fit.fun <- if (cdf) npcdist else npcdens
    family <- if (cdf) "npcdist" else "npcdens"
    bw <- bw.fun(xdat = x, ydat = y, bws = c(.22, .31),
                 bandwidth.compute = FALSE, regtype = "lc")
    args <- list(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey)
    default <- do.call(fit.fun, args)
    expect_identical(default[["se", exact = TRUE]], FALSE)
    expect_null(default[["conderr", exact = TRUE]])
    expect_null(default[["congerr", exact = TRUE]])
    expect_error(se(default), "were not computed.*without repeating bandwidth search")
    expect_error(do.call(fit.fun, c(args, list(se = NA))), "se")
    expect_error(do.call(fit.fun, c(args, list(sse = TRUE))), "unused")
    for (gradient in c(FALSE, TRUE)) {
      off <- do.call(fit.fun, c(args, list(gradients = gradient, se = FALSE)))
      on <- do.call(fit.fun, c(args, list(gradients = gradient, se = TRUE)))
      expect_identical(fitted(off), fitted(on))
      expect_identical(off$congrad, on$congrad)
      expect_null(off$conderr)
      expect_null(off$congerr)
      expect_identical(on[["se", exact = TRUE]], TRUE)
      expect_identical(se(on), on$conderr)
      if (gradient) {
        expect_identical(gradients(off), gradients(on))
        expect_identical(gradients(on, se = TRUE), on$congerr)
        expect_error(gradients(off, se = TRUE), "gradients = TRUE, se = TRUE")
      }
      legacy <- on
      legacy[["se"]] <- NULL
      expect_identical(se(legacy), se(on))
    }
    pred <- predict(default, exdat = ex, eydat = ey,
                    txdat = x, tydat = y, se.fit = TRUE)
    on <- do.call(fit.fun, c(args, list(se = TRUE)))
    expect_identical(pred$fit, fitted(on))
    expect_identical(pred$se.fit, se(on))
    expect_error(predict(default, se.fit = TRUE, se = FALSE), "conflicting")
    expect_error(predict(default, se.fit = FALSE, se = TRUE), "conflicting")

    method <- getS3method(family, if (cdf) "condbandwidth" else "conbandwidth")
    expect_identical(tail(names(formals(method)), 2L), c("...", "se"))
    # The old extra positional argument still belongs to dots, not the new SE.
    positional <- do.call(method, list(bw, x, y, ex, ey, FALSE, 1L,
                                      FALSE, if (cdf) "isotonic" else "project",
                                      list(), TRUE))
    expect_identical(positional[["se", exact = TRUE]], FALSE)
    expect_identical(fitted(positional), fitted(default))
  }
})

test_that("conditional formula routes keep SE requests out of bandwidth controls", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.05, .95, length.out = 16L),
                  y = .5 + .3*sin(seq_len(16L)))
  for (cdf in c(FALSE, TRUE)) {
    f <- if (cdf) npcdist else npcdens
    off <- f(y ~ x, data = d, bandwidth.compute = FALSE, se = FALSE)
    on <- f(y ~ x, data = d, bandwidth.compute = FALSE, se = TRUE)
    expect_identical(fitted(off), fitted(on))
    expect_null(off$conderr)
    expect_length(on$conderr, nrow(d))
  }
})

test_that("conditional private inference masks are subordinate to public SE", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(25L)
  x <- data.frame(x = .04 + .92*i/26,
                  u = factor(rep(c("a", "b", "c"), length.out = 25L)))
  y <- data.frame(y = .03 + .94*((i*7L) %% 29L)/29)
  for (cdf in c(FALSE, TRUE)) {
    bw.fun <- if (cdf) npcdistbw else npcdensbw
    fit.fun <- if (cdf) npcdist else npcdens
    b <- bw.fun(xdat = x, ydat = y, bws = c(.2, .3, .2),
                bandwidth.compute = FALSE, regtype = "lc")
    a <- list(bws = b, txdat = x, tydat = y, exdat = x[3:5,,drop=FALSE],
              eydat = y[3:5,,drop=FALSE], gradients = TRUE)
    off <- do.call(fit.fun, c(a, list(se = FALSE)))
    on <- do.call(fit.fun, c(a, list(se = TRUE)))
    expect_identical(fitted(off), fitted(on))
    expect_identical(gradients(off), gradients(on))
    expect_null(off$congerr)
    expect_true(is.matrix(on$congerr))
  }
})

test_that("beta conditional point owners accept omitted uncertainty outputs", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(16L)
  x <- data.frame(x = .05 + .9*i/17)
  y <- data.frame(y = .04 + .92*((i*7L) %% 17L)/17)
  for (cdf in c(FALSE, TRUE)) for (side in c("X", "Y", "XY")) {
    xb <- side %in% c("X", "XY")
    yb <- side %in% c("Y", "XY")
    a <- list(xdat = x, ydat = y, bws = c(.2, .3),
      bandwidth.compute = FALSE, regtype = "lc",
      cxkertype = if (xb) "beta" else "gaussian",
      cykertype = if (yb) "beta" else "gaussian")
    if (xb) a <- c(a, list(cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1))
    if (yb) a <- c(a, list(cykerbound = "fixed", cykerlb = 0, cykerub = 1))
    b <- do.call(if (cdf) npcdistbw else npcdensbw, a)
    f <- if (cdf) npcdist else npcdens
    a <- list(bws = b, txdat = x, tydat = y, gradients = TRUE,
      exdat = x[c(4, 9), , drop = FALSE], eydat = y[c(4, 9), , drop = FALSE])
    on <- do.call(f, c(a, list(se = TRUE)))
    off <- do.call(f, c(a, list(se = FALSE)))
    expect_identical(fitted(on), fitted(off))
    expect_identical(gradients(on), gradients(off))
    expect_null(off$conderr)
    expect_null(off$congerr)
  }
})
