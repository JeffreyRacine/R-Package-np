test_that("external compiled hats distinguish strict failure from typed empty rows", {
  invoke <- function(name, ...) .Call(name, ..., PACKAGE = "np")
  w <- cbind(1, seq(-1, 1, length.out = 12))
  e <- rbind(c(1, -.3), c(1, .2), c(1, .6))
  k <- cbind(rep(.6, 12), rep(0, 12), seq(.2, .9, length.out = 12))
  strict <- function(k, wt = w, ev = e)
    invoke("C_np_reghat_lp_matrix_fast", k, wt, ev)
  external <- function(k, wt = w, ev = e, norm = FALSE)
    invoke("C_np_reghat_lp_matrix_external", k, wt, ev, norm)
  expect_error(strict(k, w, e), "invalid wider-LP", fixed = TRUE)
  out <- external(k, w, e)
  expect_identical(attr(out, ".np.empty.rows"), c(0L, 1L, 0L))
  expect_true(all(is.na(out[2L, ])))
  expect_identical(unname(out[c(1L, 3L), ]),
                   strict(k[, c(1L, 3L)], w, e[c(1L, 3L), ]))
  norms <- external(k, w, e, TRUE)
  expect_identical(norms$hat, out)
  expect_true(all(is.na(norms$norm[2L, ])))
  healthy <- k
  healthy[, 2L] <- .7
  expect_identical(external(healthy, w, e), strict(healthy, w, e))
  expect_identical(external(healthy, w, e, TRUE),
                   invoke("C_np_reghat_lp_matrix_norm", healthy, w, e))
  expect_null(attr(external(healthy, w, e), ".np.empty.rows"))
  signed <- healthy
  signed[, 2L] <- rep(c(-1, 1), 6L)
  outcome <- function(expr) tryCatch(expr, error = conditionMessage)
  expect_identical(outcome(external(signed, w, e)), outcome(strict(signed, w, e)))
  k[1L, 2L] <- NaN
  expect_error(external(k, w, e), "non-finite", fixed = TRUE)
  expect_error(invoke("C_np_reghat_lp_matrix_external", healthy, w, e, NA),
               "external hat norm request", fixed = TRUE)
})

test_that("conditional external base rows preserve supported fits and required purpose", {
  opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(opts), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 80))
  y <- data.frame(y = sin(3*x$x) + .1*cos(seq_len(nrow(x))))
  ex <- data.frame(x = c(-.4, .4, 4))
  ey <- data.frame(y = c(-.3, .3, 0))
  capture <- function(expr) {
    notices <- character()
    out <- withCallingHandlers(expr, warning = function(w) {
      notices <<- c(notices, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(out = out, notices = notices)
  }
  for (cdf in c(FALSE, TRUE)) for (reg in c("lc", "lp")) {
    bwfun <- if (cdf) npcdistbw else npcdensbw
    fitfun <- if (cdf) npcdist else npcdens
    args <- list(xdat = x, ydat = y, bws = c(.5, .6),
                 bandwidth.compute = FALSE, regtype = reg,
                 cxkertype = "epanechnikov", cykertype = "epanechnikov")
    if (reg == "lp") args$degree <- 3L
    bw <- do.call(bwfun, args)
    for (inference in c(FALSE, TRUE)) {
      z <- capture(fitfun(bws = bw, txdat = x, tydat = y,
        exdat = ex, eydat = ey, gradients = inference, se = inference))
      good <- fitfun(bws = bw, txdat = x, tydat = y,
        exdat = ex[1:2, , drop = FALSE], eydat = ey[1:2, , drop = FALSE],
        gradients = inference, se = inference)
      expect_identical(attr(z$out, ".np.empty.base.rows"), c(0L, 0L, 1L))
      expect_identical(fitted(z$out)[1:2], fitted(good))
      expect_true(is.na(fitted(z$out)[3L]))
      expect_length(z$notices, 1L)
      if (inference) expect_true(all(is.na(gradients(z$out)[3L, , drop = FALSE])))
    }
    if (reg == "lp")
      expect_error(fitfun(bws = bw, txdat = x, tydat = y, exdat = ex,
                         eydat = ey, .np.require.complete = TRUE),
                   "LP solve failed", fixed = TRUE)
  }
})
