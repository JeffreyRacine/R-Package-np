test_that("quantile and mode SE controls are named-only without consuming dots", {
  for (name in c("npqreg.formula", "npqreg.default", "npqreg.condbandwidth",
                 "npconmode.formula", "npconmode.default", "npconmode.conbandwidth")) {
    fn <- getFromNamespace(name, "np")
    fml <- formals(fn)
    expect_identical(fml$se, FALSE)
    expect_gt(match("se", names(fml)), match("...", names(fml)))
    # Populate every old positional formal, then one forwarded dots argument.
    npos <- match("...", names(fml)) - 1L
    args <- rep(list(quote(positional)), npos + 1L)
    cl <- as.call(c(list(as.name(name)), args, list(se = TRUE)))
    matched <- match.call(fn, cl, expand.dots = FALSE)
    expect_identical(matched$se, TRUE)
    expect_length(matched[["..."]], 1L)
    expect_identical(matched[["..."]][[1L]], quote(positional))
  }
  strip <- getFromNamespace(".npqreg_strip_fit_controls_from_bw_call", "np")
  expect_identical(strip(quote(npcdistbw(y ~ x, se = TRUE, tau = .25))),
                   quote(npcdistbw(y ~ x)))
})

test_that("Q extractors read historical errors and preserve tau in missing-work hints", {
  old <- structure(list(quanterr = matrix(c(0, NA_real_), 2, 1),
                        quantgerr = array(0.2, c(2, 1, 2)),
                        tau = c(.2, .7)), class = "qregression")
  expect_identical(se(old), old$quanterr)
  expect_identical(gradients(old, se = TRUE), old$quantgerr)
  off <- old
  off$se <- FALSE
  expect_error(se(off), "without repeating bandwidth search", fixed = TRUE)
  expect_error(se(off), "tau = c(0.2, 0.7)", fixed = TRUE)
  expect_error(gradients(off, se = TRUE), "gradients = TRUE, se = TRUE", fixed = TRUE)
  absent <- structure(list(tau = c(.2, .7)), class = "qregression")
  expect_error(se(absent), "were not computed", fixed = TRUE)
  expect_error(predict(old, se.fit = TRUE, se = FALSE), "conflicting", fixed = TRUE)

  mode <- structure(list(conmode = factor(c("a", "b")),
                         xndim = 1L, xeval = data.frame(x = c(.2, .8)),
                         probabilities = matrix(c(.8, .2, .2, .8), 2),
                         probability.errors = matrix(.1, 2, 2)),
                    class = "conmode")
  expect_identical(predict(mode, type = "prob", se.fit = TRUE)$se.fit,
                   mode$probability.errors)
  mode$se <- FALSE
  expect_error(predict(mode, type = "prob", se.fit = TRUE),
               "probabilities = TRUE, se = TRUE", fixed = TRUE)
  expect_identical(predict(mode), mode$conmode)
  expect_error(predict(mode, type = "prob", se.fit = TRUE, se = FALSE),
               "conflicting", fixed = TRUE)
  plot.data <- getFromNamespace(".np_plot_conmode_data", "np")
  expect_error(plot.data(mode, errors = "asymptotic"),
               "probabilities = TRUE, se = TRUE", fixed = TRUE)
  expect_equal(plot.data(mode)$x$probability, mode$probabilities[, 1L])
})

test_that("mode selection does not depend on an error matrix", {
  select <- getFromNamespace(".npConmodeSelect", "np")
  probabilities <- matrix(c(NA, NA, 0, 0, .4, .4, .2, .6), 4, 2, byrow = TRUE)
  on <- select(probabilities, matrix(.1, 4, 2))
  off <- select(probabilities, NULL)
  expect_identical(off$indices, on$indices)
  expect_identical(off$condens, on$condens)
  expect_null(off$conderr)
  expect_identical(off$indices, c(0L, 0L, 1L, 2L))
})

test_that("quantile helper preserves the LSQ forwarding owner", {
  helper <- getFromNamespace(".np_plot_quantile_eval", "np")
  scope <- new.env(parent = environment(helper))
  scope$.np_plot_lsqregression_eval <- function(...) list(...)
  environment(helper) <- scope
  bw <- structure(list(tau = c(.25, .75)), class = "lsqregressionbandwidth")
  args <- list(bws = bw, txdat = data.frame(x = 1:3), tydat = 1:3,
               exdat = data.frame(x = 1:2), tau = bw$tau, gradients = TRUE,
               gradient.order = 2L)
  for (request in c(FALSE, TRUE)) {
    out <- do.call(helper, c(args, list(need.errors = request)))
    expect_identical(out$need.errors, request)
    for (name in names(args)) expect_identical(out[[name]], args[[name]], info = name)
    expect_false("se" %in% names(out))
  }
})

test_that("quantile SE demand leaves point and gradient estimates unchanged", {
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  y <- data.frame(y = sin(x$x) + .1 * cos(seq_len(nrow(x))))
  ex <- data.frame(x = c(-.3, NA_real_, .4))
  bw <- npcdistbw(xdat = x, ydat = y, bws = c(.6, .6),
                  bandwidth.compute = FALSE)
  args <- list(bws = bw, txdat = x, tydat = y, exdat = ex, tau = c(.25, .65))
  off <- do.call(npqreg, args)
  on <- do.call(npqreg, c(args, list(se = TRUE)))
  grad <- do.call(npqreg, c(args, list(gradients = TRUE)))
  both <- do.call(npqreg, c(args, list(gradients = TRUE, se = TRUE)))
  expect_false(off$se)
  expect_true(on$se)
  expect_identical(off$quantile, on$quantile)
  expect_identical(grad$quantile, both$quantile)
  expect_identical(grad$quantgrad, both$quantgrad)
  expect_null(off$quanterr)
  expect_null(off$quantgerr)
  expect_null(grad$quanterr)
  expect_null(grad$quantgerr)
  expect_equal(dim(on$quanterr), c(3L, 2L))
  expect_equal(dim(both$quantgerr), c(3L, 1L, 2L))
  expect_true(all(is.na(on$quanterr[2L, ])))
  expect_true(all(is.na(grad$quantgrad[2L, , ])))
  expect_identical(se(on), on$quanterr)
  expect_error(se(off), "without repeating bandwidth search", fixed = TRUE)
  pred <- predict(off, exdat = ex, txdat = x, tydat = y, se.fit = TRUE)
  expect_identical(pred$fit, on$quantile)
  expect_identical(pred$se.fit, on$quanterr)
  expect_error(npqreg(bws = bw, txdat = x, tydat = y, se = NA),
               "'se' must be TRUE or FALSE", fixed = TRUE)
})

test_that("mode SE demand leaves class, probability and gradient estimates unchanged", {
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  y <- data.frame(y = factor(rep(c("a", "b", "b", "a"), 6L)))
  bw <- npcdensbw(xdat = x, ydat = y, bws = c(.15, .6),
                  bandwidth.compute = FALSE, regtype = "ll")
  args <- list(bws = bw, txdat = x, tydat = y, exdat = x[c(5, 12, 20), , drop = FALSE],
               probabilities = TRUE, gradients = TRUE, level = "a")
  off <- do.call(npconmode, args)
  on <- do.call(npconmode, c(args, list(se = TRUE)))
  expect_false(off$se)
  expect_true(on$se)
  for (field in c("conmode", "condens", "probabilities", "probability.gradients",
                  "proper.info", "probability.gradient.level"))
    expect_identical(off[[field]], on[[field]], info = field)
  expect_null(off$conderr)
  expect_null(off$probability.errors)
  expect_identical(dim(on$probability.errors), dim(on$probabilities))
  expect_identical(predict(off, type = "prob"), on$probabilities)
  expect_error(predict(off, type = "prob", se.fit = TRUE),
               "without repeating bandwidth search", fixed = TRUE)
  pred <- predict(off, exdat = args$exdat, txdat = x, tydat = y,
                   type = "prob", se.fit = TRUE)
  expect_identical(pred$fit, on$probabilities)
  expect_identical(pred$se.fit, on$probability.errors)
  expect_error(npconmode(bws = bw, txdat = x, tydat = y, se = c(TRUE, FALSE)),
               "'se' must be TRUE or FALSE", fixed = TRUE)
})
