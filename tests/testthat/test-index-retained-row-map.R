test_that("single-index one-call bandwidth replay retains its training row map", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1909204)
  d <- data.frame(x = runif(36), z = runif(36), y = rnorm(36))
  d$x[c(3, 7)] <- NA_real_
  for (case in list(list(bt = "fixed", rt = "lc", h = .5),
                    list(bt = "generalized_nn", rt = "ll", h = 12),
                    list(bt = "adaptive_nn", rt = "lp", h = 12))) {
    for (policy in list(na.omit, na.exclude)) {
      args <- list(formula = y ~ x + z, data = d, bws = c(1, .3, case$h),
                   bandwidth.compute = FALSE, bwtype = case$bt,
                   regtype = case$rt, na.action = policy)
      if (case$rt == "lp") args$degree <- 2
      original <- do.call(npindex, c(args, list(se = TRUE, gradients = TRUE, residuals = TRUE)))
      replay <- npindex(unserialize(serialize(original$bws, NULL)),
                        se = TRUE, gradients = TRUE, residuals = TRUE)
      expect_identical(fitted(replay), fitted(original))
      expect_identical(se(replay), se(original))
      expect_identical(gradients(replay), gradients(original))
      expect_identical(gradients(replay, se = TRUE), gradients(original, se = TRUE))
      expect_identical(residuals(replay), residuals(original))
      expect_identical(replay$omit, original$omit)
      next.fit <- npindex(replay$bws, se = FALSE)
      expect_identical(fitted(next.fit), fitted(original))
      lazy <- do.call(npindex, c(args, list(se = FALSE)))
      expect_identical(residuals(lazy), residuals(original))
      nd <- d[c(1, 4, 9), ]
      ev <- npindex(original$bws, newdata = nd, residuals = TRUE, se = FALSE)
      native <- npindex(original$bws, txdat = na.omit(d)[c("x", "z")],
                        tydat = na.omit(d)$y, exdat = nd[c("x", "z")],
                        residuals = TRUE, se = FALSE)
      expect_identical(fitted(ev), fitted(native))
      expect_identical(residuals(ev), residuals(original))
      expect_length(fitted(ev), 3L)
      expect_length(residuals(native), 34L)
    }
  }
})

test_that("single-index retained replay does not reevaluate formula transformations", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  state <- new.env(parent = emptyenv()); state$n <- 0L
  transform.x <- function(x) { state$n <- state$n + 1L; x + .2 }
  d <- data.frame(x = seq(.1, .9, length.out = 30), y = sin(seq_len(30)))
  d$x[4] <- NA_real_
  expect_warning(fit <- npindex(y ~ transform.x(x), data = d, subset = seq_len(28),
                 bws = c(1, .4), bandwidth.compute = FALSE,
                 na.action = na.exclude, se = FALSE), "xdat has one dimension")
  expect_identical(state$n, 1L)
  replay <- npindex(fit$bws, se = FALSE, residuals = TRUE)
  expect_identical(state$n, 1L)
  expect_identical(fitted(replay), fitted(fit))
  expect_length(residuals(replay), 28L)
  # The retained subset also applies to a full explicit replacement; use the
  # same original sample size to keep this a row-map rather than subset test.
  replacement <- d; replacement$x[4] <- .3; replacement$x[6] <- NA_real_
  a <- npindex(fit$bws, data = replacement, se = FALSE, residuals = TRUE)
  expect_identical(which(is.na(fitted(a))), 6L)
  expect_identical(which(is.na(residuals(a))), 6L)
  b <- npindex(a$bws, se = FALSE, residuals = TRUE)
  expect_identical(fitted(b), fitted(a))
  expect_identical(residuals(b), residuals(a))
  # A historical object without a snapshot still uses its native-value reader;
  # no original omission map can be reconstructed from compact values.
  legacy <- fit$bws; legacy[[".np.formula.training"]] <- NULL
  old.fit <- npindex(legacy, se = FALSE)
  expect_identical(unname(fitted(old.fit)), unname(fitted(fit)[-4]))
})
