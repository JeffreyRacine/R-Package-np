test_that("evaluation exclusion restores certificates with gradient rows only", {
  package <- environmentName(environment(npreg))
  restore <- getFromNamespace(".npreg_restore_eval_exclude", package)
  cert <- matrix(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), 3L, 2L,
    dimnames = list(c("a", "c", "d"), c("u", "x")))
  omitted <- c(b = 2L)
  fixture <- list(nobs = 3L, eval.rows.omit = omitted, mean = 1:3,
    merr = rep(1, 3L), grad = matrix(1:6, 3L, 2L, dimnames = dimnames(cert)),
    gerr = matrix(1, 3L, 2L, dimnames = dimnames(cert)))
  attr(fixture, ".np.gradient.structural.zero") <- cert
  full <- restore(fixture, fields = c("mean", "merr", "grad", "gerr"))
  got <- attr(full, ".np.gradient.structural.zero", exact = TRUE)
  expected <- stats::napredict(structure(omitted, class = "exclude"), cert)
  expect_identical(got, expected)
  expect_identical(dim(got), dim(full$grad))
  expect_type(got, "logical")
  expect_true(all(is.na(got[2L, ])))
  expect_identical(restore(full, fields = c("mean", "merr", "grad", "gerr")), full)

  # predict's mean-only restore must keep its still-compact gradient and
  # certificate aligned rather than expanding the certificate alone.
  mean.only <- restore(fixture)
  expect_identical(mean.only$grad, fixture$grad)
  expect_identical(attr(mean.only, ".np.gradient.structural.zero"), cert)
  errors.only <- restore(fixture, fields = "gerr")
  expect_identical(errors.only$grad, fixture$grad)
  expect_identical(attr(errors.only, ".np.gradient.structural.zero"), cert)
  expect_identical(restore(mean.only, fields = c("grad", "gerr")), full)
  no.cert <- fixture
  attr(no.cert, ".np.gradient.structural.zero") <- NULL
  expect_null(attr(restore(no.cert, fields = c("grad", "gerr")),
                   ".np.gradient.structural.zero", exact = TRUE))
  none <- fixture
  none$eval.rows.omit <- integer()
  expect_identical(restore(none, fields = c("grad", "gerr")), none)

  bad <- fixture
  attr(bad, ".np.gradient.structural.zero") <- matrix(FALSE, 2L, 2L)
  expect_error(restore(bad, fields = c("grad", "gerr")),
               "inconsistent regression evaluation rows")
  attr(bad, ".np.gradient.structural.zero") <- matrix(0, 3L, 2L)
  expect_error(restore(bad, fields = c("grad", "gerr")),
               "invalid regression structural-zero certificate")
  attr(bad, ".np.gradient.structural.zero") <- matrix(FALSE, 3L, 1L)
  expect_error(restore(bad, fields = c("grad", "gerr")),
               "invalid regression structural-zero certificate")
})

test_that("external formula na.exclude retains unknown-donor certificate shape", {
  package <- environmentName(environment(npreg))
  run <- function(expression) {
    if (identical(package, "npRmpi"))
      getFromNamespace(".npRmpi_with_local_regression", package)(expression)
    else
      expression
  }
  run({
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  dat <- data.frame(y = c(1, 3, 2), u = factor(c("a", "b", "c")),
                    x = c(0, .5, 1))
  bw <- npregbw(y ~ u + x, data = dat, bws = c(.2, .15),
    bandwidth.compute = FALSE, regtype = "lc", ckertype = "beta",
    ckerorder = 2L, ckerbound = "fixed", ckerlb = 0, ckerub = 1)
  ev <- dat[c(3L, 2L, 1L, 2L), c("u", "x")]
  rownames(ev) <- c("last", "missing", "first", "middle")
  ev$x[2L] <- NA_real_
  compact <- suppressWarnings(npreg(bws = bw, data = dat,
    newdata = ev[-2L, , drop = FALSE], gradients = TRUE, se = TRUE))
  fit <- suppressWarnings(npreg(bws = bw, data = dat, newdata = ev,
    gradients = TRUE, se = TRUE, na.action = stats::na.exclude))
  original <- attr(compact, ".np.gradient.structural.zero", exact = TRUE)
  got <- attr(fit, ".np.gradient.structural.zero", exact = TRUE)
  expect_true(is.matrix(original))
  expect_type(original, "logical")
  expect_true(anyNA(compact$merr))
  expect_true(any(original))
  expect_identical(dim(got), dim(fit$grad))
  expect_identical(unname(got[-2L, , drop = FALSE]), unname(original))
  expect_true(all(is.na(got[2L, ])))
  expect_identical(unname(fit$grad[-2L, , drop = FALSE]), unname(compact$grad))
  expect_identical(unname(fit$gerr[-2L, , drop = FALSE]), unname(compact$gerr))
  expect_identical(as.vector(fit$mean[-2L]), as.vector(compact$mean))
  expect_identical(as.vector(fit$merr[-2L]), as.vector(compact$merr))
  expect_true(all(is.na(fit$grad[2L, ])))
  })
})
