test_that("single-index factors use retained R contrast coordinates", {
  old <- options(np.messages = FALSE, contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old), add = TRUE)
  set.seed(619)
  n <- 120L
  d <- data.frame(x = rnorm(n),
                  race = factor(rep(c("A", "B", "C"), length.out = n)),
                  score = ordered(rep(c("0", "2", "6"), each = 40)))
  d$yc <- sin(d$x + .4 * (d$race == "B") - .6 * (d$race == "C")) + rnorm(n, sd = .2)
  d$yb <- rbinom(n, 1, plogis(d$x + .4 * (d$race == "B") - .6 * (d$race == "C")))
  for (method in c("ichimura", "kleinspady")) {
    d$y <- if (method == "ichimura") d$yc else d$yb
    mm <- model.matrix(y ~ x + race + score, d)[, -1L, drop = FALSE]
    v <- c(1, .4, -.6, .2, -.1, .8)
    b <- npindexbw(y ~ x + race + score, data = d, method = method,
                   bws = v, bandwidth.compute = FALSE)
    expect_identical(names(coef(b)), colnames(mm))
    expect_identical(b$index.design$raw.names, c("x", "race", "score"))
    f <- npindex(bws = b, se = TRUE, gradients = TRUE)
    bm <- npindexbw(xdat = as.data.frame(mm), ydat = d$y, bws = v,
                    method = method, bandwidth.compute = FALSE)
    fm <- npindex(bws = bm, se = TRUE, gradients = TRUE)
    expect_equal(as.vector(f$index), as.vector(mm %*% b$beta), tolerance = 1e-12)
    expect_identical(fitted(f), fitted(fm))
    expect_equal(unname(gradients(f)), unname(gradients(fm)), tolerance = 1e-12)
    expect_identical(colnames(gradients(f)), colnames(mm))
    expect_equal(vcov(f), vcov(fm), tolerance = 1e-12)
    eval <- d[seq_len(9L), c("score", "race", "x", "y")]
    eval$race <- factor(eval$race, levels = c("C", "B", "A"))
    before <- predict(f, newdata = eval)
    options(contrasts = c("contr.sum", "contr.helmert"))
    expect_identical(predict(f, newdata = eval), before)
    options(contrasts = c("contr.treatment", "contr.poly"))
    expect_equal(before, fitted(f)[seq_len(9L)], tolerance = 1e-12)
    one.level <- droplevels(d[d$race == "B", ])
    expect_equal(predict(f, newdata = one.level), fitted(f)[d$race == "B"], tolerance = 1e-12)
    unknown <- eval
    unknown$race <- factor(rep("new", nrow(eval)))
    expect_error(predict(f, newdata = unknown), "new level")
    expect_error(npindexbw(y ~ x + race + score, data = d, bws = c(1, .4, .8),
                           bandwidth.compute = FALSE), "requires coefficients")
    expect_error(npindexbw(y ~ 0 + x + race, data = d, bws = rep(1, 5),
                           bandwidth.compute = FALSE), "unidentified constant")
  }
})

test_that("explicit contrasts, raw retention and factor-first normalization agree", {
  old <- options(np.messages = FALSE, contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old), add = TRUE)
  set.seed(719)
  d <- data.frame(x = rnorm(120), group = factor(rep(c("A", "B", "C"), 40)))
  d$y <- rnorm(120)
  contrasts(d$group) <- contr.sum(3)
  mm <- model.matrix(y ~ group + x, d)[, -1L, drop = FALSE]
  b <- npindexbw(y ~ group + x, data = d, bws = c(1, -.3, .7, .9),
                 bandwidth.compute = FALSE)
  f <- npindex(bws = b, se = FALSE, gradients = TRUE)
  expect_identical(names(coef(f)), colnames(mm))
  expect_equal(as.vector(f$index), as.vector(mm %*% b$beta), tolerance = 1e-12)
  native <- npindexbw(xdat = d[c("group", "x")], ydat = d$y,
                      bws = c(1, -.3, .7, .9), bandwidth.compute = FALSE)
  nf <- npindex(bws = native, se = FALSE)
  expect_equal(nf$index, f$index, tolerance = 1e-12)
  expect_true(is.factor(native[[".np.native.training"]]$xdat$group))
  expect_equal(predict(nf, newdata = d[1:7, ]), fitted(nf)[1:7], tolerance = 1e-12)
  expect_identical(native$index.design$contrasts$group, contrasts(d$group))
  d$x[2L] <- NA_real_
  bna <- npindexbw(y ~ group + x, data = d, na.action = na.exclude,
                   bws = c(1, -.3, .7, .9), bandwidth.compute = FALSE)
  fn <- npindex(bws = bna, se = FALSE)
  expect_length(fitted(fn), 120L)
  expect_true(is.na(fitted(fn)[2L]))
})


test_that("single-index construction freezes session contrast options", {
  old <- options(np.messages = FALSE, contrasts = c("contr.sum", "contr.helmert"))
  on.exit(options(old), add = TRUE)
  set.seed(902)
  d <- data.frame(x = rnorm(120), g = factor(rep(letters[1:3], 40)),
                  o = ordered(rep(1:3, each = 40)))
  d$y <- rbinom(120, 1, plogis(d$x))
  mm <- model.matrix(y ~ x + g + o, d)[, -1L, drop = FALSE]
  b <- npindexbw(y ~ x + g + o, data = d, method = "kleinspady",
                 bws = c(1, .3, -.2, .1, -.1, .9), bandwidth.compute = FALSE)
  f <- npindex(bws = b, se = FALSE)
  expect_identical(names(coef(b)), colnames(mm))
  expect_equal(as.vector(f$index), as.vector(mm %*% coef(b)), tolerance = 1e-12)
})
