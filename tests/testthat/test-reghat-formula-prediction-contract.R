test_that("formula hat prediction reuses trained predictor transformations", {
  old <- options(np.messages = FALSE, na.action = "na.omit")
  on.exit(options(old), add = TRUE)
  ctr <- function(x, center = mean(x)) {
    structure(x - center, center = center, class = c("numeric", "r21hatcenter"))
  }
  make.center <- function(var, call) {
    call$center <- attr(var, "center")
    call
  }
  previous <- getS3method("makepredictcall", "r21hatcenter", optional = TRUE)
  registerS3method("makepredictcall", "r21hatcenter", make.center,
                   envir = asNamespace("stats"))
  on.exit({
    if (is.null(previous)) {
      table <- get(".__S3MethodsTable__.", envir = asNamespace("stats"))
      rm("makepredictcall.r21hatcenter", envir = table)
    } else registerS3method("makepredictcall", "r21hatcenter", previous,
                            envir = asNamespace("stats"))
  }, add = TRUE)
  d <- data.frame(x = seq(.25, 3, length.out = 31))
  d$y <- sin(log(d$x)) + .2 * d$x
  d$y[5L] <- NA_real_
  e <- data.frame(x = c(.35, .7, 1.2, 2.7))
  keep <- !is.na(d$y)
  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    h <- if (bwtype == "fixed") .4 else 8
    for (shape in c("ordinary", "log", "center")) {
      fm <- switch(shape, ordinary = y ~ x, log = y ~ log(x), center = y ~ ctr(x))
      b <- npregbw(fm, data = d, bws = h, bandwidth.compute = FALSE,
                   bwtype = bwtype, regtype = "ll")
      H <- npreghat(b)
      train <- data.frame(z = switch(shape, ordinary = d$x, log = log(d$x),
                                    center = d$x - mean(d$x))[keep])
      grid <- data.frame(z = switch(shape, ordinary = e$x, log = log(e$x),
                                   center = e$x - mean(d$x)))
      names(train) <- names(grid) <- b$xnames
      for (s in 0:1) {
        native <- npreghat(b, txdat = train, exdat = grid, s = s)
        direct <- npreghat(b, newdata = e, s = s)
        replay <- predict(H, newdata = e, s = s)
        expect_equal(unname(as.vector(direct)), unname(as.vector(native)), tolerance = 1e-11)
        expect_equal(unname(as.vector(replay)), unname(as.vector(native)), tolerance = 1e-11)
        expect_equal(predict(H, newdata = e, s = s, y = d$y[keep], output = "apply"),
                     as.vector(native %*% d$y[keep]), tolerance = 1e-11)
        expect_equal(as.vector(predict(H, newdata = e[1:2, , drop = FALSE], s = s)),
                     as.vector(native[1:2, , drop = FALSE]), tolerance = 1e-11)
        expect_equal(as.vector(predict(H, exdat = grid, s = s)), as.vector(native),
                     tolerance = 1e-11)
        expect_equal(as.vector(predict(H, newdata = stop("unused raw grid"),
                                      exdat = grid, s = s)), as.vector(native), tolerance = 1e-11)
        expect_equal(as.vector(npreghat(b, newdata = data.frame(not_x = 1),
                                       exdat = grid, s = s)), as.vector(native), tolerance = 1e-11)
      }
      expect_equal(as.vector(predict(unserialize(serialize(H, NULL)), newdata = e)),
                   as.vector(npreghat(b, newdata = e)), tolerance = 1e-11)
      expect_identical(predict(H), H)
      expect_error(predict(H, newdata = data.frame(not_x = e$x)), "must contain columns")
      incomplete <- e
      incomplete$x[2L] <- NA_real_
      expect_equal(as.vector(predict(H, newdata = incomplete)),
                   as.vector(npreghat(b, newdata = incomplete)), tolerance = 1e-11)
    }
  }
})
