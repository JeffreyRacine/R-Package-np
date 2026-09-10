dc1_conditional_fixture <- function(cdf, ridge, partial = FALSE) {
  x <- c(rep(0, 7L), rep(.5, 8L), 2)
  if (!ridge) x[7L] <- .25
  tx <- data.frame(x = x)
  ex <- data.frame(x = .2)
  degree <- 2L
  h <- .35
  if (partial) {
    tx$z <- seq(.1, .9, length.out = length(x))
    ex$z <- .5
    degree <- c(2L, 0L)
    h <- c(h, 1)
  }
  ty <- data.frame(y = c(-.6, -.2, .1, .4, .7, -.4, .2, .5,
                         -.5, -.1, .3, .6, .8, -.3, .15, .45))
  ey <- data.frame(y = .1)
  hy <- .3
  b <- do.call(if (cdf) npcdistbw else npcdensbw,
    list(xdat = tx, ydat = ty, bws = c(hy, h), bandwidth.compute = FALSE,
         bwscaling = FALSE, bwtype = "fixed", regtype = "lp",
         degree = degree, basis = "glp", bernstein.basis = FALSE,
         cxkertype = "epanechnikov", cykertype = "gaussian"))

  # Independent kernel, design, moments and accepted-map oracle. The CDF
  # retains its historical 0.7071067810 literal rather than exact 1/sqrt(2).
  w <- rep(1, length(x))
  for (j in seq_along(h)) {
    u <- (tx[[j]] - ex[[j]])/h[j]
    w <- w * ifelse(abs(u) < sqrt(5),
                    3/(4*sqrt(5))*(1-u^2/5), 0)/h[j]
  }
  z <- if (cdf) pnorm(sqrt(2)*0.7071067810*(ey$y-ty$y)/hy) else
    dnorm((ey$y-ty$y)/hy)/hy
  B <- cbind(1, x, x^2)
  A <- crossprod(B, B*w)
  M2 <- crossprod(B, B*w^2)
  rhs <- drop(crossprod(B, w*z))
  r <- if (ridge) max(abs(diag(A)))/length(x) else 0
  C <- solve(A + diag(r, ncol(B)))
  D <- diag(ncol(B))
  D[1L, 1L] <- 1+r/A[1L, 1L]
  d <- cbind(level = c(1, ex$x, ex$x^2), first = c(0, 1, 2*ex$x))
  a <- D %*% t(C) %*% d
  sigma2 <- sum(w*z^2)/sum(w)-(sum(w*z)/sum(w))^2
  variance <- sigma2*diag(crossprod(a, M2 %*% a))
  list(bws = b, txdat = tx, tydat = ty, exdat = ex, eydat = ey,
       point = drop(crossprod(d, C %*% D %*% rhs)),
       error = sqrt(variance), ridge = r,
       accepted.rcond = rcond(A + diag(r, ncol(B))))
}

dc1_expected_warnings <- function(expr, partial) {
  expected <- paste0(
    c("[np] npcdens", "[np] npcdist", "[npRmpi] npcdens", "[npRmpi] npcdist"),
    ": requested derivative order would exceed polynomial degree or is unavailable for ",
    "continuous[2] (requested order 1, degree 0); returning NA for unavailable gradient component(s)")
  seen <- unexpected <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    message <- conditionMessage(w)
    if (partial && message %in% expected) {
      seen <<- c(seen, message)
      invokeRestart("muffleWarning")
    }
    unexpected <<- c(unexpected, message)
  })
  expect_length(unexpected, 0L)
  if (partial) expect_true(length(seen) > 0L)
  value
}

test_that("conditional GENERAL covariance uses its accepted full/selected map", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  for (cdf in c(FALSE, TRUE)) for (ridge in c(FALSE, TRUE))
    for (partial in c(FALSE, TRUE)) dc1_expected_warnings({
      d <- dc1_conditional_fixture(cdf, ridge, partial)
      fun <- if (cdf) npcdist else npcdens
      args <- d[c("bws", "txdat", "tydat", "exdat", "eydat")]
      fit <- do.call(fun, c(args,
        list(gradients = TRUE, gradient.order = 1L, se = TRUE)))
      off <- do.call(fun, c(args,
        list(gradients = TRUE, gradient.order = 1L, se = FALSE)))
      expect_equal(as.numeric(fitted(fit)), unname(d$point[1L]),
                   tolerance = 2e-11)
      expect_equal(as.numeric(gradients(fit)[, 1L]), unname(d$point[2L]),
                   tolerance = 2e-11)
      expect_equal(as.numeric(se(fit)), unname(d$error[1L]),
                   tolerance = 2e-11)
      expect_equal(as.numeric(gradients(fit, se = TRUE)[, 1L]),
                   unname(d$error[2L]), tolerance = 2e-11)
      expect_identical(fitted(fit), fitted(off))
      expect_identical(gradients(fit), gradients(off))
      expect_true(d$accepted.rcond > 1e-5)
      if (partial) {
        expect_true(all(is.na(gradients(fit)[, 2L])))
        expect_true(all(is.na(gradients(fit, se = TRUE)[, 2L])))
        none <- do.call(fun, c(args,
          list(gradients = TRUE, se = TRUE, .np_lp_first_se_demand = FALSE)))
        expect_identical(fitted(none), fitted(fit))
        expect_identical(se(none), se(fit))
        expect_identical(gradients(none), gradients(fit))
        expect_true(all(is.na(gradients(none, se = TRUE))))
      }
    }, partial)
})
