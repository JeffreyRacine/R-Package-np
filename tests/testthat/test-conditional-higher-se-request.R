higher_se_fixture <- function(cdf, degree, type) {
  x <- data.frame(x = seq(.08, .92, length.out = 28L))
  y <- data.frame(y = .4 + .2*sin(7*x$x) + .08*cos(seq_len(nrow(x))))
  ex <- data.frame(x = c(.72, .3, .5))
  ey <- data.frame(y = c(.42, .3, .5))
  adaptive <- identical(type, "adaptive_nn")
  bw <- do.call(if (cdf) npcdistbw else npcdensbw, list(
    xdat = x, ydat = y, bws = if (adaptive) c(15, 17) else c(.17, .2),
    bandwidth.compute = FALSE, bwscaling = FALSE, bwtype = type,
    regtype = "lp", degree = degree, bernstein.basis = FALSE))
  radius <- function(z, k) vapply(seq_along(z), function(i)
    sort(abs(z[-i] - z[i]))[k], numeric(1))
  hx <- if (adaptive) radius(x$x, 17L) else .2
  hy <- if (adaptive) radius(y$y, 15L) else .17
  B <- outer(x$x, 0:degree, `^`)
  expected <- vapply(seq_len(nrow(ex)), function(j) {
    w <- dnorm((x$x - ex$x[j])/hx)/hx
    u <- (ey$y[j] - y$y)/hy
    z <- if (cdf) pnorm(u*(sqrt(2)*0.7071067810)) else dnorm(u)/hy
    v <- max(0, sum(w*z^2)/sum(w) - (sum(w*z)/sum(w))^2)
    d <- c(rep(0, degree), factorial(degree))
    ell <- w*drop(B %*% solve(crossprod(B, B*w), d))
    sqrt(v*sum(ell^2))
  }, numeric(1))
  list(call = list(bws = bw, txdat = x, tydat = y, exdat = ex, eydat = ey,
                   gradients = TRUE, gradient.order = degree),
       expected = expected)
}

test_that("higher conditional SEs use the actual derivative row and working variance", {
  for (cdf in c(FALSE, TRUE)) for (degree in c(2L, 3L))
    for (type in c("fixed", "adaptive_nn")) {
      d <- higher_se_fixture(cdf, degree, type)
      f <- if (cdf) npcdist else npcdens
      on <- do.call(f, c(d$call, list(se = TRUE)))
      off <- do.call(f, c(d$call, list(se = FALSE)))
      expect_equal(as.numeric(on$congerr), d$expected, tolerance = 3e-10)
      expect_identical(on$congrad, off$congrad)
      expect_identical(fitted(on), fitted(off))
      expect_null(off$congerr)
    }
})

test_that("higher conditional uncertainty has a separate optional demand", {
  request <- getFromNamespace(".np_conditional_higher_se_request", "npRmpi")
  expect_null(request(FALSE, TRUE, "lp", 2L, TRUE))
  expect_null(request(TRUE, FALSE, "lp", 2L, TRUE))
  expect_null(request(TRUE, TRUE, "lp", 1L, TRUE))
  expect_null(request(TRUE, TRUE, "lp", 2L, TRUE, FALSE))
  expect_identical(request(TRUE, TRUE, "lp", c(2L, 1L, 3L),
                           c(TRUE, TRUE, FALSE)), c(TRUE, FALSE, FALSE))
})

test_that("unused higher uncertainty never reaches the supplementary hat owner", {
  d <- higher_se_fixture(FALSE, 2L, "fixed")
  point.owner <- getFromNamespace(".np_conditional_higher_hat", "npRmpi")
  testthat::local_mocked_bindings(
    .np_conditional_higher_hat = function(..., return.norm = TRUE) {
      if (return.norm) stop("unrequested higher hat norm")
      point.owner(..., return.norm = FALSE)
    },
    .package = "npRmpi")
  off <- do.call(npcdens, c(d$call, list(se = FALSE)))
  expect_null(off$congerr)
  first <- d$call
  first$gradient.order <- 1L
  expect_no_error(do.call(npcdens, c(first, list(se = TRUE))))
  levels <- d$call
  levels$gradients <- FALSE
  expect_no_error(do.call(npcdens, c(levels, list(se = TRUE))))
  selected <- getFromNamespace(".np_conditional_eval_selected", "npRmpi")
  masked <- selected(bws = d$call$bws, xdat = d$call$txdat, ydat = d$call$tydat,
    exdat = d$call$exdat, eydat = d$call$eydat, gradients = TRUE,
    gradient.order = 2L, se = TRUE, lp.higher.se.demand = FALSE)
  expect_equal(masked$congrad, off$congrad, tolerance = 0)
  expect_true(all(is.na(masked$congerr)))
})

test_that("conditional higher-SE metadata distinguishes invalid, zero and scaled rows", {
  combine <- getFromNamespace(".np_conditional_higher_se", "npRmpi")
  metadata <- cbind(variance = c(4, 0, 4, 4), log.scale = c(0, 0, 0, 2),
                    invalid = c(0, 0, 1, 0), owners = rep(1, 4))
  norms <- cbind(scale = rep(3, 4), sumsq = rep(25/9, 4), invalid = rep(0, 4))
  expect_equal(combine(metadata, norms, 4L), c(10, 0, NA, 10*exp(2)))
  metadata[1L, 4L] <- 2
  expect_error(combine(metadata, norms, 4L), "row/ownership mismatch")
  for (row in list(c(3, 4), c(3e200, 4e200), c(3e-200, 4e-200), c(0, 0))) {
    value <- .Call("C_np_reghat_row_norm", row, PACKAGE = "npRmpi")
    expect_equal(value[1L]*sqrt(value[2L]), abs(row[1L])*5/3)
    expect_identical(value[3L], 0)
  }
})
