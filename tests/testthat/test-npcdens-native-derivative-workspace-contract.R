test_that("npcdens categorical-response gradients preserve native derivative workspace", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260513)
  n <- 80L
  x <- sort(runif(n, -1, 1))
  z <- factor(sample(c("a", "b"), n, TRUE), levels = c("a", "b"))
  y <- sin(pi * x) + rnorm(n, sd = 0.18)
  cls <- factor(
    ifelse(y + 0.4 * (z == "b") > 0, "hi", "lo"),
    levels = c("lo", "hi")
  )
  dat <- data.frame(y = y, cls = cls, x = x, z = z)

  bw <- npcdensbw(
    cls ~ x + z,
    data = dat,
    bws = c(0.35, 0.25, 0.15),
    bwscaling = FALSE,
    bandwidth.compute = FALSE
  )

  fit_level <- function(level, gradients) {
    npcdens(
      txdat = dat[c("x", "z")],
      tydat = dat["cls"],
      exdat = dat[c("x", "z")],
      eydat = factor(rep(level, n), levels = levels(dat$cls)),
      bws = bw,
      gradients = gradients
    )
  }
  expect_gradient_payload <- function(fit) {
    expect_length(fit$condens, n)
    expect_identical(dim(fit$congrad), c(n, 2L))
    expect_identical(dim(fit$congerr), c(n, 2L))
    expect_true(all(is.finite(fit$condens)))
    expect_true(all(is.finite(fit$congrad)))
    expect_true(all(is.finite(fit$congerr)))
  }

  low.nograd <- fit_level("lo", gradients = FALSE)
  hi.grad.after.nograd <- fit_level("hi", gradients = TRUE)
  expect_length(low.nograd$condens, n)
  expect_gradient_payload(hi.grad.after.nograd)

  low.grad <- fit_level("lo", gradients = TRUE)
  hi.grad.after.grad <- fit_level("hi", gradients = TRUE)
  expect_gradient_payload(low.grad)
  expect_gradient_payload(hi.grad.after.grad)
})
