test_that("npqreg prediction and plot evaluation leave workers usable", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old <- options(
    np.messages = FALSE,
    np.plot.progress = FALSE,
    npRmpi.autodispatch = TRUE
  )
  on.exit(options(old), add = TRUE)

  set.seed(20260513)
  n <- 80L
  x <- sort(runif(n, -1, 1))
  z <- factor(sample(c("a", "b"), n, TRUE))
  y <- sin(pi * x) + rnorm(n, sd = 0.18)
  cls <- factor(ifelse(y + 0.4 * (z == "b") > 0, "hi", "lo"),
                levels = c("lo", "hi"))
  dat <- data.frame(y = y, cls = cls, x = x, z = z)
  nd <- data.frame(
    x = seq(-0.8, 0.8, length.out = 11L),
    z = factor(rep(c("a", "b"), length.out = 11L), levels = levels(z))
  )

  q_bw <- npcdistbw(y ~ x, data = dat, bws = c(0.30, 0.25),
                    bwscaling = FALSE, bandwidth.compute = FALSE)
  q_fit <- npqreg(bws = q_bw, tau = c(0.25, 0.5, 0.75),
                  gradients = TRUE)

  expect_silent(predict(q_fit, newdata = nd, se.fit = TRUE))
  expect_silent(plot(q_fit, tau = c(0.25, 0.5, 0.75), output = "data"))
  expect_silent(plot(q_fit, gradients = TRUE, tau = 0.5, output = "data"))

  cm_bw <- npcdensbw(cls ~ x + z, data = dat, bws = c(0.35, 0.25, 0.15),
                     bwscaling = FALSE, bandwidth.compute = FALSE)
  cm_fit <- npconmode(bws = cm_bw, probabilities = TRUE, gradients = TRUE,
                      level = "hi", proper = TRUE)

  expect_true(all(is.finite(as.vector(gradients(cm_fit)))))
})
