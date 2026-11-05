test_that("smooth coefficient plot terms identify the selected coefficient", {
  if (TRUE) {
    skip_on_cran()
    if (!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES", "1"))))
      stop("required MPI pool unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(11)
  n <- 90L
  d <- data.frame(a = rnorm(n), b = rnorm(n), z = runif(n))
  d$y <- 100 + 10*d$a - 10*d$b + rnorm(n, sd = .03)
  X <- cbind(Intercept = 1, as.matrix(d[c("a", "b")]))
  bw <- npscoefbw(xdat = d[c("a", "b")], zdat = d["z"], ydat = d$y,
                  bws = .2, regtype = "lc", bandwidth.compute = FALSE)
  fit <- npscoef(bws = bw, txdat = d[c("a", "b")], tzdat = d["z"],
                 tydat = d$y, betas = TRUE)
  expect_identical(colnames(coef(fit)), colnames(X))
  original <- serialize(fit, NULL)
  rng <- .Random.seed
  snapshots <- list()
  for (i in 1:3) {
    out <- plot(fit, neval = 5L, perspective = FALSE, coef = TRUE,
                coef_index = i, output = "data")
    expect_true(all(out$targets$term == colnames(X)[i]))
    expect_identical(.Random.seed, rng)
    expect_identical(serialize(fit, NULL), original)
    for (k in seq_along(out$panels)) {
      panel <- out$panels[[k]]
      zkey <- panel$variables$coordinate[panel$variables$name == "z"]
      ze <- panel$grid[[zkey]]
      truth <- vapply(ze, function(at) {
        w <- exp(-.5*((d$z-at)/.2)^2)
        solve(crossprod(X, X*w), crossprod(X, d$y*w))[i]
      }, numeric(1L))
      dat <- as.data.frame(out, panel = k)
      expect_equal(dat$estimate, truth, tolerance = 1e-9)
    }
    snapshots[[as.character(i)]] <- out
  }
  expect_error(plot(fit, coef = TRUE, gradients = TRUE, output = "data"),
               "cannot be combined")
  if (identical(Sys.getenv("N5_FULL"), "true")) {
    formula.fit <- npscoef(y ~ a + b | z, data = d, bws = .2,
                          regtype = "lc", bandwidth.compute = FALSE, betas = TRUE)
    for (i in 1:3) {
      out <- plot(formula.fit, neval = 5L, perspective = FALSE, coef = TRUE,
                  coef_index = i, output = "data")
      expect_identical(out$targets$term, snapshots[[as.character(i)]]$targets$term)
      for (k in seq_along(out$panels))
        expect_equal(out$panels[[k]]$values, snapshots[[as.character(i)]]$panels[[k]]$values,
                     tolerance = 1e-12)
    }
    for (gradient in c(FALSE, TRUE)) {
      out <- plot(fit, neval = 5L, perspective = FALSE, gradients = gradient,
                  coef_index = 1L, output = "data")
      expect_true(all(is.na(out$targets$term)))
      snapshots[[paste0("gradient-", gradient)]] <- out
    }
    out <- plot(fit, neval = 5L, perspective = FALSE, coef = TRUE,
                coef_index = 4L, output = "data")
    expect_true(all(out$targets$term == "b"))
    for (k in seq_along(out$panels))
      expect_identical(out$panels[[k]]$values, snapshots[["3"]]$panels[[k]]$values)
  }
  if (nzchar(Sys.getenv("N5_SNAPSHOT"))) saveRDS(snapshots, Sys.getenv("N5_SNAPSHOT"))
})
