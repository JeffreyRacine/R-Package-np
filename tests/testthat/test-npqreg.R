test_that("npqreg basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("cps71")
  cps71_sub <- cps71[1:50, ]
  # Quantile regression needs a condbandwidth object
  bw <- npcdistbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                                bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # Median regression
  model <- npqreg(bws=bw, tau=0.5)
  
  expect_s3_class(model, "qregression")
  expect_type(predict(model), "double")
  expect_equal(length(predict(model)), 50)
  expect_equal(model$tau, 0.5)
  expect_error(gradients(model), "fit the model with gradients=TRUE", fixed = TRUE)
  expect_error(gradients(model, errors = TRUE), "fit the model with gradients=TRUE", fixed = TRUE)
  
  expect_output(summary(model))
})

test_that("npqreg works with multiple taus", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("cps71")
  cps71_sub <- cps71[1:30, ]
  bw <- npcdistbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                                bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  model_q25 <- npqreg(bws=bw, tau=0.25)
  expect_equal(model_q25$tau, 0.25)
})

test_that("npqreg gradients populate qregression objects", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(70507)
  n <- 150L
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  ydat <- xdat$x1 - 0.5 * xdat$x2 + rnorm(n, sd = 0.15)
  exdat <- data.frame(
    x1 = seq(-0.40, 0.40, length.out = 5L),
    x2 = seq(0.35, -0.35, length.out = 5L)
  )
  bw <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.80, 0.80, 0.80),
    bandwidth.compute = FALSE
  )

  fit <- npqreg(bws = bw, txdat = xdat, tydat = ydat, exdat = exdat,
                tau = 0.45, gradients = TRUE, tol = 1e-6, small = 1e-7)
  ref <- getFromNamespace(".np_plot_quantile_eval", "npRmpi")(
    bws = bw,
    txdat = xdat,
    tydat = ydat,
    exdat = exdat,
    tau = 0.45,
    gradients = TRUE,
    tol = 1e-6,
    small = 1e-7
  )

  expect_s3_class(fit, "qregression")
  expect_true(fit$gradients)
  expect_equal(dim(fit$quantgrad), c(nrow(exdat), ncol(xdat)))
  expect_equal(dim(fit$quantgerr), c(nrow(exdat), ncol(xdat)))
  expect_true(all(is.finite(fit$quantgrad)))
  expect_true(all(is.finite(fit$quantgerr)))
  expect_true(all(is.finite(se(fit))))
  expect_true(all(se(fit) > 0))
  expect_equal(fit$quantgrad, ref$quantgrad, tolerance = 0)
  expect_equal(gradients(fit), fit$quantgrad, tolerance = 0)
  expect_equal(gradients(fit, errors = TRUE), fit$quantgerr, tolerance = 0)

  qdelta <- getFromNamespace(".npqreg_quantile_delta_from_conditional", "npRmpi")(
    bws = bw,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    quantile = fitted(fit),
    gradients = TRUE
  )
  expect_equal(se(fit), qdelta$quanterr, tolerance = 0)
  expect_equal(fit$quantgrad, qdelta$quantgrad, tolerance = 0)
  expect_equal(fit$quantgerr, qdelta$quantgerr, tolerance = 0)

  eps <- 1e-2
  for (j in seq_along(exdat)) {
    xp <- xm <- exdat
    xp[[j]] <- xp[[j]] + eps
    xm[[j]] <- xm[[j]] - eps
    qp <- fitted(npqreg(bws = bw, txdat = xdat, tydat = ydat, exdat = xp,
                        tau = 0.45, gradients = FALSE, tol = 1e-6, small = 1e-7))
    qm <- fitted(npqreg(bws = bw, txdat = xdat, tydat = ydat, exdat = xm,
                        tau = 0.45, gradients = FALSE, tol = 1e-6, small = 1e-7))
    fd <- (qp - qm) / (2 * eps)
    expect_lt(max(abs(fit$quantgrad[, j] - fd), na.rm = TRUE), 1e-3)
  }
})
