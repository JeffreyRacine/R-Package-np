test_that("pooled quadratic Gaussian GNN means and slopes reproduce their fitted physical rows", {
  skip_on_cran()
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  # Retained donors from the saved four-neighbour upper boundary.
  x <- data.frame(x=c(0.0034768893383443356, 0.970864646602422, 0.9729915813077241,
    0.9739157187286764, 0.9745722964871675, 0.9782712699379772,
    0.978393143741414, 0.9813846340402961, 0.9818960200063884,
    0.9849397426005453, 0.9850722078699619, 0.9864174819085747,
    0.9946712902747095, 0.9995494075119495))
  y <- sin(seq_len(nrow(x)))
  degree <- 2L
  k <- 5L
  train <- x
  train$a <- factor(rep(0:1, length.out=nrow(x)))
  train$b <- factor(rep(c(0,0,1,1), length.out=nrow(x)))
  bandwidth <- c(k,.5,.5)
  at <- train[8L,,drop=FALSE]
  h <- sort(abs(x$x-at$x))[k]
  t <- (x$x-at$x)/h
  sw <- exp(-t^2/4)
  design <- outer(t, 0:degree, "^")
  Q <- qr.solve(design*sw, diag(sw))
  L <- Q[2L,]/h
  bw <- npregbw(xdat=train, ydat=y, bws=bandwidth, bandwidth.compute=FALSE,
                regtype="lp", basis="glp", degree=degree, bernstein.basis=TRUE,
                bwtype="generalized_nn", ckertype="gaussian", ckerbound="fixed",
                ckerlb=-1, ckerub=1)
  fit <- npreg(bw, txdat=train, tydat=y, exdat=at, gradients=TRUE, se=FALSE)
  H <- npreghat(bw, txdat=train, exdat=at, s=1L, output="matrix")
  capture <- function() .np_plot_lsq_gradient_projection(bw,train,numeric(nrow(train)),at,
      gradient.order=1L,capture.gradient=TRUE)
  actual <- .npRmpi_with_local_regression(capture())
  expect_identical(dim(actual$gate), c(3L, ncol(train)))
  expect_identical(as.numeric(attr(H, "ridge.used")), as.numeric(actual$ridge))
  expect_equal(as.numeric(attr(H, "ridge.used")), 0)
  expect_lt(max(abs(H-L)), 1e-8)
  expect_lt(abs(fit$grad[1,1]-sum(L*y)), 1e-8)
  M <- npreghat(bw, txdat=train, exdat=at, s=0L, output="matrix")
  expect_lt(max(abs(M-Q[1L,])), 1e-8)
  expect_lt(abs(fit$mean[1L]-sum(Q[1L,]*y)), 1e-8)
})
