test_that("adaptive fold moments initialize row state on repeated calls", {
  old <- options(np.messages = FALSE, np.largeh = FALSE, np.tree = FALSE)
  on.exit(options(old))
  z <- data.frame(z = c(.02, .06, .13, .22, .35, .51, .64, .79, .98))
  y <- seq_len(nrow(z)) / nrow(z)
  bw <- npregbw(xdat = z, ydat = y, bws = 3,
                bandwidth.compute = FALSE, bwtype = "adaptive_nn")
  reference <- vapply(seq_len(nrow(z)), function(i) {
    npksum(txdat = z[-i, , drop = FALSE], exdat = z[i, , drop = FALSE],
           tydat = y[-i], bws = bw, bandwidth.divide = TRUE)$ksum[1L]
  }, numeric(1))
  for (iteration in seq_len(4L)) {
    actual <- .np_estimator_loo_ksum(txdat = z, tydat = y, bws = bw,
                                   leave.one.out = TRUE, bandwidth.divide = TRUE)$ksum
    expect_equal(as.double(actual), reference, tolerance = 2e-12)
  }
})

test_that("raw fold radius partitions retain parent approximation ownership", {
  old <- options(np.messages=FALSE,np.largeh=TRUE,np.largelambda=TRUE,np.tree=TRUE)
  on.exit(options(old))
  # A preceding ordinary fit refreshes the parent's default approximation policy.
  invisible(npreg(txdat=data.frame(x=c(0,1,2)),tydat=c(0,1,2),bws=.4))
  set.seed(542)
  n <- 23L
  z <- data.frame(a=sample(seq(.02,.97,length.out=n)),b=runif(n,.04,.96),
    u=factor(rep(letters[1:3],length.out=n)),o=ordered(rep(1:3,length.out=n)))
  z[2,1:2] <- z[1,1:2]
  y <- rnorm(n)
  indices <- c(21L,2L,1L,2L,5L)
  for(tree in c(FALSE,TRUE)) for(kernel in c("gaussian","epanechnikov")) {
    options(np.tree=tree)
    bw <- npregbw(xdat=z,ydat=y,bws=c(9,10,.2,.25),bwtype="adaptive_nn",
      bandwidth.compute=FALSE,ckertype=kernel,ckerbound="fixed",ckerlb=c(0,0),ckerub=c(1,1))
    all <- .np_kernel_weights_direct(bw,z,leave.one.out=TRUE)
    part <- .np_kernel_weights_direct(bw,z,z[indices,,drop=FALSE],fold.rows=indices)
    oracle <- vapply(indices,function(i) {
      out <- numeric(n)
      out[-i] <- .np_kernel_weights_direct(bw,z[-i,,drop=FALSE],z[i,,drop=FALSE])
      out
    },numeric(n))
    expect_equal(part,oracle,tolerance=2e-12,ignore_attr=TRUE)
    expect_equal(part,all[,indices],tolerance=2e-12,ignore_attr=TRUE)
  }
})

test_that("private NN fold weights and moments delete occurrences before geometry", {
  old <- options(np.messages = FALSE, np.largeh = FALSE, np.largelambda = FALSE, np.tree = TRUE)
  on.exit(options(old))
  set.seed(542)
  n <- 23L
  z <- data.frame(a = sample(seq(.02, .97, length.out = n)),
                  b = runif(n, .04, .96),
                  u = factor(rep(letters[1:3], length.out = n)),
                  o = ordered(rep(1:3, length.out = n)))
  z[2, 1:2] <- z[1, 1:2]
  y <- rnorm(n)
  w <- cbind(1, seq_len(n)/n)
  indices <- c(21L, 2L, 1L, 2L, 5L)
  for (tree in c(FALSE, TRUE)) for (type in c("generalized_nn", "adaptive_nn"))
    for (kernel in c("gaussian", "epanechnikov", "beta")) {
      options(np.tree = tree)
      bw <- npregbw(xdat = z, ydat = y, bws = c(9, 10, .2, .25),
                    bwtype = type, bandwidth.compute = FALSE, regtype = "lc",
                    ckertype = kernel, ckerbound = "fixed",
                    ckerlb = c(0, 0), ckerub = c(1, 1))
      for (power in c(1, 2)) {
        k <- .np_kernel_weights_direct(bw, z, leave.one.out = TRUE, kernel.pow = power)
        oracle <- .npRmpi_with_local_regression(vapply(seq_len(n), function(i) {
          row <- numeric(n)
          row[-i] <- .np_kernel_weights_direct(bw, z[-i, , drop = FALSE],
                                               z[i, , drop = FALSE], kernel.pow = power)
          row
        }, numeric(n)))
        expect_equal(k, oracle, tolerance = 2e-12, ignore_attr = TRUE)
        chunk <- .np_kernel_weights_direct(bw, z, z[indices, , drop = FALSE],
                                            fold.rows = indices, kernel.pow = power)
        expect_equal(chunk, oracle[, indices], tolerance = 2e-12, ignore_attr = TRUE)
        sums <- .np_estimator_loo_ksum(txdat = z, tydat = w, weights = w, bws = bw,
                     leave.one.out = TRUE, kernel.pow = power, bandwidth.divide = TRUE)$ksum
        expected <- array(0, c(2L, 2L, n))
        for (i in seq_len(n)) expected[, , i] <- .npRmpi_with_local_regression(npksum(txdat = z[-i, , drop = FALSE],
          exdat = z[i, , drop = FALSE], tydat = w[-i, , drop = FALSE],
          weights = w[-i, , drop = FALSE], bws = bw, kernel.pow = power,
          bandwidth.divide = TRUE)$ksum)
        expect_equal(sums, expected, tolerance = 2e-12, ignore_attr = TRUE)
        partial <- .np_estimator_loo_ksum(txdat = z, exdat = z[indices, , drop = FALSE],
          tydat = w, weights = w, bws = bw, kernel.pow = power,
          bandwidth.divide = TRUE, .np.internal.fold.train.index = indices)$ksum
        expect_equal(partial, expected[, , indices], tolerance = 2e-12, ignore_attr = TRUE)
      }
    }
})

test_that("smooth coefficient LOO hats fits and search share deleted raw moments", {
  old <- options(np.messages = FALSE, np.largeh = FALSE, np.largelambda = FALSE)
  on.exit(options(old))
  set.seed(920520)
  n <- 35L
  z <- data.frame(z = sort(runif(n, .04, .96)))
  x <- data.frame(x = runif(n))
  y <- sin(4*z$z)*x$x + rnorm(n, sd = .1)
  for (type in c("fixed", "generalized_nn", "adaptive_nn"))
    for (rt in c("lc", "ll", "lp")) {
      a <- list(xdat = x, zdat = z, ydat = y, bws = if(type == "fixed") .3 else 18,
                bandwidth.compute = FALSE, bwtype = type, regtype = rt)
      if (rt == "lp") a$degree <- 2L
      bw <- do.call(npscoefbw, a)
      fit <- npscoef(bw, txdat = x, tzdat = z, tydat = y, iterate = FALSE,
                     leave.one.out = TRUE, se = TRUE)
      hat <- npscoefhat(bw, txdat = x, tzdat = z, y = y, output = "apply",
                        leave.one.out = TRUE)
      literal <- .npRmpi_with_local_regression(vapply(seq_len(n), function(i) fitted(npscoef(bw,
        txdat = x[-i, , drop = FALSE], tzdat = z[-i, , drop = FALSE], tydat = y[-i],
        exdat = x[i, , drop = FALSE], ezdat = z[i, , drop = FALSE], iterate = FALSE)), numeric(1)))
      expect_equal(as.double(fitted(fit)), literal, tolerance = 2e-10)
      expect_equal(as.double(hat), literal, tolerance = 2e-10)
      H <- npscoefhat(bw, txdat=x, tzdat=z, output="matrix", leave.one.out=TRUE)
      expect_equal(as.vector(H %*% y), literal, tolerance=2e-10)
      expected.se <- sqrt(as.vector(H^2 %*% (y-literal)^2))
      expect_equal(as.double(se(fit)), expected.se, tolerance=2e-10)
      ctx <- list(n = n, W = cbind(1, x$x), ydat = y, zdat.df = z)
      objective <- .npscoefbw_nomad_eval_direct(ctx, bw)
      expect_true(objective$raw.valid)
      expect_equal(objective$objective, mean((y-literal)^2), tolerance = 2e-10,
                   info = paste(type, rt))
      idx <- c(31L, 3L, 1L, 3L)
      partial <- .npscoefbw_nomad_eval_subset(ctx, bw, idx, sqrt(.Machine$double.xmax))
      expect_identical(partial$invalid, 0L)
      expect_equal(partial$sse, sum((y[idx]-literal[idx])^2), tolerance = 2e-10,
                   info = paste(type, rt))
    }
})

test_that("raw fold counts use the deleted sample at NN boundaries", {
  old <- options(np.messages=FALSE, np.extendednn=TRUE, np.largeh=FALSE)
  on.exit(options(old))
  x <- data.frame(x=c(.02,.05,.12,.18,.4,.56,.68,.87,.98))
  n <- nrow(x)
  for (type in c("generalized_nn","adaptive_nn"))
    for (kernel in c("gaussian","beta"))
      for (k in c(2L,n-2L,n-1L,n,2L*n)) {
        bw <- npregbw(xdat=x,ydat=seq_len(n),bws=k,bwtype=type,
          bandwidth.compute=FALSE,ckertype=kernel,
          ckerbound="fixed",ckerlb=0,ckerub=1)
        actual <- .np_kernel_weights_direct(bw,x,leave.one.out=TRUE)
        oracle <- vapply(seq_len(n),function(i) {
          row <- numeric(n)
          row[-i] <- .np_kernel_weights_direct(bw,x[-i,,drop=FALSE],x[i,,drop=FALSE])
          row
        },numeric(n))
        expect_equal(actual,oracle,tolerance=2e-12,ignore_attr=TRUE,
                     info=paste(type,kernel,k))
      }
})

test_that("raw fold geometry errors leave the next kernel call usable", {
  old <- options(np.messages = FALSE)
  on.exit(options(old))
  x <- data.frame(x = c(0, 0, 0, 0, 1, 2, 3))
  bw <- npregbw(xdat = x, ydat = seq_len(nrow(x)), bws = 2,
                bwtype = "adaptive_nn", bandwidth.compute = FALSE)
  expect_error(.np_kernel_weights_direct(bw, x, leave.one.out = TRUE),
               "nearest|radius|delete-one")
  expect_error(.np_kernel_weights_direct(bw, x, x[1:2, , drop = FALSE],
                                         fold.rows = c(0L, 2L)), "occurrence map")
  expect_true(all(is.finite(npksum(txdat = x, bws = .4)$ksum)))
})

test_that("smooth coefficient CV rejects counts inadmissible after deletion", {
  old <- options(np.messages=FALSE, np.extendednn=FALSE)
  on.exit(options(old))
  n <- 9L
  x <- data.frame(x=seq_len(n)/n)
  y <- sin(x$x)
  for (type in c("generalized_nn", "adaptive_nn")) {
    bw <- npscoefbw(xdat=x, zdat=x, ydat=y, bws=n-1L,
                    bwtype=type, bandwidth.compute=FALSE)
    expect_error(.npscoef_nn_assert_training_radius(bw,x,"test CV"),
                 class="np_nn_candidate_invalid")
    result <- .npscoefbw_nomad_eval_direct(
      list(n=n,W=cbind(1,x$x),ydat=y,zdat.df=x),bw)
    expect_false(result$raw.valid)
    expect_true(is.finite(result$objective))
  }
})
