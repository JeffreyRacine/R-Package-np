test_that("bias controls preserve explicit construction and auxiliary basis", {
  control <- np_bias_control()
  expect_s3_class(control, "np_bias_control")
  expect_identical(control$iterations, 16L)
  expect_true(control$bernstein.basis)
  expect_identical(np_boot_control()$wild, "mammen")
  expect_identical(np_boot_control(wild="rademacher")$wild, "rademacher")
  expect_error(np_bias_control(degree=3, bernstein.basis=FALSE), "Bernstein|bernstein")
  expect_error(np_bias_control(iterations=1.2), "positive integer")
  expect_error(np_bias_control(p0.scale=4, p0.bandwidth=.2), "mutually exclusive")
  expect_error(np_bias_control(construction="legacy", degree=5), "do not apply")
  control$construction <- "typo"
  expect_error(npRmpi:::.np_bias_validate_control(control), "invalid np_bias_control")
})

test_that("corrected plot families nest and singleton limits coincide", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(27051)
  x <- data.frame(x=runif(80))
  y <- sin(3*x$x) + rnorm(80, sd=.2)
  bw <- npregbw(xdat=x, ydat=y, bws=.23, bandwidth.compute=FALSE,
                regtype="lp", degree=2, bernstein.basis=TRUE)
  fit <- npreg(bw, txdat=x, tydat=y)
  original <- fit
  rng <- .Random.seed
  out <- plot(fit, errors="asymptotic", center="bias-corrected", band="all",
              gaussian.draws=63L, neval=3L, output="data")[[1L]]
  expect_identical(fit, original)
  expect_identical(.Random.seed, rng)
  expect_true(out$inference$settings$bernstein.basis)
  expect_identical(out$inference$settings$main$degree, bw$degree)
  b <- out$inference$bands
  expect_true(all(b$bonferroni[,1,] <= b$simultaneous[,1,]))
  expect_true(all(b$simultaneous[,1,] <= b$pointwise[,1,]))
  expect_true(all(b$pointwise[,2,] <= b$simultaneous[,2,]))
  expect_true(all(b$simultaneous[,2,] <= b$bonferroni[,2,]))
  state <- npRmpi:::.np_bias_regression_prepare(bw, x, y)
  for (law in c("gaussian", "mammen", "rademacher", "inid", "fixed", "geom")) {
    z <- npRmpi:::.np_with_seed(42L, npRmpi:::.np_bias_inference(
      state, x[1L,,drop=FALSE], "mean", .05, 63L, law))
    expect_identical(z$bands$pointwise, z$bands$simultaneous)
    expect_identical(z$bands$pointwise, z$bands$bonferroni)
  }
})

test_that("iterated intervals describe maximum calibration without relabeling legacy bands", {
  label <- npRmpi:::.np_plot_all_band_legend_labels
  descriptor <- npRmpi:::.np_plot_variability_descriptor
  expect_match(label("bootstrap")[["simultaneous"]], "rank-based", fixed=TRUE)
  expect_match(label("bootstrap", "iterated")[["simultaneous"]], "maximum", fixed=TRUE)
  expect_match(label("bootstrap", "iterated")[["bonferroni"]], "envelope", fixed=TRUE)
  text <- descriptor("bootstrap", "simultaneous", .05, "bias-corrected",
                     construction="iterated")$candidates
  expect_true(all(grepl("maximum", text, fixed=TRUE)))
  expect_false(any(grepl("rank", text, fixed=TRUE)))
  text <- descriptor("asymptotic", "pointwise", .05, "bias-corrected",
                     construction="iterated")$candidates
  expect_true(all(grepl("bias-corrected", text, fixed=TRUE)))
})


test_that("auxiliary degrees resolve by coordinate without altering the main model", {
  resolve <- npRmpi:::.np_bias_degree_spec
  spec <- resolve(np_bias_control(), c("a", "b", "c"), c(1, 2, 6), 200)
  expect_identical(spec$effective, c(a=5L, b=5L, c=6L))
  expect_identical(spec$main, c(a=1L, b=2L, c=6L))
  spec <- resolve(np_bias_control(degree=c(b=4, a=2)), c("a", "b"), c(1, 3), 100)
  expect_identical(spec$effective, c(a=2L, b=4L))
  expect_error(resolve(np_bias_control(degree=2), c("a", "b"), c(1, 3), 100), "cover")
  expect_error(resolve(np_bias_control(), paste0("x",1:10), rep(1,10), 100), "residual degrees")
  expect_error(np_bias_control(degree=c(a=2,a=3)), "unique names")
  expect_error(np_bias_control(degree=0), "positive integers")
  old <- np_bias_control(degree=5L)
  old$version <- 1L
  expect_identical(npRmpi:::.np_bias_validate_control(old), old)
  for (version in list(1, 2, "1", "2", NA_integer_, TRUE, c(1L,2L))) {
    edited <- np_bias_control()
    edited$version <- version
    expect_error(npRmpi:::.np_bias_validate_control(edited), "invalid np_bias_control")
  }
})

test_that("target blocks cover both axes without changing target-major indices", {
  for (ng in c(1L,26L)) for (nt in c(20L,21L,501L)) {
    blocks <- npRmpi:::.np_bias_target_blocks(ng,nt)
    rows <- unlist(lapply(blocks, `[[`, "rows"))
    expect_identical(sort(as.integer(rows)), seq_len(ng*nt))
    expect_true(all(vapply(blocks,function(b) length(b$grid)<=25L && length(b$rows)<=500L,logical(1))))
  }
})

test_that("projector geometry does not depend on an unused HC3 denominator", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(42);n <- 120L;x <- runif(n);f <- factor(sample(c("a","b"),n,TRUE))
  y <- sin(4*x)+(f=="b")+rnorm(n,sd=.3);X <- data.frame(lx=log(x),f=f)
  bw <- npregbw(xdat=X,ydat=y,bws=c(.3,.3),regtype="ll",bandwidth.compute=FALSE)
  state <- npRmpi:::.np_bias_regression_prepare(bw,X,y)
  anchor <- state$x[which.min(state$x$lx),,drop=FALSE]
  blocks <- npRmpi:::.np_bias_projector(state,anchor)$blocks
  expect_length(blocks,2L)
  for(b in blocks) {
    expect_named(b,c("indices","Q","coefficient"))
    expect_equal(b$coefficient%*%b$Q,diag(ncol(b$Q)),tolerance=1e-7)
  }
})

test_that("regularized Bernstein GNN rows retain the fitted model and correction covariance", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(1501080011);x <- data.frame(x1=runif(800),x2=runif(800))
  y <- sin(3*x$x1)+x$x2^2
  bw <- npregbw(xdat=x,ydat=y,bws=c(7,4),regtype="lp",degree=c(1,6),
    bernstein.basis=TRUE,bwtype="generalized_nn",bandwidth.compute=FALSE)
  saved <- bw;grid <- data.frame(x1=c(.34,.58),x2=c(.5,.9))
  state <- npRmpi:::.np_bias_regression_prepare(bw,x,y)
  map <- expand.grid(row=1:2,target=c("mean","d1","d2"),stringsAsFactors=FALSE)
  value <- npRmpi:::.np_bias_regression_rows(state,grid,map)
  fit <- npreg(bw,txdat=x,tydat=y,exdat=grid,gradients=TRUE)
  expect_identical(bw,saved)
  expect_true(any(value$main.regularization$ridge>0))
  expect_equal(value$original_center,c(fitted(fit),as.vector(gradients(fit))),tolerance=1e-8)
  expect_equal(drop(value$W%*%(1+rowSums(as.matrix(x)^2))),
    c(1+rowSums(as.matrix(grid)^2),as.vector(2*as.matrix(grid))),tolerance=1e-7)
  covariance <- matrix(0,nrow(map),nrow(map))
  for(j in seq_len(nrow(x)))covariance <- covariance+outer(value$W[,j],value$W[,j])*state$variance[j]
  expect_equal(tcrossprod(value$factor),covariance,tolerance=1e-10)
})
