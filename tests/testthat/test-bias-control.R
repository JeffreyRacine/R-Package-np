test_that("higher derivative rows differentiate weights and preserve polynomial moments", {
  u <- seq(-1,1,length.out=24L); h <- .7
  Q <- cbind(1,u,u^2); k <- exp(-u^2/(2*h^2))
  ratios <- np:::.np_bias_gaussian_derivatives(u/h,h,2L)
  expect_equal(ratios[,3L],((u/h)^2-1)/h^2)
  V <- rbind(c(1,0,0),c(0,1,0),c(0,0,2))
  H <- np:::.np_bias_weighted_derivative(Q,k,V,ratios)
  expect_equal(drop(H %*% Q),c(0,0,2),tolerance=1e-11)
  expect_null(np:::.np_bias_derivative_request(list(degree.engine=1L),1L,1L))
  expect_identical(np:::.np_bias_derivative_request(list(degree.engine=c(3L,4L)),2L,c(2L,3L)),c(2L,3L))
  expect_error(np:::.np_bias_derivative_request(list(degree.engine=2L),1L,2L),"degree - 1",fixed=TRUE)
})

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
  expect_error(np:::.np_bias_validate_control(control), "invalid np_bias_control")
})

test_that("corrected plot families nest and singleton limits coincide", {
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
  expect_identical(out$inference$studentization$construction, "residual-pilot-v1")
  expect_identical(out$inference$settings$iterations, 16L)
  expect_true(out$inference$settings$bernstein.basis)
  expect_identical(out$inference$settings$main$degree, bw$degree)
  b <- out$inference$bands
  expect_true(all(b$bonferroni[,1,] <= b$simultaneous[,1,]))
  expect_true(all(b$simultaneous[,1,] <= b$pointwise[,1,]))
  expect_true(all(b$pointwise[,2,] <= b$simultaneous[,2,]))
  expect_true(all(b$simultaneous[,2,] <= b$bonferroni[,2,]))
  state <- np:::.np_bias_regression_prepare(bw, x, y)
  for (law in c("gaussian", "mammen", "rademacher", "inid", "fixed", "geom")) {
    z <- np:::.np_with_seed(42L, np:::.np_bias_inference(
      state, x[1L,,drop=FALSE], "mean", .05, 63L, law))
    expect_identical(z$bands$pointwise, z$bands$simultaneous)
    expect_identical(z$bands$pointwise, z$bands$bonferroni)
  }
})

test_that("iterated intervals describe maximum calibration without relabeling legacy bands", {
  label <- np:::.np_plot_all_band_legend_labels
  descriptor <- np:::.np_plot_variability_descriptor
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
  expect_true(all(grepl("studentized", text, fixed=TRUE)))
  expect_false(any(grepl("normal quantile", text, fixed=TRUE)))
  expect_match(label("asymptotic", "iterated")[["pointwise"]], "quantiles", fixed=TRUE)
  expect_match(label("asymptotic")[["pointwise"]], "SE", fixed=TRUE)
  expect_match(descriptor("asymptotic", "pmzsd", .05, "bias-corrected",
    construction="iterated")$candidates[1], "normal quantile", fixed=TRUE)
})


test_that("auxiliary degrees resolve by coordinate without altering the main model", {
  resolve <- np:::.np_bias_degree_spec
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
  expect_identical(np:::.np_bias_validate_control(old), old)
  for (version in list(1, 2, "1", "2", NA_integer_, TRUE, c(1L,2L))) {
    edited <- np_bias_control()
    edited$version <- version
    expect_error(np:::.np_bias_validate_control(edited), "invalid np_bias_control")
  }
})

test_that("target blocks cover both axes without changing target-major indices", {
  for (ng in c(1L,26L)) for (nt in c(20L,21L,501L)) {
    blocks <- np:::.np_bias_target_blocks(ng,nt)
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
  state <- np:::.np_bias_regression_prepare(bw,X,y)
  anchor <- state$x[1L,,drop=FALSE]
  anchor$lx <- seq(min(state$x$lx),max(state$x$lx),length.out=50L)[2L]
  anchor$f <- factor("a",levels=levels(f))
  blocks <- np:::.np_bias_projector(state,anchor)$blocks
  expect_length(blocks,2L)
  expect_true(any(vapply(blocks,function(b)any(1-rowSums(b$Q*t(b$coefficient))<=0),logical(1))))
  for(b in blocks) {
    expect_named(b,c("indices","Q","coefficient"))
    expect_true(all(is.finite(b$coefficient)))
  }
})

test_that("regularized Bernstein GNN rows retain the fitted model and correction covariance", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(1501080011);x <- data.frame(x1=runif(800),x2=runif(800))
  y <- sin(3*x$x1)+x$x2^2
  bw <- npregbw(xdat=x,ydat=y,bws=c(7,4),regtype="lp",degree=c(1,6),
    bernstein.basis=TRUE,bwtype="generalized_nn",bandwidth.compute=FALSE)
  saved <- bw;grid <- data.frame(x1=c(.34,.58),x2=c(.5,.9))
  state <- np:::.np_bias_regression_prepare(bw,x,y)
  map <- expand.grid(row=1:2,target=c("mean","d1","d2"),stringsAsFactors=FALSE)
  value <- np:::.np_bias_regression_rows(state,grid,map)
  fit <- npreg(bw,txdat=x,tydat=y,exdat=grid,gradients=TRUE)
  expect_identical(bw,saved)
  expect_true(any(value$main.regularization$ridge>0))
  expect_equal(value$original_center,c(fitted(fit),as.vector(gradients(fit))),tolerance=1e-8)
  expect_equal(drop(value$W%*%(1+rowSums(as.matrix(x)^2))),
    c(1+rowSums(as.matrix(grid)^2),as.vector(2*as.matrix(grid))),tolerance=1e-7)
  covariance <- matrix(0,nrow(map),nrow(map))
  for(j in seq_len(nrow(x)))covariance <- covariance+outer(value$W[,j],value$W[,j])*state$variance[j]
  expect_equal(tcrossprod(value$factor),covariance,tolerance=1e-10)
  for(law in c("mammen","rademacher","inid","fixed","geom")) {
    set.seed(61028);rng <- .Random.seed
    expect_error(np:::.np_bias_inference(state,grid,c("mean","d1","d2"),
      .05,63L,law,if(law%in%c("fixed","geom"))3L else NULL),
      "bootstrap calibration.*remains unqualified")
    expect_identical(.Random.seed,rng)
  }
})

test_that("streamed calibration uses whole-family sizes and common-draw maxima", {
  set.seed(61024);groups <- rep(c("mean","d1","d2"),each=26L)
  absolute <- abs(matrix(rnorm(78*63),78,63));alpha <- c(.01,.05,.1)
  state <- np:::.np_bias_limits_prepare(groups,alpha,FALSE,63L)
  for(b in np:::.np_bias_target_blocks(26L,3L))
    state <- np:::.np_bias_limits_update(state,absolute[b$rows,,drop=FALSE],b$rows)
  result <- np:::.np_bias_limits_finish(state)
  for(g in unique(groups)) {
    ii <- which(groups==g)
    maximum <- apply(absolute[ii,,drop=FALSE],2L,max)
    critical <- unname(quantile(maximum,1-alpha,type=1L))
    expect_identical(unname(result$raw_max_critical[,g]),critical)
    for(j in ii) {
      point <- unname(quantile(absolute[j,],1-alpha,type=1L))
      bonf <- pmax(critical,unname(quantile(absolute[j,],1-alpha/length(ii),type=1L)))
      expect_identical(unname(result$limits$pointwise[j,2,]),point)
      expect_identical(unname(result$limits$bonferroni[j,2,]),bonf)
    }
  }
})


test_that("public adaptation controls preserve frozen and explicit constructions", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(42)
  x <- data.frame(x=runif(60,-2,2)); y <- x$x+rnorm(60)
  bw <- npregbw(xdat=x,ydat=y,bws=10000,bwscaling=TRUE,
    regtype="ll",bandwidth.compute=FALSE)
  saved <- bw
  out <- plot(bw,errors="asymptotic",center="bias-corrected",band="all",
    neval=3L,gaussian.draws=63L,output="data")[[1L]]$inference
  expect_true(out$settings$adaptation$active)
  expect_identical(bw,saved)
  expect_true(all(is.finite(out$center)))
  expect_equal(out$center,out$estimate,tolerance=1e-6)
  expect_true(out$settings$bandwidth$p0 > out$settings$reference_bandwidth$p0)
  for(control in list(np_bias_control(adaptation="none"),
                      np_bias_control(p0.bandwidth=.9),np_bias_control(p0.scale=4))) {
    state <- .np_bias_regression_prepare(bw,x,y,control,adapt=TRUE)
    expect_null(state$adaptation)
    expect_false(state$settings$adaptation$active)
  }
  old.control <- np_bias_control();old.control$version <- 2L;old.control$adaptation <- NULL
  expect_identical(.np_bias_resolve(old.control,x,60L,1L)$adaptation_request,"none")
  bad <- np_bias_control();bad$adaptation <- c("auto","none")
  expect_error(.np_bias_validate_control(bad),"invalid bias adaptation")
})


test_that("studentization streams frozen draw blocks and preserves zero targets", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(27052)
  x <- data.frame(x=runif(64));y <- sin(3*x$x)+rnorm(64,sd=.2)
  bw <- npregbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE,regtype="ll")
  state <- np:::.np_bias_regression_prepare(bw,x,y)
  prepared <- np:::.np_bias_student_prepare(state)
  for (response.scale in c(1e-8,1e8)) {
    scaled <- state
    for (name in c("y","residual_raw","residual"))
      scaled[[name]] <- state[[name]]*response.scale
    scaled$variance <- state$variance*response.scale^2
    check <- np:::.np_bias_student_prepare(scaled)
    expect_equal(check$rmu/response.scale,prepared$rmu,tolerance=1e-10)
  }
  inconsistent <- state
  inconsistent$residual[1L] <- inconsistent$residual[1L]+.01*max(abs(state$y))
  expect_error(np:::.np_bias_student_prepare(inconsistent),"disagrees")
  map <- expand.grid(row=1:2,target=c("mean","d1"),stringsAsFactors=FALSE)
  value <- np:::.np_bias_regression_rows(state,x[1:2,,drop=FALSE],map)
  W <- rbind(value$W,0);se <- c(value$se,0)
  target <- np:::.np_bias_student_target(prepared,W,se)
  for(law in c("gaussian","mammen")) {
    set.seed(27053)
    expected <- cbind(np:::.np_bias_student_chunk(prepared,W,target,256L,law),
      np:::.np_bias_student_chunk(prepared,W,target,1L,law))
    rng <- .Random.seed
    set.seed(27053)
    actual <- np:::.np_bias_student_draws(prepared,W,se,257L,law)
    expect_identical(actual,expected)
    expect_identical(.Random.seed,rng)
    expect_true(all(actual[nrow(W),]==0))
  }
  bad <- target;bad$A[,] <- 0
  expect_error(np:::.np_bias_student_chunk(prepared,W,bad,2L,"gaussian"),
    "nonpositive replicate studentization scale")
})


test_that("pooled profile admission is explicit and weighted projection is reproducible", {
  x <- data.frame(u=factor(rep(1:2,each=8)), v=factor(rep(1:2,8)), z=seq_len(16)/16)
  b <- list(icon=c(FALSE,FALSE,TRUE),iuno=c(TRUE,TRUE,FALSE),iord=rep(FALSE,3))
  expect_true(np:::.np_bias_profile_layout(b,x))
  b$iuno[2] <- FALSE; b$iord[2] <- TRUE
  expect_false(np:::.np_bias_profile_layout(b,x))
  x$v <- factor(rep(1:4,4))
  expect_true(np:::.np_bias_profile_layout(b,x))
  Q <- cbind(1,x$z,x$z^2);g <- exp(-x$z^2);y <- sin(x$z)
  p <- np:::.np_bias_profile_solve(Q,g)
  expect_equal(drop(p$coefficient %*% y),unname(lm.wfit(Q,y,g)$coefficients),tolerance=1e-10)
  expect_equal(p$coefficient %*% Q,diag(ncol(Q)),tolerance=1e-10)
  expect_error(np:::.np_bias_profile_solve(cbind(Q,Q[,2]),g),"numerically unidentified")
})
