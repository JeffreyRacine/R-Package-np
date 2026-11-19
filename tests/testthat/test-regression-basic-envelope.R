test_that("basic envelopes invert independently sorted endpoint-unit roots", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(98031);x <- data.frame(x=runif(90));y <- sin(3*x$x)+rnorm(90,sd=.3)
  bw <- npregbw(xdat=x,ydat=y,bws=.25,regtype="ll",bandwidth.compute=FALSE)
  state <- .np_bias_regression_prepare(bw,x,y,np_bias_control(),adapt=TRUE,law="inid")
  grid <- data.frame(x=seq(.15,.85,length.out=26L));fit <- .np_plot_regression_eval(bws=bw,xdat=x,ydat=y,exdat=grid,need.bias.capture=TRUE)
  B <- 199L;alpha <- c(.03,.13,.24)
  z <- .np_with_seed(1901,.np_bias_inference(state,grid,B=B,law="inid",fit=fit,
    alpha=alpha,calibration="basic-envelope",resolution.warn=FALSE))
  W <- do.call(rbind,lapply(split(seq_len(nrow(grid)),ceiling(seq_len(nrow(grid))/25)),function(ii)
    .np_bias_regression_rows(state,grid[ii,,drop=FALSE],data.frame(row=seq_along(ii),target="mean"))$W))
  F <- sweep(W,2L,state$residual,"*");F <- F-rowMeans(F)
  counts <- .np_bias_stream(z$rng$start,rmultinom(B,nrow(x),rep(1/nrow(x),nrow(x)))-1)
  raw <- F %*% counts$value
  root <- z$se*(raw/apply(raw,1L,sd))
  expect_identical(counts$end,z$rng$end)
  for(a in seq_along(alpha)) {
    k <- 1L+floor(B*alpha[a]/2);kb <- 1L+floor(B*alpha[a]/(2*nrow(grid)))
    sorted <- t(apply(root,1,sort));abs.sorted <- t(apply(abs(root),1,sort))
    sim <- c(sort(apply(root,2,min))[k],sort(apply(root,2,max))[B+1L-k])
    abs.sim <- sort(apply(abs(root),2,max))[ceiling(B*(1-alpha[a]))]
    for(band in names(z$bands)) {
      v <- switch(band,pointwise=sorted[,c(k,B+1L-k)],
        simultaneous=matrix(rep(sim,each=nrow(grid)),ncol=2L),
        bonferroni=cbind(pmin(sorted[,kb],sim[1]),pmax(sorted[,B+1L-kb],sim[2])))
      q <- switch(band,pointwise=abs.sorted[,ceiling(B*(1-alpha[a]))],
        simultaneous=rep(abs.sim,nrow(grid)),
        bonferroni=pmax(abs.sorted[,ceiling(B*(1-alpha[a]/nrow(grid)))],abs.sim))
      lo <- pmin(v[,1],-q);hi <- pmax(v[,2],q)
      expect_equal(unname(z$bands[[band]][,,a]),unname(cbind(z$center-hi,z$center-lo)),tolerance=1e-12)
    }
  }
  expect_true(all(z$bands$simultaneous[,1,]<=z$bands$pointwise[,1,]))
  expect_true(all(z$bands$simultaneous[,2,]>=z$bands$pointwise[,2,]))
  expect_true(all(z$bands$bonferroni[,1,]<=z$bands$simultaneous[,1,]))
  expect_true(all(z$bands$bonferroni[,2,]>=z$bands$simultaneous[,2,]))
  # Mirroring the response must mirror basic endpoints using the same counts.
  mirror <- .np_bias_regression_prepare(bw,x,-y,np_bias_control(),adapt=TRUE,law="inid")
  mfit <- .np_plot_regression_eval(bws=bw,xdat=x,ydat=-y,exdat=grid,need.bias.capture=TRUE)
  mz <- .np_with_seed(1901,.np_bias_inference(mirror,grid,B=B,law="inid",fit=mfit,
    alpha=alpha,calibration="basic-envelope",resolution.warn=FALSE))
  for(band in names(z$bands))
    expect_equal(unname(mz$bands[[band]]),unname(-z$bands[[band]][,2:1,,drop=FALSE]),tolerance=1e-12)
  expect_match(z$root.units,"endpoint units")
  expect_error(.np_bias_inference(state,grid,law="mammen",calibration="basic-envelope"),"basic-envelope requires inid, fixed or geom")
})

test_that("public basic envelopes preserve objects, state and explicit admission", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(98032);x <- data.frame(x=runif(90));y <- x$x+rnorm(90,sd=.3)
  bw <- npregbw(xdat=x,ydat=y,bws=.25,regtype="ll",bandwidth.compute=FALSE)
  fit <- npreg(bws=bw,txdat=x,tydat=y);saved <- fit;rng <- .Random.seed
  args <- list(x=fit,errors="bootstrap",center="bias-corrected",bootstrap="inid",B=199L,
    band="all",neval=4L,output="data",boot.control=np_boot_control(calibration="basic-envelope"))
  out <- suppressWarnings(do.call(plot,args));d <- as.data.frame(out,band="simultaneous")
  expect_equal(nrow(d),4L);expect_true(all(is.finite(d$lower)));expect_true(all(d$lower<=d$upper))
  expect_identical(fit,saved);expect_identical(.Random.seed,rng)
  for(law in "wild") {
    bad <- args;bad$bootstrap <- law
    expect_error(do.call(plot,bad),"basic-envelope requires inid, fixed or geom and fixed-bandwidth regression")
  }
  bad <- args;bad$gradients <- TRUE
  derivative <- suppressWarnings(do.call(plot,bad))
  expect_true(all(is.finite(as.data.frame(derivative,band="simultaneous")$lower)))
  bad <- args;bad$band <- "pmzsd"
  expect_error(do.call(plot,bad),"empirical bands")
  bad <- args;bad$center <- "estimate"
  expect_error(do.call(plot,bad),"iterated bias-corrected regression")
  expect_identical(.Random.seed,rng)
})


test_that("basic envelopes keep zero certificates ahead of scale alignment", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(0,1,length.out=60));y <- sin(3*x$x)+.2*cos(19*x$x)
  bw <- npregbw(xdat=x,ydat=y,bws=.3,regtype="ll",bandwidth.compute=FALSE)
  state <- .np_bias_regression_prepare(bw,x,y,np_bias_control(),adapt=FALSE,law="inid")
  grid <- data.frame(x=seq(.1,.9,length.out=26L))
  fit <- .np_plot_regression_eval(bws=bw,xdat=x,ydat=y,exdat=grid,need.bias.capture=TRUE)
  owner <- .np_bias_regression_rows;inconsistent <- FALSE
  fake <- function(spec,grid,map,...) {
    z <- owner(spec,grid,map,...)
    zero <- which(grid[[1L]][map$row] %in% c(.1,.9))
    z$se[zero] <- 0
    if(!inconsistent) {z$W[zero,] <- 0;z$factor[zero,] <- 0}
    z
  }
  local_mocked_bindings(.np_bias_regression_rows=fake,
    .package=getNamespaceName(environment(.np_bias_inference)))
  z <- .np_with_seed(1981,.np_bias_inference(state,grid,B=99,law="inid",fit=fit,
    calibration="basic-envelope",resolution.warn=FALSE))
  for(band in names(z$bands))for(side in 1:2)
    expect_equal(unname(z$bands[[band]][c(1,26),side,1]),z$center[c(1,26)])
  inconsistent <- TRUE
  expect_error(.np_with_seed(1981,.np_bias_inference(state,grid,B=99,law="inid",fit=fit,
    calibration="basic-envelope",resolution.warn=FALSE)),"zero target scale with nonzero signed process draws")
})
