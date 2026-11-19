test_that("pairs correction uses the accepted fitted rows including regularization", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  p <- readRDS(test_path("fixtures", "pairs-fitted-failure.rds"))
  x <- data.frame(x=p$x)
  bw <- npregbw(xdat=x,ydat=p$y,bws=p$h,bandwidth.compute=FALSE,
    regtype="lp",degree=3L,bernstein.basis=TRUE)
  grid <- data.frame(x=seq(min(p$x),max(p$x),length.out=50L)[c(10L,25L,46L,47L,49L,50L)])
  fit <- npreg(bws=bw,txdat=x,tydat=p$y,exdat=grid,se=TRUE,gradients=TRUE,.np.fit.inference=TRUE)
  saved <- serialize(fit,NULL)
  targets <- c("mean","d1")
  map <- expand.grid(row=seq_len(nrow(grid)),target=targets,stringsAsFactors=FALSE)
  for (law in c("inid","fixed","geom")) {
    state <- npRmpi:::.np_bias_regression_prepare(bw,x,p$y,np_bias_control(adaptation="none"),law=law)
    expect_identical(state$main.capture,list(version="fixed-pairs-fitted-v1",law=law))
    expect_null(state$local.inference)
    expect_identical(state$variance,npRmpi:::.np_bias_variance(state,state$residual))
    capture <- npRmpi:::.np_bias_capture_prepare(state,fit,grid,targets)
    value <- capture$rows(seq_len(nrow(grid)),map)
    expect_true(any(value$main.regularization$ridge>0))
    expected <- c(fit$mean,fit$grad[,1L])
    expect_lt(max(abs(value$original_center-expected)),1e-13*(1+max(abs(expected))))
    set.seed(107L)
    inference <- npRmpi:::.np_bias_inference(state,grid,targets,B=39L,alpha=.2,law=law,
      blocklen=if(law=="inid")NULL else 1L,fit=fit,resolution.warn=FALSE)
    expect_true(all(is.finite(inference$center)))
    expect_equal(inference$estimate,expected,tolerance=1e-12)
    expect_identical(inference$law,law)
    other <- if(law=="inid")"fixed" else "inid"
    expect_error(npRmpi:::.np_bias_inference(state,grid,targets,B=39L,law=other,fit=fit),"match its prepared law")
    expect_error(npRmpi:::.np_bias_inference(state,grid,targets,B=39L,law=other,
      main.capture=capture),"match its prepared law")
    bad.version <- state; bad.version$main.capture$version <- "invalid"
    expect_error(npRmpi:::.np_bias_pairs_capture(bad.version),"invalid prepared pairs")
    bad.eligibility <- state; bad.eligibility$gradient.order <- 2L
    expect_error(npRmpi:::.np_bias_pairs_capture(bad.eligibility),"invalid prepared pairs")
    changed <- state;changed$main.capture$law <- "rademacher"
    expect_error(npRmpi:::.np_bias_inference(changed,grid,targets,B=39L,law="rademacher",fit=fit),"invalid prepared pairs")
    for (field in c("ckertype","type","gradient.order")) {
      unsupported <- state
      if(field=="gradient.order") unsupported$gradient.order <- 2L else
        unsupported$bws[[field]] <- if(field=="type")"generalized_nn" else "epanechnikov"
      expect_false(npRmpi:::.np_bias_main_capture_supported(unsupported,law))
    }
  }
  expect_identical(serialize(fit,NULL),saved)
})

# Protect limits at the existing operator seam; no large allocation is needed.
test_that("pairs capture preserves law-specific limits and legacy law refusals", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  i <- seq_len(40L)
  all.x <- data.frame(x=sin(i*sqrt(2)),z=cos(i*sqrt(3)))
  real.operator <- npRmpi:::.np_plot_regression_capture_operator
  seen <- numeric()
  testthat::local_mocked_bindings(.np_plot_regression_capture_operator =
    function(fit,bws,x,y,ex,deriv,categorical.index,max.bytes) {
      seen <<- c(seen,max.bytes)
      real.operator(fit,bws,x,y,ex,deriv,categorical.index,max.bytes)
    }, .package="npRmpi")
  for (q in 1:2) {
    x <- all.x[seq_len(q)]; y <- sin(x$x)+.1*cos(i)
    bw <- npregbw(xdat=x,ydat=y,bws=rep(.6,q),bandwidth.compute=FALSE,
      regtype="lp",degree=rep(1L,q),bernstein.basis=FALSE)
    grid <- x[c(10L,30L),,drop=FALSE]
    fit <- npreg(bws=bw,txdat=x,tydat=y,exdat=grid,se=TRUE,.np.fit.inference=TRUE)
    for (law in c("gaussian","mammen","inid","fixed","geom")) {
      state <- npRmpi:::.np_bias_regression_prepare(bw,x,y,
        np_bias_control(adaptation="none"),law=law)
      seen <- numeric()
      capture <- npRmpi:::.np_bias_capture_prepare(state,fit,grid,"mean")
      expect_length(seen,1L)
      if (law %in% c("gaussian","mammen")) {
        expect_identical(seen,64*1024^2)
        for (pairs in c("inid","fixed","geom"))
          expect_error(npRmpi:::.np_bias_inference(state,grid,"mean",B=39L,
            law=pairs,fit=fit),"Gaussian or Mammen")
      } else {
        width <- ncol(state$Q)
        expected <- 256*(nrow(x)*nrow(grid)+(nrow(x)+nrow(grid))*(q+width+1)+width^2)
        expect_identical(seen,as.double(expected))
        rows <- capture$rows(1:2,data.frame(row=1:2,target="mean"))
        expect_equal(rows$original_center,fit$mean,tolerance=1e-10)
      }
    }
  }
})
