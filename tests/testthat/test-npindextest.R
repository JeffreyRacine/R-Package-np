test_that("index test is exported with registered methods and validates arguments", {
  ns <- asNamespace("npRmpi")
  expect_identical(getExportedValue("npRmpi","npindextest"), ns$npindextest)
  for (generic in c("print","summary"))
    expect_identical(getS3method(generic,"indextest",envir=ns),
                     get(paste0(generic,".indextest"),ns,inherits=FALSE))
  expect_error(npindextest(NULL,B=2L),"at least 9")
  expect_error(npindextest(NULL),"object must")
  expect_error(npindextest(NULL,bw.method="cv"),"arg")
})

test_that("the native index score agrees with the ordered-pair statistic", {
  x <- cbind(c(-1, 0, 1, .5, -.5), c(1, .2, -.5, 1.2, 1.4))
  y <- c(.5, 1, 1.5, -1, 2)
  z <- x[,1]
  g <- 0.4
  h <- c(.3,.5)
  dz <- outer(z,z,"-")
  u <- rowSums(outer(y,y,"-")*exp(-.5*(dz/g)^2))
  base <- tcrossprod(u)*exp(-.5*outer(x[,2],x[,2],"-")^2)
  diag(base) <- 0
  ref <- vapply(h,function(width) {
    a <- base*exp(-.5*(dz/width)^2)
    sqrt(5/4)*sum(a)/sqrt(2*sum(a*a))
  },numeric(1))
  expect_equal(.Call(C_np_indextest_score,y,z,x,g,h),ref,tolerance=1e-12)
})

test_that("bootstrap refits match an independently refitted response stream", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  set.seed(711)
  x <- data.frame(a = rnorm(45), b = rnorm(45), c = rnorm(45))
  y <- x$a + 0.4*x$b + rnorm(45)
  bw <- npindexbw(xdat = x, ydat = y, method = "ichimura", nmulti = 1L)
  data <- .indextest_setup(bw, "published")
  index <- .indextest_index(data, data$beta, data$bw)
  data$pilot <- .indextest_pilot(data, index)
  set.seed(412)
  reference <- vapply(seq_len(2L), function(i) {
    response <- .indextest_response(data)
    fit <- npindexbw(xdat = x, ydat = response, method = "ichimura", nmulti = 1L)
    .indextest_score(response, .indextest_index(data, fit$beta, fit$bw), data)
  }, numeric(1))
  set.seed(412)
  actual <- vapply(seq_len(2L), function(i)
    .npRmpi_with_local_bootstrap(.np_plot_with_local_compiled_eval(
      .indextest_refit(.indextest_response(data), data, i))), numeric(1))
  expect_false(reference[1L] == reference[2L])
  expect_equal(actual, reference, tolerance = 1e-8)
})
