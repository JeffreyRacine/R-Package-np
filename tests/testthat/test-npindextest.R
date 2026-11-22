test_that("index test is exported with registered methods and validates arguments", {
  ns <- asNamespace("np")
  expect_identical(getExportedValue("np","npindextest"), ns$npindextest)
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
  # The reference stays in ordinary units. These binary scalings are exact
  # for this response, including the subnormal case.
  for (units in c(1, 1e80, 1e-80, 2^1022, -2^1022, 2^-1073)) {
    expect_equal(.Call(C_np_indextest_score, units*y, z, x, g, h),
                 ref, tolerance = 1e-12)
  }
  expect_error(.Call(C_np_indextest_score, rep(1e300, length(y)), z, x, g, h),
               "Degenerate single-index pairwise variance")
})

test_that("bootstrap refits match an independently refitted response stream", {
  set.seed(42)
  x <- as.data.frame(matrix(rnorm(72L * 3L), 72L, 3L))
  y <- x[[1L]] + 0.4*x[[3L]] + rnorm(72L)
  fit.response <- function(v)
    npindexbw(xdat = x, ydat = v, method = "ichimura", bws = c(1,.1,.3,.4),
              optim.method = "CG", nmulti = 2L, random.seed = 32L,
              optim.maxit = 30L, scale.factor.search.lower = .05)
  bw <- fit.response(y)
  data <- .indextest_setup(bw, "published")
  index <- .indextest_index(data, data$beta, data$bw)
  data$pilot <- .indextest_pilot(data, index)
  # The pre-streaming full-matrix pilot is an independent small-sample oracle.
  logw <- -0.5*(outer(index$z, index$z, "-")/index$g)^2
  diag(logw) <- -Inf
  weights <- exp(logw - apply(logw, 1L, max))
  oracle <- data
  oracle$pilot <- drop(weights %*% y / rowSums(weights))
  set.seed(42)
  reference <- vapply(seq_len(2L), function(i) {
    response <- .indextest_response(oracle)
    fit <- fit.response(response)
    .indextest_score(response, .indextest_index(data, fit$beta, fit$bw), data)
  }, numeric(1))
  set.seed(42)
  actual <- vapply(seq_len(2L), function(i)
    .indextest_refit(.indextest_response(data), data, i), numeric(1))
  expect_false(reference[1L] == reference[2L])
  expect_equal(actual, reference, tolerance = 1e-8)
})
