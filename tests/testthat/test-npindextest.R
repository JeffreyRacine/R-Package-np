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

test_that("guarded bootstrap CV respects search controls and the native criterion", {
  old <- options(np.tree = "auto")
  on.exit(options(old), add = TRUE)
  recipe <- list(optim.method = "Nelder-Mead")
  expect_true(.indextest_grid_eligible(160L, recipe))
  expect_false(.indextest_grid_eligible(159L, recipe))
  expect_false(.indextest_grid_eligible(160L, list(optim.method = "BFGS")))
  expect_false(.indextest_grid_eligible(160L, c(recipe, optim.reltol = 1e-12)))
  options(np.tree = TRUE)
  expect_false(.indextest_grid_eligible(160L, recipe))
  expect_false(.indextest_replay(list(fval = 1, fval.vector = c(1, 2))))
  expect_true(.indextest_replay(list(fval = 1, fval.vector = c(1, 1+1e-12))))
  expect_true(.indextest_replay(list(fval = 1.1, fval.vector = 1)))

  z <- seq(-2, 2, length.out = 160L)
  h <- .3
  kernel <- exp(-.5*(outer(z,z,"-")/h)^2)
  diag(kernel) <- 0
  counts <- new.env()
  counts$grid <- 0
  for (method in c("ichimura", "kleinspady")) {
    binary <- method == "kleinspady"
    y <- if (binary) as.double(seq_along(z) %% 3L == 0L) else sin(z)
    mean.loo <- drop(kernel %*% y)/rowSums(kernel)
    expected <- if (binary) -mean(ifelse(y == 1, log(mean.loo), log1p(-mean.loo))) else
      mean((mean.loo-y)^2)
    leaf <- .indextest_grid_leaf(function(...) stop("native trial"), method, counts)
    expect_equal(leaf(z, y, h, NULL, NULL, NULL)$objective, expected, tolerance = 1e-10)
    expect_error(leaf(z, y, h/1e6, NULL, NULL, NULL), "native trial")
    tiny <- .Machine$double.xmin * .Machine$double.eps
    if (tiny > 0) {
      near.node <- c(0, tiny, seq(.01, 2, length.out = 158L))
      expect_error(leaf(near.node, y, h, NULL, NULL, NULL), "native trial")
    }
    if (!binary) expect_error(leaf(z, 1e80*y, h, NULL, NULL, NULL), "native trial")
  }
  expect_equal(counts$grid, 2)
})

test_that("Ichimura grid loss preserves small variation around a large response level", {
  # The last row has positive donor mass below the sparse-row threshold.
  z <- c(seq(-1, 1, length.out = 156L), 1.5, 2.25, 3, 4.5)
  h <- .25
  # Dyadic responses preserve these differences exactly under either shift.
  variation <- round(256 * sin(2*z))/16384
  logw <- -.5 * (outer(z, z, "-")/h)^2
  diag(logw) <- -Inf
  last.mass <- sum(exp(logw[length(z), ]))
  expect_gt(last.mass, 1e-14)
  expect_lt(last.mass, 1e-3)
  weights <- exp(logw - apply(logw, 1L, max))
  residual <- rowSums(weights * outer(variation, variation, "-"))/rowSums(weights)
  expected <- mean(residual^2)
  counts <- new.env()
  counts$grid <- 0
  for (level in c(0, 1000, -1000)) {
    response <- level + variation
    native <- function(index, ydat, ...) {
      expect_identical(ydat, response)
      list(objective = -1)
    }
    leaf <- .indextest_grid_leaf(native, "ichimura", counts)
    value <- leaf(z, response, h, NULL, NULL, NULL)
    expect_lt(abs(value$objective - expected), 1e-13)
    # The distant last row forces late delegation after the FFT.
    sparse <- z
    sparse[length(sparse)] <- 8
    expect_identical(leaf(sparse, response, h, NULL, NULL, NULL)$objective, -1)
  }
  expect_equal(counts$grid, 3)
})
