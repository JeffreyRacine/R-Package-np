test_that("copula density SE includes the marginal covariance", {
  pkg <- getNamespaceName(environment(npcopula))
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  ns <- asNamespace(pkg)
  # Marginal construction already warns about the stored, inapplicable
  # uniform order. Do not hide any other warning, notably unavailable SEs.
  ignore.uniform.order <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (grepl("ignoring kernel order specified with uniform kernel type",
              conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
  raw.local <- get(".npcopula_density_se_fixed", ns)
  local <- function(...) ignore.uniform.order(raw.local(...))
  x <- data.frame(x = c(.1, .2, 1.2, 1.4, -1.3, -1.5),
                  z = c(.1, 1.2, .2, 1.4, .3, -1.5))
  q <- data.frame(x = 0, z = 0)
  bw <- npudensbw(dat = x, bws = c(1, 1), bandwidth.compute = FALSE,
                  ckertype = "uniform")
  # Centered influence is (2,-2,-1,1,-1,1), including donors
  # outside the joint support intersection.
  expect_equal(local(bw, x, q), sqrt(12)/6, tolerance = 2e-14)
  one <- x["x"]
  b1 <- npudensbw(dat = one, bws = 1, bandwidth.compute = FALSE,
                  ckertype = "uniform")
  expect_lt(max(abs(local(b1, one, q["x"]))), 1e-14)
  expect_warning(unavailable <- local(bw, x, data.frame(x = 100, z = 0)),
                  "marginal density is zero", fixed = TRUE)
  expect_identical(unavailable, NA_real_)

  # A public one-dimensional copula is identically one on supported rows.
  fitted <- ignore.uniform.order(npcopula(bws = b1, data = one, se = TRUE))
  expect_lt(max(abs(se(fitted))), 1e-14)
  expect_identical(predict(fitted, se.fit = TRUE)$se.fit, se(fitted))
  # Deterministic SEs must not consume random transport-tag draws in MPI.
  set.seed(417)
  saved <- .Random.seed
  invisible(local(bw, x, q))
  expect_identical(.Random.seed, saved)
  rm(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
  invisible(local(bw, x, q))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("native copula reduction is centered and validates dimensions", {
  pkg <- getNamespaceName(environment(npcopula))
  native <- function(x) .Call("C_np_copula_density_se", x, PACKAGE = pkg)
  oracle <- function(w) {
    u <- lapply(w, function(x) sweep(x, 2, colMeans(x), "/"))
    product <- Reduce("*", u)
    ratio <- colMeans(product)
    phi <- product - sweep(Reduce("+", u) - length(u) + 1, 2, ratio, "*")
    phi <- sweep(phi, 2, colMeans(phi), "-")
    sqrt(colSums(phi^2))/nrow(phi)
  }
  w <- list(matrix(c(-.1,.4,.8,.7),4,1),
            matrix(c(.3,.6,-.1,.9),4,1))
  expect_equal(native(w), oracle(w), tolerance = 2e-14)
  negative <- lapply(w, function(x) -x)
  expect_equal(native(negative), oracle(negative), tolerance = 2e-14)
  for (e in c(1e-3,1e-6,1e-9,1e-12)) {
    w <- list(matrix(1+e*sin(1:64),32,2), matrix(1+e*cos(1:64),32,2))
    expect_lt(max(abs(native(w)-oracle(w))), 16*.Machine$double.eps)
  }
  expect_error(native(list(matrix(1:4,4,1))), "double weight matrices")
  expect_error(native(list(matrix(1,4,1), matrix(1,3,1))), "dimensions differ")
  expect_identical(native(list(matrix(0,4,1))), NA_real_)
  expect_identical(native(list(matrix(Inf,4,1))), NA_real_)
})
