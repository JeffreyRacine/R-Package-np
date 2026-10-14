test_that("conmode surfaces require distinct resolved axes", {
  pkg <- getNamespaceName(environment(npconmode))
  resolve <- getFromNamespace(".np_plot_conmode_resolve_surface_vars", pkg)
  x <- data.frame(first = 1:3, second = 4:6, held = 7:9)
  object <- list(xnames = names(x), bws = list(ixcon = rep(TRUE, 3L)))
  expect_identical(resolve(x, object, c("second", "first")), c("second", "first"))
  expect_error(resolve(x, object, c("first", "first")), "two distinct conditioning variables", fixed = TRUE)
  expect_error(resolve(x, object, c(1, 1)), "character vector", fixed = TRUE)
  expect_error(resolve(x, object, c("absent", "first")), "fitted conmode object", fixed = TRUE)
  two <- x[1:2]
  object <- list(xnames = names(two), bws = list(ixcon = c(TRUE, TRUE)))
  expect_identical(resolve(two, object), names(two))
  names(two) <- object$xnames <- c("first", "first")
  expect_error(resolve(two, object), "two distinct conditioning variables", fixed = TRUE)
})

test_that("ordered-only ANN copula distribution retains empirical CDF covariance", {
  pkg <- getNamespaceName(environment(npcopula))
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  withr::local_seed(7421L)
  i <- seq_len(64L)
  x <- data.frame(a = ordered(as.character((i + floor(i/7)) %% 3 + 1), levels = as.character(1:3)),
                  b = ordered(as.character((7*i + floor(i/5)) %% 4 + 1), levels = as.character(1:4)))
  b <- npudistbw(dat = x, bws = c(.27, .43), bwtype = "adaptive_nn",
                 bandwidth.compute = FALSE, okertype = "racineliyan")
  object <- npcopula(bws = b, data = x, u = data.frame(a = c(.38, .68), b = c(.42, .73)),
                    n.quasi.inv = 24L, se = TRUE)
  ex <- object$eval[, names(x), drop = FALSE]
  # Independent normalized finite-support RLY CDF for each common donor.
  cdf <- function(train, eval, lambda) {
    support <- seq_along(levels(train))
    vapply(as.integer(train), function(donor) {
      mass <- lambda^abs(support - donor)
      mass <- mass/sum(mass)
      vapply(as.integer(eval), function(q) sum(mass[support <= q]), numeric(1))
    }, numeric(length(eval)))
  }
  K <- cdf(x$a, ex$a, .27) * cdf(x$b, ex$b, .43)
  oracle <- (K - rowMeans(K))/nrow(x)
  influence <- getFromNamespace(".np_plot_copula_influence", pkg)(object, x, ex, object$copulaerr)
  expect_identical(object$bws$type, "adaptive_nn")
  expect_equal(object$copula, rowMeans(K), tolerance = 2e-11)
  expect_equal(influence$A, oracle, tolerance = 2e-11)
  expect_equal(tcrossprod(influence$A), tcrossprod(oracle), tolerance = 2e-11)
  method <- getS3method("plot", "npcopula", envir = asNamespace(pkg))
  request <- function(band) {
    args <- list(x = object, view = "fixed", output = "data", errors = "asymptotic", band = band)
    if (band %in% c("all", "simultaneous"))
      args <- c(args, list(gaussian.draws = 9999L, random.seed = 19L))
    do.call(method, args)
  }
  pw <- request("pointwise")
  bf <- request("bonferroni")
  joint <- request("simultaneous")
  all <- request("all")
  expect_true(all(is.finite(c(joint$lower, joint$upper))))
  expect_identical(all$simultaneous.lower, joint$lower)
  expect_identical(all$simultaneous.upper, joint$upper)
  expect_identical(all$pointwise.lower, pw$lower)
  expect_identical(all$pointwise.upper, pw$upper)
  expect_identical(all$bonferroni.lower, bf$lower)
  expect_identical(all$bonferroni.upper, bf$upper)
  expect_identical(nrow(all), nrow(object$eval))
})
