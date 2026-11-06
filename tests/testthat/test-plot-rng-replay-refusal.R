test_that("unknown training classes are parked without numeric dispatch", {
  ns <- asNamespace("npRmpi")
  n <- get(".np_plot_replay_training_n", ns)
  env <- new.env(parent = environment(n))
  original <- base::is.numeric
  env$is.numeric <- function(x) {
    if (is.object(x)) stop("custom numeric dispatch must not run")
    original(x)
  }
  environment(n) <- env
  custom <- structure(1:3, class = "r13_unknown")
  expect_null(n(list(xdat = custom, ydat = 1:3)))
  expect_null(n(list(xdat = 1:3, ydat = custom)))
  expect_null(n(list(xdat = data.frame(x = c(1, NA, 3)), ydat = 1:3)))
  expect_null(n(list(xdat = data.frame(x = factor(c("a", NA, "b"))), ydat = 1:3)))
  expect_identical(n(list(xdat = data.frame(x = factor(c("a", "a", "b"))),
                         ydat = 1:3)), 3L)
})

# Unique default-lane contract: pre-seed refusal and retained cache boundaries.
# Broad public/MPI replay qualification is kept in the R13 forensic campaign.
test_that("known uncached replay is refused without consuming a pending normal", {
  ns <- asNamespace("npRmpi")
  guard <- get(".np_plot_rng_replay_preflight", ns)
  method <- get(".np_plot_rbandwidth_engine", ns)
  bws <- structure(list(), class = "rbandwidth")
  dots <- list(plot.errors.center = "bias-corrected",
    plot.errors.method = "asymptotic", gaussian.draws = 1025L,
    xdat = seq_len(2048), ydat = seq_len(2048))
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  setup <- function() {
    RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")
    set.seed(19); invisible(rnorm(1))
  }
  setup(); expected <- rnorm(7)
  setup(); state <- .Random.seed
  expect_error(guard(method, bws, dots), "unsafe RNG replay")
  expect_identical(.Random.seed, state)
  expect_identical(rnorm(7), expected)
  dots$gaussian.draws <- 1024L
  expect_silent(guard(method, bws, dots))
  dots$gaussian.draws <- 1023L
  expect_silent(guard(method, bws, dots))
  # Historical bws$nobs must not override explicit materialized data.
  bws$nobs <- 100000L
  dots$xdat <- dots$ydat <- 1:3
  expect_silent(guard(method, bws, dots))
})

test_that("CDF refusal respects target geometry and deterministic bands", {
  ns <- asNamespace("npRmpi")
  guard <- get(".np_plot_rng_replay_preflight", ns)
  method <- get(".np_plot_dbandwidth_engine", ns)
  bws <- structure(list(type = "fixed", ncon = 2L, nuno = 0L, nord = 0L),
                   class = "dbandwidth")
  dots <- list(plot.errors.center = "bias-corrected",
    plot.errors.method = "asymptotic", plot.errors.type = "simultaneous",
    neval = 6L)
  old <- RNGkind(); on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")
  expect_error(guard(method, bws, dots), "unsafe RNG replay")
  dots$neval <- 5L
  expect_silent(guard(method, bws, dots))
  dots$neval <- 26L; dots$perspective <- FALSE
  expect_error(guard(method, bws, dots), "unsafe RNG replay")
  dots$plot.errors.type <- "pointwise"
  expect_silent(guard(method, bws, dots))
  dots$plot.errors.type <- "simultaneous"
  RNGkind(normal.kind = "Inversion")
  expect_silent(guard(method, bws, dots))
})
