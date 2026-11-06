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

# Policy-table test only; actual hidden/exposed native generator fixtures are
# exercised in separate bounded subprocesses in the R13 forensic evidence.
test_that("opaque uniform refusal leaves exposed-state and legacy policy alone", {
  ns <- asNamespace("npRmpi")
  guard <- get(".np_plot_rng_replay_preflight", ns)
  env <- new.env(parent = environment(guard))
  env$RNGkind <- function() c("user-supplied", "Inversion", "Rejection")
  environment(guard) <- env
  present <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  prior <- if (present) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  on.exit(if (present) assign(".Random.seed", prior, .GlobalEnv) else
    rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
  method <- get(".np_plot_rbandwidth_engine", ns)
  bws <- structure(list(), class = "rbandwidth")
  dots <- list(plot.errors.center = "bias-corrected",
               plot.errors.method = "bootstrap")
  assign(".Random.seed", 10405L, .GlobalEnv)
  expect_error(guard(method, bws, dots), "custom uniform state absent")
  expect_identical(get(".Random.seed", .GlobalEnv), 10405L)
  dots$plot.errors.boot.calibration <- "legacy"
  expect_silent(guard(method, bws, dots))
  dots$plot.errors.boot.calibration <- NULL
  assign(".Random.seed", c(10405L, 1L, 2L), .GlobalEnv)
  expect_silent(guard(method, bws, dots))
})

# Retained values extend the same count contract, not the data resolver.
test_that("retained training classification preserves explicit suppliedness", {
  ns <- asNamespace("npRmpi")
  n <- get(".np_plot_replay_training_n", ns)
  b <- list(.np.native.training = list(xdat = data.frame(x = 1:5), ydat = 2:6))
  expect_identical(n(list(), b), 5L)
  expect_identical(n(list(xdat = 1:2, ydat = 2:3), b), 2L)
  expect_null(n(list(xdat = NULL), b))
  expect_null(n(list(ydat = NULL), b))
  frame <- model.frame(y ~ x, data.frame(x = 1:5, y = 2:6))
  b$formula <- y ~ x
  b$.np.formula.training <- list(frame = frame, na.action = NULL)
  expect_identical(n(list(), b), 5L)
  b$.np.formula.training$frame$x[1] <- NA
  expect_null(n(list(), b))
  b$.np.formula.training$frame <- structure(frame, class = c("custom", "data.frame"))
  expect_null(n(list(), b))
  b$.np.formula.training <- NULL
  expect_null(n(list(), b))
})

test_that("explicit CDF replay retains covariance and deterministic boundaries", {
  ns <- asNamespace("npRmpi")
  guard <- get(".np_plot_rng_replay_preflight", ns)
  method <- get(".np_plot_dbandwidth_engine", ns)
  b <- structure(list(type="fixed", ncon=1L, nuno=0L, nord=0L,
    ckertype="gaussian", ckerorder=2L, ckerbound="none"), class="dbandwidth")
  d <- list(plot.errors.center="bias-corrected", plot.errors.method="asymptotic",
    plot.errors.type="simultaneous", bias.control=np_udist_bias_control("cdf"))
  old <- RNGkind(); on.exit(do.call(RNGkind, as.list(old)), add=TRUE)
  RNGkind(normal.kind="Box-Muller")
  for (m in c(31L,32L,33L)) {
    d$neval <- m
    if(m>32L) expect_error(guard(method,b,d),"unsafe RNG replay") else expect_silent(guard(method,b,d))
  }
  for(band in c("pointwise","bonferroni","pmzsd")) {
    d$plot.errors.type <- band; expect_silent(guard(method,b,d))
  }
  d$plot.errors.type <- "simultaneous"; d$bias.control$version <- 99L
  expect_silent(guard(method,b,d))
})

test_that("explicit density replay retains default and covariance providers", {
  ns <- asNamespace("npRmpi")
  guard <- get(".np_plot_rng_replay_preflight", ns)
  method <- get(".np_plot_bandwidth_engine", ns)
  b <- structure(list(type="fixed",ncon=1L,nuno=0L,nord=0L,
    ckertype="gaussian",ckerorder=2L,ckerbound="none"),class="bandwidth")
  d <- list(plot.errors.center="bias-corrected",plot.errors.method="asymptotic",
    plot.errors.type="simultaneous",bias.control=np_udens_bias_control())
  old <- RNGkind(); on.exit(do.call(RNGkind,as.list(old)),add=TRUE)
  RNGkind(normal.kind="Box-Muller")
  for(m in c(511L,512L,513L)) {
    d$neval <- m
    if(m>512L) expect_error(guard(method,b,d),"unsafe RNG replay") else expect_silent(guard(method,b,d))
  }
  d$bias.control <- NULL; expect_silent(guard(method,b,d))
  d$bias.control <- np_udens_bias_control("positive-covariance")
  expect_silent(guard(method,b,d))
  d$bias.control$version <- 2L; expect_error(guard(method,b,d),"unsafe RNG replay")
  d$plot.errors.method <- "bootstrap"; expect_silent(guard(method,b,d))
})
