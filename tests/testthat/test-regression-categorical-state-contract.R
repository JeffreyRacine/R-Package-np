test_that("all-categorical fixed regression kernels do not retain stale state under MPI", {
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "npRmpi subprocess install unavailable")

  lines <- c(
    "library(npRmpi)",
    "npRmpi.init(nslaves = 1, quiet = TRUE)",
    "cat_state_design <- function(seed, n = 192L) {",
    "  set.seed(seed)",
    "  xdat <- data.frame(",
    "    z1 = factor(sample(letters[1:3], n, TRUE)),",
    "    z2 = factor(sample(letters[1:4], n, TRUE)),",
    "    z3 = factor(sample(letters[1:2], n, TRUE))",
    "  )",
    "  mu <- as.integer(xdat$z1) + 0.3 * as.integer(xdat$z2) -",
    "    0.2 * as.integer(xdat$z3)",
    "  y <- mu + rnorm(n, sd = 0.25 * sd(mu))",
    "  list(xdat = xdat, y = y, dat = data.frame(y = y, xdat))",
    "}",
    "cat_state_dense_fit <- function(lambda, xdat, y, ukertype) {",
    "  codes <- as.data.frame(lapply(xdat, as.integer))",
    "  ncat <- vapply(xdat, nlevels, integer(1))",
    "  W <- matrix(1, nrow(xdat), nrow(xdat))",
    "  for (j in seq_along(codes)) {",
    "    same <- outer(codes[[j]], codes[[j]], `==`)",
    "    Kj <- switch(ukertype,",
    "      aitchisonaitken = ifelse(same, 1 - lambda[j], lambda[j] / (ncat[j] - 1)),",
    "      liracine = ifelse(same, 1, lambda[j])",
    "    )",
    "    W <- W * Kj",
    "  }",
    "  as.numeric(W %*% y / rowSums(W))",
    "}",
    "cat_state_fit <- function(d, lambda, ukertype) {",
    "  bw <- npregbw(",
    "    y ~ z1 + z2 + z3, data = d$dat, bws = lambda,",
    "    bandwidth.compute = FALSE, bwscaling = FALSE, regtype = 'lc',",
    "    ukertype = ukertype, okertype = 'liracine')",
    "  fitted(npreg(bws = bw))",
    "}",
    "target <- cat_state_design(20260520)",
    "target_lambda <- rep(0.1, 3L)",
    "target_ref <- cat_state_dense_fit(target_lambda, target$xdat, target$y, 'liracine')",
    "stopifnot(max(abs(cat_state_fit(target, target_lambda, 'liracine') - target_ref)) < 1e-10)",
    "warm <- cat_state_design(20260519)",
    "invisible(cat_state_fit(warm, rep(0, 3L), 'aitchisonaitken'))",
    "stopifnot(max(abs(cat_state_fit(target, target_lambda, 'liracine') - target_ref)) < 1e-10)"
  )

  res <- npRmpi_run_rscript_subprocess(lines, timeout = 60L, env = env)
  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
})

