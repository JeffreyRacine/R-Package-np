extract_nomad_case_block <- function(lines, case) {
  start <- grep(sprintf("^CASE_START %s$", case), lines)
  done <- grep(sprintf("^CASE_DONE %s$", case), lines)

  if (!length(start) || !length(done)) {
    return(character())
  }

  lines[seq.int(start[[1L]], done[[1L]])]
}

expect_nomad_unknown_bound_output <- function(lines, case) {
  block <- extract_nomad_case_block(lines, case)

  expect_true(length(block) > 0L, info = case)
  expect_false(any(grepl("nomad\\+powell", block, ignore.case = TRUE)), info = case)
  expect_false(any(grepl("eval [0-9]+", block)), info = case)
  expect_false(any(grepl("fval=", block, fixed = TRUE)), info = case)
  expect_false(any(grepl("%|eta ", block)), info = case)
  expect_true(any(grepl("^\\[npRmpi\\] Selecting degree and bandwidth \\(", block)), info = case)
  expect_true(any(grepl("multistart [12]/2", block)), info = case)
  expect_true(any(grepl("iteration [0-9]+", block)), info = case)
  expect_true(any(grepl("deg \\(", block)), info = case)
  expect_true(any(grepl("best \\(", block)), info = case)
}

test_that("remaining MPI NOMAD families use unknown-bound restart progress", {
  skip_if_not_installed("crs")

  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    c(
      "options(",
      "  np.messages = TRUE,",
      "  np.tree = FALSE,",
      "  npRmpi.autodispatch = TRUE,",
      "  np.progress.start.grace.known.sec = 0,",
      "  np.progress.start.grace.unknown.sec = 0,",
      "  np.progress.interval.known.sec = 0,",
      "  np.progress.interval.unknown.sec = 0",
      ")",
      "library(npRmpi)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "ns <- asNamespace('npRmpi')",
      "set_binding <- function(name, value) {",
      "  was_locked <- bindingIsLocked(name, ns)",
      "  if (was_locked) unlockBinding(name, ns)",
      "  assign(name, value, envir = ns)",
      "  if (was_locked) lockBinding(name, ns)",
      "}",
      "old <- list(",
      "  .np_progress_is_interactive = get('.np_progress_is_interactive', envir = ns, inherits = FALSE),",
      "  .np_progress_renderer_for_surface = get('.np_progress_renderer_for_surface', envir = ns, inherits = FALSE),",
      "  .np_progress_now = get('.np_progress_now', envir = ns, inherits = FALSE)",
      ")",
      "on.exit({",
      "  for (nm in names(old)) set_binding(nm, old[[nm]])",
      "}, add = TRUE)",
      "vals <- seq(0, 60, by = 0.25)",
      "i <- 0L",
      "set_binding('.np_progress_is_interactive', function() TRUE)",
      "set_binding('.np_progress_renderer_for_surface', function(surface, capability) 'legacy')",
      "set_binding('.np_progress_now', function() {",
      "  i <<- min(i + 1L, length(vals))",
      "  vals[[i]]",
      "})",
      "set.seed(20260401)",
      "n <- 18L",
      "x <- data.frame(x = sort(runif(n)))",
      "z <- data.frame(z = sort(runif(n)))",
      "y_reg <- sin(2 * pi * x$x) + rnorm(n, sd = 0.05)",
      "y_pl <- 1 + 0.5 * x$x + cos(2 * pi * z$z) + rnorm(n, sd = 0.05)",
      "y_sc <- (1 + z$z^2) * x$x + rnorm(n, sd = 0.05)",
      "y_idx <- sin(x$x + 0.5 * z$z) + rnorm(n, sd = 0.05)",
      "cases <- list(",
      "  npcdistbw = function() npcdistbw(",
      "    y_reg ~ x,",
      "    data = data.frame(y_reg = y_reg, x = x$x),",
      "    regtype = 'lp',",
      "    degree.select = 'coordinate',",
      "    search.engine = 'nomad+powell',",
      "    degree.min = 0L,",
      "    degree.max = 1L,",
      "    degree.verify = FALSE,",
      "    bwtype = 'fixed',",
      "    bwmethod = 'cv.ls',",
      "    nmulti = 2L,",
      "    ngrid = 30L,",
      "    max.bb.eval = 8L",
      "  ),",
      "  npplregbw = function() npplregbw(",
      "    xdat = x,",
      "    zdat = z,",
      "    ydat = y_pl,",
      "    regtype = 'lp',",
      "    bernstein.basis = TRUE,",
      "    degree.select = 'coordinate',",
      "    search.engine = 'nomad+powell',",
      "    degree.min = 0L,",
      "    degree.max = 1L,",
      "    degree.verify = FALSE,",
      "    bwtype = 'fixed',",
      "    bwmethod = 'cv.ls',",
      "    nmulti = 2L,",
      "    max.bb.eval = 8L",
      "  ),",
      "  npscoefbw = function() npscoefbw(",
      "    xdat = x,",
      "    zdat = z,",
      "    ydat = y_sc,",
      "    regtype = 'lp',",
      "    bernstein.basis = TRUE,",
      "    degree.select = 'coordinate',",
      "    search.engine = 'nomad+powell',",
      "    degree.min = 0L,",
      "    degree.max = 1L,",
      "    degree.verify = FALSE,",
      "    bwtype = 'fixed',",
      "    bwmethod = 'cv.ls',",
      "    nmulti = 2L,",
      "    max.bb.eval = 8L",
      "  ),",
      "  npindexbw = function() npindexbw(",
      "    xdat = data.frame(x1 = x$x, x2 = z$z),",
      "    ydat = y_idx,",
      "    bws = c(1, 0.5, 0.35),",
      "    method = 'ichimura',",
      "    regtype = 'lp',",
      "    bernstein.basis = TRUE,",
      "    degree.select = 'coordinate',",
      "    search.engine = 'nomad+powell',",
      "    degree.min = 0L,",
      "    degree.max = 1L,",
      "    degree.verify = FALSE,",
      "    bwtype = 'fixed',",
      "    nmulti = 2L,",
      "    max.bb.eval = 8L",
      "  )",
      ")",
      "for (nm in names(cases)) {",
      "  cat(sprintf('CASE_START %s\\n', nm))",
      "  force(cases[[nm]]())",
      "  cat(sprintf('CASE_DONE %s\\n', nm))",
      "}"
    ),
    timeout = 240L,
    env = env
  )

  expect_identical(res$status, 0L, info = paste(res$output, collapse = "\n"))

  for (case in c("npcdistbw", "npplregbw", "npscoefbw", "npindexbw")) {
    expect_nomad_unknown_bound_output(res$output, case)
  }
})
