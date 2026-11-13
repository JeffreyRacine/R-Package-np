test_that("unconditional builders distinguish retained and fresh uniform metadata", {
  withr::local_options(np.messages = FALSE, np.tree = FALSE)
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  ns <- asNamespace("npRmpi")
  collect <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
  }
  builder <- get(".np_make_kbandwidth_unconditional", ns)
  old <- builder
  env <- new.env(parent = ns)
  env$.np_retained_uniform_constructor <- function(constructor, ...) constructor(...)
  environment(old) <- env
  x <- data.frame(x = seq(-2, 2, length.out = 31))
  for(kernel in c("uniform", "gaussian")) for(order in c(2L, 4L)) {
    b <- collect(npudensbw(dat = x, bws = .7, bandwidth.compute = FALSE,
      ckertype = kernel, ckerorder = order))$value
    got <- collect(builder(b, x)); before <- collect(old(b, x))
    expect_identical(got$value, before$value)
    expect_identical(got$value$ckertype, b$ckertype)
    expect_identical(got$value$ckerorder, b$ckerorder)
    expect_length(got$warnings, 0L)
    expect_length(before$warnings, as.integer(kernel == "uniform"))
  }
  fresh <- collect(get("kbandwidth.numeric", ns)(.7, ckertype = "uniform",
    ckerorder = 4L, xdati = get("untangle", ns)(x), xnames = "x"))
  expect_match(fresh$warnings, "ignoring kernel order specified", fixed = TRUE)
  b$ckertype <- "uniform"
  probe <- builder; e <- new.env(parent = ns); environment(probe) <- e
  e$kbandwidth.numeric <- function(...) {warning("unrelated sentinel"); 7L}
  expect_identical(collect(probe(b, x))$warnings, "unrelated sentinel")
  e$kbandwidth.numeric <- function(...) stop("constructor sentinel", call. = FALSE)
  expect_error(probe(b, x), "constructor sentinel", fixed = TRUE)
})
