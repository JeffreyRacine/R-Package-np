# Test-only interception for mathematical contracts about provider state that is
# deliberately absent from the public plot object (draws, covariance and
# alternative calibrations). The real public normalizer, computation and new
# builder all execute. This is not a compatibility API or a substituted fit.
capture_plot_provider <- function(expr) {
  provider <- NULL
  builder <- getFromNamespace(".np_plot_data_from_engine", "npRmpi")
  testthat::local_mocked_bindings(.np_plot_data_from_engine = function(raw, ...) {
    provider <<- raw
    builder(raw, ...)
  }, .package = "npRmpi")
  result <- force(expr)
  testthat::expect_s3_class(result, "np_plot_data")
  testthat::expect_false(inherits(result, c("npdensity", "npdistribution", "npregression")))
  testthat::expect_false(is.null(provider))
  list(result = result, provider = provider)
}

plot_target <- function(object, panel = 1L, target = 1L, band = NULL) {
  as.data.frame(object, panel = panel, target = target, band = band)
}

plot_band <- function(object, band, panel = 1L, target = 1L) {
  d <- as.data.frame(object, panel = panel, target = target, band = band)
  unname(as.matrix(d[c("lower", "upper")]))
}

plot_coordinates <- function(object, panel = 1L, role) {
  p <- object$panels[[panel]]
  v <- p$variables[p$variables$role == role, , drop = FALSE]
  d <- p$grid[v$coordinate]
  names(d) <- v$name
  d
}

# Generated subprocess contracts use the same extraction helpers as in-process
# contracts. This emits only pure table selectors, never a fitted-object shim.
plot_test_support <- function() {
  vapply(c("plot_target", "plot_band", "plot_coordinates"), function(n)
    paste(n, "<-", paste(deparse(get(n, mode = "function")), collapse = "\n")), character(1L))
}

# Geometry/calibration tests may inspect provider internals without publishing
# them. The real public builder still executes; public table tests stay separate.
capture_copula_provider <- function(expr) {
  provider <- NULL
  builder <- getFromNamespace(".np_plot_data_from_copula", "npRmpi")
  testthat::local_mocked_bindings(.np_plot_data_from_copula = function(raw, ...) {
    provider <<- raw
    builder(raw, ...)
  }, .package = "npRmpi")
  result <- force(expr)
  testthat::expect_s3_class(result, "np_plot_data")
  testthat::expect_false(is.null(provider))
  list(result = result, provider = provider)
}
