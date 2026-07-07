.np_test_plot_proto <- function(name, package = "np") {
  raw <- getFromNamespace(name, package)
  force(raw)

  function(..., errors, bootstrap, boot_control, B, center, band) {
    args <- list(...)

    if (!missing(errors)) {
      args$plot.errors.method <- errors
    }
    if (!missing(bootstrap)) {
      args$plot.errors.boot.method <- bootstrap
    }
    if (!missing(boot_control)) {
      args$plot.errors.boot.blocklen <- boot_control[["blocklen"]]
    }
    if (!missing(B)) {
      args$plot.errors.boot.num <- B
    }
    if (!missing(center)) {
      args$plot.errors.center <- center
    }
    if (!missing(band)) {
      args$plot.errors.type <- band
    }

    do.call(raw, args)
  }
}

.np_expect_plot_proto_names <- function(candidate, reference, info = NULL) {
  missing <- setdiff(names(reference), names(candidate))
  extra <- setdiff(names(candidate), names(reference))

  expect_identical(extra, character(), info = info)
  expect_true(all(missing %in% "bias.corrected"), info = info)
}

.np_expect_plot_proto_bias <- function(candidate, reference, info = NULL) {
  candidate <- as.numeric(candidate)
  reference <- as.numeric(reference)

  expect_equal(length(candidate), length(reference), info = info)
  if (length(candidate) != length(reference)) {
    return(invisible(NULL))
  }

  same <- isTRUE(all.equal(candidate, reference, check.attributes = FALSE))
  opposite <- isTRUE(all.equal(candidate, -reference, check.attributes = FALSE))
  expect_true(same || opposite, info = info)
}

.np_expect_plot_proto_error_shape <- function(candidate, reference, info = NULL) {
  if (is.null(reference)) {
    return(invisible(NULL))
  }

  expect_equal(length(candidate), length(reference), info = info)
  expect_equal(dim(candidate), dim(reference), info = info)
  expect_true(all(is.finite(candidate) | is.na(candidate)), info = info)
}
