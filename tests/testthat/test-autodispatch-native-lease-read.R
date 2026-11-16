test_that("untagging removes real dispatch tags and preserves unrelated structure", {
  tags <- c("npRmpi.dispatch.mode", "npRmpi.autodispatch.remote",
            "npRmpi.autodispatch.fingerprint", "npRmpi.autodispatch.lease")
  d <- data.frame(x = 1:8, g = factor(rep(c("a", "b"), 4)))
  attr(d$x, "label") <- "keep"
  clean <- list(training = d, empty = NULL, call = quote(f(, x, NULL)), formula = y ~ x)
  tagged <- clean
  for (name in tags) {
    attr(tagged, name) <- "remove"
    attr(tagged$training$x, name) <- "remove"
  }
  expect_identical(.npRmpi_autodispatch_untag(tagged), clean)
  expect_identical(attr(tagged$training$x, "label", exact = TRUE), "keep")
  for (name in tags) expect_identical(attr(tagged, name, exact = TRUE), "remove")
})


test_that("reading unchanged native inputs preserves their published MPI lease", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  ns <- asNamespace("npRmpi")
  resolve <- get(".npRmpi_lease_resolve", ns)
  for (family in c("npreg", "npscoef")) {
    roles <- if (family == "npreg") c("x", "y") else c("x", "y", "z")
    for (role in roles) {
      i <- seq_len(48L)
      d <- data.frame(x = sin(i * .7), y = cos(i * .3), z = cos(i * .4))
      d[[role]] <- sort(d[[role]])
      args <- list(xdat = d["x"], ydat = d$y, bws = .5,
                   bandwidth.compute = FALSE)
      if (family == "npscoef") args$zdat <- d["z"]
      b <- do.call(get(paste0(family, "bw"), ns), args)
      lease <- resolve(b, family, "bws")
      expect_false(is.null(lease))
      invisible(complete.cases(d[[role]]))
      expect_identical(resolve(b, family, "bws"), lease)
      fit <- get(family, ns)(b, se = FALSE)
      expect_true(all(is.finite(fitted(fit))))
      expect_identical(resolve(b, family, "bws"), lease)
      b[[".np.native.training"]][["ydat"]][1L] <-
        b[[".np.native.training"]][["ydat"]][1L] + .01
      expect_null(resolve(b, family, "bws"))
    }
  }
})
