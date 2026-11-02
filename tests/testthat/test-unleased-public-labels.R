test_that("unleased native labels finalize complete valid public metadata", {
  finalize <- getFromNamespace(".npRmpi_autodispatch_finalize_unleased_names", "npRmpi")
  signature <- getFromNamespace(".np_bw_procedure_signature", "npRmpi")
  old <- ".__npRmpi_autod_dat_1"
  x <- c(-1, 0, 1)
  make <- function(cls) {
    b <- structure(list(bw = .4, xnames = old, varnames = list(x = old),
      icon = setNames(TRUE, old), iuno = setNames(FALSE, old),
      iord = setNames(FALSE, old), type = "fixed",
      xdati = list(all.lev = setNames(list(old), old)),
      .np.native.training = list(dat = setNames(data.frame(x), old))), class = cls)
    b$selection.procedure <- list(signature = signature(b),
      arguments = list(unrelated.literal = old))
    b
  }
  run <- function(dat, result, external = NULL) {
    finalize(result, list(dat = dat, bws = external), run,
      owner = environment(), caller = parent.frame())
  }
  for (cls in c("bandwidth", "dbandwidth")) {
    b <- make(cls); before <- serialize(b, NULL)
    out <- run(x, b)
    expect_identical(b, unserialize(before))
    expect_identical(out$xnames, "x")
    expect_identical(out$varnames$x, "x")
    expect_identical(names(out$.np.native.training$dat), "x")
    expect_identical(out$.np.native.training$dat[[1L]], x)
    expect_identical(names(out$xdati$all.lev), "x")
    expect_identical(out$xdati$all.lev[[1L]], old)
    expect_identical(out$selection.procedure$arguments,
                     b$selection.procedure$arguments)
    expect_identical(out$selection.procedure$signature, signature(out))
    expect_identical(run(x, b, external = b), b)
    stale <- b; stale$selection.procedure$signature$bw <- .5
    expect_identical(run(x, stale), stale)
    supplied <- setNames(data.frame(x), old)
    expect_identical(run(supplied, b), b)
    fitted <- structure(list(bws = b, xnames = old, trainiseval = TRUE,
      eval = supplied), class = if (cls == "bandwidth") "npdensity" else "npdistribution")
    out <- run(x, fitted)
    expect_identical(names(out$eval), "x")
    expect_identical(out$xnames, "x")
    expect_identical(out$bws$selection.procedure$signature, signature(out$bws))
    fitted$trainiseval <- FALSE
    expect_identical(run(x, fitted)$eval, supplied)
  }
  untouched <- structure(list(xnames = old), class = "rbandwidth")
  expect_identical(run(x, untouched), untouched)
})
