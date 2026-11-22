test_that("plot entry preserves only recognized existing lease ownership", {
  ns <- asNamespace("npRmpi")
  prepare <- get(".npRmpi_plot_bandwidth", environment(get(".np_plot_call_method", ns)))
  origin <- get(".npRmpi_lease_origin", ns)
  x <- list(values = seq_len(4L), child = list(a = 1))
  handle <- new.env(parent = emptyenv())
  for (name in c("id", "session", "capability", "fingerprint"))
    assign(name, name, handle)
  handle$origin <- origin
  attr(x, "npRmpi.autodispatch.lease") <- handle
  attr(x, "npRmpi.dispatch.mode") <- "auto"
  attr(x$child, "npRmpi.autodispatch.lease") <- handle
  y <- prepare(x)
  expect_identical(attr(y, "npRmpi.autodispatch.lease", exact = TRUE), handle)
  expect_null(attr(y, "npRmpi.dispatch.mode", exact = TRUE))
  expect_null(attr(y$child, "npRmpi.autodispatch.lease", exact = TRUE))
  expect_identical(y$values, x$values)
  expect_identical(attr(x, "npRmpi.dispatch.mode", exact = TRUE), "auto")
  expect_identical(attr(x$child, "npRmpi.autodispatch.lease", exact = TRUE), handle)
  bad <- x; attr(bad, "npRmpi.autodispatch.lease") <- new.env(parent = emptyenv())
  expect_null(attr(prepare(bad), "npRmpi.autodispatch.lease", exact = TRUE))
  bad <- x; attr(bad, "npRmpi.autodispatch.lease") <- "not a handle"
  expect_null(attr(prepare(bad), "npRmpi.autodispatch.lease", exact = TRUE))
  expect_null(attr(prepare(unserialize(serialize(x, NULL))), "npRmpi.autodispatch.lease", exact = TRUE))
  expect_identical(prepare(list(a = 1)), list(a = 1))
})

test_that('plot consumption uses the current lease and rejects changed ownership', {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), 'MPI runtime not available')
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace('npRmpi')
  old <- options(l1.test.plot.hits=0L, l1.test.plot.active=FALSE)
  on.exit(options(old), add=TRUE)
  trace('.npRmpi_lease_resolve', exit=quote({
    if (isTRUE(getOption('l1.test.plot.active')) && !is.null(returnValue()))
      options(l1.test.plot.hits=getOption('l1.test.plot.hits')+1L)
  }), where=ns, print=FALSE)
  on.exit(untrace('.npRmpi_lease_resolve',where=ns),add=TRUE)
  e <- new.env(parent=baseenv()); i <- seq_len(24L)
  e$d <- data.frame(y=cos(i*.3),x=sin(i*.7))
  f <- y~x; environment(f) <- e
  b <- npregbw(f,data=e$d,bws=.5,bandwidth.compute=FALSE)
  fp <- get('.npRmpi_autodispatch_fingerprint',ns)(b)
  options(l1.test.plot.active=TRUE)
  out <- plot(b,errors='asymptotic',output='data',neval=5)
  options(l1.test.plot.active=FALSE)
  expect_gt(getOption('l1.test.plot.hits'),0L)
  expect_s3_class(out,'np_plot_data')
  joint <- suppressWarnings(plot(b, errors = "asymptotic", band = "simultaneous",
                                 gaussian.draws = 199L, output = "data", neval = 5L))
  expect_s3_class(joint, "np_plot_data")
  expect_true(all(is.finite(joint$panels[[1L]]$intervals$lower)))
  expect_identical(get('.npRmpi_autodispatch_fingerprint',ns)(b),fp)
  bad <- b; bad$bw <- bad$bw*1.1
  options(l1.test.plot.active=TRUE,l1.test.plot.hits=0L)
  changed <- plot(bad,errors='asymptotic',output='data',neval=5)
  options(l1.test.plot.active=FALSE)
  expect_identical(getOption('l1.test.plot.hits'),0L)
  expect_s3_class(changed,'np_plot_data')
  restored <- unserialize(serialize(b,NULL))
  expect_null(get('.npRmpi_lease_resolve',ns)(restored,'npreg','bws'))
})

test_that("dispatch cleanup preserves byte identity of clean formula frames", {
  d <- data.frame(x = factor(c("b", "a")), y = 1:2, row.names = c("one", "two"))
  frame <- model.frame(y ~ x, data = d)
  nested <- list(clean = frame, tagged = structure(1:2, npRmpi.dispatch.mode = "auto"), empty = NULL)
  cleaned <- .npRmpi_autodispatch_untag(nested)
  expect_identical(serialize(cleaned$clean, NULL), serialize(frame, NULL))
  expect_identical(cleaned$tagged, 1:2)
  expect_identical(names(cleaned), names(nested))
  expect_identical(attr(nested$tagged, "npRmpi.dispatch.mode"), "auto")
})

test_that("a captured regression bandwidth validates its own formula provenance", {
  f <- y ~ x
  owner <- new.env(parent = baseenv())
  environment(f) <- owner
  original <- structure(list(bw = .4, formula = f, call = quote(npregbw(y ~ x))), class = "rbandwidth")
  remote <- .npRmpi_bw_metadata_wire(original)
  captured <- remote
  captured$bw <- .35
  result <- structure(list(bws = remote, fit.request = list(bws = captured)), class = "npregression")
  out <- .npRmpi_bw_restore_public_owner(result, list(original))
  expect_identical(environment(out$fit.request$bws$formula), owner)
  expect_identical(out$fit.request$bws$bw, .35)
  result$fit.request$bws$formula <- y ~ z
  expect_error(.npRmpi_bw_restore_public_owner(result, list(original)), "metadata changed")
})
