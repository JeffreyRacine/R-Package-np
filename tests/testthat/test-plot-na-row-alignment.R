test_that("single-band render call sites preserve row identity with missing limits", {
  # Execute the actual argument expressions, then compare each retained point to
  # its finite-data counterpart at the same coordinate. This distinguishes a
  # common mask from independently dropping NAs without relying on source text.
  engines <- list(.np_plot_rbandwidth_engine, .np_plot_plbandwidth_engine,
                  .np_plot_scbandwidth_engine, .np_plot_dbandwidth_engine,
                  .np_plot_sibandwidth_engine)
  sites <- list()
  walk <- function(node) {
    if (!is.call(node) && !is.expression(node)) return(invisible(NULL))
    parts <- as.list(node)
    if (all(c("ex", "ely", "ehy") %in% names(parts)))
      sites[[length(sites) + 1L]] <<- parts[c("ex", "ely", "ehy")]
    for (i in seq_along(parts)) {
      if (identical(parts[[i]], quote(expr = ))) next
      walk(parts[[i]])
    }
  }
  for (engine in engines) walk(body(engine))
  expect_length(sites, 14L)
  fixture <- function(missing, corrected, negative) {
    e <- new.env(parent = environment())
    e$ei <- c(2, 5, 9, 14, 20, 27)
    e$temp.mean <- e$temp.dens <- e$ei^2
    e$temp.err <- cbind(seq_len(6)/10, seq_len(6)/5, e$ei^2 + 3)
    e$temp.err[missing, 1] <- NA_real_
    if (length(missing)) e$temp.err[3, 2] <- NA_real_
    e$allei <- matrix(e$ei, ncol = 1L)
    e$data.eval <- matrix(e$temp.mean, ncol = 1L)
    e$data.err <- e$temp.err
    e$tobj <- list(index = e$ei)
    e$i <- e$plot.index <- 1L
    e$i.sort <- 6:1
    e$plotOnEstimate <- !corrected
    ctr <- if (corrected) e$temp.err[, 3] else e$temp.mean
    beta <- if (negative) -2 else 2
    e$lo.i <- beta * (ctr[e$i.sort] - e$temp.err[e$i.sort, 1])
    e$hi.i <- beta * (ctr[e$i.sort] + e$temp.err[e$i.sort, 2])
    e
  }
  for (corrected in c(FALSE, TRUE)) for (negative in c(FALSE, TRUE)) {
    for (site in sites) {
      full <- lapply(site, eval, envir = fixture(integer(), corrected, negative))
      for (missing in list(2L, c(1L, 4L), 6L)) {
        partial <- lapply(site, eval, envir = fixture(missing, corrected, negative))
        expect_equal(lengths(partial), rep(6L, 3L), ignore_attr = TRUE)
        keep <- Reduce(`&`, lapply(partial, is.finite))
        expect_true(any(keep) && any(!keep))
        pos <- match(partial$ex[keep], full$ex)
        expect_equal(partial$ely[keep], full$ely[pos])
        expect_equal(partial$ehy[keep], full$ehy[pos])
      }
    }
  }
})
