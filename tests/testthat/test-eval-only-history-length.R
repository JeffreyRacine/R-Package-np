test_that("evaluation-only histories describe one actual evaluation", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(12L)
  dat <- data.frame(x = i / 12, z = sin(i))
  for (family in c("npudens", "npudist")) {
    bw <- do.call(get(paste0(family, "bw")), list(
      dat = dat, bws = c(0.5, 0.5), bwscaling = FALSE,
      bandwidth.compute = FALSE))
    method <- getFromNamespace(paste0(family, "bw.",
      if (family == "npudens") "bandwidth" else "dbandwidth"), "np")
    raw <- method(dat = dat, bws = bw, eval.only = TRUE, nmulti = 3L)
    fields <- raw[c("fval.history", "eval.history", "invalid.history")]
    expect_identical(lengths(fields), setNames(rep(1L, 3L), names(fields)))
    expect_true(all(is.finite(unlist(fields))))
    # A baseline penalty may also evaluate its anchor; it is still one cell.
    expect_gte(as.double(raw$eval.history), 1)
    expect_identical(as.double(raw$invalid.history), 0)
  }
})
