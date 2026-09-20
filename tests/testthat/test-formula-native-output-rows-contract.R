test_that("formula restoration belongs to the actual evaluation rows", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(920344)
  d <- data.frame(x = runif(40), z = rnorm(40), y = rnorm(40))
  d$x[3] <- NA_real_
  ee <- data.frame(x = c(.1, .3, .6, .8), y = c(-1, -.3, .3, 1),
                   z = c(-1, 1, .2, -.2))
  for (family in c("npreg", "npudens", "npudist", "npcdens", "npcdist",
                   "npscoef", "npindex", "npplreg")) {
    form <- switch(family, npudens = ~x, npudist = ~x,
                   npscoef = y~x|z, npplreg = y~x|z, npindex = y~x+z, y~x)
    widths <- switch(family, npcdens = c(.4,.3), npcdist = c(.4,.3),
                     npindex = c(1,.5,.3), npplreg = matrix(.3,2,1), .3)
    bw <- do.call(get(paste0(family, "bw")), list(form, data = d,
       bws = widths, bandwidth.compute = FALSE, na.action = na.exclude))
    args <- list(bws = bw, se = family != "npindex")
    if (family == "npscoef") args <- c(args, list(iterate = FALSE, betas = TRUE))
    if (family %in% c("npreg", "npcdens", "npcdist")) args$gradients <- TRUE
    ff <- get(family)
    training <- do.call(ff, args)
    expect_equal(NROW(fitted(training)), nrow(d), info = family)
    expect_true(is.na(fitted(training)[3]), info = family)
    for (missing.eval in c(FALSE, TRUE)) {
      e <- ee
      if (missing.eval) e$x[2] <- NA_real_
      native <- switch(family, npudens = list(edat = e["x"]),
        npudist = list(edat = e["x"]),
        npcdens = list(exdat = e["x"], eydat = e["y"]),
        npcdist = list(exdat = e["x"], eydat = e["y"]),
        npscoef = list(exdat = e["x"], ezdat = e["z"]),
        npplreg = list(exdat = e["x"], ezdat = e["z"]),
        npindex = list(exdat = e[c("x", "z")]), list(exdat = e["x"]))
      a <- do.call(ff, c(args, native))
      b <- do.call(ff, c(args, list(newdata = e)))
      both <- do.call(ff, c(args, native, list(newdata = ee[1:2, ])))
      # Partially linear native evaluation already pads its own omitted rows.
      compared <- if (family == "npplreg" && missing.eval)
        fitted(a)[-2L] else fitted(a)
      expect_equal(compared, fitted(b), tolerance = 1e-12, info = family)
      expect_equal(fitted(a), fitted(both), tolerance = 0, info = family)
      expect_equal(NROW(fitted(a)),
        if (family == "npplreg") 4L else 4L - missing.eval, info = family)
      if (!family %in% c("npindex", "npplreg")) {
        sa <- se(a)
        sb <- se(b)
        expect_equal(sa, sb, tolerance = 1e-12, info = family)
        expect_equal(NROW(sa), 4L - missing.eval, info = family)
      }
      if (isTRUE(args$gradients)) {
        expect_equal(gradients(a), gradients(b), tolerance = 1e-12, info = family)
        expect_equal(NROW(gradients(a)), 4L - missing.eval, info = family)
      }
      if (family == "npscoef")
        expect_equal(a$beta, b$beta, tolerance = 1e-12)
    }
  }
})
