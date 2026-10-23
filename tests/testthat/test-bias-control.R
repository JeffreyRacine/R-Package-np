test_that("bias controls preserve explicit construction and auxiliary basis", {
  control <- np_bias_control()
  expect_s3_class(control, "np_bias_control")
  expect_identical(control$iterations, 16L)
  expect_true(control$bernstein.basis)
  expect_identical(np_boot_control()$wild, "mammen")
  expect_identical(np_boot_control(wild="rademacher")$wild, "rademacher")
  expect_error(np_bias_control(degree=3, bernstein.basis=FALSE), "Bernstein|bernstein")
  expect_error(np_bias_control(iterations=1.2), "positive integer")
  expect_error(np_bias_control(p0.scale=4, p0.bandwidth=.2), "mutually exclusive")
  expect_error(np_bias_control(construction="legacy", degree=5), "do not apply")
  control$construction <- "typo"
  expect_error(np:::.np_bias_validate_control(control), "invalid np_bias_control")
})

test_that("corrected plot families nest and singleton limits coincide", {
  old <- options(np.messages=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(27051)
  x <- data.frame(x=runif(80))
  y <- sin(3*x$x) + rnorm(80, sd=.2)
  bw <- npregbw(xdat=x, ydat=y, bws=.23, bandwidth.compute=FALSE,
                regtype="lp", degree=2, bernstein.basis=TRUE)
  fit <- npreg(bw, txdat=x, tydat=y)
  original <- fit
  rng <- .Random.seed
  out <- plot(fit, errors="asymptotic", center="bias-corrected", band="all",
              gaussian.draws=63L, neval=3L, output="data")[[1L]]
  expect_identical(fit, original)
  expect_identical(.Random.seed, rng)
  expect_true(out$inference$settings$bernstein.basis)
  expect_identical(out$inference$settings$main$degree, bw$degree)
  b <- out$inference$bands
  expect_true(all(b$bonferroni[,1,] <= b$simultaneous[,1,]))
  expect_true(all(b$simultaneous[,1,] <= b$pointwise[,1,]))
  expect_true(all(b$pointwise[,2,] <= b$simultaneous[,2,]))
  expect_true(all(b$simultaneous[,2,] <= b$bonferroni[,2,]))
  state <- np:::.np_bias_regression_prepare(bw, x, y)
  for (law in c("gaussian", "mammen", "rademacher", "inid", "fixed", "geom")) {
    z <- np:::.np_with_seed(42L, np:::.np_bias_inference(
      state, x[1L,,drop=FALSE], "mean", .05, 63L, law))
    expect_identical(z$bands$pointwise, z$bands$simultaneous)
    expect_identical(z$bands$pointwise, z$bands$bonferroni)
  }
})

test_that("iterated intervals describe maximum calibration without relabeling legacy bands", {
  label <- np:::.np_plot_all_band_legend_labels
  descriptor <- np:::.np_plot_variability_descriptor
  expect_match(label("bootstrap")[["simultaneous"]], "rank-based", fixed=TRUE)
  expect_match(label("bootstrap", "iterated")[["simultaneous"]], "maximum", fixed=TRUE)
  expect_match(label("bootstrap", "iterated")[["bonferroni"]], "envelope", fixed=TRUE)
  text <- descriptor("bootstrap", "simultaneous", .05, "bias-corrected",
                     construction="iterated")$candidates
  expect_true(all(grepl("maximum", text, fixed=TRUE)))
  expect_false(any(grepl("rank", text, fixed=TRUE)))
  text <- descriptor("asymptotic", "pointwise", .05, "bias-corrected",
                     construction="iterated")$candidates
  expect_true(all(grepl("bias-corrected", text, fixed=TRUE)))
})


test_that("auxiliary degrees resolve by coordinate without altering the main model", {
  resolve <- np:::.np_bias_degree_spec
  spec <- resolve(np_bias_control(), c("a", "b", "c"), c(1, 2, 6), 200)
  expect_identical(spec$effective, c(a=5L, b=5L, c=6L))
  expect_identical(spec$main, c(a=1L, b=2L, c=6L))
  spec <- resolve(np_bias_control(degree=c(b=4, a=2)), c("a", "b"), c(1, 3), 100)
  expect_identical(spec$effective, c(a=2L, b=4L))
  expect_error(resolve(np_bias_control(degree=2), c("a", "b"), c(1, 3), 100), "cover")
  expect_error(resolve(np_bias_control(), paste0("x",1:10), rep(1,10), 100), "residual degrees")
  expect_error(np_bias_control(degree=c(a=2,a=3)), "unique names")
  expect_error(np_bias_control(degree=0), "positive integers")
  old <- np_bias_control(degree=5L)
  old$version <- 1L
  expect_identical(np:::.np_bias_validate_control(old), old)
})

test_that("target blocks cover both axes without changing target-major indices", {
  for (ng in c(1L,26L)) for (nt in c(20L,21L,501L)) {
    blocks <- np:::.np_bias_target_blocks(ng,nt)
    rows <- unlist(lapply(blocks, `[[`, "rows"))
    expect_identical(sort(as.integer(rows)), seq_len(ng*nt))
    expect_true(all(vapply(blocks,function(b) length(b$grid)<=25L && length(b$rows)<=500L,logical(1))))
  }
})
