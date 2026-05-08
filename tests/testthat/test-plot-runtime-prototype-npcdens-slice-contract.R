test_that("npcdens LC fixed no-error plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(123)

  n <- 80L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lc",
    bwtype = "fixed"
  )

  old <- plot(
    bw,
    xdat = x,
    ydat = y,
    output = "data",
    errors = "none",
    view = "fixed",
    neval = 9L,
    perspective = TRUE
  )
  candidate <- proto(bw, xdat = x, ydat = y, neval = 9L)
  stages <- proto(bw, xdat = x, ydat = y, neval = 9L, return.stages = TRUE)

  expect_named(candidate, names(old))
  expect_s3_class(candidate$cd1, "condensity")
  expect_named(candidate$cd1, names(old$cd1))
  expect_equal(candidate$cd1$xeval, old$cd1$xeval)
  expect_equal(candidate$cd1$yeval, old$cd1$yeval)
  expect_equal(candidate$cd1$condens, old$cd1$condens)
  expect_equal(candidate$cd1$conderr, old$cd1$conderr)
  expect_equal(candidate$cd1$bias, old$cd1$bias)
  expect_identical(candidate$cd1$proper.requested, old$cd1$proper.requested)
  expect_identical(candidate$cd1$proper.applied, old$cd1$proper.applied)
  expect_named(stages, c("state", "target_grid", "evaluator", "intervals", "bootstrap", "proper_projection", "plot_data"))
  expect_null(stages$intervals)
  expect_null(stages$bootstrap)
  expect_null(stages$proper_projection)
  expect_equal(stages$plot_data, candidate)
  expect_equal(stages$evaluator$condens, old$cd1$condens)
})

test_that("npcdens plot prototype fails early outside its vertical slice", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(124)

  n <- 60L
  x <- data.frame(x = factor(rbinom(n, 1L, 0.5)))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lc",
    bwtype = "fixed"
  )

  expect_error(
    proto(bw, xdat = x, ydat = y, neval = 9L),
    "continuous/ordered surface variables",
    fixed = TRUE
  )

  expect_error(
    proto(bw, neval = 9L),
    "explicit xdat and ydat",
    fixed = TRUE
  )
})

test_that("npcdens LC fixed asymptotic plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_asymptotic_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(125)

  n <- 80L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lc",
    bwtype = "fixed"
  )

  for (band in c("pmzsd", "pointwise", "bonferroni", "simultaneous", "all")) {
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "asymptotic",
      band = band,
      view = "fixed",
      neval = 9L,
      perspective = TRUE
    ))
    candidate <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 9L,
      band = band
    )
    stages <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 9L,
      band = band,
      return.stages = TRUE
    )

    expect_equal(candidate$cd1$xeval, old$cd1$xeval, info = band)
    expect_equal(candidate$cd1$yeval, old$cd1$yeval, info = band)
    expect_equal(candidate$cd1$condens, old$cd1$condens, info = band)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = band)
    expect_named(stages, c("state", "target_grid", "evaluator", "intervals", "bootstrap", "proper_projection", "plot_data"))
    expect_equal(stages$plot_data, candidate, info = band)
    expect_identical(stages$intervals$method, "asymptotic")
    expect_identical(stages$intervals$type, band)
    expect_null(stages$bootstrap)
    expect_null(stages$proper_projection)
  }
})

test_that("npcdens LC fixed inid bootstrap plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_bootstrap_inid_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(126)

  n <- 70L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lc",
    bwtype = "fixed"
  )

  cases <- expand.grid(
    band = c("pmzsd", "pointwise", "bonferroni"),
    center = c("estimate", "bias-corrected"),
    stringsAsFactors = FALSE
  )
  for (ii in seq_len(nrow(cases))) {
    band <- cases$band[ii]
    center <- cases$center[ii]
    boot.seed <- 1000 + ii
    set.seed(boot.seed)
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "bootstrap",
      bootstrap = "inid",
      B = 19L,
      center = center,
      band = band,
      view = "fixed",
      neval = 7L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 7L,
      B = 19L,
      center = center,
      band = band
    ))
    set.seed(boot.seed)
    stages <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 7L,
      B = 19L,
      center = center,
      band = band,
      return.stages = TRUE
    ))

    expect_equal(candidate$cd1$xeval, old$cd1$xeval, info = paste(band, center))
    expect_equal(candidate$cd1$yeval, old$cd1$yeval, info = paste(band, center))
    expect_equal(candidate$cd1$condens, old$cd1$condens, info = paste(band, center))
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = paste(band, center))
    expect_named(stages, c("state", "target_grid", "evaluator", "intervals", "bootstrap", "proper_projection", "plot_data"))
    expect_equal(stages$plot_data, candidate, info = paste(band, center))
    expect_identical(stages$intervals$method, "bootstrap")
    expect_identical(stages$intervals$type, band)
    expect_identical(stages$bootstrap$method, "inid")
    expect_identical(stages$bootstrap$B, 19L)
    expect_equal(stages$bootstrap$boot.err[, 1:2, drop = FALSE],
                 candidate$cd1$conderr,
                 info = paste(band, center))
    expect_null(stages$proper_projection)
  }
})

test_that("npcdens LC fixed proper plot-data prototype matches current route", {
  none_proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_none_data", "np")
  boot_proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_bootstrap_inid_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(127)

  n <- 75L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lc",
    bwtype = "fixed"
  )

  old_none <- plot(
    bw,
    xdat = x,
    ydat = y,
    output = "data",
    errors = "none",
    proper = TRUE,
    view = "fixed",
    neval = 7L,
    perspective = TRUE
  )
  cand_none <- none_proto(
    bw,
    xdat = x,
    ydat = y,
    neval = 7L,
    proper = TRUE
  )
  stages_none <- none_proto(
    bw,
    xdat = x,
    ydat = y,
    neval = 7L,
    proper = TRUE,
    return.stages = TRUE
  )

  expect_equal(cand_none$cd1$condens, old_none$cd1$condens)
  expect_equal(cand_none$cd1$conderr, old_none$cd1$conderr)
  expect_identical(cand_none$cd1$proper.requested, old_none$cd1$proper.requested)
  expect_identical(cand_none$cd1$proper.applied, old_none$cd1$proper.applied)
  expect_identical(stages_none$proper_projection$requested, TRUE)
  expect_identical(stages_none$proper_projection$applied, cand_none$cd1$proper.applied)

  boot.seed <- 9001L
  set.seed(boot.seed)
  old_boot <- suppressWarnings(plot(
    bw,
    xdat = x,
    ydat = y,
    output = "data",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 19L,
    band = "pmzsd",
    proper = TRUE,
    view = "fixed",
    neval = 7L,
    perspective = TRUE,
    random.seed = boot.seed
  ))
  set.seed(boot.seed)
  cand_boot <- suppressWarnings(boot_proto(
    bw,
    xdat = x,
    ydat = y,
    neval = 7L,
    B = 19L,
    band = "pmzsd",
    proper = TRUE
  ))
  set.seed(boot.seed)
  stages_boot <- suppressWarnings(boot_proto(
    bw,
    xdat = x,
    ydat = y,
    neval = 7L,
    B = 19L,
    band = "pmzsd",
    proper = TRUE,
    return.stages = TRUE
  ))

  expect_equal(cand_boot$cd1$condens, old_boot$cd1$condens)
  expect_equal(cand_boot$cd1$conderr, old_boot$cd1$conderr)
  expect_identical(cand_boot$cd1$proper.requested, old_boot$cd1$proper.requested)
  expect_identical(cand_boot$cd1$proper.applied, old_boot$cd1$proper.applied)
  expect_identical(stages_boot$proper_projection$requested, TRUE)
  expect_identical(stages_boot$proper_projection$applied, cand_boot$cd1$proper.applied)
  expect_equal(stages_boot$bootstrap$boot.err[, 1:2, drop = FALSE],
               cand_boot$cd1$conderr)
})

test_that("npcdens LC fixed block bootstrap plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_bootstrap_block_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(128)

  n <- 70L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lc",
    bwtype = "fixed"
  )

  for (method in c("fixed", "geom")) {
    boot.seed <- if (identical(method, "fixed")) 9101L else 9102L
    set.seed(boot.seed)
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "bootstrap",
      bootstrap = method,
      boot_control = np_boot_control(blocklen = 5L),
      B = 19L,
      band = "pmzsd",
      view = "fixed",
      neval = 7L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 7L,
      bootstrap = method,
      boot_control = np_boot_control(blocklen = 5L),
      B = 19L,
      band = "pmzsd"
    ))
    set.seed(boot.seed)
    stages <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 7L,
      bootstrap = method,
      boot_control = np_boot_control(blocklen = 5L),
      B = 19L,
      band = "pmzsd",
      return.stages = TRUE
    ))

    expect_equal(candidate$cd1$condens, old$cd1$condens, info = method)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = method)
    expect_identical(stages$bootstrap$method, method)
    expect_identical(stages$bootstrap$blocklen, 5L)
    expect_equal(stages$bootstrap$boot.err[, 1:2, drop = FALSE],
                 candidate$cd1$conderr,
                 info = method)
  }
})

test_that("npcdens fixed LL/LP no-error plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(129)

  n <- 75L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  cases <- list(
    ll = list(regtype = "ll"),
    lp1 = list(regtype = "lp", degree = 1L),
    lp2 = list(regtype = "lp", degree = 2L)
  )

  for (nm in names(cases)) {
    bw <- do.call(
      npcdensbw,
      c(list(xdat = x, ydat = y, nmulti = 1L, bwtype = "fixed"), cases[[nm]])
    )
    old <- plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "none",
      view = "fixed",
      neval = 7L,
      perspective = TRUE
    )
    candidate <- proto(bw, xdat = x, ydat = y, neval = 7L)
    stages <- proto(bw, xdat = x, ydat = y, neval = 7L, return.stages = TRUE)

    expect_equal(candidate$cd1$xeval, old$cd1$xeval, info = nm)
    expect_equal(candidate$cd1$yeval, old$cd1$yeval, info = nm)
    expect_equal(candidate$cd1$condens, old$cd1$condens, info = nm)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = nm)
    expect_equal(stages$plot_data, candidate, info = nm)
  }
})

test_that("npcdens fixed LL/LP inid bootstrap plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_fixed_bootstrap_inid_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(130)

  n <- 70L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  cases <- list(
    ll = list(regtype = "ll"),
    lp1 = list(regtype = "lp", degree = 1L),
    lp2 = list(regtype = "lp", degree = 2L)
  )

  for (ii in seq_along(cases)) {
    nm <- names(cases)[ii]
    bw <- do.call(
      npcdensbw,
      c(list(xdat = x, ydat = y, nmulti = 1L, bwtype = "fixed"), cases[[ii]])
    )
    boot.seed <- 9200L + ii
    set.seed(boot.seed)
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "bootstrap",
      bootstrap = "inid",
      B = 17L,
      band = "pmzsd",
      view = "fixed",
      neval = 7L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 7L,
      B = 17L,
      band = "pmzsd"
    ))
    set.seed(boot.seed)
    stages <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 7L,
      B = 17L,
      band = "pmzsd",
      return.stages = TRUE
    ))

    expect_equal(candidate$cd1$condens, old$cd1$condens, info = nm)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = nm)
    expect_identical(stages$bootstrap$method, "inid")
    expect_equal(stages$bootstrap$boot.err[, 1:2, drop = FALSE],
                 candidate$cd1$conderr,
                 info = nm)
  }
})

test_that("npcdens generalized/adaptive NN no-error plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(131)

  n <- 70L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  cases <- c("generalized_nn", "adaptive_nn")

  for (bwtype in cases) {
    bw <- npcdensbw(
      xdat = x,
      ydat = y,
      nmulti = 1L,
      regtype = "lc",
      bwtype = bwtype
    )
    old <- plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "none",
      view = "fixed",
      neval = 7L,
      perspective = TRUE
    )
    candidate <- proto(bw, xdat = x, ydat = y, neval = 7L)
    stages <- proto(bw, xdat = x, ydat = y, neval = 7L, return.stages = TRUE)

    expect_equal(candidate$cd1$xeval, old$cd1$xeval, info = bwtype)
    expect_equal(candidate$cd1$yeval, old$cd1$yeval, info = bwtype)
    expect_equal(candidate$cd1$condens, old$cd1$condens, info = bwtype)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = bwtype)
    expect_equal(stages$plot_data, candidate, info = bwtype)
  }
})

test_that("npcdens generalized/adaptive NN inid bootstrap plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_fixed_bootstrap_inid_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(132)

  n <- 55L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  cases <- c("generalized_nn", "adaptive_nn")

  for (ii in seq_along(cases)) {
    bwtype <- cases[ii]
    bw <- npcdensbw(
      xdat = x,
      ydat = y,
      nmulti = 1L,
      regtype = "lc",
      bwtype = bwtype
    )
    boot.seed <- 9300L + ii
    set.seed(boot.seed)
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "bootstrap",
      bootstrap = "inid",
      B = 11L,
      band = "pmzsd",
      view = "fixed",
      neval = 5L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 5L,
      B = 11L,
      band = "pmzsd"
    ))
    set.seed(boot.seed)
    stages <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 5L,
      B = 11L,
      band = "pmzsd",
      return.stages = TRUE
    ))

    expect_equal(candidate$cd1$condens, old$cd1$condens, info = bwtype)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = bwtype)
    expect_identical(stages$bootstrap$method, "inid")
    expect_equal(stages$bootstrap$boot.err[, 1:2, drop = FALSE],
                 candidate$cd1$conderr,
                 info = bwtype)
  }
})

test_that("npcdist fixed no-error plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdist_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(134)

  n <- 75L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  cases <- list(
    lc = list(regtype = "lc"),
    ll = list(regtype = "ll"),
    lp1 = list(regtype = "lp", degree = 1L),
    lp2 = list(regtype = "lp", degree = 2L)
  )

  for (nm in names(cases)) {
    bw <- do.call(
      npcdistbw,
      c(list(xdat = x, ydat = y, nmulti = 1L, bwtype = "fixed"), cases[[nm]])
    )
    old <- plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "none",
      view = "fixed",
      neval = 7L,
      perspective = TRUE
    )
    candidate <- proto(bw, xdat = x, ydat = y, neval = 7L)
    stages <- proto(bw, xdat = x, ydat = y, neval = 7L, return.stages = TRUE)

    expect_s3_class(candidate$cd1, "condistribution")
    expect_equal(candidate$cd1$xeval, old$cd1$xeval, info = nm)
    expect_equal(candidate$cd1$yeval, old$cd1$yeval, info = nm)
    expect_equal(candidate$cd1$condist, old$cd1$condist, info = nm)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = nm)
    expect_identical(stages$state$family, "npcdist")
    expect_true(stages$state$cdf)
    expect_equal(stages$plot_data, candidate, info = nm)
  }
})

test_that("npcdist fixed inid bootstrap plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdist_fixed_bootstrap_inid_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(135)

  n <- 65L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  cases <- list(
    lc = list(regtype = "lc"),
    ll = list(regtype = "ll"),
    lp1 = list(regtype = "lp", degree = 1L),
    lp2 = list(regtype = "lp", degree = 2L)
  )

  for (ii in seq_along(cases)) {
    nm <- names(cases)[ii]
    bw <- do.call(
      npcdistbw,
      c(list(xdat = x, ydat = y, nmulti = 1L, bwtype = "fixed"), cases[[ii]])
    )
    boot.seed <- 9400L + ii
    set.seed(boot.seed)
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "bootstrap",
      bootstrap = "inid",
      B = 13L,
      band = "pmzsd",
      view = "fixed",
      neval = 6L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 6L,
      B = 13L,
      band = "pmzsd"
    ))
    set.seed(boot.seed)
    stages <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 6L,
      B = 13L,
      band = "pmzsd",
      return.stages = TRUE
    ))

    expect_equal(candidate$cd1$condist, old$cd1$condist, info = nm)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = nm)
    expect_identical(stages$bootstrap$method, "inid")
    expect_equal(stages$bootstrap$boot.err[, 1:2, drop = FALSE],
                 candidate$cd1$conderr,
                 info = nm)
  }
})

test_that("npcdist generalized/adaptive NN plot-data prototype matches current route", {
  none_proto <- getFromNamespace(".np_plot_proto_npcdist_fixed_none_data", "np")
  boot_proto <- getFromNamespace(".np_plot_proto_npcdist_fixed_bootstrap_inid_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(136)

  n <- 60L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))

  for (ii in seq_along(c("generalized_nn", "adaptive_nn"))) {
    bwtype <- c("generalized_nn", "adaptive_nn")[ii]
    bw <- npcdistbw(
      xdat = x,
      ydat = y,
      nmulti = 1L,
      regtype = "lc",
      bwtype = bwtype
    )
    old <- plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "none",
      view = "fixed",
      neval = 6L,
      perspective = TRUE
    )
    candidate <- none_proto(bw, xdat = x, ydat = y, neval = 6L)
    expect_equal(candidate$cd1$xeval, old$cd1$xeval, info = bwtype)
    expect_equal(candidate$cd1$yeval, old$cd1$yeval, info = bwtype)
    expect_equal(candidate$cd1$condist, old$cd1$condist, info = bwtype)
    expect_equal(candidate$cd1$conderr, old$cd1$conderr, info = bwtype)

    boot.seed <- 9500L + ii
    set.seed(boot.seed)
    old.boot <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "bootstrap",
      bootstrap = "inid",
      B = 11L,
      band = "pmzsd",
      view = "fixed",
      neval = 5L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate.boot <- suppressWarnings(boot_proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 5L,
      B = 11L,
      band = "pmzsd"
    ))
    expect_equal(candidate.boot$cd1$condist, old.boot$cd1$condist, info = bwtype)
    expect_equal(candidate.boot$cd1$conderr, old.boot$cd1$conderr, info = bwtype)
  }
})

test_that("npcdist proper projection plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdist_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(137)

  n <- 70L
  x <- data.frame(x = runif(n, -1, 1))
  y <- data.frame(y = sin(2 * pi * x$x) + rnorm(n, sd = 0.2))
  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L
  )
  old <- suppressWarnings(plot(
    bw,
    xdat = x,
    ydat = y,
    output = "data",
    errors = "none",
    view = "fixed",
    neval = 8L,
    perspective = TRUE,
    proper = TRUE
  ))
  candidate <- suppressWarnings(proto(
    bw,
    xdat = x,
    ydat = y,
    neval = 8L,
    proper = TRUE
  ))
  stages <- suppressWarnings(proto(
    bw,
    xdat = x,
    ydat = y,
    neval = 8L,
    proper = TRUE,
    return.stages = TRUE
  ))

  expect_true(isTRUE(candidate$cd1$proper.requested))
  expect_equal(candidate$cd1$proper.applied, old$cd1$proper.applied)
  expect_equal(candidate$cd1$condist, old$cd1$condist)
  expect_equal(candidate$cd1$conderr, old$cd1$conderr)
  expect_named(stages$proper_projection, c("requested", "applied", "method", "info"))
})

test_that("npcdens staged plot-data object renders through base renderer smoke", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_fixed_none_data", "np")
  render <- getFromNamespace(".np_plot_proto_npcdens_surface_base_render", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(133)

  n <- 65L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed"
  )
  candidate <- proto(bw, xdat = x, ydat = y, neval = 7L)

  out.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(out.file)
  on.exit({
    if (grDevices::dev.cur() > 1L)
      grDevices::dev.off()
  }, add = TRUE)
  expect_identical(render(candidate, perspective = FALSE), candidate)
  expect_identical(render(candidate, perspective = TRUE), candidate)
  grDevices::dev.off()

  expect_true(file.exists(out.file))
  expect_gt(file.info(out.file)$size, 0)
})

test_that("npcdist staged plot-data object renders through base renderer smoke", {
  proto <- getFromNamespace(".np_plot_proto_npcdist_fixed_none_data", "np")
  render <- getFromNamespace(".np_plot_proto_npcdens_surface_base_render", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(138)

  n <- 65L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed"
  )
  candidate <- proto(bw, xdat = x, ydat = y, neval = 7L)

  out.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(out.file)
  on.exit({
    if (grDevices::dev.cur() > 1L)
      grDevices::dev.off()
  }, add = TRUE)
  expect_identical(render(candidate, perspective = FALSE), candidate)
  expect_identical(render(candidate, perspective = TRUE), candidate)
  grDevices::dev.off()

  expect_true(file.exists(out.file))
  expect_gt(file.info(out.file)$size, 0)
})
