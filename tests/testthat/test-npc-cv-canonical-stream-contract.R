locate_npc_cv_source <- function(filename) {
  candidates <- c(
    testthat::test_path("..", "..", "src", filename),
    testthat::test_path("..", "..", "..", "src", filename),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", filename),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", filename),
    file.path(getwd(), "src", filename),
    file.path(getwd(), "..", "src", filename)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) {
    return(NULL)
  }
  hits[[1L]]
}

test_that("conditional LP objective streams retain bounded memory topology", {
  src_file <- locate_npc_cv_source("jksum.c")
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)

  owners <- c(
    "np_conditional_density_cvml_lp_stream",
    "np_conditional_density_cvls_lp_stream",
    "np_conditional_density_cvls_lp_stream_impl",
    "np_conditional_distribution_cvls_lp_stream",
    "np_conditional_lp_all_large_ctx_prepare_core"
  )
  for (owner in owners) {
    body <- np_test_extract_c_function(lines, owner)
    expect_false(
      grepl(
        paste0(
          "(malloc|calloc|alloc_matd|alloc_tmatd)\\(",
          "[^\\)]*(num_obs|num_obs_train_extern)",
          "[^\\)]*(num_obs|num_obs_train_extern)"
        ),
        body
      ),
      info = owner
    )
  }
})

test_that("conditional CVLS production uses canonical block cores", {
  src_file <- locate_npc_cv_source("jksum.c")
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  wrapper <- np_test_extract_c_function(
    lines, "np_conditional_density_cvls_lp_stream"
  )
  body <- np_test_extract_c_function(
    lines, "np_conditional_density_cvls_lp_stream_impl"
  )

  expect_match(
    wrapper, "np_conditional_density_cvls_lp_stream_impl(", fixed = TRUE
  )
  expect_true(grepl(
    "np_conditional_x_weight_block_stream_core(",
    body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "np_conditional_y_block_stream_op_core(",
    body,
    fixed = TRUE
  ))
})

test_that("conditional density LL remains canonical LP1 across tree modes", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026073111L)
  n <- 42L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y = x$x1 - x$x2 + rnorm(n, sd = 0.08))
  common <- list(
    xdat = x,
    ydat = y,
    bws = c(0.31, 0.36, 0.42),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls"
  )
  bw_ll <- do.call(npcdensbw, c(common, list(regtype = "ll")))
  bw_lp <- do.call(npcdensbw, c(
    common,
    list(
      regtype = "lp",
      basis = "glp",
      degree = c(1L, 1L),
      bernstein.basis = FALSE
    )
  ))

  options(np.tree = FALSE)
  ll_plain <- np:::.npcdensbw_eval_only(x, y, bw_ll)$objective
  lp_plain <- np:::.npcdensbw_eval_only(x, y, bw_lp)$objective
  options(np.tree = TRUE)
  ll_tree <- np:::.npcdensbw_eval_only(x, y, bw_ll)$objective
  lp_tree <- np:::.npcdensbw_eval_only(x, y, bw_lp)$objective

  expect_equal(ll_plain, lp_plain, tolerance = 1e-12)
  expect_equal(ll_tree, lp_tree, tolerance = 1e-12)
  expect_equal(lp_plain, lp_tree, tolerance = 1e-12)
})

test_that("conditional generalized-NN objectives retain LL-to-LP1 ownership", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026073112L)
  n <- 38L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(
    y = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.12)
  )
  common <- list(
    xdat = x,
    ydat = y,
    bws = c(4L, 7L, 5L),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE
  )
  density_ll <- do.call(npcdensbw, c(
    common,
    list(regtype = "ll", bwmethod = "cv.ml")
  ))
  density_lp <- do.call(npcdensbw, c(
    common,
    list(
      regtype = "lp",
      basis = "glp",
      degree = c(1L, 1L),
      bernstein.basis = FALSE,
      bwmethod = "cv.ml"
    )
  ))
  distribution_ll <- do.call(npcdistbw, c(
    common,
    list(regtype = "ll", bwmethod = "cv.ls")
  ))
  distribution_lp <- do.call(npcdistbw, c(
    common,
    list(
      regtype = "lp",
      basis = "glp",
      degree = c(1L, 1L),
      bernstein.basis = FALSE,
      bwmethod = "cv.ls"
    )
  ))

  density_ll_value <- np:::.npcdensbw_eval_only(
    x, y, density_ll
  )$objective
  density_lp_value <- np:::.npcdensbw_eval_only(
    x, y, density_lp
  )$objective
  distribution_ll_value <- np:::.npcdistbw_eval_only(
    x, y, bws = distribution_ll, do.full.integral = TRUE
  )$objective
  distribution_lp_value <- np:::.npcdistbw_eval_only(
    x, y, bws = distribution_lp, do.full.integral = TRUE
  )$objective

  expect_true(all(is.finite(c(
    density_ll_value,
    density_lp_value,
    distribution_ll_value,
    distribution_lp_value
  ))))
  expect_equal(density_ll_value, density_lp_value, tolerance = 1e-10)
  expect_equal(distribution_ll_value, distribution_lp_value, tolerance = 1e-10)
})

test_that("conditional LP objectives remain degree-sensitive", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026073113L)
  n <- 35L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y = x$x1^2 - x$x2 + rnorm(n, sd = 0.05))

  make_density <- function(degree) {
    npcdensbw(
      xdat = x,
      ydat = y,
      bws = c(0.24, 0.31, 0.42),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ml",
      regtype = "lp",
      basis = "glp",
      degree = rep.int(degree, ncol(x))
    )
  }
  make_distribution <- function(degree) {
    npcdistbw(
      xdat = x,
      ydat = y,
      bws = c(0.24, 0.31, 0.42),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      regtype = "lp",
      basis = "glp",
      degree = rep.int(degree, ncol(x))
    )
  }

  density_d1 <- np:::.npcdensbw_eval_only(
    x, y, make_density(1L)
  )$objective
  density_d2 <- np:::.npcdensbw_eval_only(
    x, y, make_density(2L)
  )$objective
  distribution_d1 <- np:::.npcdistbw_eval_only(
    x, y, bws = make_distribution(1L), do.full.integral = TRUE
  )$objective
  distribution_d2 <- np:::.npcdistbw_eval_only(
    x, y, bws = make_distribution(2L), do.full.integral = TRUE
  )$objective

  expect_gt(abs(density_d2 - density_d1), 1e-6)
  expect_gt(abs(distribution_d2 - distribution_d1), 1e-8)
})

test_that("all-large conditional LP routes retain canonical fast ownership", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026073114L)
  n <- 48L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rbeta(n, 1 + x$x, 2 - x$x))

  density_common <- list(
    xdat = x,
    ydat = y,
    bws = c(0.17, 1e6),
    bandwidth.compute = FALSE
  )
  for (method in c("cv.ml", "cv.ls")) {
    bw_ll <- do.call(npcdensbw, c(
      density_common,
      list(regtype = "ll", bwmethod = method)
    ))
    bw_lp <- do.call(npcdensbw, c(
      density_common,
      list(
        regtype = "lp",
        basis = "glp",
        degree = 1L,
        bernstein.basis = FALSE,
        bwmethod = method
      )
    ))
    eval_ll <- np:::.npcdensbw_eval_only(x, y, bw_ll)
    eval_lp <- np:::.npcdensbw_eval_only(x, y, bw_lp)
    expect_equal(eval_ll$objective, eval_lp$objective, tolerance = 1e-10)
    expect_gt(eval_ll$num.feval.fast, 0)
    expect_gt(eval_lp$num.feval.fast, 0)
  }

  distribution_common <- c(
    density_common,
    list(bwmethod = "cv.ls")
  )
  dist_ll <- do.call(
    npcdistbw,
    c(distribution_common, list(regtype = "ll"))
  )
  dist_lp <- do.call(npcdistbw, c(
    distribution_common,
    list(
      regtype = "lp",
      basis = "glp",
      degree = 1L,
      bernstein.basis = FALSE
    )
  ))
  eval_ll <- np:::.npcdistbw_eval_only(
    x, y, bws = dist_ll, do.full.integral = TRUE
  )
  eval_lp <- np:::.npcdistbw_eval_only(
    x, y, bws = dist_lp, do.full.integral = TRUE
  )
  expect_equal(eval_ll$objective, eval_lp$objective, tolerance = 1e-10)
  expect_gt(eval_ll$num.feval.fast, 0)
  expect_gt(eval_lp$num.feval.fast, 0)
})
