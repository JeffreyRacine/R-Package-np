library(np)

empty_shadow_matrix <- function(n) {
  matrix(numeric(0), nrow = n, ncol = 0)
}

shadow_cker_code <- function(kernel) {
  switch(kernel,
    gaussian = 0L,
    epanechnikov = 4L,
    uniform = 8L,
    truncated = 9L,
    stop("unsupported continuous kernel for shadow test")
  )
}

shadow_uker_code <- function(kernel) {
  switch(kernel,
    aitchisonaitken = 0L,
    liracine = 1L,
    stop("unsupported unordered kernel for shadow test")
  )
}

shadow_oker_code <- function(kernel) {
  switch(kernel,
    wangvanryzin = 0L,
    liracine = 2L,
    racineliyan = 3L,
    stop("unsupported ordered kernel for shadow test")
  )
}

shadow_basis_code <- function(basis) {
  switch(basis,
    additive = 0L,
    glp = 1L,
    tensor = 2L,
    stop("unsupported LP basis for shadow test")
  )
}

shadow_rbw <- function(bw) {
  c(
    bw$xbw[bw$ixcon],
    bw$ybw[bw$iycon],
    bw$ybw[bw$iyuno],
    bw$ybw[bw$iyord],
    bw$xbw[bw$ixuno],
    bw$xbw[bw$ixord]
  )
}

shadow_bwtype_code <- function(bw) {
  switch(bw$type,
    fixed = 0L,
    generalized_nn = 1L,
    adaptive_nn = 2L,
    stop("unsupported bandwidth type for shadow test")
  )
}

shadow_reg_code <- function(bw) {
  if (identical(bw$regtype.engine, "lp")) 2L else 0L
}

shadow_degree <- function(bw) {
  if (identical(bw$regtype.engine, "lp")) as.integer(bw$degree.engine) else integer(0)
}

shadow_safe_call <- function(name, ...) {
  on.exit(
    tryCatch(.Call("C_np_shadow_reset_state", PACKAGE = "np"),
             error = function(e) NULL),
    add = TRUE
  )
  .Call(name, ..., PACKAGE = "np")
}

call_shadow_density <- function(bw, x, y, criterion = c("cv.ml", "cv.ls"),
                                use_tree = FALSE, compare_old = TRUE) {
  criterion <- match.arg(criterion)
  n <- nrow(x)
  shadow_safe_call(
    "C_np_shadow_cv_density_conditional",
    empty_shadow_matrix(n), empty_shadow_matrix(n), as.matrix(y),
    empty_shadow_matrix(n), empty_shadow_matrix(n), as.matrix(x),
    as.double(shadow_rbw(bw)),
    shadow_bwtype_code(bw),
    shadow_cker_code(bw$cykertype),
    shadow_uker_code(bw$uykertype),
    shadow_oker_code(bw$oykertype),
    shadow_cker_code(bw$cxkertype),
    shadow_uker_code(bw$uxkertype),
    shadow_oker_code(bw$oxkertype),
    use_tree,
    if (identical(criterion, "cv.ml")) 0L else 1L,
    shadow_reg_code(bw),
    shadow_degree(bw),
    isTRUE(bw$bernstein.basis.engine),
    shadow_basis_code(bw$basis.engine),
    compare_old
  )
}

call_shadow_distribution <- function(bw, x, ytrain, yeval = ytrain,
                                     use_tree = FALSE, cdfontrain = FALSE,
                                     compare_old = TRUE) {
  n <- nrow(x)
  ne <- nrow(yeval)
  shadow_safe_call(
    "C_np_shadow_cv_distribution_conditional",
    empty_shadow_matrix(n), empty_shadow_matrix(n), as.matrix(ytrain),
    empty_shadow_matrix(ne), empty_shadow_matrix(ne), as.matrix(yeval),
    empty_shadow_matrix(n), empty_shadow_matrix(n), as.matrix(x),
    as.double(shadow_rbw(bw)),
    shadow_bwtype_code(bw),
    shadow_cker_code(bw$cykertype),
    shadow_uker_code(bw$uykertype),
    shadow_oker_code(bw$oykertype),
    shadow_cker_code(bw$cxkertype),
    shadow_uker_code(bw$uxkertype),
    shadow_oker_code(bw$oxkertype),
    use_tree,
    shadow_reg_code(bw),
    shadow_degree(bw),
    isTRUE(bw$bernstein.basis.engine),
    shadow_basis_code(bw$basis.engine),
    cdfontrain,
    compare_old
  )
}

call_shadow_xweights_row <- function(bw, x, y, row_index, use_tree = FALSE) {
  n <- nrow(x)
  shadow_safe_call(
    "C_np_shadow_cv_xweights_conditional",
    empty_shadow_matrix(n), empty_shadow_matrix(n), as.matrix(y),
    empty_shadow_matrix(n), empty_shadow_matrix(n), as.matrix(x),
    as.double(shadow_rbw(bw)),
    shadow_bwtype_code(bw),
    shadow_cker_code(bw$cxkertype),
    shadow_uker_code(bw$uxkertype),
    shadow_oker_code(bw$oxkertype),
    use_tree,
    shadow_reg_code(bw),
    shadow_degree(bw),
    isTRUE(bw$bernstein.basis.engine),
    shadow_basis_code(bw$basis.engine),
    as.integer(row_index)
  )
}

test_that("shadow reset state is a harmless no-op when inactive", {
  expect_null(.Call("C_np_shadow_reset_state", PACKAGE = "np"))
  expect_null(.Call("C_np_shadow_reset_state", PACKAGE = "np"))
})

test_that("shadow density lc matches legacy cv objectives", {
  set.seed(42)
  n <- 40L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.35, 0.45, 0.55),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  res_ml <- call_shadow_density(bw, x, y, criterion = "cv.ml", compare_old = TRUE)
  res_ls <- call_shadow_density(bw, x, y, criterion = "cv.ls", compare_old = TRUE)

  expect_equal(res_ml$new, res_ml$old, tolerance = 1e-12)
  expect_equal(res_ls$new, res_ls$old, tolerance = 2e-3)
})

test_that("shadow density lp preserves ll canonicalization and tree parity", {
  set.seed(33)
  n <- 42L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.15))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.31, 0.43, 0.57),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.31, 0.43, 0.57),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree
  )

  ll_ml <- call_shadow_density(bw.ll, x, y, criterion = "cv.ml", compare_old = FALSE)
  lp_ml <- call_shadow_density(bw.lp, x, y, criterion = "cv.ml", compare_old = FALSE)
  ll_ls <- call_shadow_density(bw.ll, x, y, criterion = "cv.ls", compare_old = FALSE)
  lp_ls <- call_shadow_density(bw.lp, x, y, criterion = "cv.ls", compare_old = FALSE)
  lp_ml_tree <- call_shadow_density(bw.lp, x, y, criterion = "cv.ml", use_tree = TRUE, compare_old = FALSE)
  lp_ls_tree <- call_shadow_density(bw.lp, x, y, criterion = "cv.ls", use_tree = TRUE, compare_old = FALSE)

  expect_equal(ll_ml$new, lp_ml$new, tolerance = 1e-12)
  expect_equal(ll_ls$new, lp_ls$new, tolerance = 1e-12)
  expect_equal(lp_ml_tree$new, lp_ml$new, tolerance = 1e-10)
  expect_equal(lp_ls_tree$new, lp_ls$new, tolerance = 1e-12)
  expect_equal(ll_ml$prod, ll_ml$new, tolerance = 1e-12)
  expect_equal(lp_ml$prod, lp_ml$new, tolerance = 1e-12)
  expect_equal(lp_ml_tree$prod, lp_ml_tree$new, tolerance = 1e-10)
  expect_equal(ll_ls$prod, ll_ls$new, tolerance = 1e-12)
  expect_equal(lp_ls$prod, lp_ls$new, tolerance = 1e-12)
  expect_equal(lp_ls_tree$prod, lp_ls_tree$new, tolerance = 1e-12)
})

test_that("shadow density cvml preserves the large-kernel X collapse", {
  set.seed(11)
  n <- 30L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))
  big <- 1e6
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(big, big, big),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  res <- call_shadow_density(bw, x, y, criterion = "cv.ml", compare_old = FALSE)
  collapsed <- -n * log(mean(dnorm(outer(y$y1, y$y1, "-"), sd = big)))

  expect_equal(res$new, collapsed, tolerance = 1e-9)
})

test_that("shadow density generalized-nn cvml preserves ll canonicalization", {
  set.seed(202)
  n <- 38L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.12))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(4, 7, 5),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(4, 7, 5),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree
  )

  ll_ml <- call_shadow_density(bw.ll, x, y, criterion = "cv.ml", compare_old = FALSE)
  lp_ml <- call_shadow_density(bw.lp, x, y, criterion = "cv.ml", compare_old = FALSE)

  expect_true(is.finite(ll_ml$new))
  expect_true(is.finite(lp_ml$new))
  expect_equal(ll_ml$new, lp_ml$new, tolerance = 1e-10)
  expect_true(is.finite(ll_ml$prod))
  expect_true(is.finite(lp_ml$prod))
  expect_equal(ll_ml$prod, lp_ml$prod, tolerance = 1e-10)
})

test_that("shadow density generalized-nn LP cures the legacy penalty collapse", {
  set.seed(203)
  n <- 38L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.12))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(4, 7, 5),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(4, 7, 5),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree
  )

  ll_ls <- call_shadow_density(bw.ll, x, y, criterion = "cv.ls", compare_old = FALSE)
  lp_ls <- call_shadow_density(bw.lp, x, y, criterion = "cv.ls", compare_old = FALSE)

  expect_true(is.finite(ll_ls$new))
  expect_true(is.finite(lp_ls$new))
  expect_equal(ll_ls$new, lp_ls$new, tolerance = 1e-10)
  expect_equal(ll_ls$prod, ll_ls$new, tolerance = 1e-10)
  expect_equal(lp_ls$prod, lp_ls$new, tolerance = 1e-10)
  expect_gt(ll_ls$new, -1e6)
})

test_that("shadow fixed-bandwidth X-side row helper matches dense oracle and leaves self out", {
  set.seed(512)
  n <- 18L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.08))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.28, 0.37, 0.41),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  for (row in seq_len(n)) {
    res <- call_shadow_xweights_row(bw, x, y, row_index = row)
    expect_equal(res$streamed, res$dense, tolerance = 1e-12)
    expect_equal(res$streamed[row], 0, tolerance = 1e-12)
    expect_equal(sum(res$streamed), 1, tolerance = 1e-10)
  }
})

test_that("shadow fixed-bandwidth X-side row helper preserves ll == lp degree-1 and tree parity", {
  set.seed(513)
  n <- 20L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 - x$x2 + rnorm(n, sd = 0.06))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.31, 0.36, 0.42),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.31, 0.36, 0.42),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree
  )

  res.ll <- call_shadow_xweights_row(bw.ll, x, y, row_index = 7L)
  res.lp <- call_shadow_xweights_row(bw.lp, x, y, row_index = 7L)
  res.tree <- call_shadow_xweights_row(bw.lp, x, y, row_index = 7L, use_tree = TRUE)

  expect_equal(res.ll$streamed, res.lp$streamed, tolerance = 1e-12)
  expect_equal(res.lp$streamed, res.tree$streamed, tolerance = 1e-12)
})

test_that("shadow generalized-nn X-side row helper matches dense oracle and preserves ll == lp", {
  set.seed(514)
  n <- 22L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + x$x2 + rnorm(n, sd = 0.07))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(5, 7, 6),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(5, 7, 6),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree
  )

  res.ll <- call_shadow_xweights_row(bw.ll, x, y, row_index = 9L)
  res.lp <- call_shadow_xweights_row(bw.lp, x, y, row_index = 9L)

  expect_equal(res.ll$streamed, res.ll$dense, tolerance = 1e-12)
  expect_equal(res.lp$streamed, res.lp$dense, tolerance = 1e-12)
  expect_equal(res.ll$streamed, res.lp$streamed, tolerance = 1e-10)
  expect_equal(sum(res.ll$streamed), 1, tolerance = 1e-10)
  expect_equal(res.ll$streamed[9L], 0, tolerance = 1e-12)
})

locate_shadow_proof_src <- function(filename) {
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
  if (length(hits) == 0L) {
    return(NULL)
  }
  hits[[1L]]
}

read_shadow_proof_src <- function(filename) {
  src_file <- locate_shadow_proof_src(filename)
  skip_if(
    is.null(src_file),
    sprintf("source file src/%s unavailable in this test context", filename)
  )
  readLines(src_file, warn = FALSE)
}

test_that("shadow fixed-bandwidth X-side row helper stays row-streamed", {
  lines <- read_shadow_proof_src("jksum.c")
  start <- grep("^int np_shadow_proof_conditional_x_weight_row_fixed\\(", lines)
  stop <- grep("^static int np_shadow_conditional_build_y_matrix\\(", lines)

  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_lt(start, stop)

  body <- paste(lines[(start + 1L):(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_train|num_obs_train_extern)[^\\)]*(num_train|num_obs_train_extern)", body))
})

test_that("shadow generalized-nn X-side row helper stays row-streamed", {
  lines <- read_shadow_proof_src("jksum.c")
  start <- grep("^int np_shadow_proof_conditional_x_weight_row_stream\\(", lines)
  stop <- grep("^static int np_shadow_conditional_build_y_matrix\\(", lines)

  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_lt(start, stop)

  body <- paste(lines[(start + 1L):(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_train|num_obs_train_extern)[^\\)]*(num_train|num_obs_train_extern)", body))
})

test_that("fixed-bandwidth cvml LP stream avoids dense row storage", {
  lines <- read_shadow_proof_src("jksum.c")
  start <- grep("^int np_conditional_density_cvml_lp_stream\\(", lines)
  stop <- grep("^static int np_shadow_conditional_build_y_matrix\\(", lines)

  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_lt(start, stop)

  body <- paste(lines[start:(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_obs|num_obs_train_extern)[^\\)]*(num_obs|num_obs_train_extern)", body))
})

test_that("generalized-nn cvml LP stream avoids dense row storage", {
  lines <- read_shadow_proof_src("jksum.c")
  start <- grep("^int np_conditional_density_cvml_lp_stream\\(", lines)
  stop <- grep("^static int np_shadow_conditional_build_y_matrix\\(", lines)

  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_lt(start, stop)

  body <- paste(lines[start:(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_obs|num_obs_train_extern)[^\\)]*(num_obs|num_obs_train_extern)", body))
})

test_that("generalized-nn cvls LP stream avoids dense row storage", {
  lines <- read_shadow_proof_src("jksum.c")
  start <- grep("^int np_conditional_density_cvls_lp_stream\\(", lines)
  stop <- grep("^static int np_shadow_conditional_build_y_matrix\\(", lines)

  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_lt(start, stop)

  body <- paste(lines[start:(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_obs|num_obs_train_extern)[^\\)]*(num_obs|num_obs_train_extern)", body))
})

test_that("conditional density cvls production bypasses shadow block helpers", {
  lines <- read_shadow_proof_src("jksum.c")
  start <- grep("^int np_conditional_density_cvls_lp_stream\\(", lines)
  stop <- grep("^int np_conditional_distribution_cvls_lp_stream\\(", lines)

  expect_length(start, 1L)
  expect_gte(length(stop), 1L)
  expect_lt(start, stop[1L])

  body <- paste(lines[start:(stop[1L] - 1L)], collapse = "\n")

  expect_false(grepl("np_shadow_proof_conditional_x_weight_block_stream\\s*\\(", body))
  expect_false(grepl("np_shadow_conditional_y_block_stream_op\\s*\\(", body))
  expect_true(grepl("np_conditional_x_weight_block_stream_core\\s*\\(", body))
  expect_true(grepl("np_conditional_y_block_stream_op_core\\s*\\(", body))
})

test_that("fixed-bandwidth cvls LP stream avoids dense row storage", {
  lines <- read_shadow_proof_src("jksum.c")
  start_matches <- grep("^int np_conditional_density_cvls_lp_stream\\(", lines)
  stop <- grep("^static int np_shadow_conditional_build_y_matrix\\(", lines)

  expect_gte(length(start_matches), 1L)
  expect_length(stop, 1L)
  start <- tail(start_matches, 1L)
  expect_lt(start, stop)

  body <- paste(lines[start:(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_obs|num_obs_train_extern)[^\\)]*(num_obs|num_obs_train_extern)", body))
})

test_that("shadow LP objectives are sensitive to degree on oracle cells", {
  set.seed(123)
  n <- 35L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1^2 - x$x2 + rnorm(n, sd = 0.05))

  bw.d1 <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.24, 0.31, 0.42),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = rep.int(1L, ncol(x))
  )
  bw.d2 <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.24, 0.31, 0.42),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = rep.int(2L, ncol(x))
  )
  dens.d1 <- call_shadow_density(bw.d1, x, y, criterion = "cv.ml", compare_old = FALSE)
  dens.d2 <- call_shadow_density(bw.d2, x, y, criterion = "cv.ml", compare_old = FALSE)

  bw.c1 <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.24, 0.31, 0.42),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = rep.int(1L, ncol(x))
  )
  bw.c2 <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.24, 0.31, 0.42),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = rep.int(2L, ncol(x))
  )
  dist.d1 <- call_shadow_distribution(bw.c1, x, y, compare_old = FALSE)
  dist.d2 <- call_shadow_distribution(bw.c2, x, y, compare_old = FALSE)

  expect_gt(abs(dens.d2$new - dens.d1$new), 1e-6)
  expect_gt(abs(dist.d2$new - dist.d1$new), 1e-8)
})

test_that("shadow distribution lc matches legacy cvls for both cdfontrain modes", {
  set.seed(55)
  n <- 28L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))
  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.33, 0.44, 0.51),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  res_false <- call_shadow_distribution(bw, x, y, cdfontrain = FALSE, compare_old = TRUE)
  res_true <- call_shadow_distribution(bw, x, y, cdfontrain = TRUE, compare_old = TRUE)

  expect_equal(res_false$new, res_false$old, tolerance = 1e-12)
  expect_equal(res_true$new, res_true$old, tolerance = 1e-12)
})

test_that("shadow distribution generalized-nn LP preserves ll canonicalization", {
  set.seed(263)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(4, 6, 5),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(4, 6, 5),
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree
  )

  dist.ll <- call_shadow_distribution(bw.ll, x, y, compare_old = FALSE)
  dist.lp <- call_shadow_distribution(bw.lp, x, y, compare_old = FALSE)

  expect_true(is.finite(dist.ll$new))
  expect_true(is.finite(dist.lp$new))
  expect_equal(dist.ll$new, dist.lp$new, tolerance = 1e-10)
  expect_equal(dist.ll$prod, dist.ll$new, tolerance = 1e-10)
  expect_equal(dist.lp$prod, dist.lp$new, tolerance = 1e-10)
})

test_that("shadow distribution lp preserves ll canonicalization and tree parity", {
  set.seed(88)
  n <- 34L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 - x$x2 + rnorm(n, sd = 0.2))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.28, 0.39, 0.52),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.28, 0.39, 0.52),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree
  )

  ll_res <- call_shadow_distribution(bw.ll, x, y, compare_old = FALSE)
  lp_res <- call_shadow_distribution(bw.lp, x, y, compare_old = FALSE)
  lp_tree <- call_shadow_distribution(bw.lp, x, y, use_tree = TRUE, compare_old = FALSE)

  expect_equal(ll_res$new, lp_res$new, tolerance = 1e-12)
  expect_equal(lp_tree$new, lp_res$new, tolerance = 1e-12)
  expect_equal(ll_res$prod, ll_res$new, tolerance = 1e-12)
  expect_equal(lp_res$prod, lp_res$new, tolerance = 1e-12)
})

test_that("fixed-bandwidth cdist cvls LP stream avoids dense row storage", {
  lines <- read_shadow_proof_src("jksum.c")
  start_matches <- grep("^int np_conditional_distribution_cvls_lp_stream\\(", lines)
  stop <- grep("^static int np_shadow_conditional_build_y_matrix\\(", lines)

  expect_gte(length(start_matches), 1L)
  expect_length(stop, 1L)
  start <- tail(start_matches, 1L)
  expect_lt(start, stop)

  body <- paste(lines[start:(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_obs|num_obs_train_extern)[^\\)]*(num_obs|num_obs_train_extern)", body))
})

test_that("generalized-nn cdist cvls LP stream avoids dense row storage", {
  lines <- read_shadow_proof_src("jksum.c")
  start_matches <- grep("^int np_conditional_distribution_cvls_lp_stream\\(", lines)
  stop <- grep("^static int np_shadow_conditional_build_y_matrix\\(", lines)

  expect_gte(length(start_matches), 1L)
  expect_length(stop, 1L)
  start <- tail(start_matches, 1L)
  expect_lt(start, stop)

  body <- paste(lines[start:(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_obs|num_obs_train_extern)[^\\)]*(num_obs|num_obs_train_extern)", body))
})

test_that("LP all-large conditional fast helpers stay memory-lean", {
  lines <- read_shadow_proof_src("jksum.c")
  start <- grep("^static int np_conditional_lp_all_large_ctx_prepare\\(", lines)
  stop <- grep("^int np_conditional_density_cvml_lp_stream\\(", lines)

  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_lt(start, stop)

  body <- paste(lines[start:(stop - 1L)], collapse = "\n")

  expect_false(grepl("(malloc|calloc|alloc_matd|alloc_tmatd)\\([^\\)]*(num_train|num_obs_train_extern)[^\\)]*(num_train|num_obs_train_extern)", body))
  expect_false(grepl("\\bdiag\\s*\\(", body))
})

test_that("shadow density LP all-large fixed objectives stay exact and count fast hits", {
  set.seed(615)
  n <- 48L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rbeta(n, 1 + x$x, 2 - x$x))
  degree <- 1L

  bw.ll.ml <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.17, 1e6),
    bandwidth.compute = FALSE,
    regtype = "ll",
    bwmethod = "cv.ml"
  )
  bw.lp.ml <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.17, 1e6),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ml"
  )
  bw.ll.ls <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.17, 1e6),
    bandwidth.compute = FALSE,
    regtype = "ll",
    bwmethod = "cv.ls"
  )
  bw.lp.ls <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.17, 1e6),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls"
  )

  sh.ll.ml <- call_shadow_density(bw.ll.ml, x, y, criterion = "cv.ml", compare_old = FALSE)
  sh.lp.ml <- call_shadow_density(bw.lp.ml, x, y, criterion = "cv.ml", compare_old = FALSE)
  sh.ll.ls <- call_shadow_density(bw.ll.ls, x, y, criterion = "cv.ls", compare_old = FALSE)
  sh.lp.ls <- call_shadow_density(bw.lp.ls, x, y, criterion = "cv.ls", compare_old = FALSE)

  ev.ll.ml <- np:::.npcdensbw_eval_only(x, y, bw.ll.ml)
  ev.lp.ml <- np:::.npcdensbw_eval_only(x, y, bw.lp.ml)
  ev.ll.ls <- np:::.npcdensbw_eval_only(x, y, bw.ll.ls)
  ev.lp.ls <- np:::.npcdensbw_eval_only(x, y, bw.lp.ls)

  expect_equal(-ev.ll.ml$objective, sh.ll.ml$new, tolerance = 1e-10)
  expect_equal(-ev.lp.ml$objective, sh.lp.ml$new, tolerance = 1e-10)
  expect_equal(-ev.ll.ls$objective, sh.ll.ls$new, tolerance = 1e-10)
  expect_equal(-ev.lp.ls$objective, sh.lp.ls$new, tolerance = 1e-10)
  expect_equal(ev.ll.ml$objective, ev.lp.ml$objective, tolerance = 1e-10)
  expect_equal(ev.ll.ls$objective, ev.lp.ls$objective, tolerance = 1e-10)
  expect_gt(ev.ll.ml$num.feval.fast, 0)
  expect_gt(ev.lp.ml$num.feval.fast, 0)
  expect_gt(ev.ll.ls$num.feval.fast, 0)
  expect_gt(ev.lp.ls$num.feval.fast, 0)
})

test_that("shadow distribution LP all-large fixed objective stays exact and counts fast hits", {
  set.seed(616)
  n <- 44L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rbeta(n, 1 + x$x, 2 - x$x))
  degree <- 1L

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.17, 1e6),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.17, 1e6),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = degree
  )

  sh.ll <- call_shadow_distribution(
    bw.ll,
    x,
    y,
    yeval = y,
    cdfontrain = TRUE,
    compare_old = FALSE
  )
  sh.lp <- call_shadow_distribution(
    bw.lp,
    x,
    y,
    yeval = y,
    cdfontrain = TRUE,
    compare_old = FALSE
  )

  ev.ll <- np:::.npcdistbw_eval_only(x, y, bws = bw.ll, do.full.integral = TRUE)
  ev.lp <- np:::.npcdistbw_eval_only(x, y, bws = bw.lp, do.full.integral = TRUE)

  expect_equal(ev.ll$objective, sh.ll$new, tolerance = 1e-10)
  expect_equal(ev.lp$objective, sh.lp$new, tolerance = 1e-10)
  expect_equal(ev.ll$objective, ev.lp$objective, tolerance = 1e-10)
  expect_gt(ev.ll$num.feval.fast, 0)
  expect_gt(ev.lp$num.feval.fast, 0)
})

test_that("kernelcv no longer references dense shadow proof helpers", {
  lines <- read_shadow_proof_src("kernelcv.c")

  expect_false(any(grepl("np_shadow_cv_con_density_ml\\s*\\(", lines)))
  expect_false(any(grepl("np_shadow_cv_con_density_ls\\s*\\(", lines)))
  expect_false(any(grepl("np_shadow_cv_con_distribution_ls\\s*\\(", lines)))
})
