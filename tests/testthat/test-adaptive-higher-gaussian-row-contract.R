library(np)

locate_adaptive_higher_gaussian_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum.c"),
    test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) NULL else hits[[1L]]
}

adaptive_higher_gaussian_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("adaptive higher-order Gaussian row fusion is narrow and bounded", {
  src_file <- locate_adaptive_higher_gaussian_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  text <- paste(lines, collapse = "\n")
  body <- adaptive_higher_gaussian_body(
    lines,
    "^static int NP_NOINLINE np_accel_gauss_adaptive_higher_row_try\\(",
    "^/\\*"
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(
    text, "#define NP_ADAPTIVE_HIGHER_GAUSS_MIN_ROWS 1", fixed = TRUE
  )
  expect_match(compact, "(out == NULL) || (ndim < 2)", fixed = TRUE)
  expect_match(
    compact, "n < NP_ADAPTIVE_HIGHER_GAUSS_MIN_ROWS", fixed = TRUE
  )
  expect_match(
    compact,
    "(kernel_c[d] < 1) || (kernel_c[d] > 3) || (operator[d] != OP_NORMAL)",
    fixed = TRUE
  )
  expect_match(body, "np_accel_gauss_scratch_ensure(n)", fixed = TRUE)
  expect_match(
    body, "np_accel_gauss_polynomial_vector(", fixed = TRUE
  )
  expect_match(body, "vvexp(np_accel_gauss_tmp", fixed = TRUE)
  expect_match(body, "if(!isfinite(out[i]))", fixed = TRUE)
  expect_match(body, "memset(out, 0, (size_t)n*sizeof(double));", fixed = TRUE)
  expect_false(grepl("n\\s*\\*\\s*n", body))
  expect_false(grepl("malloc|calloc|realloc", body))

  calls <- gregexpr(
    "np_accel_gauss_adaptive_higher_row_try\\(", text, perl = TRUE
  )[[1L]]
  expect_length(calls[calls > 0L], 2L)

  conditional <- adaptive_higher_gaussian_body(
    lines,
    "^static int .*np_conditional_xrow_from_ctx_impl\\(",
    "^static int np_conditional_xrow_from_ctx\\("
  )
  expect_match(
    conditional, "BANDWIDTH_den_extern == BW_ADAP_NN", fixed = TRUE
  )
  expect_match(
    conditional, "num_reg_unordered_extern == 0", fixed = TRUE
  )
  expect_match(
    conditional, "num_reg_ordered_extern == 0", fixed = TRUE
  )
  expect_match(conditional, "!int_cxker_bound_extern", fixed = TRUE)
  expect_match(
    conditional, "int_TREE_X != NP_TREE_TRUE", fixed = TRUE
  )
  expect_match(
    conditional,
    "np_shadow_conditional_kernel_row_raw(ctx->kernel_cx,",
    fixed = TRUE
  )
})

test_that("adaptive higher-order conditional objectives retain scalar oracle", {
  old_options <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.largeh = FALSE,
    matprod = "default"
  )
  on.exit(options(old_options), add = TRUE)

  set.seed(2026073317L)
  n <- 96L
  xdat <- data.frame(
    x1 = runif(n, -0.85, 0.9),
    x2 = runif(n, -0.8, 0.88)
  )
  ydat <- data.frame(
    y = sin(1.4 * xdat$x1) + 0.3 * xdat$x2^2 + rnorm(n, sd = 0.2)
  )
  bws <- c(31L, 29L, 27L)

  objective <- function(route, order, degree, accelerate) {
    options(np.macMseries.accelerate = accelerate)
    common <- list(
      xdat = xdat,
      ydat = ydat,
      bws = bws,
      bandwidth.compute = FALSE,
      bwtype = "adaptive_nn",
      regtype = "lp",
      basis = "glp",
      degree = degree,
      bernstein.basis = FALSE,
      cxkertype = "gaussian",
      cxkerorder = order,
      cykertype = "gaussian",
      cykerorder = 2L
    )
    if (route == "cdens_cvml") {
      state <- do.call(npcdensbw, c(common, list(bwmethod = "cv.ml")))
      return(as.numeric(
        np:::.npcdensbw_eval_only(xdat, ydat, state)$objective
      ))
    }
    if (route == "cdens_cvls") {
      state <- do.call(npcdensbw, c(common, list(bwmethod = "cv.ls")))
      return(as.numeric(
        np:::.npcdensbw_eval_only(xdat, ydat, state)$objective
      ))
    }
    state <- do.call(npcdistbw, c(common, list(bwmethod = "cv.ls")))
    as.numeric(np:::.npcdistbw_eval_only(
      xdat, ydat, bws = state, do.full.integral = TRUE
    )$objective)
  }

  for (route in c("cdens_cvml", "cdens_cvls", "cdist_cvls")) {
    for (order in c(4L, 6L, 8L)) {
      for (degree in list(c(0L, 0L), c(2L, 2L))) {
        scalar <- objective(route, order, degree, FALSE)
        vector <- objective(route, order, degree, TRUE)
        expect_equal(vector, scalar, tolerance = 3e-9)
      }
    }
  }
})
