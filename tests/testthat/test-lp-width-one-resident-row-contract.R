locate_lp_row_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum_lp_row.c"),
    test_path("..", "..", "..", "src", "jksum_lp_row.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum_lp_row.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum_lp_row.c"),
    file.path(getwd(), "src", "jksum_lp_row.c"),
    file.path(getwd(), "..", "src", "jksum_lp_row.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits))
    return(NULL)
  hits[[1L]]
}

test_that("fixed-objective width-one rows use an implicit scalar basis", {
  row_file <- locate_lp_row_source()
  skip_if(is.null(row_file), "source file src/jksum_lp_row.c unavailable")
  source <- paste(readLines(row_file, warn = FALSE), collapse = "\n")

  start <- regexpr(
    "static void np_lp_accumulate_dense_row_1(",
    source,
    fixed = TRUE
  )
  stop <- regexpr(
    "#define NP_LP_DEFINE_RESIDENT_WIDTH(WIDTH)",
    source,
    fixed = TRUE
  )
  expect_gt(start, 0L)
  expect_gt(stop, start)
  scalar <- substr(source, start, stop - 1L)

  expect_true(grepl("if(!ctx->use_tree)", scalar, fixed = TRUE))
  expect_true(grepl("fixed_rhs += weight*ctx->response[orig_ii]", scalar,
                    fixed = TRUE))
  expect_true(grepl("ctx->rhs[orig_ii] += weight*eval_response", scalar,
                    fixed = TRUE))
  expect_true(grepl("fixed_moment += weight", scalar, fixed = TRUE))
  expect_true(grepl("ctx->moments[orig_ii] += weight", scalar,
                    fixed = TRUE))
  expect_true(grepl("np_lp_dense_support_add(", scalar, fixed = TRUE))
  expect_false(grepl("ctx->basis", scalar, fixed = TRUE))
  expect_false(grepl("F77_", scalar, fixed = TRUE))
  expect_false(grepl("LL_LC", scalar, fixed = TRUE))
  expect_false(grepl("NP_LP_DEFINE_RESIDENT_WIDTH(1)", source, fixed = TRUE))
  expect_true(grepl("NP_LP_DEFINE_RESIDENT_WIDTH(2)", source, fixed = TRUE))
})

test_that("sparse-tree width-one pairs use the same implicit scalar algebra", {
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
  skip_if(!length(hits), "source file src/jksum.c unavailable")
  source <- paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")

  start <- regexpr("static inline void np_lp_accumulate_pair(", source,
                   fixed = TRUE)
  stop <- regexpr("static inline void np_lp_cvls_support_add(", source,
                  fixed = TRUE)
  expect_gt(start, 0L)
  expect_gt(stop, start)
  pair <- substr(source, start, stop - 1L)

  scalar_start <- regexpr("if(nterms == 1)", pair, fixed = TRUE)
  wider_start <- regexpr("if(nterms == 2)", pair, fixed = TRUE)
  expect_gt(scalar_start, 0L)
  expect_gt(wider_start, scalar_start)
  scalar <- substr(pair, scalar_start, wider_start - 1L)
  expect_true(grepl("tj[0] += w*yi", scalar, fixed = TRUE))
  expect_true(grepl("ti[0] += w*eval_ybasis[0]", scalar, fixed = TRUE))
  expect_true(grepl("sj[0] += w", scalar, fixed = TRUE))
  expect_true(grepl("si[0] += w", scalar, fixed = TRUE))
  expect_false(grepl("(^|[^[:alnum:]_])basis\\[", scalar, perl = TRUE))
  expect_false(grepl("F77_", scalar, fixed = TRUE))
})

test_that("raw and Bernstein degree-zero objectives retain one scalar result", {
  package_owner <- environmentName(environment(npregbw))
  skip_if(
    identical(package_owner, "npRmpi"),
    "npRmpi numerical equality is covered by active-pool parity tests"
  )
  eval_only <- get(
    ".npregbw_eval_only",
    envir = asNamespace(package_owner),
    inherits = FALSE
  )

  set.seed(20260725)
  n <- 96L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(2 * pi * x$x1) * cos(2 * pi * x$x2) +
    rnorm(n, sd = 0.2)

  make_bw <- function(bernstein) {
    npregbw(
      xdat = x,
      ydat = y,
      bws = c(0.22, 0.24),
      regtype = "lp",
      basis = "glp",
      degree = c(0L, 0L),
      bernstein.basis = bernstein,
      bwmethod = "cv.ls",
      bandwidth.compute = FALSE
    )
  }

  raw <- eval_only(x, y, make_bw(FALSE))$objective
  bernstein <- eval_only(x, y, make_bw(TRUE))$objective
  expect_identical(raw, bernstein)
})
