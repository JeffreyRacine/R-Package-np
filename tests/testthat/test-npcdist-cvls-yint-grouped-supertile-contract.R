library(npRmpi)

locate_mpi_cdist_grouped_supertile_source <- function() {
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

mpi_cdist_grouped_supertile_body <- function(lines) {
  start <- grep(
    "^static int np_conditional_distribution_cvls_lp_supertile\\(",
    lines
  )
  stop <- grep("^#if defined\\(__clang__\\) \\|\\| defined\\(__GNUC__\\)", lines)
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("MPI conditional-distribution grouped supertile is bounded and progressive", {
  src_file <- locate_mpi_cdist_grouped_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- mpi_cdist_grouped_supertile_body(readLines(src_file, warn = FALSE))
  flat <- gsub("[[:space:]]+", " ", body)

  expect_match(body, "double **xblocks[4]", fixed = TRUE)
  expect_match(flat, "MIN(4, (nblocks / owned_stride) +", fixed = TRUE)
  expect_match(body, "int local_group_width = 1;", fixed = TRUE)
  expect_match(body, "group_width = local_group_width;", fixed = TRUE)
  expect_equal(lengths(regmatches(
    body,
    gregexpr("np_optional_tmatd\\(num_train, block_size\\)", body, perl = TRUE)
  )), 3L)
  expect_match(body, "if(requested_group_width >= 3)", fixed = TRUE)
  expect_match(body, "if(requested_group_width >= 4)", fixed = TRUE)
  expect_false(grepl("num_train\\s*\\*\\s*num_train", body, perl = TRUE))
})

test_that("MPI grouped supertile agrees on one allocation width across ranks", {
  src_file <- locate_mpi_cdist_grouped_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- mpi_cdist_grouped_supertile_body(readLines(src_file, warn = FALSE))
  flat <- gsub("[[:space:]]+", " ", body)

  expect_match(
    flat,
    "MPI_Allreduce(&local_group_width, &group_width, 1, MPI_INT, MPI_MIN, comm[1]);",
    fixed = TRUE
  )
  expect_match(body, "if(group_width < 2)", fixed = TRUE)
  expect_match(body, "status = 2;", fixed = TRUE)
  expect_match(body, "first_block + g*owned_stride", fixed = TRUE)
  expect_match(body, "group_width*owned_stride", fixed = TRUE)
})

test_that("MPI grouped supertile preserves rank-owned block order", {
  src_file <- locate_mpi_cdist_grouped_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- mpi_cdist_grouped_supertile_body(readLines(src_file, warn = FALSE))

  expect_equal(lengths(regmatches(
    body,
    gregexpr("np_conditional_y_eval_block_stream_op_core\\(", body, perl = TRUE)
  )), 1L)
  expect_equal(lengths(regmatches(
    body,
    gregexpr("np_blas_dgemm_tn_int\\(", body, perl = TRUE)
  )), 1L)
  expect_match(body, "fit_cross[ii + jj*ib]", fixed = TRUE)
  expect_match(body, "block_sum[g] += tvd*tvd;", fixed = TRUE)
  expect_match(body, "block_terms[block_id[g]] = block_sum[g];", fixed = TRUE)
  expect_match(body, "*cv += block_terms[ii];", fixed = TRUE)

  markers <- c(
    "np_conditional_y_eval_block_stream_op_core",
    "np_blas_dgemm_tn_int(ib",
    "block_sum[g] += tvd*tvd;",
    "block_terms[block_id[g]] = block_sum[g];",
    "*cv += block_terms[ii];"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))
})

test_that("MPI conditional-distribution grouped supertile remains repeatable", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess MPI contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages=FALSE, np.tree=FALSE, np.largeh=FALSE)",
      "npRmpi.init(nslaves=1L, quiet=TRUE)",
      "set.seed(20260727)",
      "n <- 1025L",
      "x <- data.frame(x1=runif(n), x2=runif(n))",
      "y <- data.frame(y=pnorm(sin(2*x$x1-x$x2)+rnorm(n,sd=0.2)))",
      "bw <- npcdistbw(xdat=x, ydat=y, bws=c(0.28,0.24,0.27),",
      "  bandwidth.compute=FALSE, bwmethod='cv.ls', bwtype='fixed',",
      "  regtype='lp', degree=c(2L,2L), bernstein.basis=TRUE)",
      "first <- npRmpi:::.npcdistbw_eval_only(xdat=x,ydat=y,bws=bw,do.full.integral=TRUE)",
      "second <- npRmpi:::.npcdistbw_eval_only(xdat=x,ydat=y,bws=bw,do.full.integral=TRUE)",
      "stopifnot(identical(writeBin(first$objective,raw(),size=8L),",
      "                    writeBin(second$objective,raw(),size=8L)),",
      "          identical(first$num.feval,1L),",
      "          identical(second$num.feval,1L))",
      "cat('MPI_CDIST_GROUPED_SUPERTILE_REPEATABLE_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_true(res$status %in% c(0L, 137L), info = info)
  expect_true(
    any(grepl("MPI_CDIST_GROUPED_SUPERTILE_REPEATABLE_OK",
              res$output, fixed = TRUE)),
    info = info
  )
})
