npRmpi_unconditional_gnn_literal_lines <- function() {
  c(
    "gnn_literal_radius <- function(train, evaluation, k, exclude=integer()) {",
    "  keep <- if (length(exclude)) setdiff(seq_along(train), exclude) else seq_along(train)",
    "  sort(abs(evaluation-train[keep]), method='radix')[[k]]",
    "}",
    "gnn_literal_density <- function(train, evaluation, k, mapped=FALSE, leave_one_out=FALSE) {",
    "  vapply(seq_along(evaluation), function(i) {",
    "    exclude <- if (mapped || leave_one_out) i else integer()",
    "    h <- gnn_literal_radius(train, evaluation[[i]], k, exclude)",
    "    keep <- if (leave_one_out) setdiff(seq_along(train), i) else seq_along(train)",
    "    mean(dnorm((evaluation[[i]]-train[keep])/h)/h)",
    "  }, numeric(1L))",
    "}",
    "gnn_literal_distribution <- function(train, evaluation, k, mapped=FALSE) {",
    "  vapply(seq_along(evaluation), function(i) {",
    "    exclude <- if (mapped) i else integer()",
    "    h <- gnn_literal_radius(train, evaluation[[i]], k, exclude)",
    "    mean(pnorm((evaluation[[i]]-train)/h))",
    "  }, numeric(1L))",
    "}"
  )
}

npRmpi_run_unconditional_gnn_geometry_case <- function(case) {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  if (is.null(env))
    return(NULL)

  common <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1L, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    npRmpi_unconditional_gnn_literal_lines()
  )

  body <- switch(
    case,
    owners = c(
      "x <- c(-1.4,-0.6,-0.1,0.2,0.9,1.7)",
      "dat <- data.frame(x=x)",
      "density_bw <- npudensbw(dat=dat, bwtype='generalized_nn', bwmethod='cv.ml', bws=2L, bandwidth.compute=FALSE)",
      "distribution_bw <- npudistbw(dat=dat, bwtype='generalized_nn', bwmethod='cv.cdf', bws=2L, bandwidth.compute=FALSE)",
      "expected_density <- gnn_literal_density(x,x,2L,mapped=TRUE)",
      "expected_distribution <- gnn_literal_distribution(x,x,2L,mapped=TRUE)",
      "for (tree in c(FALSE,TRUE)) {",
      "  options(np.tree=tree)",
      "  stopifnot(isTRUE(all.equal(npudens(bws=density_bw, tdat=dat)$dens, expected_density, tolerance=2e-10)))",
      "  stopifnot(isTRUE(all.equal(npudist(bws=distribution_bw, tdat=dat)$dist, expected_distribution, tolerance=2e-10)))",
      "}",
      "options(np.tree=FALSE)",
      "external <- dat[c(2L,4L),,drop=FALSE]",
      "stopifnot(isTRUE(all.equal(npudens(bws=density_bw, tdat=dat, edat=external)$dens, gnn_literal_density(x,external$x,2L), tolerance=2e-10)))",
      "eval_only <- npRmpi:::npudensbw.bandwidth(dat=dat, bws=density_bw, bandwidth.compute=TRUE, eval.only=TRUE, nmulti=1L)",
      "expected_cvml <- sum(log(gnn_literal_density(x,x,2L,mapped=TRUE,leave_one_out=TRUE)))",
      "stopifnot(isTRUE(all.equal(eval_only$fval[[1L]], expected_cvml, tolerance=2e-10)))",
      "cat('NPRMPI_UNCONDITIONAL_GNN_OWNERS_OK\\n')"
    ),
    zero_radius = c(
      "dat <- data.frame(x=c(0,0,0,1,2))",
      "bw <- npudensbw(dat=dat, bwtype='generalized_nn', bwmethod='cv.ml', bws=2L, bandwidth.compute=FALSE)",
      "density_error <- tryCatch({npudens(bws=bw, tdat=dat); ''}, error=conditionMessage)",
      "stopifnot(grepl('zero literal radius after occurrence exclusion', density_error, fixed=TRUE))",
      "dist_bw <- npudistbw(dat=dat, bwtype='generalized_nn', bwmethod='cv.cdf', bws=2L, bandwidth.compute=FALSE)",
      "distribution_error <- tryCatch({npudist(bws=dist_bw, tdat=dat); ''}, error=conditionMessage)",
      "stopifnot(grepl('zero literal radius after occurrence exclusion', distribution_error, fixed=TRUE))",
      "invalid <- npRmpi:::npudensbw.bandwidth(dat=dat, bws=bw, bandwidth.compute=TRUE, eval.only=TRUE, nmulti=1L)",
      "stopifnot(sum(invalid$invalid.history) > 0)",
      "cat('NPRMPI_UNCONDITIONAL_GNN_ZERO_RADIUS_OK\\n')"
    ),
    stop("unknown case")
  )

  npRmpi_run_rscript_subprocess(
    lines = c(common, body), timeout = 90L, env = env
  )
}

test_that("ordinary generalized-NN unconditional owners delete mapped occurrences", {
  result <- npRmpi_run_unconditional_gnn_geometry_case("owners")
  skip_if(is.null(result), "installed npRmpi unavailable for subprocess proof")
  info <- paste(result$output, collapse = "\n")
  expect_equal(result$status, 0L, info = info)
  expect_true(any(grepl("NPRMPI_UNCONDITIONAL_GNN_OWNERS_OK",
                        result$output, fixed = TRUE)), info = info)
})

test_that("ordinary generalized-NN unconditional fit reports zero literal radii", {
  result <- npRmpi_run_unconditional_gnn_geometry_case("zero_radius")
  skip_if(is.null(result), "installed npRmpi unavailable for subprocess proof")
  info <- paste(result$output, collapse = "\n")
  expect_equal(result$status, 0L, info = info)
  expect_true(any(grepl("NPRMPI_UNCONDITIONAL_GNN_ZERO_RADIUS_OK",
                        result$output, fixed = TRUE)), info = info)
})
