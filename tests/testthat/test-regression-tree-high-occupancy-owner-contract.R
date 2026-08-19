locate_high_occupancy_source <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  paths <- file.path(roots[nzchar(roots)], "src", "jksum.c")
  paths <- paths[file.exists(paths)]
  if (length(paths)) paths[[1L]] else NULL
}

test_that("MPI high-occupancy certificate adds no allocation or collective", {
  path <- locate_high_occupancy_source()
  skip_if(is.null(path), "package sources unavailable")
  source <- readLines(path, warn = FALSE)
  start <- grep(
    "^static int np_reg_fixed_tree_dense_high_occupancy_admitted\\($",
    source
  )
  stop <- grep("^/\\*$", source)
  stop <- stop[stop > start][1L]
  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_gt(stop, start)
  region <- paste(source[seq.int(start, stop - 1L)], collapse = "\n")
  all.source <- paste(source, collapse = "\n")

  expect_match(region, "if(terms[l] != 0)", fixed = TRUE)
  expect_match(region, "sorted = basis[0];", fixed = TRUE)
  expect_match(
    region,
    "sorted[i] = matrix_X_continuous[partial_dimension][i];",
    fixed = TRUE
  )
  expect_match(region, "sorted[i] = 1.0;", fixed = TRUE)
  expect_match(all.source, "F77_CALL(dlasrt)", fixed = TRUE)
  expect_match(region,
               "certificate_ok = np_reg_cv_support_sort_increasing(sorted, num_obs);",
               fixed = TRUE)
  expect_match(all.source, "single exit", fixed = TRUE)
  expect_false(grepl("malloc(", region, fixed = TRUE))
  expect_false(grepl("calloc(", region, fixed = TRUE))
  expect_false(grepl("realloc(", region, fixed = TRUE))
  expect_false(grepl("MPI_", region, fixed = TRUE))
  expect_false(grepl("R_CheckUserInterrupt", region, fixed = TRUE))
  expect_false(grepl("num_obs*(size_t)num_obs", region, fixed = TRUE))
  expect_false(grepl("qsort(sorted", all.source, fixed = TRUE))
  expect_match(all.source, "NP_REG_CV_LP_DENSE_MIN_TERMS = 10", fixed = TRUE)
  expect_match(all.source, "NP_REG_CV_LP_DENSE_OCCUPANCY_NUMERATOR = 9",
               fixed = TRUE)
  expect_match(all.source, "NP_REG_CV_LP_DENSE_OCCUPANCY_DENOMINATOR = 10",
               fixed = TRUE)
  expect_false(grepl("NP_REG_CV_SUPPORT_ORDER_", all.source, fixed = TRUE))
  expect_false(grepl("A10", all.source, fixed = TRUE))
})

test_that("MPI fixed CVLS honors the exact high-occupancy owner boundary", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  for (nslaves in c(1L, 3L)) {
    marker <- paste0("REGRESSION_HIGH_OCCUPANCY_OWNER_OK_", nslaves)
    result <- npRmpi_run_rscript_subprocess(
      lines = c(
        "suppressPackageStartupMessages(library(npRmpi))",
        sprintf("npRmpi.init(nslaves=%dL, quiet=TRUE)", nslaves),
        "options(np.messages=FALSE, np.tree=TRUE, np.macMseries.accelerate=FALSE)",
        "set.seed(20260815L)",
        "n <- 40L",
        "x <- data.frame(x1=sort(runif(n,-1,1)), x2=runif(n,-0.2,0.2))",
        "y <- 0.7 + 0.5*x$x1 - 0.3*x$x2 + 0.25*x$x1*x$x2 + 0.08*x$x1^2 + rnorm(n,sd=0.04)",
        "distances <- sort(outer(x$x1,x$x1,function(a,b) abs(a-b))[lower.tri(diag(n))])",
        "total <- n*(n-1L)/2L",
        "threshold <- ceiling(9*total/10)",
        "radius <- function(count) mean(distances[c(count,count+1L)])",
        "full2 <- diff(range(x$x2))*1.01",
        "core <- getFromNamespace('.npregbw_call_fixed_degree_core','npRmpi')",
        "evaluate <- function(degree, bandwidth, method='cv.ls', objective='ls', kernel='uniform', order=2L, bernstein=TRUE) {",
        "  args <- list(xdat=x,ydat=y,bws=bandwidth,bandwidth.compute=FALSE,bwtype='fixed',bwmethod=method,ckertype=kernel,regtype='lp',degree=degree,bernstein.basis=bernstein)",
        "  if (kernel == 'epanechnikov') args$ckerorder <- order",
        "  bw <- do.call(npregbw,args)",
        "  options(np.tree=TRUE); Sys.setenv(NP_LP_TREE_ORACLE='1')",
        "  messages <- capture.output(tree <- core(xdat=x,ydat=y,bws=bw,eval.only=TRUE,objective=objective),type='message')",
        "  Sys.unsetenv('NP_LP_TREE_ORACLE'); options(np.tree=FALSE)",
        "  dense <- core(xdat=x,ydat=y,bws=bw,eval.only=TRUE,objective=objective)",
        "  stopifnot(isTRUE(all.equal(tree$objective,dense$objective,tolerance=2e-11)))",
        "  stopifnot(identical(tree$num.feval,dense$num.feval))",
        "  any(grepl('NP_LP_TREE_ORACLE',messages,fixed=TRUE))",
        "}",
        "stopifnot(evaluate(c(3L,3L),c(radius(threshold-1L),full2)))",
        "stopifnot(!evaluate(c(3L,3L),c(radius(threshold),full2)))",
        "stopifnot(!evaluate(c(3L,3L),c(radius(threshold+1L),full2)))",
        "stopifnot(evaluate(c(3L,2L),c(radius(threshold),full2)))",
        "stopifnot(evaluate(c(3L,3L),c(radius(threshold),full2*0.6)))",
        "stopifnot(!evaluate(c(3L,3L),c(radius(threshold),full2)/sqrt(5),kernel='epanechnikov'))",
        "stopifnot(!evaluate(c(3L,3L),c(radius(threshold),full2),bernstein=FALSE))",
        "stopifnot(evaluate(c(3L,3L),c(radius(threshold),full2),method='cv.aic'))",
        "stopifnot(evaluate(c(3L,3L),c(radius(threshold),full2),objective='ks'))",
        "npRmpi.quit(force=TRUE)",
        sprintf("cat('%s\\n')", marker)
      ),
      timeout = 180L,
      env = env
    )
    expect_equal(result$status, 0L,
                 info = paste(result$output, collapse = "\n"))
    expect_true(any(grepl(marker, result$output, fixed = TRUE)),
                info = paste(result$output, collapse = "\n"))
  }
})
