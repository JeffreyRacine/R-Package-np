n1_selected_fixture <- function(regtype = "lc", bwtype = "fixed") {
  set.seed(120929)
  n <- 36L
  x <- data.frame(o = ordered(rep(1:4, 9L)), x = runif(n, .05, .95),
    u = factor(rep(letters[1:3], 12L)), v = factor(rep(c("a", "b"), 18L)))
  y <- sin(x$x * 3) + .3 * (x$u == "b") + rnorm(n, sd = .25)
  response <- vapply(seq_len(8L), function(k) y + rnorm(n, sd = .2), numeric(n))
  bw <- npregbw(xdat = x, ydat = y,
    bws = c(.3, if (bwtype == "fixed") .45 else 18, .25, .2),
    bandwidth.compute = FALSE, regtype = regtype, bwtype = bwtype,
    degree = if (regtype == "lp") 2L else NULL)
  list(x = x, y = y, response = response, bw = bw)
}

test_that("selected-coordinate attribute is protected across allocating native calls", {
  source <- test_path("..", "..", "src", "np.c")
  if (!file.exists(source)) skip("native source is not available in installed context")
  lines <- readLines(source, warn = FALSE)
  first <- grep("^SEXP C_np_regression\\(", lines)
  last <- grep("^SEXP C_np_density\\(", lines)
  expect_length(first, 1L)
  expect_length(last, 1L)
  wrapper <- lines[seq.int(first, last - 1L)]
  protected <- grep('SEXP gradient_coordinate = PROTECT(getAttrib(output_request, install(".np.gradient.coordinate")));',
                    wrapper, fixed = TRUE)
  expect_length(protected, 1L)
  expect_lt(protected, min(grep("coerceVector(", wrapper, fixed = TRUE)))
  # Count actual PROTECT tokens, including inline assignments but excluding
  # UNPROTECT. There are 27 unconditional protections and two guarded outputs.
  expect_length(grep("\\bPROTECT\\(", wrapper, perl = TRUE), 29L)
  empty.start <- grep("PROTECT(empty_flags = allocVector(INTSXP, en));",
                      wrapper, fixed = TRUE)
  certificate.start <- grep("SEXP certificate = PROTECT(allocVector(LGLSXP, gsize));",
                            wrapper, fixed = TRUE)
  expect_length(empty.start, 1L)
  expect_length(certificate.start, 1L)
  expect_identical(trimws(wrapper[empty.start + 1L]), "++extra_protect;")
  expect_identical(trimws(wrapper[certificate.start + 1L]), "++extra_protect;")
  expect_match(paste(wrapper[seq.int(empty.start - 2L, empty.start)], collapse = " "),
               "INTEGER(output_request)[1] == 1 && !train_is_eval)", fixed = TRUE)
  expect_identical(trimws(wrapper[certificate.start - 1L]),
                   "if(gradient_zero_out != NULL) {")
  expect_lt(empty.start, grep("empty_rows.flags = INTEGER(empty_flags);",
                              wrapper, fixed = TRUE))
  expect_lt(certificate.start,
    grep('setAttrib(out, install(".np.gradient.structural.zero"), certificate);',
         wrapper, fixed = TRUE))
  expect_length(grep("++extra_protect;", wrapper, fixed = TRUE), 2L)
  expect_length(grep("UNPROTECT(27 + extra_protect);", wrapper, fixed = TRUE), 1L)
  expect_length(grep("return out;", wrapper, fixed = TRUE), 1L)
})

test_that("private selected coordinates retain canonical full layout and arithmetic", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  direct.owner <- getFromNamespace(".np_regression_direct", "npRmpi")
  local.regression <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  direct <- function(...) local.regression(direct.owner(...))
  complete <- getFromNamespace(".npRmpi_npsig_npreg_local", "npRmpi")
  statistic <- getFromNamespace(".np_npsig_statistic", "npRmpi")
  structure <- getFromNamespace(".np_npsig_structure", "npRmpi")
  for (type in c("lc", "ll", "lp")) {
    z <- n1_selected_fixture(type)
    seed <- .Random.seed
    full <- complete(z$bw, txdat = z$x, tydat = z$y, gradients = TRUE, se = TRUE)
    plain <- direct(z$bw, z$x, z$y, gradients = TRUE, se = TRUE)
    for (field in c("mean", "merr", "grad", "gerr"))
      expect_identical(plain[[field]], full[[field]])
    for (coordinate in c(1L, 3L, 4L)) {
      selected <- direct(z$bw, z$x, z$y, gradients = TRUE, se = TRUE,
        gradient.coordinate = coordinate)
      expect_identical(selected$mean, full$mean)
      expect_identical(selected$merr, full$merr)
      for (field in c("grad", "gerr")) {
        expect_identical(dim(selected[[field]]), dim(full[[field]]))
        expect_identical(selected[[field]][, coordinate], full[[field]][, coordinate])
        expect_true(all(is.na(selected[[field]][, -coordinate, drop = FALSE])))
      }
      mask <- structure(z$bw, z$x, coordinate)
      expect_identical(statistic(selected, coordinate, TRUE, mask),
                       statistic(full, coordinate, TRUE, mask))
    }
    expect_identical(.Random.seed, seed)
  }
})

test_that("selected tiles preserve responses donors and the complete joint cache", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  z <- n1_selected_fixture("ll")
  tile <- getFromNamespace(".np_npsig_streamed_iid_tile", "npRmpi")
  joint <- getFromNamespace(".np_npsig_streamed_response_statistic", "npRmpi")
  complete <- getFromNamespace(".npRmpi_npsig_npreg_local", "npRmpi")
  direct.owner <- getFromNamespace(".np_regression_direct", "npRmpi")
  local.regression <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  direct <- function(...) local.regression(direct.owner(...))
  statistic <- getFromNamespace(".np_npsig_statistic", "npRmpi")
  mask <- getFromNamespace(".np_npsig_structure", "npRmpi")(z$bw, z$x, 3L)
  expected <- vapply(seq_len(8L), function(k) statistic(
    complete(z$bw, txdat = z$x, tydat = z$response[, k], gradients = TRUE, se = TRUE),
    3L, TRUE, mask), numeric(1L))
  count <- new.env(parent = emptyenv())
  count$full <- count$direct <- count$progress <- 0L
  local_mocked_bindings(
    .npRmpi_npsig_npreg_local = function(...) {
      count$full <- count$full + 1L
      complete(...)
    },
    .np_regression_direct = function(..., gradient.coordinate = NULL) {
      count$direct <- count$direct + 1L
      expect_identical(gradient.coordinate, 3L)
      direct(..., gradient.coordinate = gradient.coordinate)
    },
    .np_with_compiled_fit_progress = function(label, total, handoff = FALSE,
                                             handoff.detail = NULL, expr) {
      count$progress <- count$progress + 1L
      expect_identical(label, "Fitting regression")
      expect_identical(total, as.integer(2 * nrow(z$x)))
      force(expr)
    }, .package = "npRmpi")
  for (width in c(1L, 7L, 8L, 1L)) {
    count$full <- count$direct <- count$progress <- 0L
    seed <- .Random.seed
    got <- tile(z$bw, z$x, 3L, response.matrix = z$response[, seq_len(width), drop = FALSE],
      null.mean = z$y, residual.pool = z$y, pivotal = TRUE)
    expect_identical(got, expected[seq_len(width)])
    expect_identical(count$full, 0L)
    expect_identical(count$direct, width)
    expect_identical(count$progress, width)
    expect_identical(.Random.seed, seed)
  }
  count$full <- count$direct <- 0L
  joint(z$bw, z$x, c(1L, 3L, 4L), z$response, TRUE)
  expect_identical(count$full, 8L)
  expect_identical(count$direct, 0L)
  donors <- matrix(rep(seq_len(nrow(z$x)), 2L), nrow(z$x))
  expect_identical(
    tile(z$bw, z$x, 3L, donor.index = donors, null.mean = z$y,
      residual.pool = z$y, pivotal = TRUE),
    tile(z$bw, z$x, 3L, response.matrix = z$y + matrix(z$y[donors], nrow(z$x)),
      null.mean = z$y, residual.pool = z$y, pivotal = TRUE))
})

test_that("selected output validation and zero inference remain fail closed", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  z <- n1_selected_fixture()
  direct.owner <- getFromNamespace(".np_regression_direct", "npRmpi")
  local.regression <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  direct <- function(...) local.regression(direct.owner(...))
  statistic <- getFromNamespace(".np_npsig_statistic", "npRmpi")
  for (coordinate in list(0L, 2L, 5L, NA_integer_, Inf, 1.5, c(1L, 3L), TRUE, "1"))
    expect_error(direct(z$bw, z$x, z$y, gradients = TRUE, se = TRUE,
      gradient.coordinate = coordinate), "requires one categorical column")
  expect_error(direct(z$bw, z$x, z$y, gradients = TRUE,
    gradient.coordinate = 3L), "requires one categorical column")
  zero <- direct(z$bw, z$x, numeric(nrow(z$x)), gradients = TRUE, se = TRUE,
    gradient.coordinate = 3L)
  expect_identical(statistic(zero, 3L, TRUE), 0)
  zero$grad[1L, 3L] <- 1e-200
  expect_error(statistic(zero, 3L, TRUE), "zero standard error")
  zero$grad[1L, 3L] <- Inf
  expect_error(statistic(zero, 3L, TRUE), "non-finite gradient estimates")
})
