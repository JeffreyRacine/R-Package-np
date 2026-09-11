mode_empty_namespace <- function() asNamespace("npRmpi")

mode_empty_level_owner <- function(native) {
  ns <- mode_empty_namespace()
  env <- new.env(parent = ns)
  owner <- get(".npConmodeEvaluateLevels", ns)
  environment(owner) <- env
  env$npcdens <- native
  env$.npConmodeLevelBlockWidth <- function(...) 2L
  owner
}

test_that("mode level blocks preserve purpose and fold typed notices by query row", {
  calls <- list()
  flagged <- FALSE
  native <- function(txdat, tydat, exdat, eydat, bws, gradients, se, ...) {
    dots <- list(...)
    calls[[length(calls) + 1L]] <<- dots
    out <- list(condens = 10*as.integer(eydat) + exdat$x,
                conderr = if (se) rep(.1, nrow(exdat)) else NULL,
                congrad = if (gradients) matrix(exdat$x, ncol = 1L) else NULL)
    if (flagged) attr(out, ".np.empty.rows") <-
      as.integer(exdat$x == 2 & as.integer(eydat) == 2L)
    out
  }
  owner <- mode_empty_level_owner(native)
  x <- data.frame(x = c(1, 2, 1))
  lev <- factor(c("a", "b", "c"))
  args <- list(bws = list(xndim = 1L, xnames = "x"), txdat = x,
    tydat = data.frame(y = lev), xeval = x, efac = lev,
    gradients = TRUE, gradient.level.index = 2L, se = TRUE)
  healthy <- do.call(owner, args)
  expect_identical(names(healthy), c("probabilities", "errors", "gradients"))
  expect_null(attr(healthy, ".np.empty.rows"))
  expect_identical(unname(healthy$probabilities), outer(x$x, c(10, 20, 30), "+"))
  expect_true(all(vapply(calls, function(a) isTRUE(a[[".np.require.complete"]]), logical(1))))
  expect_true(all(vapply(calls, function(a) isTRUE(a[[".np.defer.empty.rows"]]), logical(1))))

  calls <- list()
  flagged <- TRUE
  mixed <- do.call(owner, c(args, list(allow.external = TRUE)))
  expect_identical(attr(mixed, ".np.empty.rows"), c(0L, 1L, 0L))
  expect_identical(mixed$probabilities, healthy$probabilities)
  expect_identical(mixed$errors, healthy$errors)
  expect_identical(mixed$gradients, healthy$gradients)
  expect_true(all(vapply(calls, function(a) identical(a[[".np.require.complete"]], FALSE), logical(1))))
  expect_true(all(vapply(calls, function(a) isTRUE(a[[".np.defer.empty.rows"]]), logical(1))))
})

test_that("mode does not reinterpret unrelated nonfinite values as typed support loss", {
  native <- function(exdat, ...) list(condens = rep(NA_real_, nrow(exdat)))
  owner <- mode_empty_level_owner(native)
  x <- data.frame(x = c(1, 2))
  args <- list(bws = list(xndim = 1L, xnames = "x"), txdat = x,
    tydat = data.frame(y = factor(c("a", "b"))), xeval = x,
    efac = factor(c("a", "b")), gradients = FALSE,
    gradient.level.index = 1L, se = FALSE)
  out <- do.call(owner, args)
  expect_true(all(is.na(out$probabilities)))
  expect_null(attr(out, ".np.empty.rows"))
  bad <- mode_empty_level_owner(function(exdat, ...) {
    structure(list(condens = rep(1, nrow(exdat))),
              .np.empty.rows = rep(1, nrow(exdat)))
  })
  expect_error(do.call(bad, args), "empty-row metadata is invalid", fixed = TRUE)
  bad <- mode_empty_level_owner(function(exdat, ...) {
    structure(list(condens = rep(1, nrow(exdat))), .np.empty.rows = 1L)
  })
  expect_error(do.call(bad, args), "empty-row metadata is invalid", fixed = TRUE)
})

test_that("mode public owner keeps purpose and one final notice explicit", {
  owner <- get("npconmode.conbandwidth", mode_empty_namespace())
  code <- paste(deparse(body(owner)), collapse = "\n")
  expect_equal(sum(gregexpr("allow.external = !no.ex", code, fixed = TRUE)[[1L]] > 0L), 2L)
  expect_equal(sum(gregexpr(".npreg_finish_empty_rows(", code, fixed = TRUE)[[1L]] > 0L), 1L)
  expect_match(code, 'owner = "npconmode"', fixed = TRUE)
  expect_match(code, "row.labels = rownames(xeval)", fixed = TRUE)
  expect_match(code, "attr(endpoint.fit, \".np.empty.rows\", exact = TRUE)", fixed = TRUE)
})
