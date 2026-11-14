test_that("significance selector controls use the original promises once", {
  pkg <- environmentName(environment(npsigtest))
  if (pkg == "npRmpi") {
    if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(-1, 1, length.out = 30L), y = sin(seq_len(30L)/4))
  values <- function(z) z[c("In", "In.bootstrap", "P", "ixvar", "boot.num",
    "pivot.effective", "bootstrap.executed", "rejectNum")]
  for (route in c("formula", "named", "native")) {
    hits <- c(nmulti = 0L, kernel = 0L)
    starts <- function() {hits["nmulti"] <<- hits["nmulti"] + 1L; 1L}
    kernel <- function() {hits["kernel"] <<- hits["kernel"] + 1L; "gaussian"}
    set.seed(754)
    a <- switch(route,
      formula = npsigtest(y ~ x, data = d, nmulti = starts(), ckertype = kernel(), B = 9L, pivot = FALSE),
      named = npsigtest(bws = y ~ x, data = d, nmulti = starts(), ckertype = kernel(), B = 9L, pivot = FALSE),
      native = npsigtest(xdat = d["x"], ydat = d$y, nmulti = starts(), ckertype = kernel(), B = 9L, pivot = FALSE))
    rng <- .Random.seed
    expect_identical(unname(hits), c(1L, 1L), info = route)
    set.seed(754)
    b <- npsigtest(xdat = d["x"], ydat = d$y, nmulti = 1L, ckertype = "gaussian", B = 9L, pivot = FALSE)
    expect_identical(.Random.seed, rng, info = route)
    expect_equal(values(a), values(b), tolerance = 0, info = route)
    expect_identical(a$bws$bw, b$bws$bw, info = route)
    expect_false(grepl("evalq|<environment|selector\\$", paste(deparse(a$bws$call), collapse = " ")), info = route)
    if (route != "native") {
      # Automatic MPI selection already omits this stored selector call.
      expect_null(a$bws$call)
    }
  }
})

test_that("selector promise binding is lazy and restores partial and null syntax", {
  pkg <- environmentName(environment(npsigtest))
  bind <- getFromNamespace(".np_npsig_selector_promises", pkg)
  restore <- getFromNamespace(".np_npsig_restore_selector_call", pkg)
  count <- 0L
  value <- function() { count <<- count + 1L; NULL }
  owner <- function(bws, xdat, ydat, ...) {
    call <- match.call(expand.dots = TRUE)
    handoff <- bind(call, environment(), parent.frame())
    expect_identical(count, 0L)
    expect_null(eval(handoff$call$degree))
    expect_null(...elt(1L))
    expect_identical(count, 1L)
    expect_identical(restore(handoff$call, handoff$replacements), call)
    handoff
  }
  x <- owner(bws = .4, degree = value(), nmul = 1L, B = 9L, subset = x > 0)
  expect_identical(x$call$B, 9L)
  expect_identical(x$call$subset, quote(x > 0))
})


test_that("selector restoration preserves pairlists and forwarded dot syntax", {
  pkg <- environmentName(environment(npsigtest))
  bind <- getFromNamespace(".np_npsig_selector_promises", pkg)
  restore <- getFromNamespace(".np_npsig_restore_selector_call", pkg)
  owner <- function(bws, xdat, ydat, ...) {
    call <- match.call(expand.dots = TRUE)
    handoff <- bind(call, environment(), parent.frame())
    expect_identical(restore(handoff$call, handoff$replacements), call)
    nested <- as.call(list(quote(npregbw), ... = as.pairlist(as.list(handoff$call)[-1L])))
    expected <- as.call(list(quote(npregbw), ... = as.pairlist(as.list(call)[-1L])))
    expect_identical(restore(nested, handoff$replacements), expected)
    expect_true(is.pairlist(restore(nested, handoff$replacements)[["..."]]))
  }
  outer <- function(...) owner(bws = .4, ...)
  outer(degree = NULL, nmulti = 1L)
  owner(bws = .4, degree = NULL, nmulti = 1L)
})


test_that("default handoff preserves absent and NULL call fields", {
  pkg <- environmentName(environment(npsigtest))
  if (pkg == "npRmpi") {
    if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  method <- getS3method("npsigtest", "default")
  selected <- NULL
  local_mocked_bindings(
    .np_eval_bw_call = function(...) selected,
    npsigtest = function(bws, ...) list(bws = bws),
    .package = pkg
  )
  for (explicit.null in c(FALSE, TRUE)) {
    selected <- structure(list(marker = "bw"), class = "rbandwidth")
    if (explicit.null) selected["call"] <- list(NULL)
    result <- method(bws = .4, xdat = data.frame(x = 1:4),
                     ydat = 1:4, nmulti = 1L)
    expect_identical(result$bws, selected)
    expect_identical("call" %in% names(result$bws), explicit.null)
  }
})
