test_that("partially linear aliases retain independent values and original owners", {
  ns <- asNamespace("npRmpi")
  restore <- get(".npRmpi_bw_restore_public_owner", ns)
  project <- get(".npRmpi_bw_metadata_wire", ns)
  without.alias <- function(result, originals) {
    alias <- result[["bw", exact = TRUE]]
    result[["bw"]] <- NULL
    out <- restore(result, originals)
    if (!is.null(alias)) out[["bw"]] <- alias
    out
  }
  owner <- new.env(parent = baseenv())
  d <- data.frame(x = 1:12, y = sin(1:12))
  f <- y ~ x
  environment(f) <- owner
  frame <- model.frame(f, d)
  subset <- list(expression = quote(x < 10), envir = owner)
  attr(frame, ".np.subset") <- subset
  original <- structure(list(
    bw = .4, call = quote(npplregbw(y ~ x, data = d, subset = x < 10)),
    formula = f, terms = terms(frame),
    .np.formula.training = list(frame = frame, na.action = stats::na.omit,
                               subset = subset)), class = "plbandwidth")
  environment(original$call) <- owner
  worker <- project(original)
  worker$.np.formula.training$subset$envir <- new.env(parent = emptyenv())
  alias <- worker
  alias$bw <- .7
  alias$.np.formula.training$frame <- alias$.np.formula.training$frame[1:5, , drop = FALSE]
  fit <- structure(list(bws = worker, bw = alias, mean = 1:12), class = "plregression")
  before <- serialize(fit, NULL)
  withr::local_seed(987)
  seed <- .Random.seed
  out <- restore(fit, list(original))
  expect_identical(serialize(fit, NULL), before)
  expect_identical(.Random.seed, seed)
  expect_identical(out$bws$bw, .4)
  expect_identical(out$bw$bw, .7)
  expect_identical(nrow(out$bws$.np.formula.training$frame), 12L)
  expect_identical(nrow(out$bw$.np.formula.training$frame), 5L)
  for (name in c("bws", "bw")) {
    expect_identical(environment(out[[name]]$formula), owner, info = name)
    expect_identical(environment(out[[name]]$terms), owner, info = name)
    expect_identical(out[[name]]$.np.formula.training$subset$envir, owner, info = name)
  }
  for (kind in c("absent", "numeric", "wrong_result_class", "wrong_alias_class", "no_original")) {
    value <- fit
    refs <- list(original)
    if (kind == "absent") value$bw <- NULL
    if (kind == "numeric") value$bw <- .5
    if (kind == "wrong_result_class") class(value) <- "density"
    if (kind == "wrong_alias_class") class(value$bw) <- "bandwidth"
    if (kind == "no_original") refs <- list()
    expect_identical(restore(value, refs)[["bw"]],
                     without.alias(value, refs)[["bw"]], info = kind)
  }
  for (field in c("bws", "bw")) {
    bad <- fit
    bad[[field]]$call <- quote(changed_call())
    expect_error(restore(bad, list(original)), "metadata changed", info = field)
  }
  for (kind in c("wrong_original", "formula_free", "missing_bws")) {
    value <- fit
    ref <- original
    if (kind == "wrong_original") class(ref) <- "rbandwidth"
    if (kind == "formula_free") ref$formula <- NULL
    if (kind == "missing_bws") value$bws <- NULL
    expect_identical(restore(value, list(ref))[["bw"]],
                     without.alias(value, list(ref))[["bw"]], info = kind)
  }
  for (literal in c(FALSE, TRUE)) {
    ref <- original
    if (literal) ref$call[["data"]] <- d
    wire <- project(ref)
    value <- structure(list(bws = wire, bw = wire), class = "plregression")
    out <- restore(value, list(ref))
    expect_identical(environment(out$bw$formula), owner)
    expect_identical(out$bw$.np.formula.training$subset$envir, owner)
  }
})

test_that("nested bandwidths validate their own metadata before owner restoration", {
  ns <- asNamespace("npRmpi")
  restore <- get(".npRmpi_bw_restore_public_owner", ns)
  project <- get(".npRmpi_bw_metadata_wire", ns)
  owner <- new.env(parent = baseenv())
  original <- structure(list(bw = .4, call = quote(nplsqregbw(y ~ x)),
                             formula = y ~ x), class = "lsqregressionbandwidth")
  environment(original$formula) <- owner
  children <- lapply(c(.25, .75), function(tau) {
    child <- original
    child$call[["tau"]] <- tau
    child$bw <- tau
    child
  })
  original$tau.bws <- children
  wire <- project(original)
  out <- restore(list(bws = wire), list(original))
  for (i in seq_along(children)) {
    expect_identical(out$bws$tau.bws[[i]], children[[i]])
  }
  bad <- wire
  bad$tau.bws[[2L]]$call[["tau"]] <- .9
  expect_error(restore(list(bws = bad), list(original)), "metadata changed")
  expect_identical(wire$tau.bws[[2L]]$call[["tau"]], .75)
})
