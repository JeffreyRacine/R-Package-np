procedure_fixture <- function(method = 'cv.ml', manual = FALSE, scaled = FALSE) {
  list(bw = .4, method = method, pmethod = if (manual) 'Manual' else 'Cross-Validation',
    scaling = scaled, type = 'fixed', ckertype = 'gaussian', ckerorder = 2L,
    ckerbound = 'none', ckerlb = NULL, ckerub = NULL,
    ukertype = 'aitchisonaitken', okertype = 'liracine',
    bandwidth = list(x = .4), sfactor = list(x = 1.2),
    xnames = 'x', icon = TRUE, iuno = FALSE, iord = FALSE)
}

test_that('density default selection forwards explicit transform.bounds', {
  calls <- list()
  selector <- function(dat, bws, ...) {
    calls[[length(calls) + 1L]] <<- list(...)
    bws
  }
  # Adapter argument assembly only; real worker ownership is covered separately.
  local_mocked_bindings(npudensbw.bandwidth = selector,
    .npRmpi_require_active_slave_pool = function(...) invisible(NULL))
  d <- data.frame(x = seq(-1, 1, length.out = 8))
  npudensbw.default(dat = d, bws = .4, bandwidth.compute = FALSE)
  npudensbw.default(dat = d, bws = .4, bandwidth.compute = FALSE, transform.bounds = FALSE)
  npudensbw.default(dat = d, bws = .4, bandwidth.compute = FALSE, transform.bounds = TRUE)
  expect_null(calls[[1L]][['transform.bounds']])
  expect_identical(calls[[2L]][['transform.bounds']], FALSE)
  expect_identical(calls[[3L]][['transform.bounds']], TRUE)
})

test_that('complete replay arguments agree with retained invariant metadata', {
  replacements <- list(bwscaling = TRUE, bwtype = 'generalized_nn',
    ckertype = 'epanechnikov', ckerorder = 4L, ckerbound = 'range',
    ckerlb = -2, ckerub = 2, ukertype = 'liracine', okertype = 'wangvanryzin',
    bws = .8)
  for (family in c('density', 'distribution')) {
    b <- procedure_fixture(if (family == 'density') 'cv.ml' else 'cv.cdf', manual = TRUE)
    frame <- list2env(list(bandwidth.compute = FALSE, eval.only = FALSE), parent = emptyenv())
    valid <- .np_bw_procedure_record(b, b, family, frame)
    for (field in names(replacements)) {
      bad <- valid
      bad$selection.procedure$arguments[field] <- replacements[field]
      expect_error(.np_bw_outer_selection(bad, family), 'invalid.*arguments', info = field)
    }
    for (start in list(c(.4, .4), Inf, NA_real_, '0.4')) {
      bad <- valid; bad$selection.procedure$arguments$bws <- start
      expect_error(.np_bw_outer_selection(bad, family), 'invalid.*arguments')
    }
    bad <- valid; bad$selection.procedure$arguments$ckertype <- NULL
    expect_error(.np_bw_outer_selection(bad, family), 'invalid.*arguments')
    expect_identical(.np_bw_outer_selection(valid, family)$diagnostics$rule, 'manual')
    cv <- valid; cv$selection.procedure$rule <- 'cv'
    cv$selection.procedure$arguments$bandwidth.compute <- TRUE
    cv$selection.procedure$arguments$bws <- .25
    expect_identical(.np_bw_outer_selection(cv, family)$diagnostics$rule, 'cv')
  }
})

test_that('evaluation-only does not erase an explicit manual replay rule', {
  for (family in c('density', 'distribution')) {
    b <- procedure_fixture(if (family == 'density') 'cv.ml' else 'cv.cdf', manual = TRUE)
    frame <- list2env(list(bandwidth.compute = FALSE, eval.only = TRUE), parent = emptyenv())
    value <- .np_bw_procedure_record(b, b, family, frame)
    expect_true(value$selection.procedure$complete)
    expect_null(value$selection.procedure$reason)
    expect_identical(.np_bw_outer_selection(value, family)$diagnostics$rule, 'manual')
    frame$bandwidth.compute <- TRUE
    automatic <- .np_bw_procedure_record(b, b, family, frame)
    expect_false(automatic$selection.procedure$complete)
    expect_identical(automatic$selection.procedure$reason, 'evaluation-only object')
  }
})


test_that('legacy outer samples are passed as data without evaluating term labels', {
  calls <- list()
  estimator <- function(tdat, ...) {
    calls[[length(calls) + 1L]] <<- list(data = tdat, dots = list(...))
    tdat
  }
  local_mocked_bindings(npudens = estimator, npudist = estimator)
  outer <- setNames(data.frame(c(.1, .3, .6)), 'sqrt(missing_outer_source)')
  for (family in c('density', 'distribution')) {
    b <- procedure_fixture(if (family == 'density') 'cv.ml' else 'cv.cdf')
    b$xnames <- names(outer)
    before <- serialize(b, NULL)
    p <- .np_bw_outer_selection(b, family)
    expect_identical(p$diagnostics$status, 'legacy-unresolved')
    expect_length(p$messages, 1L)
    expect_identical(p$fit(outer), outer)
    expect_identical(calls[[length(calls)]]$data, outer)
    expect_length(calls[[length(calls)]]$dots, 0L)
    expect_identical(serialize(b, NULL), before)
  }
  expect_length(calls, 2L)
})

test_that('outer bandwidth provenance is value-only and never evaluates a call', {
  b <- procedure_fixture(manual = TRUE, scaled = TRUE)
  b$call <- quote(stop('must not be evaluated'))
  env <- list2env(list(bandwidth.compute = FALSE, eval.only = FALSE), parent = emptyenv())
  before <- serialize(b, NULL)
  v <- .np_bw_procedure_record(b, b, 'density', env)
  expect_identical(serialize(b, NULL), before)
  expect_true(v$selection.procedure$complete)
  expect_true(.np_bw_procedure_values(v$selection.procedure))
  expect_identical(v, unserialize(serialize(v, NULL)))
  p <- .np_bw_outer_selection(v, 'density')
  expect_identical(p$diagnostics$rule, 'manual')
  expect_identical(p$diagnostics$bandwidth.units, 'scale-factor')
  expect_true(p$diagnostics$original.rule.replayed)
  expect_length(p$messages, 0L)
  expect_false(.np_bw_procedure_values(list(x = new.env())))
  expect_false(.np_bw_procedure_values(list(x = quote(runif(1)))))
})

test_that('old objects retain explicit recovery and compatibility distinctions', {
  for (family in c('density', 'distribution')) {
    b <- procedure_fixture(if (family == 'density') 'cv.ml' else 'cv.cdf')
    set.seed(713); entry <- .Random.seed
    p <- .np_bw_outer_selection(b, family)
    expect_identical(.Random.seed, entry)
    expect_identical(p$diagnostics$status, 'legacy-unresolved')
    expect_false(p$diagnostics$original.rule.replayed)
    expect_length(p$messages, 1L)
    b$pmethod <- 'Manual'
    expect_identical(.np_bw_outer_selection(b, family)$diagnostics$status, 'legacy-rule-recovered')
    b$pmethod <- 'Normal Reference'; b$method <- 'normal-reference'
    expect_identical(.np_bw_outer_selection(b, family)$diagnostics$rule, 'normal-reference')
  }
})

test_that('stale, unsafe and malformed replay descriptors reject before fitting', {
  b <- procedure_fixture(manual = TRUE)
  env <- list2env(list(bandwidth.compute = FALSE, eval.only = FALSE), parent = emptyenv())
  v <- .np_bw_procedure_record(b, b, 'density', env)
  bad <- v; bad$bw <- .5
  expect_error(.np_bw_outer_selection(bad, 'density'), 'invalid or stale')
  bad <- v; bad$selection.procedure$arguments$data <- data.frame(x = 1)
  expect_error(.np_bw_outer_selection(bad, 'density'), 'invalid.*arguments')
  bad <- v; bad$selection.procedure$arguments$bws <- quote(runif(1))
  expect_error(.np_bw_outer_selection(bad, 'density'), 'invalid or stale')
  expect_error(.np_bw_outer_selection(v, 'distribution'), 'invalid or stale')
})

test_that('outer fitting preserves CV starts and avoids unnecessary public entries', {
  calls <- list()
  estimator <- function(tdat, ...) { calls[[length(calls) + 1L]] <<- c(list(owner = 'fit'), list(...)); tdat }
  selector <- function(dat, ...) { calls[[length(calls) + 1L]] <<- c(list(owner = 'select'), list(...)); 'selected' }
  local_mocked_bindings(npudens = estimator, npudist = estimator,
    npudensbw = selector, npudistbw = selector)
  for (family in c('density', 'distribution')) {
    for (case in c('manual', 'normal-reference', 'cv-zero', 'cv-start')) {
      b <- procedure_fixture(if (case == 'normal-reference') case else
        if (family == 'density') 'cv.ml' else 'cv.cdf', manual = case == 'manual')
      args <- .np_bw_procedure_kernel(b)
      args$bws <- if (case == 'cv-zero') 0 else .4
      args$bandwidth.compute <- case != 'manual'
      b$selection.procedure <- list(version = 1L, family = family,
        complete = TRUE, rule = if (grepl('^cv', case)) 'cv' else case,
        arguments = args, signature = .np_bw_procedure_signature(b))
      calls <- list(); data <- data.frame(x = 1:3)
      expect_identical(.np_bw_outer_selection(b, family)$fit(data), data)
      if (case == 'cv-start') {
        expect_identical(vapply(calls, `[[`, '', 'owner'), c('select', 'fit'))
        expect_identical(calls[[1L]][['bws', exact = TRUE]], .4)
        expect_true(calls[[1L]]$bandwidth.compute)
        expect_identical(calls[[2L]][['bws', exact = TRUE]], 'selected')
      } else {
        expect_length(calls, 1L)
        expect_identical(calls[[1L]]$owner, 'fit')
        expect_identical(calls[[1L]][['bws', exact = TRUE]], if (case == 'manual') .4 else NULL)
      }
    }
  }
})
