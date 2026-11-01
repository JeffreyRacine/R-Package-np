procedure_fixture <- function(method = 'cv.ml', manual = FALSE, scaled = FALSE) {
  list(bw = .4, method = method, pmethod = if (manual) 'Manual' else 'Cross-Validation',
    scaling = scaled, type = 'fixed', ckertype = 'gaussian', ckerorder = 2L,
    ckerbound = 'none', ckerlb = NULL, ckerub = NULL,
    ukertype = 'aitchisonaitken', okertype = 'liracine',
    bandwidth = list(x = .4), sfactor = list(x = 1.2),
    xnames = 'x', icon = TRUE, iuno = FALSE, iord = FALSE)
}

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
