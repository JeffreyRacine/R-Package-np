test_that("CMS call dots are values while subset and genuine call formals stay syntax", {
  method <- function(formula, xdat, ydat, data, subset, ...) {
    .npRmpi_autodispatch_bind_data_promises(match.call(), environment(),
      sys.call(), sys.function(), parent.frame())
  }
  environment(method) <- asNamespace("npRmpi")
  marker <- "master-local-call"
  actual <- method(y ~ x, data = data.frame(x = 1:3, y = 3:1),
                   subset = stop("subset must stay lazy"), call = marker)
  expect_identical(actual$call, marker)
  expect_identical(actual$subset, quote(stop("subset must stay lazy")))
  formal <- function(formula, xdat, ydat, data, subset, call, ...) {
    .npRmpi_autodispatch_bind_data_promises(match.call(), environment(),
      sys.call(), sys.function(), parent.frame())
  }
  environment(formal) <- asNamespace("npRmpi")
  actual <- formal(y ~ x, call = stop("call formal must stay lazy"))
  expect_identical(actual$call, quote(stop("call formal must stay lazy")))
})

test_that("pooled CMS formula call promises match rank-local statistical results", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(np.messages = FALSE)
  set.seed(3861)
  d <- data.frame(x = runif(50)); d$y <- sin(d$x) + rnorm(50, sd = .2)
  models <- list(npcmstest = lm(y ~ x, d, x = TRUE, y = TRUE),
                 npqcmstest = quantreg::rq(y ~ x, data = d, tau = .5, model = TRUE))
  for (family in names(models)) for (searched in c(FALSE, TRUE)) {
    model <- models[[family]]
    controls <- if (searched) list(nmulti = 1L, itmax = 5L) else
      list(bws = .3, bandwidth.compute = FALSE)
    c_m <- paste("master-only", family, searched)
    args <- c(list(formula = y ~ x, data = d, model = model,
                    B = 9L, random.seed = 19L), controls)
    expected <- .npRmpi_with_local_regression(do.call(get(family),
      c(args, list(call = c_m))))
    actual <- eval(as.call(c(list(as.name(family)), args, list(call = quote(c_m)))))
    expect_equal(actual[c("Jn", "In", "P")], expected[c("Jn", "In", "P")],
                 tolerance = 1e-12, info = paste(family, searched))
  }
})
