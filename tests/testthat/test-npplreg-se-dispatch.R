test_that("partial-linear default dispatch pins the validated uncertainty request", {
  owner <- getS3method("npplreg", "default")
  resolve.arg <- get(".npRmpi_autodispatch_resolve_owned_arg",
                     envir = environment(owner), inherits = FALSE)
  private <- new.env(parent = environment(owner))
  environment(owner) <- private
  private$.npRmpi_require_active_slave_pool <- function(...) invisible(TRUE)
  private$.npRmpi_autodispatch_active <- function() TRUE
  private$.npRmpi_autodispatch_called_from_bcast <- function() FALSE
  dispatched <- 0L
  private$.npRmpi_autodispatch_call <- function(mc, caller_env) {
    dispatched <<- dispatched + 1L
    list(call = mc,
         request = resolve.arg(mc[["se"]], "se", caller_env, character(0)))
  }

  plain <- owner()
  expect_identical(plain$call[["se"]], FALSE)
  expect_identical(plain$request, FALSE)
  for (request in list(FALSE, TRUE, 0L, 1L)) {
    captured <- owner(se = request)
    expect_identical(captured$call[["se"]], as.logical(request))
    expect_identical(captured$request, as.logical(request))
  }

  evaluations <- 0L
  captured <- owner(se = {
    evaluations <- evaluations + 1L
    evaluations == 1L
  })
  expect_identical(evaluations, 1L)
  expect_identical(captured$request, TRUE)

  x <- 1:3
  y <- 2:4
  z <- 3:5
  captured <- owner(txdat = x, tydat = y, tzdat = z,
                    se = TRUE, residuals = TRUE, degree.select = "manual")
  expect_identical(captured$call[["txdat"]], quote(x))
  expect_identical(captured$call[["tydat"]], quote(y))
  expect_identical(captured$call[["tzdat"]], quote(z))
  expect_identical(captured$call[["residuals"]], TRUE)
  expect_identical(captured$call[["degree.select"]], "manual")
  before <- dispatched
  expect_error(owner(se = NA), "'se' must be TRUE or FALSE", fixed = TRUE)
  expect_identical(dispatched, before)
})
