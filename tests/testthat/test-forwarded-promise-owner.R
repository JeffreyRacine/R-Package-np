test_that("only an exact activation's formals own matched forwarded arguments", {
  resolve <- getFromNamespace(".npRmpi_autodispatch_resolve_owned_arg", "npRmpi")
  dots.only <- function(...) {
    ydat <- rep(11, 3)
    resolve(quote(..1), "ydat", environment(), "")
  }
  formal <- function(ydat, ...) {
    # ..1 has been consumed by S3 matching; the original formal now owns it.
    resolve(quote(..1), "ydat", environment(), "")
  }
  environment(formal) <- asNamespace("npRmpi")
  body(formal) <- quote(.npRmpi_autodispatch_resolve_owned_arg(
    quote(..1), "ydat", environment(), ""))
  named <- function(...) {
    ydat <- 11
    resolve(quote(..1), "ydat", environment(), "ydat")
  }
  expect_identical(dots.only(1:3), 1:3)
  expect_identical(formal(1:3, 88), 1:3)
  expect_identical(named(ydat = 1:3), 1:3)
  unrelated.formal <- function(ydat, ...) {
    resolve(quote(..1), "ydat", environment(), "")
  }
  expect_identical(unrelated.formal(11, 1:3), 1:3)
  explicit <- function(...) {
    resolve(quote(FALSE), "se", environment(), "se")
  }
  expect_identical(explicit(se = TRUE), FALSE)
  n <- 0L
  expect_identical(dots.only({ n <- n + 1L; 1:3 }), 1:3)
  expect_identical(n, 1L)
  expect_error(dots.only({ n <- n + 1L; stop("original promise") }), "original promise")
  expect_identical(n, 2L)
  targets <- getFromNamespace(".npRmpi_autodispatch_target_args", "npRmpi")()
  sweep <- function(name, ...) {
    assign(name, "unrelated local", envir = environment())
    resolve(quote(..1), name, environment(), "")
  }
  # name is a control formal in this test, not a materialized package argument.
  for (target in setdiff(targets, "name"))
    expect_identical(sweep(target, 1:3), 1:3, info = target)
})

test_that("native dispatch consumes unchanged actual method data promises once", {
  bind <- getFromNamespace(".npRmpi_autodispatch_bind_data_promises", "npRmpi")
  method <- function(xdat, ydat, unused = stop("unused"), ...) {
    xdat <- as.data.frame(xdat)
    mc <- match.call(expand.dots = FALSE)
    mc[[1L]] <- quote(npregbw)
    .npRmpi_autodispatch_bind_data_promises(mc, environment(), sys.call(),
      sys.function(), parent.frame())
  }
  environment(method) <- asNamespace("npRmpi")
  nx <- ny <- 0L
  d <- data.frame(x = 1:6)
  out <- method(xdat = { nx <- nx + 1L; d }, ydat = { ny <- ny + 1L; 6:1 })
  expect_identical(c(nx, ny), c(1L, 1L))
  expect_identical(out$xdat, d)
  expect_identical(out$ydat, 6:1)
  wrapper <- function(...) method(...)
  out <- wrapper(xdat = { nx <- nx + 1L; d }, ydat = { ny <- ny + 1L; 6:1 })
  expect_identical(c(nx, ny), c(2L, 2L))
  expect_identical(out$xdat, d)
  expect_identical(out$ydat, 6:1)
  # A rewritten leaf expression must win over the original formal.
  original <- quote(method(xdat = original.x, ydat = original.y))
  changed <- quote(npregbw(xdat = replacement.x, ydat = original.y))
  owner <- list2env(list(xdat = d, ydat = 6:1))
  out <- bind(changed, owner, original, method, environment())
  expect_identical(out$xdat, quote(replacement.x))
  expect_identical(out$ydat, 6:1)
  unrelated <- function(xdat, ydat) NULL
  expect_identical(bind(changed, owner, original, unrelated, environment()), changed)
})

test_that("public native constructors materialize stateful data once before transport", {
  withr::local_options(np.messages = FALSE)
  captured <- NULL
  testthat::local_mocked_bindings(
    .npRmpi_require_active_slave_pool = function(...) invisible(TRUE),
    .npRmpi_autodispatch_active = function(...) TRUE,
    .npRmpi_distributed_call_impl = function(mc, ...) {
      captured <<- mc
      list(diagnostic = TRUE)
    }, .package = "npRmpi")
  d <- data.frame(x = seq_len(24) / 24, y = sin(1:24), z = cos(1:24))
  for (family in c("npregbw", "npcdensbw", "npcdistbw", "npudensbw", "npudistbw")) {
    nx <- ny <- 0L
    x <- function() { nx <<- nx + 1L; d["x"] + 100 * (nx - 1L) }
    y <- function() { ny <<- ny + 1L; d$y + 100 * (ny - 1L) }
    f <- get(family)
    if (family %in% c("npudensbw", "npudistbw")) {
      f(dat = x(), bws = .4, bandwidth.compute = FALSE)
      expect_identical(nx, 1L, info = family)
      expect_identical(captured$dat, d["x"], info = family)
    } else {
      h <- if (family == "npregbw") .4 else c(.4, .6)
      f(xdat = x(), ydat = y(), bws = h, bandwidth.compute = FALSE)
      expect_identical(c(nx, ny), c(1L, 1L), info = family)
      expect_identical(captured$xdat, d["x"], info = family)
      expect_identical(as.double(unlist(captured$ydat)), d$y, info = family)
    }
  }
})
