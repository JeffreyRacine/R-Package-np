test_that("private regression dispatch retains its declared package leaf", {
  make <- getFromNamespace(".npregbw_owned_dispatch_call", "npRmpi")
  generic <- getFromNamespace(".npRmpi_autodispatch_as_generic_call", "npRmpi")
  for (owner in c("npregbw.NULL", "npregbw.default", "npregbw.rbandwidth")) {
    mc <- quote(npregbw.rbandwidth(xdat = stop("forced"), ydat = stop("forced"),
                                 bws = stop("forced"), degree = 2L))
    mc[[1L]] <- as.name(owner)
    context <- list(expanded = mc, public = FALSE)
    z <- make(context, owner)
    expect_identical(z[[1L]], as.name(owner))
    expect_identical(as.list(z)[-1L], as.list(mc)[-1L])
    expect_true(attr(z, ".npRmpi.bound.call.owner", exact = TRUE))
    context$public <- TRUE
    expect_identical(make(context, owner), generic("npregbw", mc))
    matched <- quote(npregbw.NULL(xdat = xdat, ydat = ydat, ... = pairlist(degree = 2L)))
    expect_identical(make(context, owner, matched), generic("npregbw", matched))
  }
  bad <- list(expanded = quote(foreign(xdat = stop("forced"))), public = FALSE)
  expect_error(make(bad, "npregbw.rbandwidth"), "private regression bandwidth dispatch owner")
  expect_error(make(bad, "foreign"), "private regression bandwidth dispatch owner")
  npregbw <- function(...) stop("ambient owner must not be used")
  expect_true(attr(make(list(expanded = quote(npregbw.default(xdat = x)), public = FALSE),
                        "npregbw.default"), ".npRmpi.bound.call.owner", exact = TRUE))
})
