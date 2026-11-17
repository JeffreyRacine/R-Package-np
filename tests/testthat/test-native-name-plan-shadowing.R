test_that("native name plans omit only the selected identity-owned role", {
  withr::local_options(list(scipen = 0))
  plan <- getFromNamespace(".npRmpi_native_name_plan", "npRmpi")
  make <- function(xdat, ydat, zdat, identity.roles = character()) {
    plan(sys.function(), environment(), parent.frame(),
         list(xdat = xdat, ydat = ydat, zdat = zdat), identity.roles)
  }
  x <- seq_len(10L)
  expected <- list(ynames = deparse(quote(x * 1e-4)),
                   xnames = paste(deparse(quote(x * 1e5)), collapse = ""),
                   znames = paste(deparse(quote(x * 1e6)), collapse = ""))
  expect_identical(make(x * 1e5, x * 1e-4, x * 1e6), expected)
  for (roles in list(character(), "ydat", "xdat", "zdat", c("xdat", "ydat", "zdat"))) {
    removed <- c(ydat = "ynames", xdat = "xnames", zdat = "znames")[roles]
    keep <- !names(expected) %in% removed
    expect_identical(make(x * 1e5, x * 1e-4, x * 1e6, roles),
                     if (any(keep)) expected[keep] else list())
  }
  expect_identical(make(x * 1e5, x * 1e-4, x * 1e6, "tydat"), expected)
  for (role in c("dat", "tdat")) {
    definition <- function(dat) NULL
    names(formals(definition)) <- role
    owner <- list2env(setNames(list(quote(x)), role), parent = baseenv())
    values <- setNames(list(x), role)
    expect_identical(plan(definition, owner, baseenv(), values), list(xnames = "x"))
    expect_identical(plan(definition, owner, baseenv(), values, role), list())
    expect_identical(plan(definition, owner, baseenv(), values, "xdat"), list(xnames = "x"))
  }
  expect_identical(plan(NULL, environment(), parent.frame(), list()), list())
})
