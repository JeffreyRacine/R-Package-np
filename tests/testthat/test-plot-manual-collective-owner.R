test_that("regression data plot local fallback excludes symmetric broadcasts", {
  ns <- asNamespace("npRmpi")
  for (broadcast in c(FALSE,TRUE)) for (pool in c(FALSE,TRUE)) for (mode in c("data","plot","plot-data")) {
    e <- new.env(parent=ns); e$owner <- "collective-or-ordinary"
    e$.np_plot_validate_public_dots <- function(...) NULL
    e$.np_plot_normalize_public_dots <- function(dots,...) dots
    e$.npRmpi_autodispatch_called_from_bcast <- function() broadcast
    e$.npRmpi_has_active_slave_pool <- function(...) pool
    e$.npRmpi_with_local_regression <- function(expr) {e$owner <- "local";force(expr)}
    e$.np_plot_from_slot <- function(...) e$owner
    f <- get(".np_plot_npregression",ns);environment(f)<-e
    actual <- f(list(bws=list()),xdat=1,ydat=1,plot.behavior=mode)
    expect_identical(actual,if(!broadcast&&!pool&&mode=="data")"local" else "collective-or-ordinary")
  }
})
