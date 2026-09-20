test_that("npsigtest preparation retains dispatched and collective owners", {
  prepare <- getFromNamespace(".npRmpi_npsig_npreg_prepare", "npRmpi")
  env <- new.env(parent = environment(prepare))
  environment(prepare) <- env
  env$.npRmpi_npsig_do_leaf <- function(fun, ...) "leaf"
  env$.npRmpi_npsig_do_local <- function(...) "local"
  for (auto in c(FALSE, TRUE)) for (collective in c(FALSE, TRUE)) {
    env$.npRmpi_autodispatch_active <- local({ x <- auto; function() x })
    env$.npRmpi_npsig_collective_context <- local({ x <- collective; function() x })
    expect_identical(prepare(), if (auto || collective) "leaf" else "local")
  }
})
