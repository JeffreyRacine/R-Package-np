# Exercise coordinator-only lease operations without disturbing the live pool.
local_test_lease_state <- function(.env = parent.frame()) {
  state <- list2env(list(status = "open", session = NULL,
                         store = new.env(parent = emptyenv()),
                         registry = new.env(parent = emptyenv()),
                         pending = new.env(parent = emptyenv())),
                    parent = emptyenv())
  testthat::local_mocked_bindings(.npRmpi_lease_state = state,
                                  .package = "npRmpi", .env = .env)
}
