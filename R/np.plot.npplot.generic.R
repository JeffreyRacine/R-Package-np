npplot <- function(bws = stop("'bws' has not been set"), ..., random.seed = 42){
  .npRmpi_require_active_slave_pool(where = "npplot()")
  .npRmpi_guard_no_auto_object_in_manual_bcast(bws, where = "npplot()")
  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state), add = TRUE)
  UseMethod("npplot", bws)
}
