npplot <- function(bws = stop("'bws' has not been set"), ..., random.seed = 42){
  .npRmpi_require_active_slave_pool(where = "npplot()")
  .npRmpi_guard_no_auto_object_in_manual_bcast(bws, where = "npplot()")
  .np_with_seed(random.seed, UseMethod("npplot", bws))
}
