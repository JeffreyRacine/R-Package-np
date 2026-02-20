npplot <- function(bws = stop("'bws' has not been set"), ..., random.seed = 42){
  .npRmpi_require_active_slave_pool(where = "npplot()")
  .npRmpi_guard_no_auto_object_in_manual_bcast(bws, where = "npplot()")
  ## Save seed prior to setting
  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  ## Restore seed
  on.exit(if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv))

  UseMethod("npplot",bws)
}

