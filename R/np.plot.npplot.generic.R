npplot <- function(bws = stop("'bws' has not been set"), ..., random.seed = 42){
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

