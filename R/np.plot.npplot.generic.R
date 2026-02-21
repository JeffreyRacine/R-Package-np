npplot <- function(bws = stop("'bws' has not been set"), ..., random.seed = 42){
  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state), add = TRUE)
  UseMethod("npplot", bws)
}
