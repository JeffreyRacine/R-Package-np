npplot <- function(bws = stop("'bws' has not been set"), ..., random.seed = 42){
  .np_with_seed(random.seed, UseMethod("npplot", bws))
}
