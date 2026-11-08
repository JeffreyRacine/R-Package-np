# Invocation-local ownership of temporary R RNG state. No generator is replaced.
.np_rng_scope_check <- function(state, where = "caller", uniform.only = FALSE) {
  bad <- function(reason) stop(paste0("cannot safely restore temporary RNG state (", where, "): ", reason,
    "; explicitly select a supported built-in RNG and call set.seed() before this operation"), call. = FALSE)
  if (!is.integer(state) || is.object(state) || length(state) < 2L ||
      is.na(state[1L]) || state[1L] < 0L) bad("uncertifiable seed")
  u <- state[1L] %% 100L
  n <- (state[1L] %/% 100L) %% 100L
  s <- state[1L] %/% 10000L
  if (u == 5L || (!uniform.only && n == 3L)) bad("user-supplied generator state is not fully represented")
  if (!uniform.only && n == 2L) bad("Box-Muller state is not fully represented")
  if (!u %in% c(0L, 1L, 2L, 3L, 4L, 6L, 7L) ||
      !n %in% c(0L, 1L, 2L, 3L, 4L, 5L) || !s %in% c(0L, 1L) ||
      length(state) != c(4L, 3L, 3L, 626L, 102L, 0L, 102L, 7L)[u + 1L]) bad("unrecognized seed representation")
  invisible(NULL)
}
.np_rng_scope_capture <- function(uniform.only = FALSE) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    stop("cannot safely restore temporary RNG state: caller seed is absent; explicitly select a supported built-in RNG and call set.seed() first", call. = FALSE)
  if (bindingIsActive(".Random.seed", .GlobalEnv) || bindingIsLocked(".Random.seed", .GlobalEnv))
    stop("cannot safely restore temporary RNG state: caller seed binding is active or locked", call. = FALSE)
  previous <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  .np_rng_scope_check(previous, uniform.only = uniform.only)
  list(seed = previous, uniform.only = uniform.only)
}
.np_rng_scope_restore <- function(scope) {
  assign(".Random.seed", scope$seed, envir = .GlobalEnv)
  invisible(NULL)
}
# A uniform-only transport record never certifies a seed reset or state swap.
.np_rng_scope_seed <- function(scope, random.seed) {
  if (isTRUE(scope$uniform.only))
    stop("uniform-only RNG scope cannot reset the generator", call. = FALSE)
  set.seed(random.seed)
  invisible(NULL)
}
.np_rng_scope_state <- function(scope, state) {
  if (isTRUE(scope$uniform.only))
    stop("uniform-only RNG scope cannot replace the generator state", call. = FALSE)
  .np_rng_scope_check(state, "target")
  assign(".Random.seed", state, envir = .GlobalEnv)
  invisible(NULL)
}
