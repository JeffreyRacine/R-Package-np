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
.np_rng_scope_capture_local <- function(uniform.only = FALSE) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    # RNGkind() would initialize an absent seed and can discard hidden state.
    # Inspect the current representation without a draw or generator reset.
    on.exit(if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    .Call(C_np_rng_scope_materialize)
    previous <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    # Even uniform-only draws initialize an absent binding and may reset an
    # otherwise unused hidden normal cache. Preserve the absent-state refusal.
    .np_rng_scope_check(previous, uniform.only = FALSE)
    return(list(seed = previous, uniform.only = uniform.only, absent = TRUE))
  }
  if (bindingIsActive(".Random.seed", .GlobalEnv) || bindingIsLocked(".Random.seed", .GlobalEnv))
    stop("cannot safely restore temporary RNG state: caller seed binding is active or locked", call. = FALSE)
  previous <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  .np_rng_scope_check(previous, uniform.only = uniform.only)
  list(seed = previous, uniform.only = uniform.only)
}
.np_rng_scope_restore <- function(scope) {
  if (isTRUE(scope$absent)) {
    on.exit(if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    assign(".Random.seed", scope$seed, envir = .GlobalEnv)
    # A state swap can change the internal kind. Load the saved built-in
    # representation before removing its binding, without drawing a variate.
    .Call(C_np_rng_scope_restore_absent)
  } else assign(".Random.seed", scope$seed, envir = .GlobalEnv)
  invisible(NULL)
}
# A uniform-only transport record never certifies a seed reset or state swap.
.np_rng_scope_seed_local <- function(scope, random.seed) {
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

# Called only at the enumerated numerical entries, never by generic local scopes.
.npRmpi_rng_scope_collective <- function() {
  !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
    (.npRmpi_autodispatch_in_context() ||
       .npRmpi_autodispatch_called_from_bcast())
}
.npRmpi_rng_scope_agree <- function(failure, phase, comm = 1L, deferred = FALSE) {
  count <- mpi.allreduce(as.integer(!is.null(failure)), type = 1L,
                         op = "sum", comm = comm)
  if (count == 0L) return(invisible(NULL))
  rank <- mpi.comm.rank(comm)
  owner <- mpi.allreduce(if (!is.null(failure)) rank else mpi.comm.size(comm),
                         type = 1L, op = "min", comm = comm)
  message <- if (rank == owner) conditionMessage(failure) else NULL
  received <- mpi.bcast.Robj(message, rank = owner, comm = comm)
  if (rank != owner) message <- received
  refusal <- structure(list(message = paste0("MPI temporary RNG scope ", phase,
    " refused: ", message), call = NULL),
    class = c("np_rng_scope_refusal", "error", "condition"))
  if (isTRUE(deferred)) return(refusal)
  stop(refusal)
}

.np_rng_scope_capture <- function(uniform.only = FALSE, collective = FALSE, deferred = FALSE,
                                  same.kind = FALSE) {
  if (!isTRUE(collective)) return(.np_rng_scope_capture_local(uniform.only))
  failure <- NULL
  scope <- tryCatch(.np_rng_scope_capture_local(uniform.only),
    error = function(e) { failure <<- e; NULL },
    interrupt = function(e) { failure <<- e; NULL })
  refusal <- .npRmpi_rng_scope_agree(failure, "admission", deferred = deferred)
  if (inherits(refusal, "np_rng_scope_refusal")) return(refusal)
  # Rank-owned resampling plans require the same generator tuple. Inspect the
  # certified representation without resetting any rank or drawing a variate.
  if (isTRUE(same.kind)) {
    kind <- scope$seed[1L]
    bounds <- mpi.allreduce(c(kind, -kind), type = 1L, op = "max")
    if (bounds[1L] != -bounds[2L]) {
      refusal <- .npRmpi_rng_scope_agree(simpleError(
        "participating ranks use different RNG kinds; select the same uniform, normal and sample kinds before this test"),
        "kind agreement", deferred = deferred)
      if (inherits(refusal, "np_rng_scope_refusal")) return(refusal)
    }
  }
  scope$collective <- TRUE
  scope$deferred <- deferred
  scope
}
.np_rng_scope_seed <- function(scope, random.seed) {
  if (!isTRUE(scope$collective))
    return(.np_rng_scope_seed_local(scope, random.seed))
  failure <- NULL
  tryCatch(.np_rng_scope_seed_local(scope, random.seed),
    error = function(e) { failure <<- e; NULL },
    interrupt = function(e) { failure <<- e; NULL })
  .npRmpi_rng_scope_agree(failure, "activation", deferred = scope$deferred)
}
