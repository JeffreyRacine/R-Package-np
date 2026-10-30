# Density transport uses the incumbent literal block laws and F4 reference units.
.np_udens_block_request <- function(x, law, blocklen = NULL) {
  automatic <- is.null(blocklen)
  selected <- if (automatic) b.star(x, round = TRUE) else NULL
  if (automatic) blocklen <- unname(selected[1L, 1L])
  request <- .np_udens_b2_boot_request(law, blocklen, nrow(x))
  request$automatic <- automatic
  request$selected <- selected
  request$source <- if (automatic) "b.star(xdat, round=TRUE)[1,1]" else "supplied block length"
  request
}

# The cache primitives are factor-generic; their CDF owner and selector stay
# unchanged. A single local cache serves both location and log-scale factors.
.np_udens_block_cache <- function(n, B, m, tile, law, blocklen, budget, max.bytes) {
  start <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  plan <- .np_udist_cache_plan(n, B, m, tile, blocklen, budget, max.bytes, length(start))
  if (is.null(plan)) return(NULL)
  stream <- .np_bias_stream(start, .np_udist_cache_prepare(n, B, law, blocklen))
  list(tasks = stream$value, plan = plan, n = n, B = B, law = law,
       blocklen = blocklen, start = start, end = stream$end)
}

.np_udens_block_draws <- function(F, B, law, blocklen, cache) {
  if (is.null(cache)) return(.np_bias_process_draws(F, B, law, blocklen))
  if (ncol(F) != cache$n || B != cache$B || !identical(law, cache$law) ||
      blocklen != cache$blocklen ||
      !identical(get(".Random.seed", .GlobalEnv, inherits = FALSE), cache$start))
    stop("density block cache does not match the requested process", call. = FALSE)
  .np_udist_cache_draws(F, cache$tasks, B, law, blocklen)
}
