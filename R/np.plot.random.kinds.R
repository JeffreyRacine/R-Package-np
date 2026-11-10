# Validate explicit sampler request metadata without querying or changing RNG state.
# Retain every supplied kind; children replay it on the exact same R runtime.
.np_plot_rng_kind_shape <- function(kind, minimum) {
  size <- length(kind)
  is.character(kind) && !anyNA(kind) &&
    ((size >= minimum && size <= 3L) ||
     (size == 4L && .np_rng_kind_count() == 4L)) &&
    (size < 3L || kind[3L] %in% c("Rounding", "Rejection")) &&
    (size < 4L || kind[4L] %in% c("Buggy BTPE", "BTPE"))
}
