# Exact pushforward identity for a categorical pilot, not a numerical cutoff.
.np_udens_b2_empirical_pilot_identity <- function(state) {
  isTRUE(state$q==0L)&&all(vapply(state$category,function(z)
    z$lambda==0||(!z$ordered&&(z$k==1L||z$pooling)),logical(1L)))
}
