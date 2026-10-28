# Pilot-anchored finite-sample selected-block scale; no change to literal draws.
.np_bias_pilot_block_scale <- function(F, pilot.se, blocklen, law) {
  F <- F-rowMeans(F)
  independent <- rowSums(F^2)
  dependent <- .np_bias_score_variance(F,blocklen,law)
  active <- pilot.se>0
  if (length(pilot.se)!=nrow(F) || any(!is.finite(pilot.se)) || any(pilot.se<0) ||
      any(!is.finite(independent)) || any(!is.finite(dependent)) ||
      any(independent[active]<=0) || any(dependent[active]<=0) ||
      any(independent[!active]!=0) || any(dependent[!active]!=0))
    stop("undefined pilot-anchored selected-block variance",call.=FALSE)
  scale <- numeric(length(pilot.se));ratio <- numeric(length(pilot.se))
  ratio[active] <- dependent[active]/independent[active]
  scale[active] <- pilot.se[active]*sqrt(ratio[active])
  if (any(!is.finite(scale)) || any(scale[active]<=0))
    stop("nonpositive pilot-anchored block scale",call.=FALSE)
  list(se=scale,pilot.se=pilot.se,independent=independent,
       dependent=dependent,variance.ratio=ratio)
}
