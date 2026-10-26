# Charge replicated resident processes before starting new private calibration.
.np_udens_b2_budget <- function(state,m,B,tile,covariance) {
  .np_udens_bc_budget(state,m,B,tile,covariance,processes=1L)
}
