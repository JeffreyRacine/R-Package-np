# Positive separable representation of the same frozen joint finite quadrature.
.np_udens_beta_series_rule <- function(theta,order) {
  base<-c((theta[1:4]+1)*.5-1,0)
  axes<-list(.np_udens_beta_marginal_rule(base[1:2],order),.np_udens_beta_marginal_rule(base[3:4],order))
  list(base=base,axes=axes,A=axes[[1]]$A+axes[[2]]$A)
}




