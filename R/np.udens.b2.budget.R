# Charge replicated resident processes before starting new private calibration.
.np_udens_b2_budget <- function(state,m,B,tile,covariance) {
  budget<-.np_udens_bc_budget(state,m,B,tile,covariance,processes=1L)
  nn<-isTRUE(state$q>0L)&&isTRUE(state$type%in%c("generalized_nn","adaptive_nn"))
  extra<-if(nn)8*as.double(state$n)*min(tile,m)*
    (8*as.integer(state$type=="adaptive_nn")+4*(state$r+state$s)*as.integer(isTRUE(state$cache.enabled))) else 0
  budget$per.process<-budget$per.process+extra
  budget$bytes<-budget$per.process*budget$processes;budget$nn.extra.per.process<-extra
  if(!is.finite(budget$bytes)||budget$bytes>state$max.bytes)
    stop("density correction live-byte forecast exceeds the requested budget",call.=FALSE)
  budget
}
