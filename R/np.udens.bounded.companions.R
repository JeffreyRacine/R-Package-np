# Analytic companions of the retained guide-region/residual construction.
# The guide region already covers the complete family at risk alpha/2.
.np_udens_bounded_analytic <- function(rhat,region,scales,alpha,m) {
  if(length(m)!=1L||!is.finite(m)||m<1||m!=length(rhat)||
     length(region$lower)!=m||length(region$upper)!=m||
     !length(scales)||any(vapply(scales,length,integer(1L))!=m)||
     any(!is.finite(c(rhat,region$lower,region$upper,unlist(scales))))||
     any(unlist(scales)<0)||any(region$lower>region$upper)||
     length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1)
    stop('invalid bounded density analytic companion request',call.=FALSE)
  pointwise<-qnorm(1-alpha/4)
  critical<-c(pointwise=pointwise,
    bonferroni=if(m==1L)pointwise else qnorm(alpha/(4*m),lower.tail=FALSE))
  lapply(critical,function(z) {
    parts<-lapply(scales,function(s)list(lower=rhat+region$lower-z*s,
      upper=rhat+region$upper+z*s))
    list(lower=Reduce(pmin,lapply(parts,`[[`,'lower')),
      upper=Reduce(pmax,lapply(parts,`[[`,'upper')),critical=z)
  })
}

# Assemble companions in original density units; offsets stay signed.
.np_udens_bounded_public_companions <- function(out,result,band,jacobian,alpha) {
  out$inference$settings$family.unique.targets<-result$family_size
  out$inference$resolution$bonferroni.signed.tail.probability<-alpha/(4*result$family_size)
  if(!is.null(result$sd_companion)) {
    out$inference$resolution$bonferroni<-result$sd_companion
    out$inference$settings$bonferroni.construction<-result$sd_companion$construction
    out$inference$settings$resolution.limited<-result$sd_companion$resolution.limited
  }
  if(band!='all')return(out)
  intervals<-lapply(result$bands,function(x)array(c(x$lower,x$upper)/jacobian,
    c(length(result$center),2L,1L),dimnames=list(NULL,c('lower','upper'),as.character(alpha))))
  offsets<-lapply(intervals,function(x)cbind(result$center-x[,1L,1L],x[,2L,1L]-result$center))
  for(key in c('bands','intervals','absolute.intervals'))out$inference[[key]]<-intervals
  out$boot.all.err<-offsets
  out$boot.err<-cbind(offsets$pointwise,result$center)
  out$inference$settings$available.bands<-names(intervals)
  out$inference$settings$companion.critical<-lapply(result$bands,`[[`,'critical')
  if(is.null(result$sd_companion))out$inference$resolution$gaussian.quantile.type<-7L
  out
}

# Empirical Bonferroni companions with an unconditional simultaneous hull.
# The hull protects the incumbent finite-family coverage without extrapolating
# unobserved bootstrap tails or increasing the caller's draw count.
.np_udens_bounded_sd <- function(root,rhat,scales,region,center,alpha,band,calibration,
                                progress=function(...)NULL) {
  if(!band%in%c('bonferroni','all')||!calibration%in%c('symmetric-sd','asymmetric-sd')||
     !length(scales)||length(rhat)!=root$m||length(center)!=root$m||
     any(vapply(scales,length,integer(1L))!=root$m))
    stop('invalid bounded density SD companion request',call.=FALSE)
  m<-root$m;B<-root$B
  combine<-function(parts)list(lower=Reduce(pmin,lapply(parts,`[[`,'lower')),
    upper=Reduce(pmax,lapply(parts,`[[`,'upper')))
  sim<-point<-raw<-bonf<-vector('list',length(scales))
  for(j in seq_along(scales)) {
    progress('SD band reduction',j-1L,length(scales))
    sim[[j]]<-.np_udens_beta_stored_reduce(root,rhat,scales[[j]],region,center,
      alpha,'simultaneous',calibration)
    if(band=='all')point[[j]]<-.np_udens_beta_stored_reduce(root,rhat,scales[[j]],region,center,
      alpha,'pointwise',calibration)
    raw[[j]]<-.np_udens_beta_stored_reduce(root,rhat,scales[[j]],region,center,
      alpha/m,'pointwise',calibration)
    bonf[[j]]<-combine(list(raw[[j]],sim[[j]]))
    progress('SD band reduction',j,length(scales))
  }
  s<-combine(sim);b<-combine(bonf);r<-combine(raw)
  result<-if(band=='all')s else b
  if(band=='all')result$bands<-list(pointwise=combine(point),simultaneous=s,bonferroni=b)
  result$B<-B
  primary.alpha<-if(band=='all')alpha else alpha/m
  result$k_signed<-floor(B*primary.alpha/4)+1L
  result$k_absolute<-ceiling(B*(1-primary.alpha/2))
  tail<-alpha/(if(calibration=='asymmetric-sd')4*m else 2*m)
  result$sd_companion<-list(construction='empirical Bonferroni with simultaneous hull',
    family.unique.targets=m,signed.lower.index=floor(B*(alpha/m)/4)+1L,
    signed.upper.index=B-floor(B*(alpha/m)/4),
    absolute.index=ceiling(B*(1-(alpha/m)/2)),signed.tail.probability=alpha/(4*m),
    absolute.tail.probability=alpha/(2*m),expected.tail.draws=B*tail,
    resolution.limited=floor(B*tail)<20L,
    hull.expanded.targets=sum(b$lower<r$lower|b$upper>r$upper),
    simultaneous.containment=TRUE)
  result
}
