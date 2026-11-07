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
  if(band!='all')return(out)
  intervals<-lapply(result$bands,function(x)array(c(x$lower,x$upper)/jacobian,
    c(length(result$center),2L,1L),dimnames=list(NULL,c('lower','upper'),as.character(alpha))))
  offsets<-lapply(intervals,function(x)cbind(result$center-x[,1L,1L],x[,2L,1L]-result$center))
  for(key in c('bands','intervals','absolute.intervals'))out$inference[[key]]<-intervals
  out$boot.all.err<-offsets
  out$boot.err<-cbind(offsets$pointwise,result$center)
  out$inference$settings$available.bands<-names(intervals)
  out$inference$settings$companion.critical<-lapply(result$bands,`[[`,'critical')
  out$inference$resolution$gaussian.quantile.type<-7L
  out
}
