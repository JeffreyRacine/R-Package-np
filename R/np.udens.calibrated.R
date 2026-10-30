# Private MX1 capability entry. Existing methods never call this module.
.np_udens_mx1 <- local({
  legacy <- local({
# External WFC1 research reducer; contains no DGP or Gaussian truth.
wfc_levels <- c(seq.int(501L,999L)/1000,1)
wfc_prepare <- function(T, q=wfc_levels) {
 if(!is.matrix(T)||!is.numeric(T)||!nrow(T)||!ncol(T)||any(!is.finite(T))) stop('invalid finite curve matrix')
 if(!is.numeric(q)||!length(q)||any(!is.finite(q))||any(q<=.5|q>1)||is.unsorted(q,strictly=TRUE)) stop('invalid calibration levels')
 m<-nrow(T);B<-ncol(T);V<-matrix(0,m,B);depth<-rep(B,B)
 for(j in seq_len(m)) {
  v<-T[j,];V[j,]<-sort(v,method='radix')
  depth<-pmin(depth,rank(v,ties.method='max'),rank(-v,ties.method='max'))
 }
 kp<-1+floor(B*(1-q)/2);kb<-1+floor(B*(1-q)/(2*m))
 ks<-pmin(kp,sort(depth,method='radix')[B-ceiling(B*q)+1])
 stopifnot(all(kp>=ks),all(ks>=kb),all(diff(ks)<=0))
 list(V=V,B=B,m=m,q=q,k=list(pointwise=kp,simultaneous=ks,bonferroni=kb),depth=depth,
  resolution=list(pointwise=B*(1-q)/2,bonferroni=B*(1-q)/(2*m),
   limited=any(B*(1-q)/(2*m)<1),q1='finite retained extrema, not an infinite band'))
}
wfc_depth <- function(z,truth) {
 if(!is.numeric(truth)||length(truth)!=z$m||any(!is.finite(truth)))stop('invalid finite target')
 min(vapply(seq_len(z$m),function(j)min(sum(z$V[j,]<=truth[j]),sum(z$V[j,]>=truth[j])),0))
}
wfc_required <- function(z,truth) {
 ok<-which(z$k$simultaneous<=wfc_depth(z,truth))
 if(length(ok))z$q[ok[1]] else Inf
}
wfc_select <- function(scores,coverage=.95) {
 if(!is.numeric(scores)||!length(scores)||anyNA(scores)||any(scores<=.5)||any(scores>1&is.finite(scores))||
    length(coverage)!=1L||!is.finite(coverage)||coverage<=0||coverage>=1)stop('invalid calibration scores')
 k<-ceiling((length(scores)+1)*coverage)
 if(k>length(scores))return(list(q=Inf,k=k,status='unattainable outer resolution'))
 q<-sort(scores)[k]
 list(q=q,k=k,status=if(is.finite(q))'finite pilot calibration' else 'unattainable retained envelope')
}
wfc_bands <- function(z,q) {
 a<-match(q,z$q);if(is.na(a))stop('level absent from frozen set')
 lapply(z$k,function(k)cbind(lower=z$V[cbind(seq_len(z$m),k[a])],upper=z$V[cbind(seq_len(z$m),z$B+1-k[a])]))
}
wfc_score <- function(z,q,truth) {
 b<-wfc_bands(z,q)
 do.call(rbind,lapply(names(b),function(band) {
  v<-b[[band]];below<-truth<v[,1];above<-truth>v[,2]
  data.frame(band=band,q=q,covered=!any(below|above),below=any(below),above=any(above),
   pointwise.mean=mean(!below&!above),mean.width=mean(v[,2]-v[,1]),
   lower.negative=mean(v[,1]<0),extreme.miss=any((below|above)[c(1,z$m)]),
   below.targets=paste(which(below),collapse=';'),above.targets=paste(which(above),collapse=';'))
 }))
}
# Only the new WFC2 external candidate consumes this policy.
wfc_companion <- function(z,nominal=.95) {
 a<-match(nominal,z$q);if(is.na(a))stop('nominal content absent')
 kp<-z$k$pointwise[a];kb<-z$k$bonferroni[a]
 z$raw.k<-z$k;z$nominal<-nominal
 z$k$simultaneous<-pmin(kp,z$raw.k$simultaneous)
 z$k$pointwise<-if(z$m==1L)z$k$simultaneous else rep(kp,length(z$q))
 z$k$bonferroni<-pmin(kb,z$k$simultaneous)
 stopifnot(all(z$k$pointwise>=z$k$simultaneous),all(z$k$simultaneous>=z$k$bonferroni))
 z
}

 as.list(environment())
 })
sc1 <- local({
 old <- legacy
  certificate <- list(family='npudens',dimension=1L,kernel='gaussian2',
                      bwtype='fixed',law='inid',smoothing='positive Gaussian')
  support <- function(h,contract=certificate) {
    if(!identical(contract,certificate))stop('SC1 support certificate does not cover this method',call.=FALSE)
    if(!is.numeric(h)||length(h)!=1L||!is.finite(h)||h<=0)
      stop('SC1 requires one positive finite fixed bandwidth',call.=FALSE)
    bounds <- c(lower=-stats::dnorm(sqrt(5))/h,upper=1.5*stats::dnorm(0)/h)
    if(any(!is.finite(bounds))||bounds[1]>=bounds[2])
      stop('SC1 support arithmetic is undefined',call.=FALSE)
    bounds
  }
  prepare <- function(T,h,contract=certificate) {
    bound<-support(h,contract)
    raw<-old$wfc_prepare(T)
    tol<-64*.Machine$double.eps*max(1,abs(bound))
    if(min(T)<bound[1]-tol||max(T)>bound[2]+tol)
      stop('SC1 retained curves contradict the certified support',call.=FALSE)
    # Outward floating-point hull preserves literal finite-array nesting.
    z<-old$wfc_companion(raw)
    z$support<-cbind(lower=pmin(rep(bound[1],z$m),z$V[,1]),
                     upper=pmax(rep(bound[2],z$m),z$V[,z$B]))
    z$analytic.support<-bound
    z
  }
  complete_required <- function(required,truth,h,contract=certificate) {
    bounds<-support(h,contract)
    if(!is.numeric(required)||length(required)!=1L||is.na(required)||required<=.5||
       (is.finite(required)&&required>1))stop('SC1 invalid retained required score',call.=FALSE)
    if(!is.numeric(truth)||!length(truth)||any(!is.finite(truth)))
      stop('SC1 invalid pilot target',call.=FALSE)
    if(is.finite(required))return(required)
    if(all(bounds[1]<=truth&truth<=bounds[2]))1 else Inf
  }
  required <- function(z,truth) {
    r<-old$wfc_required(z,truth)
    if(is.finite(r))return(r)
    if(all(z$support[,1]<=truth&truth<=z$support[,2]))1 else Inf
  }
  bands <- function(z,q) {
    out<-old$wfc_bands(z,q)
    if(q==1) {
      out$simultaneous<-out$bonferroni<-z$support
      if(z$m==1L)out$pointwise<-z$support
    }
    out
  }
  reduce <- function(T,h,outer,warn=TRUE,contract=certificate) {
    if(!is.list(outer)||length(outer)!=39L)
      stop('SC1 requires exactly 39 retained outer records',call.=FALSE)
    if(!is.logical(warn)||length(warn)!=1L||is.na(warn))stop('SC1 invalid warning control',call.=FALSE)
    scores<-vapply(outer,function(a)complete_required(a$row$q.required,a$truth,a$row$h,contract),0)
    retained<-vapply(outer,function(a)a$row$q.required,0)
    selection<-old$wfc_select(scores)
    z<-prepare(T,h,contract) # Exactly one original preparation.
    attainable<-is.finite(selection$q)
    q<-if(attainable)selection$q else 1
    out<-bands(z,q);boundary<-q==1
    at<-match(q,z$q);indices<-vapply(z$k,function(k)k[at],0)
    messages<-character()
    if(boundary)messages<-c(messages,'Full bootstrap-support endpoint used; simultaneous bounds may be very wide')
    if(any(indices==1))messages<-c(messages,'Empirical tail resolution is limited; every requested band is retained')
    if(!attainable)messages<-c(messages,'Calibration remains unattainable; returned endpoint is diagnostic only')
    if(warn&&length(messages))warning(paste(messages,collapse='; '),call.=FALSE)
    list(bands=out,selection=selection,scores=scores,retained.scores=retained,
         attainable=attainable,diagnostic.only=!attainable,boundary.used=boundary,
         outer.completions=sum(!is.finite(retained)&is.finite(scores)),
         order.index=indices,analytic.support=z$analytic.support,
         messages=messages,calibration='OP1 with SC1 support completion',
         endpoint.source=if(boundary)'analytic conditional bootstrap support' else 'retained empirical ranks')
  }
  list(support=support,prepare=prepare,complete_required=complete_required,
       required=required,bands=bands,reduce=reduce,certificate=certificate,old=old)
})

  # MX1 joint pilot and support. No draws, truth, or coverage outcomes enter selection.
  numeric_mask <- function(state) if (is.null(state$numeric)) rep(TRUE,ncol(state$x)) else state$numeric
  layout <- function(state) if (state$r+state$s) list(numeric=state$numeric,category=state$category) else NULL
  legacy_domain <- function(state) state$q==1L && state$r+state$s==0L
  lap_zero <- function(q,r) (-1)^r*exp(r*log(2)+lgamma(q/2+r)-lgamma(q/2)-q/2*log(2*pi))
  lap_gaussian <- function(radius,q,r,scale) {
    v<-radius/(2*scale^2);alpha<-q/2-1
    previous<-rep(1,length(v));current<-1+alpha-v
    if(r==0L)current<-previous
    if(r>1L)for(k in 2:r) {
      next_value<-((2*k-1+alpha-v)*current-(k-1+alpha)*previous)/k
      previous<-current;current<-next_value
    }
    (-2)^r*factorial(r)*current*exp(-v-q/2*log(2*pi)-(q+2*r)*log(scale))
  }
  cat_overlap <- function(category,a,b) {
    d<-abs(outer(a,b,'-'));l<-category$lambda;k<-category$k
    if(category$ordered) {
      if(l==0)return(1*(d==0))
      return(((1-l)/(1+l))^2*l^d*(d+1+2*l^2/(1-l^2)))
    }
    if(k==1L)return(matrix(1,length(a),length(b)))
    same<-(1-l)^2+l^2/(k-1)
    different<-2*(1-l)*l/(k-1)+(k-2)*(l/(k-1))^2
    ifelse(d==0,same,different)
  }
  joint_functional <- function(state,z,r=0L,scale=1,block=32L) {
    n<-nrow(z);total<-0
    for(a in seq.int(1L,n,by=block))for(b in seq.int(1L,n,by=block)) {
      ia<-a:min(n,a+block-1L);ib<-b:min(n,b+block-1L)
      w<-matrix(1,length(ia),length(ib));radius<-matrix(0,length(ia),length(ib))
      for(j in seq_len(ncol(z)))radius<-radius+outer(z[ia,j],z[ib,j],'-')^2
      for(c in state$category)w<-w*cat_overlap(c,c$training[ia],c$training[ib])
      value<-if(r==0L)rep(1,length(radius)) else lap_gaussian(as.vector(radius),ncol(z),r,scale)
      total<-total+sum(w*value)
    }
    total/n^2
  }
  pilot <- function(state) {
    if(legacy_domain(state))return(list(g=ks::hpi(state$x[[1]],nstage=2,binned=TRUE,deriv.order=4),
      method='unchanged OP1 univariate hpi derivative 4',eta=numeric()))
    nm<-numeric_mask(state);q<-sum(nm)
    if(q<1L)stop('MX1 requires at least one continuous coordinate',call.=FALSE)
    z<-sweep(as.matrix(state$x[nm]),2,state$h[nm],'/')
    sigma<-sqrt(mean(apply(z,2,stats::var)))
    diagonal<-1
    for(c in state$category)diagonal<-diagonal*cat_overlap(c,c$training[1],c$training[1])[1]
    average<-joint_functional(state,z)
    psi8<-average*lap_zero(q,8L)/(sqrt(2)*sigma)^(q+16)
    a7<-(-2*diagonal*lap_zero(q,7L)/(state$n*psi8))^(1/(q+16))
    psi7<-joint_functional(state,z,7L,a7)
    a6<-(-2*diagonal*lap_zero(q,6L)/(state$n*psi7))^(1/(q+14))
    psi6<-joint_functional(state,z,6L,a6)
    roughness<-lap_zero(q,4L)/sqrt(2)^(q+8)
    t<-((q+8)*diagonal*roughness/(state$n*psi6))^(1/(q+12))
    if(any(!is.finite(c(sigma,psi8,a7,psi7,a6,psi6,t)))||
       any(c(sigma,psi8,a7,a6,psi6,t)<=0)||psi7>=0)
      stop('MX1 joint derivative pilot arithmetic is undefined',call.=FALSE)
    g<-state$h;g[nm]<-t*g[nm]
    list(g=g,eta=g[!nm],method='joint scaled Laplacian-fourth plug-in; retained categorical smoothing',
      diagnostics=list(sigma=sigma,D=diagonal,Cbar=average,psi8=psi8,a7=a7,psi7=psi7,a6=a6,psi6=psi6,t=t))
  }
  support <- function(state) {
    if(legacy_domain(state))return(sc1$support(state$h))
    nm<-numeric_mask(state);q<-sum(nm)
    c0<-exp(-q/2*log(2*pi)-sum(log(state$h[nm])))
    bounds<-if(q) c(-c0*exp(-(q+4)/2),c0*(q+2)/2) else c(1,1)
    for(c in state$category) {
      l<-c$lambda;k<-c$k
      if(c$ordered) {
        if(l==0)z<-c(0,1) else {
          a<-(1-l)/(1+l);A<-2-a-2*l^2/(1+l)^2
          stationary<-A/a-1/log(l)
          d<-unique(c(0,max(0,floor(stationary)),max(0,ceiling(stationary))))
          z<-c(0,a*l^d*(A-a*d))
        }
      } else if(k==1L)z<-1 else if(!is.null(c$constant))z<-2*c$constant-k*c$constant^2 else {
        e<-k*l^2/(k-1);z<-c(1-e,e/(k-1))
      }
      bounds<-range(outer(bounds,range(z),'*'))
    }
    names(bounds)<-c('lower','upper')
    if(any(!is.finite(bounds))||bounds[1]>=bounds[2])stop('MX1 certified support arithmetic is undefined',call.=FALSE)
    bounds
  }
  complete <- function(score,truth,bounds) {
    if(length(score)!=1L||!is.numeric(score)||is.na(score)||score<=.5||
       (is.finite(score)&&score>1)||!is.numeric(truth)||!length(truth)||any(!is.finite(truth)))
      stop('invalid MX1 calibration score or pilot target',call.=FALSE)
    if(is.finite(score))return(score)
    if(all(bounds[1]<=truth & truth<=bounds[2]))1 else Inf
  }
  # Positive generating density, evaluated on physical mixed coordinates.
  pilot_value <- function(state,g,grid) {
    if(legacy_domain(state))return(vapply(grid[[1]],function(t)mean(stats::dnorm((t-state$x[[1]])/g)/g),0))
    nm<-numeric_mask(state);x<-toMatrix(state$x);target<-toMatrix(grid)
    out<-numeric(nrow(grid))
    for(i in seq_len(nrow(grid))) {
      K<-rep(1,state$n)
      for(j in which(nm))K<-K*stats::dnorm((target[i,j]-x[,j])/g[j])/g[j]
      for(c in state$category) {
        j<-c$index;d<-abs(target[i,j]-x[,j]);l<-g[j];k<-c$k
        S<-if(c$ordered)(1-l)/(1+l)*l^d else if(k==1L)rep(1,state$n) else ifelse(d==0,1-l,l/(k-1))
        K<-K*S
      }
      out[i]<-mean(K)
    }
    out
  }

  # MPI transport may use coordinator RNG for its own bookkeeping. Each private
  # tile replay owns explicit task streams and must preserve the caller state.
  preserve_rng <- function(start,expr) {
    existed<-exists('.Random.seed',.GlobalEnv,inherits=FALSE)
    previous<-if(existed)get('.Random.seed',.GlobalEnv,inherits=FALSE) else NULL
    kind<-RNGkind()
    on.exit({
      do.call(RNGkind,as.list(kind))
      if(existed)assign('.Random.seed',previous,.GlobalEnv) else
        if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm(list='.Random.seed',envir=.GlobalEnv)
    },add=TRUE)
    if(!existed)assign('.Random.seed',start,.GlobalEnv)
    force(expr)
  }
  provider <- function(state,grid,B,streams,law='inid',blocklen=NULL,tile=25L,chunk=32L) preserve_rng(streams$count,{
    .np_udens_bc_grid(state,grid)
    B<-.np_plot_gaussian_integer(B,'MX1 draws',9999L)
    tile<-.np_plot_gaussian_integer(tile,'MX1 target tile',25L)
    chunk<-.np_plot_gaussian_integer(chunk,'MX1 draw chunk',32L)
    law<-match.arg(law,c('inid','fixed','geom'))
    if(law!='inid' && !identical(as.numeric(blocklen),1))
      stop('MX1 nonunit block outer calibration is not yet admitted',call.=FALSE)
    if(!is.list(streams)||!identical(names(streams),c('count','jitter'))||
       !all(vapply(streams,is.integer,logical(1))))stop('invalid MX1 explicit streams',call.=FALSE)
    bound<-support(state);m<-nrow(grid);lay<-layout(state)
    x<-if(is.null(lay))as.matrix(state$x) else toMatrix(state$x)
    budget<-.np_udens_b2_budget(state,m,B,tile,FALSE)
    extra<-8*(8*min(m,tile)*state$n*chunk+6*state$n*ncol(x)*chunk+12*B+12*m+8*32^2)+
      32*ceiling(B/chunk)*(length(streams$count)+length(streams$jitter)+32)
    budget$per.process<-budget$per.process+extra
    budget$bytes<-budget$per.process*budget$processes
    if(!is.finite(budget$bytes)||budget$bytes>state$max.bytes)
      stop('MX1 live-byte forecast exceeds the requested budget',call.=FALSE)
    tasks<-.np_udens_smoothed_plan(x,state$h,B,law,blocklen,streams,chunk,FALSE,lay)
    end<-attr(tasks,'streams_end',exact=TRUE)
    get<-function(ii) {
      if(!is.numeric(ii)||!length(ii)||length(ii)>tile||anyNA(ii)||any(ii!=as.integer(ii))||any(ii<1L|ii>m))
        stop('invalid MX1 target tile',call.=FALSE)
      preserve_rng(streams$count,.np_udens_smoothed_draws(tasks,x,state$h,grid[ii,,drop=FALSE],B,law,blocklen,lay))
    }
    list(get=get,m=m,B=B,tile=tile,bound=bound,start=streams,end=end,budget=budget)
  })
  scan <- function(p,truth=NULL,nominal=.95) {
    m<-p$m;B<-p$B;depth<-rep(B,B);target.depth<-B
    empirical<-matrix(0,m,2L,dimnames=list(NULL,c('lower','upper')))
    if(!is.null(truth)&&(!is.numeric(truth)||length(truth)!=m||any(!is.finite(truth))))
      stop('invalid MX1 pilot target',call.=FALSE)
    tolerance<-64*.Machine$double.eps*max(1,abs(p$bound))
    for(a in seq.int(1L,m,by=p$tile)) {
      ii<-a:min(m,a+p$tile-1L);T<-p$get(ii)
      if(!is.matrix(T)||!identical(dim(T),c(length(ii),as.integer(B)))||any(!is.finite(T)))
        stop('invalid MX1 finite curve tile',call.=FALSE)
      if(min(T)<p$bound[1]-tolerance||max(T)>p$bound[2]+tolerance)
        stop('MX1 curves contradict certified support',call.=FALSE)
      for(j in seq_along(ii)) {
        v<-T[j,];depth<-pmin(depth,rank(v,ties.method='max'),rank(-v,ties.method='max'))
        empirical[ii[j],]<-range(v)
        if(!is.null(truth))target.depth<-min(target.depth,sum(v<=truth[ii[j]]),sum(v>=truth[ii[j]]))
      }
    }
    q<-sort(unique(c(legacy$wfc_levels,nominal)))
    kp<-1+floor(B*(1-q)/2);kb<-1+floor(B*(1-q)/(2*m))
    ks<-pmin(kp,sort(depth,method='radix')[B-ceiling(B*q)+1])
    z<-legacy$wfc_companion(list(B=B,m=m,q=q,k=list(pointwise=kp,simultaneous=ks,bonferroni=kb),depth=depth),nominal)
    z$support<-cbind(lower=pmin(rep(p$bound[1],m),empirical[,1]),upper=pmax(rep(p$bound[2],m),empirical[,2]))
    z$analytic.support<-p$bound
    z$target.depth<-target.depth
    ok<-which(z$k$simultaneous<=target.depth)
    z$required<-if(length(ok))q[ok[1]] else Inf
    z
  }
  reduce <- function(p,outer,warn=TRUE,nominal=.95) {
    if(!is.list(outer)||length(outer)!=39L)stop('MX1 requires exactly 39 outer records',call.=FALSE)
    if(!is.logical(warn)||length(warn)!=1L||is.na(warn))stop('invalid MX1 warning control',call.=FALSE)
    retained<-vapply(outer,function(a)a$row$q.required,0)
    scores<-vapply(outer,function(a)complete(a$row$q.required,a$truth,
      if(!is.null(a$support))a$support else sc1$support(a$row$h)),0)
    selection<-legacy$wfc_select(scores,coverage=nominal)
    z<-scan(p,nominal=nominal);attainable<-is.finite(selection$q);q<-if(attainable)selection$q else 1
    at<-match(q,z$q);indices<-vapply(z$k,function(k)k[at],0)
    out<-lapply(indices,function(k)matrix(0,p$m,2L,dimnames=list(NULL,c('lower','upper'))))
    for(a in seq.int(1L,p$m,by=p$tile)) {
      ii<-a:min(p$m,a+p$tile-1L);T<-p$get(ii)
      for(j in seq_along(ii)) {
        v<-sort(T[j,],method='radix')
        for(b in names(out)) {k<-indices[b];out[[b]][ii[j],]<-v[c(k,p$B+1-k)]}
      }
    }
    boundary<-q==1
    if(boundary) {
      out$simultaneous<-out$bonferroni<-z$support
      if(p$m==1L)out$pointwise<-z$support
    }
    messages<-character()
    if(boundary)messages<-c(messages,'Full bootstrap-support endpoint used; simultaneous bounds may be very wide')
    if(any(indices==1))messages<-c(messages,'Empirical tail resolution is limited; every requested band is retained')
    if(!attainable)messages<-c(messages,'Calibration remains unattainable; returned endpoint is diagnostic only')
    if(warn&&length(messages))warning(paste(messages,collapse='; '),call.=FALSE)
    list(bands=out,selection=selection,scores=scores,retained.scores=retained,
      attainable=attainable,diagnostic.only=!attainable,boundary.used=boundary,
      outer.completions=sum(!is.finite(retained)&is.finite(scores)),order.index=indices,
      analytic.support=z$analytic.support,messages=messages,calibration='OP1 with SC1 support completion',
      endpoint.source=if(boundary)'analytic conditional bootstrap support' else 'retained empirical ranks')
  }

  # Reconstruct an outer fit frame without mutating original level registries.
  frame <- function(state,y) {
    nm<-numeric_mask(state);out<-as.data.frame(y);names(out)<-names(state$x)
    for(c in state$category) {
      j<-c$index
      if(c$ordered) {
        values<-sort(unique(c(c$support,y[,j])))
        out[[j]]<-factor(as.character(y[,j]),levels=as.character(values),ordered=TRUE)
      } else {
        at<-match(y[,j],c$support)
        if(anyNA(at))stop('MX1 unordered draw left declared support',call.=FALSE)
        out[[j]]<-factor(c$levels[at],levels=c$levels,ordered=FALSE)
      }
    }
    if(anyNA(out)||!identical(unname(toMatrix(out)),unname(y)))
      stop('MX1 outer frame changed physical coordinates',call.=FALSE)
    out
  }
  outer_sample <- function(state,g,seed) {
    x<-if(legacy_domain(state))state$x[[1]] else toMatrix(state$x)
    .np_bias_stream(seed,{
      donor<-sort(sample.int(state$n,state$n,replace=TRUE))
      if(legacy_domain(state))return_value<-data.frame(X1=x[donor]+g*stats::rnorm(state$n)) else {
        y<-x[donor,,drop=FALSE];nm<-numeric_mask(state)
        for(j in seq_len(ncol(y))) {
          if(nm[j])y[,j]<-y[,j]+g[j]*stats::rnorm(state$n) else {
            c<-state$category[[match(j,vapply(state$category,function(a)a$index,integer(1)))]]
            l<-g[j]
            if(c$ordered)y[,j]<-y[,j]+stats::rgeom(state$n,1-l)-stats::rgeom(state$n,1-l) else {
              u<-stats::runif(state$n);selected<-match(y[,j],c$support);change<-u>=1-l
              if(l>0&&any(change)) {
                other<-floor((u[change]-(1-l))/l*(c$k-1))+1L
                selected[change]<-other+(other>=selected[change])
              }
              y[,j]<-c$support[selected]
            }
          }
        }
        return_value<-frame(state,y)
      }
      names(return_value)<-names(state$x)
      return_value
    })
  }
  family <- function(state,points,profiles=NULL,max.targets=100000L) {
    nm<-numeric_mask(state);q<-sum(nm)
    if(length(points)!=q||anyNA(points)||any(!is.finite(points))||any(points<1|points!=as.integer(points)))
      stop('MX1 requires positive integer target counts per numeric coordinate',call.=FALSE)
    if(!is.null(profiles)&&(!is.matrix(profiles)||ncol(profiles)!=sum(!nm)||!nrow(profiles)||any(!is.finite(profiles))))
      stop('invalid MX1 requested categorical profiles',call.=FALSE)
    if(sum(!nm)&&is.null(profiles))stop('MX1 requires explicit categorical profiles',call.=FALSE)
    size<-prod(as.double(points))*if(is.null(profiles))1 else nrow(profiles)
    if(!is.finite(size)||size>max.targets)stop('MX1 requested family exceeds declared target budget',call.=FALSE)
    axes<-lapply(which(nm),function(j)seq(min(state$x[[j]]),max(state$x[[j]]),length.out=points[match(j,which(nm))]))
    numeric.grid<-expand.grid(axes,KEEP.OUT.ATTRS=FALSE);names(numeric.grid)<-names(state$x)[nm]
    if(is.null(profiles))return(numeric.grid)
    y<-matrix(0,size,ncol(state$x));m<-nrow(numeric.grid)
    for(a in seq_len(nrow(profiles))) {
      ii<-(a-1L)*m+seq_len(m);y[ii,nm]<-as.matrix(numeric.grid)
      y[ii,!nm]<-matrix(rep(profiles[a,],each=m),m,sum(!nm))
    }
    out<-as.data.frame(y);names(out)<-names(state$x)
    for(c in state$category) {
      j<-c$index;at<-match(y[,j],c$support)
      if(anyNA(at))stop('MX1 requested profile absent from fitted category registry',call.=FALSE)
      out[[j]]<-factor(c$levels[at],levels=c$levels,ordered=c$ordered)
    }
    out
  }
  capture <- function(expr) {
    warnings<-character()
    value<-withCallingHandlers(force(expr),warning=function(e){warnings<<-c(warnings,conditionMessage(e));invokeRestart('muffleWarning')})
    list(value=value,warnings=warnings)
  }
  run <- function(fit,x,grid,streams,outer,B=9999L,points=NULL,profiles=NULL,
                  law='inid',blocklen=NULL,max.bytes=256*1024^2,tile=25L,chunk=32L,warn=TRUE,progress=NULL,
                  alpha=.05,family.builder=NULL,prepared=NULL) {
    # Explicit streams are sufficient; leave ambient RNG presence/kind untouched.
    existed<-exists('.Random.seed',.GlobalEnv,inherits=FALSE)
    previous<-if(existed)get('.Random.seed',.GlobalEnv,inherits=FALSE) else NULL
    kind<-RNGkind()
    on.exit({
      do.call(RNGkind,as.list(kind))
      if(existed)assign('.Random.seed',previous,.GlobalEnv) else
        if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm(list='.Random.seed',envir=.GlobalEnv)
    },add=TRUE)
    if(!existed) {
      if(!is.list(streams)||!is.integer(streams$count))stop('invalid MX1 explicit streams',call.=FALSE)
      assign('.Random.seed',streams$count,.GlobalEnv)
    }
    if(!is.numeric(alpha)||length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=.5)
      stop('smoothed-percentile requires one alpha in (0, .5)',call.=FALSE)
    if(!is.null(family.builder)&&!is.function(family.builder))
      stop('invalid calibrated family builder',call.=FALSE)
    state<-if(is.null(prepared)) .np_udens_bc_prepare(fit$bws,x,max.bytes,admit.effective.constant=TRUE,admit.ordered.literal=TRUE) else prepared
    if(state$q<1L)stop('MX1 pure categorical calibration is not yet admitted',call.=FALSE)
    if(!is.list(outer)||length(outer)!=39L||!all(vapply(outer,function(z)
      is.list(z)&&all(c('data','cv','count','jitter')%in%names(z))&&
      all(vapply(z[c('data','cv','count','jitter')],is.integer,logical(1))),logical(1))))
      stop('MX1 requires 39 explicit independent outer stream records',call.=FALSE)
    if(!is.null(progress)&&!is.function(progress))stop('invalid MX1 progress callback',call.=FALSE)
    if(is.null(points))points<-if(legacy_domain(state))50L else
      stop('MX1 requires explicit outer numeric family design',call.=FALSE)
    nm<-numeric_mask(state)
    if(state$r+state$s && is.null(profiles))profiles<-unique(toMatrix(grid)[,!nm,drop=FALSE])
    # Validate target design before selecting a pilot or consuming outer work.
    design<-if(is.null(family.builder))family(state,points,profiles) else family.builder(state)
    sorted_coordinates<-function(z) {
      mat<-toMatrix(z);mat[do.call(order,as.data.frame(mat)),,drop=FALSE]
    }
    if(!identical(unname(sorted_coordinates(design)),unname(sorted_coordinates(grid))))
      stop('MX1 original and outer family designs differ',call.=FALSE)
    original<-provider(state,grid,B,streams,law,blocklen,tile,chunk)
    retained.bytes<-8*39*nrow(grid)*(2*ncol(x)+12)+8*39*(500+4*length(streams$count))
    if(state$r>0L)retained.bytes<-retained.bytes+4096*39*state$r
    if(state$s>0L)retained.bytes<-retained.bytes+4096*39*state$s
    if(!is.finite(retained.bytes)||original$budget$bytes+retained.bytes>max.bytes)
      stop('MX1 retained outer diagnostics exceed the requested budget',call.=FALSE)
    original$budget$retained.outer.bytes<-retained.bytes
    original$budget$bytes<-original$budget$bytes+retained.bytes
    pilot.capture<-capture(pilot(state));selected<-pilot.capture$value
    records<-vector('list',39L)
    for(j in seq_len(39L)) {
      os<-outer[[j]]
      context<-list(index=j,stage='outer sample',seeds=os,pilot=selected)
      records[[j]]<-tryCatch({
      sample<-outer_sample(state,selected$g,os$data);y<-sample$value
      context$data<-y;context$stage<-'outer CV'
      cv<-capture(.np_bias_stream(os$cv,{
        if(legacy_domain(state)&&identical(names(y),'X1'))npudens(~X1,data=y) else
          npudens(stats::reformulate(names(y)),data=y)
      }))
      fy<-cv$value$value;context$bandwidth<-fy$bws;context$stage<-'corrected provider preparation'
      sy<-.np_udens_bc_prepare(fy$bws,y,max.bytes,admit.effective.constant=TRUE,admit.ordered.literal=TRUE)
      context$stage<-'outer family and pilot target'
      gy<-if(is.null(family.builder))family(sy,points,profiles) else family.builder(sy)
      truth<-pilot_value(state,selected$g,gy)
      context$grid<-gy;context$stage<-'inner corrected draws and ranks'
      py<-provider(sy,gy,B,list(count=os$count,jitter=os$jitter),law,blocklen,tile,chunk)
      scan.y<-scan(py,truth,nominal=1-alpha)
      record<-list(row=list(index=j,h=sy$h,q.required=scan.y$required),truth=truth,
        support=py$bound,depth=scan.y$target.depth,kS=scan.y$k$simultaneous,
        data.end=sample$end,cv.end=cv$value$end,draw.end=py$end,warnings=cv$warnings,
        grid=gy,bandwidth=sy$h)
      constant.owners<-Filter(function(c)!is.null(c$effective.constant),sy$category)
      if(length(constant.owners))record$effective.categorical<-lapply(constant.owners,
        function(c)c[c('index','k','lambda','constant','owner','effective.constant')])
      literal.owners<-Filter(function(c)!is.null(c$literal.lattice),sy$category)
      if(length(literal.owners))record$literal.ordered<-lapply(literal.owners,
        function(c)c[c('index','lambda','owner','literal.lattice')])
      record
      },error=function(e){
        attr(e,'mx1.context')<-context
        attr(e,'mx1.completed.outer')<-records[seq_len(j-1L)]
        stop(e)
      })
      if(!is.null(progress))progress(j,39L)
    }
    answer<-reduce(original,records,warn,nominal=1-alpha)
    list(inference=answer,pilot=selected,outer=records,pilot.warnings=pilot.capture$warnings,
      rng=list(start=original$start,end=original$end),grid=grid,
      budget=original$budget,
      resolution=list(B=B,nominal=1-alpha,family.size=nrow(grid),order.index=answer$order.index,
        nominal.tail.draws=c(pointwise=B*alpha/2,bonferroni=B*alpha/(2*nrow(grid))),
        limited=any(answer$order.index==1)),settings=list(method='MX1 refined calibrated unshifted percentile',
        B=B,outer=39L,rank=as.integer(answer$selection$k),q=state$q,r=state$r,s=state$s,law=law,blocklen=blocklen,
        coverage.status='capability prototype; new-domain coverage unassessed'))
  }
  list(pilot=pilot,support=support,complete=complete,pilot_value=pilot_value,
       lap_zero=lap_zero,lap_gaussian=lap_gaussian,cat_overlap=cat_overlap,
       joint_functional=joint_functional,provider=provider,scan=scan,reduce=reduce,
       frame=frame,outer_sample=outer_sample,family=family,run=run,legacy=legacy,sc1=sc1)

})
.np_udens_calibrated_percentile <- function(...) .np_udens_mx1$run(...)
