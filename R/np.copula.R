## This function accepts a joint distribution bandwidth object,
## optionally a list of marginal distribution bandwidth objects, and
## optionally a matrix with columns being vectors of taus (u) for
## which the copula is to be computed. It returns a data frame with
## the copula `copula', the copula density `copula density' and u
## (columns for the grid constructed from the columns of the u matrix)

npcopula <- function(bws.joint,bws.univariate=NULL,u=NULL) {

  ## We make a copy of the joint bandwidth object then copy marginal
  ## bandwidths into each position if supplied. All attributes of the
  ## joint object are preserved. Note that for computational
  ## efficiency we call npudist directly using tdat arguments so here
  ## we must copy over all attributes.
  
  bws.marginal <- bws.joint  

  ## Basic error checking

  if(missing(bws.joint)) stop("You must provide a joint distribution bandwidth object")
  if(!is.null(u)) if(any(u>1) || any(u<0)) stop("u must lie in [0,1]")
  num.var <- length(bws.joint$xnames)
  
  ## Test for compatible quantile vector/matrix if provided
  if(!is.null(u)) {
    if(is.vector(u)) u <- matrix(u,1,length(u))
    if(ncol(u) != num.var) stop(paste("matrix u must have ", num.var," columns",sep=""))
  }

  ## If marginal bandwidth object list is provided copy each into
  ## appropriate place in joint object

  if(!is.null(bws.univariate)) {
    if(!is.list(bws.univariate)) stop("bws.univariate must be a list containing bandwidth objects")
    if(length(bws.univariate)!=num.var) stop(paste("bws.univariate must contain ", num.var, " bandwidth objects",sep=""))
    for(j in 1:num.var) bws.marginal$bw[j] <- bws.univariate[[j]]$bw[1]
  }

  console <- newLineConsole()

  if(is.null(u)) {
    ## Compute the copula distribution and density for the sample
    ## realizations (joint CDF)
    copula <- fitted(npudist(bws=bws.joint))
    copula.density <- fitted(npudens(bws=bws.joint))    
    ## Compute the marginal quantiles from the marginal CDFs (u_i=\hat
    ## F(x_i)) and divide copula density by its marginals
    u <- matrix(NA,bws.joint$nobs,num.var)
    for(j in 1:num.var) {
      bws.F <- npudensbw(formula(paste("~",bws.joint$xnames[j])),
                         bws=bws.marginal$bw[j],
                         bandwidth.compute=FALSE,
                         bwmethod=bws.marginal$method,
                         bwtype=bws.marginal$type,
                         ckerorder=bws.marginal$ckerorder,
                         ckertype=bws.marginal$ckertype,
                         okertype=bws.marginal$okertype,
                         ukertype=bws.marginal$ukertype)
      
      u[,j] <- fitted(npudist(bws=bws.F))
      ## For copula density we require marginal densities. Desirable
      ## to have the same bws in numerator and denominator, so use
      ## those from the joint (mirror regression, conditional density
      ## estimation etc.)
      bws.f <- npudensbw(formula(paste("~",bws.joint$xnames[j])),
                         bws=bws.joint$bw[j],
                         bandwidth.compute=FALSE,
                         bwmethod=bws.marginal$method,
                         bwtype=bws.marginal$type,
                         ckerorder=bws.marginal$ckerorder,
                         ckertype=bws.marginal$ckertype,
                         okertype=bws.marginal$okertype,
                         ukertype=bws.marginal$ukertype)

      copula.density <- copula.density/NZD(fitted(npudens(bws=bws.f)))
    }
  } else {
    ## User wishes to compute copula for inputted u matrix. To
    ## economize on computation we obtain the quantiles from the
    ## marginals then expand to compute the copula
    x.u <- matrix(NA,nrow(u),num.var)
    for(j in 1:num.var) {
      for(i in 1:nrow(u)) {
        ## Compute the quasi inverse (Definition 2.3.6, Nelson (2006)).
        ## Here we take pains to span a sufficiently rich set of
        ## evaluation points to cover a range of possibilities. In
        ## particular, we extend the range of the variable by a
        ## fraction 1 on min/max and create an equally spaced grid on
        ## this range to try to cover long-tailed distributions but
        ## provide a sufficiently fine grid on this extended range
        ## (uniform spacing, add and subtract the range of the data
        ## to/from max/min). We also use a sequence of equi-quantile
        ## spaced points from the quantiles of the raw data again to
        ## provide a sufficiently fine grid. Finally, we use the
        ## datapoints themselves. We then concatenate and sort the
        ## equally space extended grid, the equi-quantile grid, and
        ## the data points themselves.
        x.er <- eval(parse(text=paste("extendrange(", bws.joint$xnames[j],",f=1)",sep="")))
        x.q <- eval(parse(text=paste("quantile(", bws.joint$xnames[j],",seq(0,1,length=500))",sep="")))
        x.eval <- sort(c(seq(x.er[1],x.er[2],length=500),x.q,eval(parse(text=bws.joint$xnames[j]))))
        F <- eval(parse(text=paste("fitted(npudist(tdat=",bws.joint$xnames[j],",edat=x.eval,bws=bws.marginal$bw[j],bwmethod=bws.marginal$method,bwtype=bws.marginal$type,ckerorder=bws.marginal$ckerorder,ckertype=bws.marginal$ckertype,okertype=bws.marginal$okertype,ukertype=bws.marginal$ukertype))",sep="")))
        x.u[i,j] <- ifelse(u[i,j]>=0.5, max(x.eval[F<=u[i,j]]), min(x.eval[F>=u[i,j]]))
      }
    }
    ## Convert to evaluation data frame and name so that we can use
    ## predict(,newdata=...). To compute the copula we expand the grid
    ## of marginal quantiles so that every combination of the columns
    ## of u is constructed
    console <- printPop(console)
    console <- printPush(msg = "Expanding the u matrix and computing the copula and copula density", console)
    x.u <- expand.grid(data.frame(x.u))
    names(x.u) <- bws.joint$xnames
    copula <- predict(npudist(bws=bws.joint),newdata=x.u)
    copula.density <- predict(npudens(bws=bws.joint),newdata=x.u)
    ## For the copula density require marginal densities. Desirable to
    ## have the same bws in numerator and denominator, so use those
    ## from the joint (mirror regression, conditional density
    ## estimation etc.)
    for(j in 1:num.var) {
      bws.f <- npudensbw(formula(paste("~",bws.joint$xnames[j])),
                         bws=bws.joint$bw[j],
                         bandwidth.compute=FALSE,
                         bwmethod=bws.marginal$method,
                         bwtype=bws.marginal$type,
                         ckerorder=bws.marginal$ckerorder,
                         ckertype=bws.marginal$ckertype,
                         okertype=bws.marginal$okertype,
                         ukertype=bws.marginal$ukertype)
      xeval <- data.frame(x.u[,j])
      names(xeval) <- bws.joint$xnames[j]
      ## Divide copula density by its marginals
      copula.density <- copula.density/NZD(predict(npudens(bws=bws.f),newdata=xeval))
    }
  }
  ## Convert to data frame and name
  console <- printPop(console)
  console <- printClear(console)
  u <- data.frame(u)
  names(u) <- paste("u",1:num.var,sep="")
  return(data.frame(copula,copula.density,u))

}

