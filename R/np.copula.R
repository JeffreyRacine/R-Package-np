## This function accepts a joint distribution bandwidth object,
## optionally a list of marginal distribution bandwidth objects, and
## optionally a matrix with columns being vectors of taus (u) for
## which the copula is to be computed. It returns a data frame with
## the copula `copula', the copula density `copula density' and u
## (columns for the grid constructed from the columns of the u matrix)

npcopula <- function(bws,data,u=NULL,density=FALSE) {

  ## Basic error checking

  if(missing(data)) stop("You must provide a data frame")
  if(!is.data.frame(data)) stop("Object `data' must be a data frame")
  if(!density&&missing(bws)) stop("You must provide a joint distribution bandwidth object")
  if(density&&missing(bws)) stop("You must provide a joint density bandwidth object")
  if(!is.null(u)) if(any(u>1) || any(u<0)) stop("u must lie in [0,1]")
  num.var <- length(bws$xnames)
  
  ## Test for compatible quantile vector/matrix if provided
  if(!is.null(u)) {
    if(is.vector(u)) u <- matrix(u,1,length(u))
    if(ncol(u) != num.var) stop(paste("matrix u must have ", num.var," columns",sep=""))
  }

  console <- newLineConsole()

  if(is.null(u)) {
    ## Compute the copula distribution and density for the sample
    ## realizations (joint CDF)
    console <- printClear(console)
    if(!density) {
      console <- printPush(msg = "Computing the copula for the sample realizations...", console)
      copula <- fitted(npudist(bws=bws,data=data))
    } else {
      console <- printPush(msg = "Computing the copula density for the sample realizations...", console)
      copula <- fitted(npudens(bws=bws,data=data))
    }
    ## Compute the marginal quantiles from the marginal CDFs (u_i=\hat
    ## F(x_i)) and divide copula density by its marginals
    u <- matrix(NA,bws$nobs,num.var)
    for(j in 1:num.var) {
      console <- printClear(console)
      console <- printPush(msg = paste("Computing the marginal of ",bws$xnames[j]," for the sample realizations...",sep=""), console)
      bws.F <- npudensbw(formula(paste("~",bws$xnames[j])),
                         bws=bws$bw[j],
                         bandwidth.compute=FALSE,
                         bwmethod=bws$method,
                         bwtype=bws$type,
                         ckerorder=bws$ckerorder,
                         ckertype=bws$ckertype,
                         okertype=bws$okertype,
                         ukertype=bws$ukertype,
                         data=data)

      u[,j] <- fitted(npudist(bws=bws.F,data=data))
      ## For copula density we require marginal densities. Desirable
      ## to have the same bws in numerator and denominator, so use
      ## those from the joint (mirror regression, conditional density
      ## estimation etc.)
      if(density) {
        bws.f <- npudensbw(formula(paste("~",bws$xnames[j])),
                           bws=bws$bw[j],
                           bandwidth.compute=FALSE,
                           bwmethod=bws$method,
                           bwtype=bws$type,
                           ckerorder=bws$ckerorder,
                           ckertype=bws$ckertype,
                           okertype=bws$okertype,
                           ukertype=bws$ukertype,
                           data=data)
        copula <- copula/NZD(fitted(npudens(bws=bws.f,data=data)))
      }
    }
  } else {
    ## User wishes to compute copula for inputted u matrix. To
    ## economize on computation we obtain the quantiles from the
    ## marginals then expand to compute the copula
    n.u <- nrow(u)
    x.u <- matrix(NA,n.u,num.var)
    for(j in 1:num.var) {
      console <- printClear(console)
      console <- printPush(msg = paste("Computing the pseudo-inverse for the marginal of ",bws$xnames[j],"...",sep=""), console)
      ## Compute the quasi inverse (Definition 2.3.6, Nelson
      ## (2006)).  Here we take pains to span a sufficiently rich
      ## set of evaluation points to cover a range of
      ## possibilities. In particular, we extend the range of the
      ## variable by a fraction 1 on min/max and create an equally
      ## spaced grid on this range to try to cover long-tailed
      ## distributions but provide a sufficiently fine grid on this
      ## extended range (uniform spacing, add and subtract the range
      ## of the data to/from max/min). We also use a sequence of
      ## equi-quantile spaced points from the quantiles of the raw
      ## data again to provide a sufficiently fine grid.  We then
      ## concatenate and sort the equally space extended grid and
      ## the equi-quantile grid.
      x.marginal <- eval(parse(text=paste("data$",bws$xnames[j],sep="")))
      x.er <- extendrange(x.marginal,f=1)
      x.q <- quantile(x.marginal,seq(0,1,length=500))
      x.eval <- sort(c(seq(x.er[1],x.er[2],length=500),x.q))
      ## Compute the CDF at this set of evaluation points
      F <- fitted(npudist(tdat=x.marginal,
                          edat=x.eval,
                          bws=bws$bw[j],
                          bwmethod=bws$method,
                          bwtype=bws$type,
                          ckerorder=bws$ckerorder,
                          ckertype=bws$ckertype,
                          okertype=bws$okertype,
                          ukertype=bws$ukertype,data=data))
      ## Now compute the psuedo inverse from the estimated F for the
      ## evaluation points
      for(i in 1:n.u) {
        x.u[i,j] <- ifelse(u[i,j]>=0.5, max(x.eval[F<=u[i,j]]), min(x.eval[F>=u[i,j]]))
      }
    }
    ## Convert to evaluation data frame and name so that we can use
    ## predict(,newdata=...). To compute the copula we expand the grid
    ## of marginal quantiles so that every combination of the columns
    ## of u is constructed
    console <- printClear(console)
    console <- printPush(msg = "Expanding the u matrix...", console)
    x.u <- expand.grid(data.frame(x.u))
    names(x.u) <- bws$xnames
    console <- printClear(console)
    if(!density) {
      console <- printPush(msg = "Computing the copula for the expanded grid...", console)
      copula <- predict(npudist(bws=bws),data=data,newdata=x.u)
    } else {
      console <- printPush(msg = "Computing the copula density for the expanded grid...", console)
      copula <- predict(npudens(bws=bws),data=data,newdata=x.u)
      ## For the copula density require marginal densities. Desirable to
      ## have the same bws in numerator and denominator, so use those
      ## from the joint (mirror regression, conditional density
      ## estimation etc.)
      for(j in 1:num.var) {
        console <- printClear(console)
        console <- printPush(msg = paste("Computing the marginal of ",bws$xnames[j]," for the expanded grid...",sep=""), console)
        bws.f <- npudensbw(formula(paste("~",bws$xnames[j])),
                           bws=bws$bw[j],
                           bandwidth.compute=FALSE,
                           bwmethod=bws$method,
                           bwtype=bws$type,
                           ckerorder=bws$ckerorder,
                           ckertype=bws$ckertype,
                           okertype=bws$okertype,
                           ukertype=bws$ukertype,
                           data=data)
        xeval <- data.frame(x.u[,j])
        names(xeval) <- bws$xnames[j]
        ## Divide copula density by its marginals
        copula <- copula/NZD(predict(npudens(bws=bws.f,data=data),newdata=xeval))
      }
    }
  }
  ## Convert to data frame and name
  console <- printClear(console)
  u <- data.frame(u)
  names(u) <- paste("u",1:num.var,sep="")
  return(data.frame(copula,u))

}

