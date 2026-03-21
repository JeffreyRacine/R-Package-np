## This function accepts a joint distribution bandwidth object,
## optionally a list of marginal distribution bandwidth objects, and
## optionally a matrix with columns being vectors of taus (u) for
## which the copula is to be computed. It returns a data frame with
## the copula `copula', the copula density `copula density' and u
## (columns for the grid constructed from the columns of the u matrix)

npcopula <- function(bws,
                     data,
                     u=NULL,
                     n.quasi.inv=1000,
                     er.quasi.inv=1) {
  .npRmpi_require_active_slave_pool(where = "npcopula()")
  # Keep npcopula as a local orchestrator in session mode.
  # Dispatching the full routine can deadlock because it nests many estimator
  # calls; those inner calls are already MPI-aware and dispatch independently.

  ## Basic error checking

  if(missing(data)) stop("You must provide a data frame")
  if(!is.data.frame(data)) stop("Object `data' must be a data frame")
  if(missing(bws)) stop("You must provide a bandwidth object")
#  if(class(bws)!="dbandwidth"&class(bws)!="bandwidth") stop("you must provide a density (npudensbw) or distribution (npudistbw) object")
  if(!isa(bws,"dbandwidth") && !isa(bws,"bandwidth")) stop("you must provide a density (npudensbw) or distribution (npudistbw) object")  
  density <- FALSE
#  if(!missing(bws)&&class(bws)=="bandwidth") density <- TRUE
  if(!missing(bws) && isa(bws,"bandwidth")) density <- TRUE
  if(!is.null(u)) {
    if(anyNA(u)) stop("u must not contain missing values")
    if(any(u>1 | u<0, na.rm = TRUE)) stop("u must lie in [0,1]")
  }
  num.var <- length(bws$xnames)
  if(!is.null(u) && (ncol(u)!=num.var)) stop("u and bws are incompatible")
  if(n.quasi.inv < 1) stop("n.quasi.inv must be a positive integer")
  if(any(is.na(data))) stop("NA values present in data - recompute bw object with na.omit on data")

  ## Check for unordered factors

  if(bws$nuno>0) stop("unordered factors not suitable for copula estimation")

  u.provided <- !is.null(u)
  bw.type <- bws$type
  bw.ckerorder <- bws$ckerorder
  bw.ckertype <- bws$ckertype
  bw.okertype <- bws$okertype
  bw.ukertype <- bws$ukertype
  
  ## Test for compatible quantile vector/matrix if provided
  if(!is.null(u)) {
    if(is.vector(u)) u <- matrix(u,1,length(u))
    if(ncol(u) != num.var) stop(paste("matrix u must have ", num.var," columns",sep=""))
  }

  if(is.null(u)) {
    ## Compute the copula distribution or copula density for the
    ## sample realizations (joint CDF)
    if(!density) {
      .np_progress_note("Computing the copula for the sample realizations")
      copula <- fitted(npudist(bws=bws,data=data))
    } else {
      .np_progress_note("Computing the copula density for the sample realizations")
      copula <- fitted(npudens(bws=bws,data=data))
    }
    ## Compute the marginal quantiles from the marginal CDFs (u_i=\hat
    ## F(x_i))
    u <- matrix(NA,bws$nobs,num.var)
    for (j in seq_len(num.var)) {
      .np_progress_note(
        sprintf(
          "Computing the marginal of %s for the sample realizations",
          bws$xnames[j]
        )
      )
      bw.j <- bws$bw[j]
      bws.F.marginal <- do.call(npudistbw, list(
                         dat = data[, bws$xnames[j], drop = FALSE],
                         bws = bw.j,
                         bandwidth.compute = FALSE,
                         bwtype = bw.type,
                         ckerorder = bw.ckerorder,
                         ckertype = bw.ckertype,
                         okertype = bw.okertype))

      u[,j] <- fitted(npudist(bws=bws.F.marginal,data=data))
      ## For the copula density we require marginal densities.
      if(density) {
        bws.f.marginal <- do.call(npudensbw, list(
                           dat = data[, bws$xnames[j], drop = FALSE],
                           bws = bw.j,
                           bandwidth.compute = FALSE,
                           bwtype = bw.type,
                           ckerorder = bw.ckerorder,
                           ckertype = bw.ckertype,
                           okertype = bw.okertype,
                           ukertype = bw.ukertype))
        ## Divide the copula density by its marginals
        copula <- copula/NZD(fitted(npudens(bws=bws.f.marginal,data=data)))

      }
    }
  } else {
    ## Here the user wishes to compute copula for the input u
    ## matrix. To economize on computation we obtain the quantiles
    ## from the marginals then call expand.grid to compute the copula
    n.u <- nrow(u)
    x.u <- data.frame(matrix(NA,n.u,num.var))
    names(x.u) <- bws$xnames
    for (j in seq_len(num.var)) {
      .np_progress_note(
        sprintf(
          "Computing the quasi-inverse for the marginal of %s",
          bws$xnames[j]
        )
      )
      ## Compute the quasi-inverse (Definition 2.3.6, Nelson
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
      x.marginal <- data[[bws$xnames[j]]]
      quantile.seq <- seq(0,1,length=round(n.quasi.inv/2))
      if(is.numeric(x.marginal)) {
        x.er <- extendrange(x.marginal,f=er.quasi.inv)
        x.q <- quantile(x.marginal,quantile.seq)
        x.eval <- sort(c(seq(x.er[1],x.er[2],length=round(n.quasi.inv/2)),x.q))
      } else {
        x.u[,j] <- ordered(x.u[,j],levels=levels(x.marginal))
        x.q <- sapply(seq_len(round(n.quasi.inv/2)), function(i) { uocquantile(x.marginal, quantile.seq[i]) })
        x.eval <- sort(ordered(c(as.character(x.q),as.character(x.q)),levels=levels(x.marginal)))
      }
      ## Compute the CDF at this set of evaluation points.
      bw.j <- bws$bw[j]
      bws.F.marginal <- do.call(npudistbw, list(
                                  dat = x.marginal,
                                  bws = bw.j,
                                  bandwidth.compute = FALSE,
                                  bwtype = bw.type,
                                  ckerorder = bw.ckerorder,
                                  ckertype = bw.ckertype,
                                  okertype = bw.okertype,
                                  ukertype = bw.ukertype))
      F <- fitted(npudist(tdat = x.marginal,
                          edat = x.eval,
                          bws = bws.F.marginal))
      ## Now compute the quasi-inverse from the estimated F for the
      ## evaluation points. If u is input and any value lies beyond
      ## the CDF values for the evaluation points, reset them to the
      ## min/max CDF values for the evaluation data (otherwise the
      ## quantiles are undefined).
      for (i in seq_len(n.u)) {
        u[u[,j]<min(F),j] <- min(F)
        u[u[,j]>max(F),j] <- max(F)        
        x.u[i,j] <-  min(x.eval[F>=u[i,j]])
      }
    }
    ## To compute the copula we expand the grid of marginal quantiles
    ## so that every combination of the columns of u is constructed.
    .np_progress_note("Expanding the u matrix")
    x.u <- expand.grid(data.frame(x.u))
    for (k in seq_len(ncol(x.u))) {
      if(is.ordered(data[,k])) x.u[,k] <- ordered(x.u[,k],levels=levels(data[,k]))
    }
    unit.weights <- rep_len(1, nrow(data))
    if(!density) {
      .np_progress_note("Computing the copula for the expanded grid")
      copula <- npudisthat(
        bws = bws,
        tdat = data,
        edat = x.u,
        y = unit.weights,
        output = "apply"
      )
    } else {
      .np_progress_note("Computing the copula density for the expanded grid")
      copula <- npudenshat(
        bws = bws,
        tdat = data,
        edat = x.u,
        y = unit.weights,
        output = "apply"
      )
      ## For the copula density require marginal densities. Desirable to
      ## have the same bws in numerator and denominator, so use those
      ## from the joint (mirror regression, conditional density
      ## estimation etc.)
      for (j in seq_len(num.var)) {
        .np_progress_note(
          sprintf(
            "Computing the marginal of %s for the expanded grid",
            bws$xnames[j]
          )
        )
        bw.j <- bws$bw[j]
        bws.f.marginal <- do.call(npudensbw, list(
                           dat = data[, bws$xnames[j], drop = FALSE],
                           bws = bw.j,
                           bandwidth.compute = FALSE,
                           bwtype = bw.type,
                           ckerorder = bw.ckerorder,
                           ckertype = bw.ckertype,
                           okertype = bw.okertype,
                           ukertype = bw.ukertype))
        xeval <- data.frame(x.u[,j])
        names(xeval) <- bws$xnames[j]
        ## Divide copula density by its marginals
        copula <- copula/NZD(npudenshat(
          bws = bws.f.marginal,
          tdat = data[, bws$xnames[j], drop = FALSE],
          edat = xeval,
          y = unit.weights,
          output = "apply"
        ))
      }
    }
  }

  if(!u.provided) {
    u <- data.frame(u)
    names(u) <- paste("u", seq_len(num.var), sep = "")
    return(data.frame(copula,u))
  } else {
    ## If u was provided we expand its grid as was done for the
    ## marginals
    u <- expand.grid(data.frame(u))
    names(u) <- paste("u", seq_len(num.var), sep = "")
    return(data.frame(copula,u,x.u))
  }

}
