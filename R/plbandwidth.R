plbandwidth <-
  function(bws = stop("plbandwidth: bws missing"),
           regtype = c("lc","ll"),
           bwmethod = c("cv.ls","cv.aic"),
           bwscaling = FALSE,
           bwtype = c("fixed","generalized_nn","adaptive_nn"),
           ckertype = c("gaussian","truncated gaussian","epanechnikov","uniform"), 
           ckerorder = c(2,4,6,8),
           ukertype = c("aitchisonaitken", "liracine"),
           okertype = c("liracine","wangvanryzin"),
           xdati, ydati, zdati,
           xnames, ynames, znames,
           nobs = NA,
           rows.omit = NA,
           bandwidth.compute = TRUE,
           total.time = NA,...){

    regtype = match.arg(regtype)
    bwmethod = match.arg(bwmethod)
    bwtype = match.arg(bwtype)
    ckertype = match.arg(ckertype)

    if(missing(ckerorder))
      ckerorder = 2
    else if (ckertype == "uniform")
      warning("ignoring kernel order specified with uniform kernel type")
    else {
      kord = eval(formals()$ckerorder) 
      if (!any(kord == ckerorder))
        stop("ckerorder must be one of ", paste(kord,collapse=" "))
    }

    if (ckertype == "truncated gaussian" && ckerorder != 2)
      warning("using truncated gaussian of order 2, higher orders not yet implemented")

    ukertype = match.arg(ukertype)
    okertype = match.arg(okertype)

    porder = switch( ckerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order",
      "Eighth-Order" )
    ## chug chug
    sfactor <- lapply(1:length(bws), function(i) { unlist(bws[[i]]$sfactor) })
    bandwidth <- lapply(1:length(bws), function(i) { unlist(bws[[i]]$bandwidth) })
    sumNum <- lapply(1:length(bws), function(i) { unlist(bws[[i]]$sumNum) })

    names(sfactor) <- names(bandwidth) <- names(sumNum) <- rep("z", length(bws))

    if (length(rows.omit) == 0)
      rows.omit <- NA

    mybw = list(
      bw=bws,
      regtype = regtype,
      pregtype = switch(regtype,
        lc = "Local-Constant",
        ll = "Local-Linear"),
      method = bwmethod,
      pmethod = bwmToPrint(bwmethod),
      scaling = bwscaling,
      pscaling = ifelse(bwscaling, "Scale Factor(s)", "Bandwidth(s)"),
      type = bwtype,
      ptype = bwtToPrint(bwtype),
      ckertype = ckertype,    
      ckerorder = ckerorder,
      pckertype = cktToPrint(ckertype, order = porder),
      ukertype = ukertype,
      pukertype = uktToPrint(ukertype),
      okertype = okertype,
      pokertype = oktToPrint(okertype),
      nobs = nobs,
      zndim = bws$yzbw$ndim,
      xndim = length(bws)-1,
      xdati = xdati,
      ydati = ydati,
      zdati = zdati,
      xnames = xnames,
      ynames = ynames,
      znames = znames,
      fval = NA,
      sfactor = sfactor,
      bandwidth = bandwidth,
      sumNum = sumNum,
      varnames = list(x = xnames, y = ynames, z = znames),
      vartitle = list(x = "Explanatory", y = "Dependent", z = "Explanatory"),
      vartitleabb = list(x = "Exp.", y = "Dep.", z = "Exp."),
      dati = list(x = xdati, y = ydati, z = zdati),
      rows.omit = rows.omit,
      nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)),
      total.time = total.time)

    mybw$klist = list(z = list(ckertype = ckertype,
                        pckertype = mybw$pckertype,
                        ukertype = ukertype,
                        pukertype = mybw$pukertype,
                        okertype = okertype,
                        pokertype = mybw$pokertype))

    if(!bandwidth.compute)
      mybw$pmethod <- "Manual"

    class(mybw) = "plbandwidth"
    if(!any(is.na(mybw$bandwidth)))
      validateBandwidth(mybw)
    mybw
    
  }

print.plbandwidth <- function(x, digits=NULL, ...){
  cat("\nPartially Linear Model",
      "\nRegression data (",x$nobs," observations,\n",
      x$xndim," linear parametric variable(s) ",
      x$zndim," nonparametric variable(s):\n",sep="")

  bwmat = matrix(data = 0, nrow = x$xndim+1, ncol = x$bw$yzbw$ndim)
  
  for (i in 1:length(x$bw))
    bwmat[i,] = x$bw[[i]]$bw
  
  ## perhaps add a column for objective function value?
  print(matrix(bwmat[1,], ncol=x$zndim,
               dimnames=list(paste(x$pscaling,":",sep=""),
                 c("y(z)", replicate(x$zndim-1,"")))))
  cat("\n")
  print(matrix(bwmat[2:(1+x$xndim),], ncol=x$zndim,
               dimnames=list(c(paste(x$pscaling,":",sep=""), replicate(x$xndim-1,"")),
                 c("x(z)", replicate(x$zndim-1,"")))))

  cat(genBwSelStr(x))
  cat(genBwKerStrs(x))

  cat("\n\n")
  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

predict.plbandwidth <- function(...) { eval(npplreg(...), envir = parent.frame()) }
plot.plbandwidth <- function(...) { npplot(...) }

summary.plbandwidth <- function(object, ...){
  cat("\nPartially Linear Model",
      "\nRegression data (",object$nobs," observations,\n",
      object$xndim," linear parametric variable(s) ",
      object$zndim," nonparametric variable(s):\n",sep="")

  cat(genOmitStr(object))
  cat(genBwSelStr(object))

  scale <- genBwScaleStrs(object)
  nm <- unlist(object$varnames[c("y","x")])

  cat(paste("\n\n", nm, " on z:", scale, sep=""))
  
  cat(genBwKerStrs(object))
  cat(genTimingStr(object))
  cat("\n\n")  
}
