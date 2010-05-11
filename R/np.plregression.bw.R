npplregbw <-
  function(...){
    args = list(...)
    if (is(args[[1]],"formula"))
      UseMethod("npplregbw",args[[1]])
    else if (!is.null(args$formula))
      UseMethod("npplregbw",args$formula)
    else
      UseMethod("npplregbw",args[[which(names(args)=="bws")[1]]])
  }

npplregbw.formula <-
  function(formula, data, subset, na.action, call, ...){

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    if(!missing(call) && is.call(call)){
      ## rummage about in the call for the original formula
      for(i in 1:length(call)){
        if(tryCatch(class(eval(call[[i]])) == "formula",
                    error = function(e) FALSE))
          break;
      }
      mf[[2]] <- call[[i]]
      
      formula.args <- c("data", "subset", "na.action")
      mc.call <- match(formula.args, names(call), nomatch = 0)
      mc.mf <- match(formula.args, names(mf), nomatch = 0)
      if(any(mc.mf > 0))
        mf[mc.mf] <- call[mc.call]
    }

    mf[[1]] <- as.name("model.frame")

    ## mangle formula ...
    chromoly <- explodePipe(mf[["formula"]])

    if (length(chromoly) != 3) ## stop if malformed formula
      stop("invoked with improper formula, please see npplregbw documentation for proper use")

    ## make formula evaluable, then eval
    bronze <- lapply(chromoly, paste, collapse = " + ")
    mf[["formula"]] <- as.formula(paste(bronze[[1]]," ~ ", bronze[[2]], " + ", bronze[[3]]),
                                  env = environment(formula))
    mf <- eval(mf, parent.frame())
    
    ydat <- model.response(mf)
    xdat <- mf[, chromoly[[2]], drop = FALSE]
    zdat <- mf[, chromoly[[3]], drop = FALSE]
    
    tbw <- npplregbw(xdat = xdat, ydat = ydat, zdat = zdat, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")
    tbw$chromoly <- chromoly

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)

    tbw
  }


npplregbw.NULL =
  function(xdat = stop("invoked without data `xdat'"),
           ydat = stop("invoked without data `ydat'"),
           zdat = stop("invoked without data `zdat'"),
           bws, ...){

    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain z names and 'toFrame'
    zdat <- toFrame(zdat)

    ## bandwidths

    bws = matrix(data = 0, nrow = 1+ncol(xdat), ncol = ncol(zdat))

    tbw <- npplregbw.default(xdat = xdat, ydat = ydat, zdat = zdat,
                             bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <-
      updateBwNameMetadata(nameList = list(ynames = deparse(substitute(ydat))),
                           bws = tbw)

    tbw
  }

npplregbw.plbandwidth =
  function(xdat = stop("invoked without data `xdat'"),
           ydat = stop("invoked without data `ydat'"),
           zdat = stop("invoked without data `zdat'"),
           bws, nmulti, ...){

    ## n observations
    ## rows are observations, columns are variables
    ## x is n x k
    ## y is n x 1
    ## z is n x p

    if (missing(nmulti)){
      nmulti <- min(5,dim(zdat)[2])
    }
    
    xdat = toFrame(xdat)
    zdat = toFrame(zdat)

    ## catch and destroy NA's
    goodrows = 1:dim(xdat)[1]
    rows.omit = na.action(na.omit(data.frame(xdat,ydat,zdat)))
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Training data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows]
    zdat = zdat[goodrows,,drop = FALSE]

    ## y on z
    bws$bw$yzbw = npregbw(xdat = zdat, ydat = ydat,
      bws = bws$bw$yzbw, nmulti = nmulti, ...)
    
    ## x on z

    for (i in 1:ncol(xdat)){
      bws$bw[[i+1]] = npregbw(xdat=zdat, ydat=xdat[,i],
              bws = bws$bw[[i+1]], nmulti = nmulti, ...)
    }

    bws <- plbandwidth(bw = bws$bw,
                       regtype = bws$regtype,
                       bwmethod = bws$method,
                       bwscaling = bws$scaling,
                       bwtype = bws$type,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       ukertype = bws$ukertype,
                       okertype = bws$okertype,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       zdati = bws$zdati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       znames = bws$znames,
                       nobs = bws$nobs,
                       rows.omit = rows.omit)

  }

npplregbw.default = 
  function(xdat = stop("invoked without data `xdat'"),
           ydat = stop("invoked without data `ydat'"),
           zdat = stop("invoked without data `zdat'"),
           bws, ...,
           bandwidth.compute = TRUE,
           nmulti, remin, itmax,
           ftol, tol, small){

    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain z names and 'toFrame'
    zdat <- toFrame(zdat)

    plband = list()
    plband$yzbw = npregbw(xdat = zdat, ydat = ydat,
      bws = bws[1,], ..., bandwidth.compute = FALSE)

    for (i in 1:dim(xdat)[2])
      plband[[i+1]] = npregbw(xdat = zdat, ydat = xdat[,i],
              bws = bws[i+1,], ..., bandwidth.compute = FALSE)

    tbw <- plbandwidth(bws = plband,
                       nobs = dim(xdat)[1],
                       ...,
                       xdati = untangle(xdat),
                       ydati = untangle(data.frame(ydat)),
                       zdati = untangle(zdat),
                       xnames = names(xdat),
                       ynames = deparse(substitute(ydat)),
                       znames = names(zdat),
                       bandwidth.compute = bandwidth.compute)


    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("nmulti", "remin", "itmax", "ftol", "tol", "small")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (bandwidth.compute)
      tbw <- eval(parse(text=paste("npplregbw.plbandwidth(xdat=xdat, ydat=ydat, zdat = zdat, bws=tbw",
                          ifelse(any.m, ",",""),
                          paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                          ")")))

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }
