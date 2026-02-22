npplregbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    dots <- mc$...

    if (length(dots) == 0L)
      stop("invoked without arguments")

    dot.names <- names(dots)

    if (!is.null(dot.names) && any(dot.names == "formula")) {
      formula.val <- eval(dots[[which(dot.names == "formula")[1L]]], envir = parent.frame())
      return(UseMethod("npplregbw", formula.val))
    }

    first.val <- eval(dots[[1L]], envir = parent.frame())
    if (inherits(first.val, "formula"))
      return(UseMethod("npplregbw", first.val))

    if (!is.null(dot.names) && any(dot.names == "bws")) {
      bws.val <- eval(dots[[which(dot.names == "bws")[1L]]], envir = parent.frame())
      return(UseMethod("npplregbw", bws.val))
    }

    if (!is.null(dot.names) && any(dot.names %in% c("xdat", "ydat", "zdat")))
      return(UseMethod("npplregbw", NULL))

    UseMethod("npplregbw", first.val)
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
        if(tryCatch(inherits(eval(call[[i]]), "formula"),
                    error = function(e) FALSE))
          break;
      }
      mf[[2]] <- call[[i]]
      
    }

    mf.xf <- mf
    
    mf[[1]] <- as.name("model.frame")
    mf.xf[[1]] <- as.name("model.frame")
    
    ## mangle formula ...
    formula_to_explode <- eval(mf[["formula"]], parent.frame())
    chromoly <- explodePipe(formula_to_explode, env = environment(formula))

    if (length(chromoly) != 3) ## stop if malformed formula
      stop("invoked with improper formula, please see npplregbw documentation for proper use")

    ## make formula evaluable, then eval
    bronze <- lapply(chromoly, paste, collapse = " + ")

    mf.xf[["formula"]] <- as.formula(paste(" ~ ", bronze[[2]]),
                                    env = environment(formula))

    mf[["formula"]] <- as.formula(paste(bronze[[1]]," ~ ", bronze[[3]]),
                                  env = environment(formula))

    formula.all <- if(missing(data)) {
        terms(as.formula(paste(" ~ ",bronze[[1]]," + ",bronze[[2]], " + ",bronze[[3]]),
                                  env = environment(formula)))
    } else {
        terms(as.formula(paste(" ~ ",bronze[[1]]," + ",bronze[[2]], " + ",bronze[[3]]),
                                  env = environment(formula)), data = data)
    }

    orig.ts <- if (missing(data))
      sapply(eval(attr(formula.all, "variables"), environment(formula.all)), inherits, "ts")
    else sapply(eval(attr(formula.all, "variables"), data, environment(formula.all)), inherits, "ts")

    arguments.mfx <- chromoly[[2]]
    arguments.mf <- c(chromoly[[1]],chromoly[[3]])

    mf[["formula"]] <- terms(mf[["formula"]])
    mf.xf[["formula"]] <- terms(mf.xf[["formula"]])
    
    if(all(orig.ts)){
      arguments <- (as.list(attr(formula.all, "variables"))[-1])
      attr(mf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments)))))[,.(match(arguments.mf,arguments)),drop = FALSE])
      attr(mf.xf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments)))))[,.(match(arguments.mfx,arguments)),drop = FALSE])
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(formula.all, "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      attr(mf[["formula"]], "predvars") <- bquote((.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])[,.(match(arguments.mf,arguments)),drop = FALSE])
      attr(mf.xf[["formula"]], "predvars") <- bquote((.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])[,.(match(arguments.mfx,arguments)),drop = FALSE])
    }
    
    mf <- eval(mf, parent.frame())
    mf.xf <- eval(mf.xf,parent.frame())

    ydat <- model.response(mf)
    xdat <- mf.xf
    zdat <- mf[, chromoly[[3]], drop = FALSE]

    tbw <- npplregbw(xdat = xdat, ydat = ydat, zdat = zdat, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")
    tbw$xterms <- attr(mf.xf,"terms")
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
    total.time <-
      system.time({
        bws$bw$yzbw  <- npregbw(xdat = zdat, ydat = ydat,
                                bws = bws$bw$yzbw, nmulti = nmulti, ...)
        
        ## x on z

        for (i in 1:ncol(xdat)){
          bws$bw[[i+1]] <- npregbw(xdat=zdat, ydat=xdat[,i],
                  bws = bws$bw[[i+1]], nmulti = nmulti, ...)
        }
      })[1]
    num.feval <- sum(sapply(bws$bw, function(bwi) {
      if (is.null(bwi$num.feval) || identical(bwi$num.feval, NA)) 0 else bwi$num.feval
    }))
    fval <- {
      fv <- unlist(lapply(bws$bw, function(bwi) bwi$fval))
      if (length(fv) == 0 || all(!is.finite(fv))) NA else sum(fv[is.finite(fv)])
    }
    bws <- plbandwidth(bws = bws$bw,
                       regtype = bws$regtype,
                       bwmethod = bws$method,
                       bwscaling = bws$scaling,
                       bwtype = bws$type,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       ckerbound = bws$ckerbound,
                       ckerlb = bws$ckerlb,
                       ckerub = bws$ckerub,
                       ukertype = bws$ukertype,
                       okertype = bws$okertype,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       zdati = bws$zdati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       znames = bws$znames,
                       nobs = bws$nobs,
                       fval = fval,
                       num.feval = num.feval,
                       rows.omit = rows.omit,
                       total.time = total.time)

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
                       fval = {
                         fv <- unlist(lapply(plband, function(bwi) bwi$fval))
                         if (length(fv) == 0 || all(!is.finite(fv))) NA else sum(fv[is.finite(fv)])
                       },
                       num.feval = sum(sapply(plband, function(bwi) {
                         if (is.null(bwi$num.feval) || identical(bwi$num.feval, NA)) 0 else bwi$num.feval
                       })),
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

    if (bandwidth.compute) {
      bwsel.args <- list(xdat = xdat, ydat = ydat, zdat = zdat, bws = tbw)
      if (any.m) {
        nms <- mc.names[m]
        bwsel.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
      }
      tbw <- do.call(npplregbw.plbandwidth, bwsel.args)
    }

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }
