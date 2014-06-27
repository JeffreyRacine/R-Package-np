npcdensbw <-
  function(...){
    args = list(...)
    if (is(args[[1]],"formula"))
      UseMethod("npcdensbw",args[[1]])
    else if (!is.null(args$formula))
      UseMethod("npcdensbw",args$formula)
    else
      UseMethod("npcdensbw",args[[which(names(args)=="bws")[1]]])
  }

npcdensbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    orig.class <- if (missing(data))
      sapply(eval(attr(terms(formula), "variables"), environment(formula)),class)
    else sapply(eval(attr(terms(formula), "variables"), data, environment(formula)),class)

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
      
    }
                     

    mf[[1]] <- as.name("model.frame")

    variableNames <- explodeFormula(mf[["formula"]])
    
    ## make formula evaluable, then eval
    varsPlus <- lapply(variableNames, paste, collapse=" + ")
    mf[["formula"]] <- as.formula(paste(" ~ ", varsPlus[[1]]," + ",
                                        varsPlus[[2]]),
                                  env = environment(formula))
    mf[["formula"]] <- terms(mf[["formula"]])
    if(all(orig.class == "ts")){
      args <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      attr(mf[["formula"]], "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
    }else if(any(orig.class == "ts")){
      arguments <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      arguments.normal <- arguments[which(orig.class != "ts")]
      arguments.timeseries <- arguments[which(orig.class == "ts")]

      ix <- sort(c(which(orig.class == "ts"),which(orig.class != "ts")),index.return = TRUE)$ix
      attr(mf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
    }

    mf <- eval(mf, parent.frame())

    ydat <- mf[, variableNames[[1]], drop = FALSE]
    xdat <- mf[, variableNames[[2]], drop = FALSE]
    
    tbw = npcdensbw(xdat = xdat, ydat = ydat, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")
    tbw$variableNames <- variableNames

    tbw
  }

npcdensbw.conbandwidth <- 
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           bws, bandwidth.compute = TRUE,
           nmulti, remin = TRUE, itmax = 10000, 
           ftol = 1.490116e-07, tol = 1.490116e-04, small = 1.490116e-05,
           memfac = 500.0, lbc.dir = 0.5, dfc.dir = 3, cfac.dir = 2.5*(3.0-sqrt(5)), initc.dir = 1.0, 
           lbd.dir = 0.1, hbd.dir = 1, dfac.dir = 0.25*(3.0-sqrt(5)), initd.dir = 1.0, 
           lbc.init = 0.1, hbc.init = 2.0, cfac.init = 0.5, 
           lbd.init = 0.1, hbd.init = 0.9, dfac.init = 0.375, 
           scale.init.categorical.sample=FALSE,
           ...){

    ydat = toFrame(ydat)
    xdat = toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- min(5,(dim(ydat)[2]+dim(xdat)[2]))
    }

    if (length(bws$ybw) != dim(ydat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of", "'ydat'"))

    if (length(bws$xbw) != dim(xdat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of", "'xdat'"))

    if (dim(ydat)[1] != dim(xdat)[1])
      stop(paste("number of rows of", "'ydat'", "does not match", "'xdat'"))

    yccon = unlist(lapply(as.data.frame(ydat[,bws$iycon]),class))
    if ((any(bws$iycon) && !all((yccon == class(integer(0))) | (yccon == class(numeric(0))))) ||
        (any(bws$iyord) && !all(unlist(lapply(as.data.frame(ydat[,bws$iyord]),class)) ==
                               class(ordered(0)))) ||
        (any(bws$iyuno) && !all(unlist(lapply(as.data.frame(ydat[,bws$iyuno]),class)) ==
                               class(factor(0)))))
      stop(paste("supplied bandwidths do not match", "'ydat'", "in type"))

    xccon = unlist(lapply(as.data.frame(xdat[,bws$ixcon]),class))
    if ((any(bws$ixcon) && !all((xccon == class(integer(0))) | (xccon == class(numeric(0))))) ||
        (any(bws$ixord) && !all(unlist(lapply(as.data.frame(xdat[,bws$ixord]),class)) ==
                               class(ordered(0)))) ||
        (any(bws$ixuno) && !all(unlist(lapply(as.data.frame(xdat[,bws$ixuno]),class)) ==
                               class(factor(0)))))
      stop(paste("supplied bandwidths do not match", "'xdat'", "in type"))

    ## catch and destroy NA's
    goodrows <- 1:dim(xdat)[1]
    rows.omit <- unclass(na.action(na.omit(data.frame(xdat,ydat))))
    goodrows[rows.omit] <- 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows,,drop = FALSE]

    
    nrow = nrow(ydat)
    yncol = ncol(ydat)
    xncol = ncol(xdat)

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.
    
    ydat = toMatrix(ydat)

    yuno = ydat[, bws$iyuno, drop = FALSE]
    ycon = ydat[, bws$iycon, drop = FALSE]
    yord = ydat[, bws$iyord, drop = FALSE]


    xdat = toMatrix(xdat)

    xuno = xdat[, bws$ixuno, drop = FALSE]
    xcon = xdat[, bws$ixcon, drop = FALSE]
    xord = xdat[, bws$ixord, drop = FALSE]

    tbw <- bws

    mysd <- EssDee(data.frame(xcon,ycon))
    nconfac <- nrow^(-1.0/(2.0*bws$cxkerorder+bws$ncon))
    ncatfac <- nrow^(-2.0/(2.0*bws$cxkerorder+bws$ncon))

    if (bandwidth.compute){
      myopti = list(num_obs_train = nrow,
        iMultistart = ifelse(nmulti==0,IMULTI_FALSE,IMULTI_TRUE),
        iNum_Multistart = nmulti,
        int_use_starting_values = ifelse(all(bws$ybw==0) && all(bws$xbw==0),
          USE_START_NO, USE_START_YES),
        int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
        BANDWIDTH_den_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=ifelse(remin,RE_MIN_TRUE,RE_MIN_FALSE), 
        int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE), 
        bwmethod = switch(bws$method,
          cv.ml = CBWM_CVML,
          cv.ls = CBWM_CVLS,
          cv.ls.np = CBWM_NPLS),        
        xkerneval = switch(bws$cxkertype,
          gaussian = CKER_GAUSS + bws$cxkerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$cxkerorder/2 - 1,
          uniform = CKER_UNI,
          "truncated gaussian" = CKER_TGAUSS),
        ykerneval = switch(bws$cykertype,
          gaussian = CKER_GAUSS + bws$cykerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$cykerorder/2 - 1,
          uniform = CKER_UNI,
          "truncated gaussian" = CKER_TGAUSS),
        uxkerneval = switch(bws$uxkertype,
          aitchisonaitken = UKER_AIT,
          liracine = UKER_LR),
        uykerneval = switch(bws$uykertype,
          aitchisonaitken = UKER_AIT,
          liracine = UKER_LR),
        oxkerneval = switch(bws$oxkertype,
          wangvanryzin = OKER_WANG,
          liracine = OKER_LR),
        oykerneval = switch(bws$oykertype,
          wangvanryzin = OKER_WANG,
          liracine = OKER_NLR),
        ynuno = dim(yuno)[2],
        ynord = dim(yord)[2],
        yncon = dim(ycon)[2],
        xnuno = dim(xuno)[2],
        xnord = dim(xord)[2],
        xncon = dim(xcon)[2],
        fast = FALSE,
        old.cdens = FALSE,
        int_do_tree = ifelse(options('np.tree'), DO_TREE_YES, DO_TREE_NO),
        scale.init.categorical.sample = scale.init.categorical.sample,
        dfc.dir = dfc.dir)
      
      myoptd = list(ftol=ftol, tol=tol, small=small, memfac = memfac,
        lbc.dir = lbc.dir, cfac.dir = cfac.dir, initc.dir = initc.dir, 
        lbd.dir = lbd.dir, hbd.dir = hbd.dir, dfac.dir = dfac.dir, initd.dir = initd.dir, 
        lbc.init = lbc.init, hbc.init = hbc.init, cfac.init = cfac.init, 
        lbd.init = lbd.init, hbd.init = hbd.init, dfac.init = dfac.init, 
        nconfac = nconfac, ncatfac = ncatfac)

      if (bws$method != "normal-reference"){
        total.time <-
          system.time(myout <- 
          .C("np_density_conditional_bw", as.double(yuno), as.double(yord), as.double(ycon),
             as.double(xuno), as.double(xord), as.double(xcon),
             as.double(mysd),
             as.integer(myopti), as.double(myoptd), 
             bw = c(bws$xbw[bws$ixcon],bws$ybw[bws$iycon],
               bws$ybw[bws$iyuno],bws$ybw[bws$iyord],
               bws$xbw[bws$ixuno],bws$xbw[bws$ixord]),
             fval = double(2), fval.history = double(max(1,nmulti)),
             timing = double(1),
             PACKAGE="np" )[c("bw","fval","fval.history","timing")])[1]
      } else {
        nbw = double(yncol+xncol)
        gbw = bws$yncon+bws$xncon
        if (gbw > 0){
          nbw[1:gbw] = 1.059224
          if(!bws$scaling)
            nbw[1:gbw]=nbw[1:gbw]*mysd*nconfac
        }
        myout= list( bw = nbw, fval = c(NA,NA) )
        total.time <- NA
      }

      yr = 1:yncol
      xr = 1:xncol
      rorder = numeric(yncol + xncol)

      ## bandwidths are passed back from the C routine in an unusual order
      ## xc, y[cuo], x[uo]
      
      rxcon = xr[bws$ixcon]
      rxuno = xr[bws$ixuno] 
      rxord = xr[bws$ixord] 

      rycon = yr[bws$iycon] 
      ryuno = yr[bws$iyuno] 
      ryord = yr[bws$iyord] 


      ## rorder[c(rxcon,rycon,ryuno,ryord,rxuno,rxord)]=1:(yncol+xncol)

      tbw <- bws
      tbw$ybw[c(rycon,ryuno,ryord)] <- myout$bw[yr+bws$xncon]
      tbw$xbw[c(rxcon,rxuno,rxord)] <- myout$bw[setdiff(1:(yncol+xncol),yr+bws$xncon)]

      tbw$fval = myout$fval[1]
      tbw$ifval = myout$fval[2]
      tbw$fval.history <- myout$fval.history
      tbw$timing <- myout$timing
      tbw$total.time <- total.time
    }
    
    ## bandwidth metadata
    tbw$sfactor <- tbw$bandwidth <- list(x = tbw$xbw, y = tbw$ybw)

    bwf <- function(i){
      tbw$bandwidth[[i]][tl[[i]]] <<- (tbw$bandwidth[[i]])[tl[[i]]]*dfactor[[i]]
    }

    sff <- function(i){
      tbw$sfactor[[i]][tl[[i]]] <<- (tbw$sfactor[[i]])[tl[[i]]]/dfactor[[i]]
    }

    myf <- if(tbw$scaling) bwf else sff
    
    if ((tbw$xnuno+tbw$ynuno) > 0){
      dfactor <- ncatfac
      dfactor <- list(x = dfactor, y = dfactor)

      tl <- list(x = tbw$xdati$iuno, y = tbw$ydati$iuno)

      lapply(1:length(tl), myf)
    }

    if ((tbw$xnord+tbw$ynord) > 0){
      dfactor <- ncatfac
      dfactor <- list(x = dfactor, y = dfactor)

      tl <- list(x = tbw$xdati$iord, y = tbw$ydati$iord)

      lapply(1:length(tl), myf)
    }

      
    if (tbw$ncon > 0){
      dfactor <- nconfac
      dfactor <- list(x = EssDee(xcon)*dfactor, y = EssDee(ycon)*dfactor)

      tl <- list(x = tbw$xdati$icon, y = tbw$ydati$icon)

      lapply(1:length(tl), myf)
    }
  
    tbw <- conbandwidth(xbw = tbw$xbw,
                        ybw = tbw$ybw,
                        bwmethod = tbw$method,
                        bwscaling = tbw$scaling,
                        bwtype = tbw$type,
                        cxkertype = tbw$cxkertype,
                        cxkerorder = tbw$cxkerorder,
                        uxkertype = tbw$uxkertype,
                        oxkertype = tbw$oxkertype,
                        cykertype = tbw$cykertype,
                        cykerorder = tbw$cykerorder,
                        uykertype = tbw$uykertype,
                        oykertype = tbw$oykertype,
                        fval = tbw$fval,
                        ifval = tbw$ifval,
                        fval.history = tbw$fval.history,
                        nobs = tbw$nobs,
                        xdati = tbw$xdati,
                        ydati = tbw$ydati,      
                        xnames = tbw$xnames,
                        ynames = tbw$ynames,
                        sfactor = tbw$sfactor,
                        bandwidth = tbw$bandwidth,
                        rows.omit = rows.omit,
                        nconfac = nconfac,
                        ncatfac = ncatfac,
                        sdev = mysd,
                        bandwidth.compute = bandwidth.compute,
                        timing = tbw$timing,
                        total.time = tbw$total.time)
           
    tbw
  }

npcdensbw.NULL <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           bws, ...){

    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain y names and 'toFrame'
    ydat <- toFrame(ydat)

    ## do bandwidths
    
    bws = double(ncol(ydat)+ncol(xdat))

    tbw <- npcdensbw.default(xdat = xdat, ydat = ydat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw
  }

npcdensbw.default <-
  function(xdat = stop("data 'xdat' missing"),
           ydat = stop("data 'ydat' missing"),
           bws, 
           bandwidth.compute = TRUE,
           nmulti, remin, itmax, 
           ftol, tol, small,memfac,
           lbc.dir, dfc.dir, cfac.dir,initc.dir, 
           lbd.dir, hbd.dir, dfac.dir, initd.dir, 
           lbc.init, hbc.init, cfac.init, 
           lbd.init, hbd.init, dfac.init, 
           scale.init.categorical.sample,
           ## dummy arguments for conbandwidth() function call
           bwmethod, bwscaling, bwtype,
           cxkertype, cxkerorder,
           cykertype, cykerorder,
           uxkertype, uykertype,
           oxkertype, oykertype,
           ...){

    ## maintain x names and 'toFrame'
    xdat <- toFrame(xdat)

    ## maintain y names and 'toFrame'
    ydat <- toFrame(ydat)

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("bwmethod", "bwscaling", "bwtype", "cxkertype", "cxkerorder",
               "cykertype", "cykerorder", "uxkertype", "uykertype", "oxkertype",
               "oykertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    tbw <- eval(parse(text=paste("conbandwidth(",
                        "xbw = bws[length(ydat)+1:length(xdat)],",
                        "ybw = bws[1:length(ydat)],",
                        paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                        ifelse(any.m, ",",""),
                        "nobs = nrow(xdat),",
                        "xdati = untangle(xdat),",
                        "ydati = untangle(ydat),",
                        "xnames = names(xdat),",
                        "ynames = names(ydat),",
                        "bandwidth.compute = bandwidth.compute)")))
                        
    ## next grab dummies for actual bandwidth selection and perform call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("bandwidth.compute", "nmulti", "remin", "itmax", "ftol",
               "tol", "small", "memfac",
               "lbc.dir", "dfc.dir", "cfac.dir","initc.dir", 
               "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir", 
               "lbc.init", "hbc.init", "cfac.init", 
               "lbd.init", "hbd.init", "dfac.init", 
               "scale.init.categorical.sample")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    tbw <- eval(parse(text=paste("npcdensbw.conbandwidth(xdat=xdat, ydat=ydat, bws=tbw",
                        ifelse(any.m, ",",""),
                        paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                        ")")))

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }

