npregbw <-
  function(...){
    args = list(...)
    if (is(args[[1]],"formula"))
      UseMethod("npregbw",args[[1]])
    else if (!is.null(args$formula))
      UseMethod("npregbw",args$formula)
    else
      UseMethod("npregbw",args[[which(names(args)=="bws")[1]]])
  }

npregbw.formula <-
  function(formula, data, subset, na.action, call, ...){

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]
    
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]
    
    tbw <- npregbw(xdat = xdat, ydat = ydat, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)
    
    tbw
  }

npregbw.NULL <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws, ...){

    xdat <- toFrame(xdat)

    bws = double(dim(xdat)[2])
    
    tbw <- npregbw.default(xdat = xdat, ydat = ydat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <- updateBwNameMetadata(nameList =
                                list(ynames = deparse(substitute(ydat))),
                                bws = tbw)
    
    tbw
  }

npregbw.rbandwidth <- 
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws, bandwidth.compute = TRUE,
           nmulti, remin = TRUE, itmax = 10000,
           ftol = 1.19209e-07, tol = 1.49012e-08, small = 2.22045e-16, ...){

    xdat <- toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- min(5,dim(xdat)[2])
    }

    if(!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (length(bws$bw) != dim(xdat)[2])
      stop("length of bandwidth vector does not match number of columns of 'xdat'")

    ccon = unlist(lapply(xdat[,bws$icon, drop = FALSE],class))
    if ((any(bws$icon) && !all((ccon == class(integer(0))) | (ccon == class(numeric(0))))) ||
        (any(bws$iord) && !all(unlist(lapply(xdat[,bws$iord, drop = FALSE],class)) ==
                               class(ordered(0)))) ||
        (any(bws$iuno) && !all(unlist(lapply(xdat[,bws$iuno, drop = FALSE],class)) ==
                               class(factor(0)))))
      stop("supplied bandwidths do not match 'xdat' in type")

    if (dim(xdat)[1] != length(ydat))
      stop("number of regression data and response data do not match")

    ## catch and destroy NA's
    goodrows = 1:dim(xdat)[1]
    rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows]
    
    nrow = dim(xdat)[1]
    ncol = dim(xdat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat = toMatrix(xdat)

    runo = xdat[, bws$iuno, drop = FALSE]
    rcon = xdat[, bws$icon, drop = FALSE]
    rord = xdat[, bws$iord, drop = FALSE]

    tbw <- bws

    if (bandwidth.compute){
      myopti = list(num_obs_train = dim(xdat)[1], 
        iMultistart = ifelse(nmulti==0,IMULTI_FALSE,IMULTI_TRUE),
        iNum_Multistart = nmulti,
        int_use_starting_values = ifelse(all(bws$bw==0),USE_START_NO, USE_START_YES),
        int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
        BANDWIDTH_reg_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=ifelse(remin,RE_MIN_TRUE,RE_MIN_FALSE), 
        int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE), 
        bwmethod = switch(bws$method,
          cv.aic = BWM_CVAIC,
          cv.ls = BWM_CVLS),
        kerneval = switch(bws$ckertype,
          gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
          uniform = CKER_UNI),
        ukerneval = switch(bws$ukertype,
          aitchisonaitken = UKER_AIT,
          liracine = UKER_LR),
        okerneval = switch(bws$okertype,
          wangvanryzin = OKER_WANG,
          liracine = OKER_LR),
        nuno = bws$nuno,
        nord = bws$nord,
        ncon = bws$ncon,
        regtype = switch(bws$regtype,
          lc = REGTYPE_LC,
          ll = REGTYPE_LL),
        int_do_tree = ifelse(options('np.tree'), DO_TREE_YES, DO_TREE_NO))
      
      myoptd = list(ftol=ftol, tol=tol, small=small)

      myout=
        .C("np_regression_bw",
           as.double(runo), as.double(rord), as.double(rcon), as.double(ydat),
           as.integer(myopti), as.double(myoptd), 
           bw = c(bws$bw[bws$icon],bws$bw[bws$iuno],bws$bw[bws$iord]),
           fval = double(2),
           PACKAGE="npRmpi" )[c("bw","fval")]
      

      rorder = numeric(ncol)
      rorder[c((1:ncol)[bws$icon], (1:ncol)[bws$iuno], (1:ncol)[bws$iord])]=1:ncol

      tbw$bw <- myout$bw[rorder]
      tbw$fval <- myout$fval[1]
      tbw$ifval <- myout$fval[2]
    }

    tbw$sfactor <- tbw$bandwidth <- tbw$bw

    nfactor <- nrow^(-2.0/(2.0*tbw$ckerorder+tbw$ncon))

    if (tbw$nuno > 0){
      if(tbw$scaling){ 
        tbw$bandwidth[tbw$xdati$iuno] <- tbw$bandwidth[tbw$xdati$iuno]*nfactor
      } else {
        tbw$sfactor[tbw$xdati$iuno] <- tbw$sfactor[tbw$xdati$iuno]/nfactor
      }
    }
    
    if (tbw$nord > 0){
      if(tbw$scaling){
        tbw$bandwidth[tbw$xdati$iord] <- tbw$bandwidth[tbw$xdati$iord]*nfactor
      } else {
        tbw$sfactor[tbw$xdati$iord] <- tbw$sfactor[tbw$xdati$iord]/nfactor
      }
    }

    if (tbw$ncon > 0){
      dfactor <- EssDee(rcon)*nrow^(-1.0/(2.0*tbw$ckerorder+tbw$ncon))

      if (tbw$scaling) {
        tbw$bandwidth[tbw$xdati$icon] <- tbw$bandwidth[tbw$xdati$icon]*dfactor
      } else {
        tbw$sfactor[tbw$xdati$icon] <- tbw$sfactor[tbw$xdati$icon]/dfactor
      }
    }
    


    tbw <- rbandwidth(bw = tbw$bw,
                      regtype = tbw$regtype,
                      bwmethod = tbw$method,
                      bwscaling = tbw$scaling,
                      bwtype = tbw$type,
                      ckertype = tbw$ckertype,
                      ckerorder = tbw$ckerorder,
                      ukertype = tbw$ukertype,
                      okertype = tbw$okertype,
                      fval = tbw$fval,
                      ifval = tbw$ifval,
                      nobs = tbw$nobs,
                      xdati = tbw$xdati,
                      ydati = tbw$ydati,
                      xnames = tbw$xnames,
                      ynames = tbw$ynames,
                      sfactor = tbw$sfactor,
                      bandwidth = tbw$bandwidth,
                      rows.omit = rows.omit,
                      bandwidth.compute = bandwidth.compute)
    tbw
  }

npregbw.default <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws,
           bandwidth.compute = TRUE, nmulti,
           remin, itmax, ftol, tol, small,
           ## dummy arguments for later passing into rbandwidth()
           regtype, bwmethod, bwscaling, bwtype,
           ckertype, ckerorder, ukertype, okertype,
           ...){

    xdat <- toFrame(xdat)

    if(!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector")

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("regtype", "bwmethod", "bwscaling", "bwtype",
               "ckertype", "ckerorder", "ukertype", "okertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    tbw <- eval(parse(text=paste("rbandwidth(bws",
                        ifelse(any.m, ",",""),
                        paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                        ", nobs = dim(xdat)[1],",
                        "xdati = untangle(xdat),",
                        "ydati = untangle(data.frame(ydat)),",
                        "xnames = names(xdat),",
                        "ynames = deparse(substitute(ydat)),",
                        "bandwidth.compute = bandwidth.compute)")))

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("bandwidth.compute", "nmulti", "remin", "itmax", "ftol", "tol",
               "small")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    tbw <- eval(parse(text=paste("npregbw.rbandwidth(xdat=xdat, ydat=ydat, bws=tbw",
                        ifelse(any.m, ",",""),
                        paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                        ")")))

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
    
  }

