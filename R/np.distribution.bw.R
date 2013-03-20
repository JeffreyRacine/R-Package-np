npudistbw <- function(...){
  args = list(...)
  if (is(args[[1]],"formula"))
    UseMethod("npudistbw",args[[1]])
  else if (!is.null(args$formula))
    UseMethod("npudistbw",args$formula)
  else
    UseMethod("npudistbw",args[[which(names(args)=="bws")[1]]])
}

npudistbw.formula <-
  function(formula, data, subset, na.action, call, gdata = NULL, ...){
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]
                     
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, envir = parent.frame())

    if (attr(attr(mf, "terms"), "response") != 0)
      stop("invalid distribution formula")
    
    dat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    if ((has.gval <- !is.null(gdata))) {
      gmf <- match.call(expand.dots = FALSE)
      gm <- match(c("formula", "gdata"), names(gmf), nomatch = 0)
      gmf <- gmf[c(1,gm)]
                     
      gmf[[1]] <- as.name("model.frame")
      names(gmf)[3] <- "data"
      gmf <- eval(gmf, envir = parent.frame())

      gdat <- gmf[, attr(attr(gmf, "terms"),"term.labels"), drop = FALSE]

    }

    tbw <- eval(parse(text=paste("npudistbw(dat = dat,",
                        ifelse(has.gval,"gdat = gdat,",""),
                        "...)")))
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$terms <- attr(mf,"terms")
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw
  }


npudistbw.NULL <-
  function(dat = stop("invoked without input data 'dat'"),
           bws, ...){

    t.names <- NULL
    if(!is.data.frame(dat) && !is.matrix(dat))
      t.names <- deparse(substitute(dat))

    dat = toFrame(dat)
    
    if(!is.null(t.names))
      names(dat) <- t.names

    bws = double(dim(dat)[2])

    tbw <- npudistbw.default(dat = dat, bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw
  }

npudistbw.dbandwidth <- 
  function(dat = stop("invoked without input data 'dat'"),
           bws, gdat = NULL, bandwidth.compute = TRUE, nmulti, remin = TRUE, itmax = 10000,
           fast.cdf = TRUE, do.full.integral = FALSE, ngrid = 100,
           ftol=1.19209e-07, tol=1.49012e-08, small=2.22045e-16, ...){

    dat = toFrame(dat)

    nofi <- missing(do.full.integral)
    nogi <- missing(ngrid)
    
    if (missing(nmulti)){
      nmulti <- min(5,dim(dat)[2])
    }

    if (length(bws$bw) != dim(dat)[2])
      stop(paste("length of bandwidth vector does not match number of columns of",
           "'dat'"))

    ccon = unlist(lapply(as.data.frame(dat[,bws$icon]),class))
    if ((any(bws$icon) && !all((ccon == class(integer(0))) | (ccon == class(numeric(0))))) ||
        (any(bws$iord) && !all(unlist(lapply(as.data.frame(dat[,bws$iord]),class)) ==
                               class(ordered(0)))) ||
        (any(bws$iuno) && !all(unlist(lapply(as.data.frame(dat[,bws$iuno]),class)) ==
                               class(factor(0)))))
      stop(paste("supplied bandwidths do not match", "'dat'", "in type"))

    if(any(bws$iuno))
      stop("distribution bandwidth selection does not support unordered data types")

    dat <- na.omit(dat)
    rows.omit <- unclass(na.action(dat))

    nrow = dim(dat)[1]
    ncol = dim(dat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    odat <- dat
    dat = toMatrix(dat)

    duno = dat[, bws$iuno, drop = FALSE]
    dcon = dat[, bws$icon, drop = FALSE]
    dord = dat[, bws$iord, drop = FALSE]

    tbw <- bws

    ## these are the points where we evaluate the CDF
    ## for now we default to the training points (hence cdf_on_train = TRUE below)
    if(!is.null(gdat)){
      gdat <- toFrame(gdat)
      if(any(is.na(gdat)))
        stop("na's not allowed to be present in cdf gdata")
      gdat <- toMatrix(gdat)

      guno = gdat[, bws$iuno, drop = FALSE]
      gord = gdat[, bws$iord, drop = FALSE]
      gcon = gdat[, bws$icon, drop = FALSE]
      cdf_on_train = FALSE
      nog = nrow(gdat)
      
    } else {
      if(do.full.integral) {
        cdf_on_train = TRUE
        nog = 0
        guno = data.frame()
        gord = data.frame()
        gcon = data.frame()
      } else {
        cdf_on_train = FALSE
        nog = ngrid
        probs <- seq(0,1,length.out = nog)
        ev <- odat[1:nog,,drop = FALSE]
        for(i in 1:ncol(ev)){
          ev[,i] <- uocquantile(odat[,i], probs)
        }

        ev <- toMatrix(ev)
        
        guno = ev[, bws$iuno, drop = FALSE]
        gord = ev[, bws$iord, drop = FALSE]
        gcon = ev[, bws$icon, drop = FALSE]
      }
    }

    if (bandwidth.compute){
      myopti = list(num_obs_train = dim(dat)[1],
        num_obs_eval = nog,
        iMultistart = ifelse(nmulti==0,IMULTI_FALSE,IMULTI_TRUE),
        iNum_Multistart = nmulti,
        int_use_starting_values = ifelse(all(bws$bw==0),USE_START_NO, USE_START_YES),
        int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
        BANDWIDTH_den_extern = switch(bws$type,
          fixed = BW_FIXED,
          generalized_nn = BW_GEN_NN,
          adaptive_nn = BW_ADAP_NN),
        itmax=itmax, int_RESTART_FROM_MIN=ifelse(remin,RE_MIN_TRUE,RE_MIN_FALSE), 
        int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE), 
        bwmethod = switch(bws$method,
          cv.cdf = DBWM_CVLS),
        kerneval = switch(bws$ckertype,
          gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
          epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
          uniform = CKER_UNI),
        cdf_on_train = cdf_on_train,
        nuno = dim(duno)[2],
        nord = dim(dord)[2],
        ncon = dim(dcon)[2],
        fast.cdf = fast.cdf)
      
      myoptd = list(ftol=ftol, tol=tol, small=small)

      if (bws$method != "normal-reference"){
        myout=
          .C("np_distribution_bw", as.double(duno), as.double(dord), as.double(dcon),
             as.double(guno), as.double(gord), as.double(gcon),
             as.integer(myopti), as.double(myoptd), 
             bw = c(bws$bw[bws$icon],bws$bw[bws$iuno],bws$bw[bws$iord]),
             fval = double(2),
             PACKAGE="npRmpi" )[c("bw","fval")]
      } else {
        nbw = double(ncol)
        gbw = sum(bws$icon)
        if (gbw > 0){
          nbw[1:gbw] = 1.587
          if(!bws$scaling)
            nbw[1:gbw]=nbw[1:gbw]*EssDee(dcon)*nrow^(-1.0/3.0)
        }
        myout= list( bw = nbw, fval = c(NA,NA) )
      }

      rorder = numeric(ncol)
      rorder[c((1:ncol)[bws$icon], (1:ncol)[bws$iuno], (1:ncol)[bws$iord])]=1:ncol

      tbw$bw <- myout$bw[rorder]

      tbw$fval = myout$fval[1]
      tbw$ifval = myout$fval[2]
    }
    
    tbw$sfactor <- tbw$bandwidth <- tbw$bw
    
    nfactor <- nrow^(-2.0/(1.0 + tbw$ckerorder))

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
      dfactor <- EssDee(dcon)*nrow^(-1.0/(1.0+tbw$ckerorder))

      if (tbw$scaling) {
        tbw$bandwidth[tbw$xdati$icon] <- tbw$bandwidth[tbw$xdati$icon]*dfactor
      } else {
        tbw$sfactor[tbw$xdati$icon] <- tbw$sfactor[tbw$xdati$icon]/dfactor
      }
    }

    tbw <- dbandwidth(bw = tbw$bw,
                      bwmethod = tbw$method,
                      bwscaling = tbw$scaling,
                      bwtype = tbw$type,
                      ckertype = tbw$ckertype,
                      ckerorder = tbw$ckerorder,
                      ukertype = c("aitchisonaitken"),
                      okertype = tbw$okertype,
                      fval = tbw$fval,
                      ifval = tbw$ifval,
                      nobs = tbw$nobs,
                      xdati = tbw$xdati,
                      xnames = tbw$xnames,
                      sfactor = tbw$sfactor,
                      bandwidth = tbw$bandwidth,
                      rows.omit = rows.omit,
                      bandwidth.compute)
    
    tbw
  }

npudistbw.default <-
  function(dat = stop("invoked without input data 'dat'"),
           bws, gdat, bandwidth.compute = TRUE,
           ## dummy arguments for later passing into npudistbw.bandwidth
           nmulti, remin, itmax, fast.cdf, do.full.integral, ngrid, ftol, tol, small,
           ## dummy arguments for later passing into bandwidth()
           bwmethod, bwscaling, bwtype,
           ckertype, ckerorder, okertype,
           ...){

    t.names <- NULL
    if(!is.data.frame(dat) && !is.matrix(dat))
      t.names <- deparse(substitute(dat))

    dat <- toFrame(dat)
    
    if(!is.null(t.names))
      names(dat) <- t.names

    ## first grab dummy args for bandwidth() and perform 'bootstrap'
    ## bandwidth() call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
               "okertype")


    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    tbw <- eval(parse(text=paste("dbandwidth(bws",
                        ifelse(any.m, ",",""),
                        paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                        ", nobs = dim(dat)[1], xdati = untangle(dat),",
                        "xnames = names(dat),",
                        "bandwidth.compute = bandwidth.compute)")))


    ## next grab dummies for actual bandwidth selection and perform call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("gdat","bandwidth.compute", "nmulti", "remin", "itmax", "fast.cdf", "do.full.integral", "ngrid", "ftol", "tol",
               "small")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    tbw <- eval(parse(text=paste("npudistbw.dbandwidth(dat=dat, bws=tbw",
                        ifelse(any.m, ",",""),
                        paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                        ")")))

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }

