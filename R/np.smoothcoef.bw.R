npscoefbw <-
  function(...){
    args = list(...)
    if (is(args[[1]],"formula"))
      UseMethod("npscoefbw",args[[1]])
    else if (!is.null(args$formula))
      UseMethod("npscoefbw",args$formula)
    else
      UseMethod("npscoefbw",args[[which(names(args)=="bws")[1]]])
  }

npscoefbw.formula <-
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

    chromoly <- explodePipe(mf[["formula"]])

    bronze <- sapply(chromoly, paste, collapse = " + ")
    mf[["formula"]] <-
      as.formula(paste(bronze[1]," ~ ",
                       paste(bronze[2:length(bronze)],
                             collapse =" + ")),
                 env = environment(formula))
    
    mf <- eval(mf, parent.frame())
    
    ydat <- model.response(mf)
    xdat <- mf[, chromoly[[2]], drop = FALSE]
    if (!(miss.z <- !(length(chromoly) == 3)))
      zdat <- mf[, chromoly[[3]], drop = FALSE]
    
    tbw <- eval(parse(text = paste('npscoefbw(xdat = xdat, ydat = ydat,',
                        ifelse(miss.z,'','zdat = zdat,'), '...)')))

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

npscoefbw.NULL <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws, ...){

    miss.z <- missing(zdat)
    ##class(xdat)

    xdat <- toFrame(xdat)

    if(!miss.z)
      zdat <- toFrame(zdat)

    bws <- double(ifelse(miss.z, ncol(xdat), ncol(zdat)))

    tbw <- eval(parse(text = paste("npscoefbw.default(xdat = xdat, ydat = ydat, bws = bws,",
                        ifelse(miss.z,"", "zdat = zdat,"), "...)")))

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <-
      updateBwNameMetadata(nameList = list(ynames = deparse(substitute(ydat))),
                           bws = tbw)

    tbw
  }

npscoefbw.scbandwidth <- 
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws,
           nmulti,
           random.seed = 42,
           cv.iterate = FALSE,
           cv.num.iterations = 1,
           backfit.iterate = FALSE,
           backfit.maxiter = 100,
           backfit.tol = .Machine$double.eps,
           bandwidth.compute = TRUE,
           optim.method = c("Nelder-Mead", "BFGS", "CG"),
           optim.maxattempts = 10,
           optim.reltol = sqrt(.Machine$double.eps),
           optim.abstol = .Machine$double.eps,
           optim.maxit = 500,
           ...){
    
    ## Save seed prior to setting

    if(exists(".Random.seed", .GlobalEnv)) {
      save.seed <- get(".Random.seed", .GlobalEnv)
      exists.seed = TRUE
    } else {
      exists.seed = FALSE
    }
    
    set.seed(random.seed)

    miss.z <- missing(zdat)
    
    xdat <- toFrame(xdat)

    if (!miss.z)
      zdat <- toFrame(zdat)
    
    if (missing(nmulti)){
      nmulti <- min(5,length(bws$bw))
    }

    if (!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector")

    eval(parse(text = paste("bwMatch(",
                 ifelse(miss.z, "xdat, bws$xdati", "zdat, bws$zdati"),")")))
    
    if (dim(xdat)[1] != length(ydat))
      stop("number of regression data and response data do not match")

    if (ncol(xdat) == 1 && missing(cv.iterate))
      cv.iterate = FALSE

    if (cv.num.iterations < 1 & cv.iterate){
      stop("invalid number of iterations specified")
    }

    if (!all(bws$xdati$icon))
      stop("Only continuous 'x' regressors are supported in this version.")

    optim.method <- match.arg(optim.method)
      
    ## catch and destroy NA's
    goodrows = 1:dim(xdat)[1]
    rows.omit =
      eval(parse(text = paste("attr(na.omit(data.frame(xdat,ydat",
                   ifelse(miss.z,'',',zdat'),')), "na.action")')))
    
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows]

    if(!miss.z)
      zdat <- zdat[goodrows,, drop = FALSE]
    
    nrow = dim(xdat)[1]
    ncol = dim(xdat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat <- toMatrix(xdat)

    ## if (!miss.z)
    ##  zdat <- toMatrix(zdat)
    
    ## bad data
    if (qr(xdat)$rank < ncol(xdat)){
      stop("columns of the independent variable (xdat) are linearly dependent") 
    }

    n <- nrow(xdat)
    
    ## ... do bandwidth selection
    
    ## construct 'W' matrix
    ## in the future one will be able to use a switch to npksum
    ## to emulate W

    W <- as.matrix(data.frame(1,xdat))
    yW <- as.matrix(data.frame(ydat,1,xdat))
    
    if (miss.z){
      zdat <- xdat
      dati <- bws$xdati
    }
    else
      dati <- bws$zdati

    if (bandwidth.compute){
      maxPenalty <- sqrt(.Machine$double.xmax)

      overall.cv.ls <- function(param) {
        if (any(param < 0) || ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
          return(maxPenalty)
        
        bws$bw <- param
        
        ##tyw <- npksum(txdat = zdat, tydat = ydat, weights = W, bws = bws,
        ##              leave.one.out = TRUE)$ksum
        
        tww <- npksum(txdat = zdat, tydat = yW, weights = yW, bws = bws,
                      leave.one.out = TRUE)$ksum

        mean.loo <- rep(maxPenalty,n)
        epsilon <- 1.0/n
        ridge <- double(n)
        doridge <- !logical(n)

        nc <- ncol(tww[-1,-1,1])

        ridger <- function(i) {
          doridge[i] <<- FALSE
          ridge.val <- ridge[i]*tww[-1,1,i][1]/
            (ifelse(tww[-1,-1,i][1,1]>=0, 1, -1)*max(.Machine$double.eps,abs(tww[-1,-1,i][1,1])))
          W[i,, drop = FALSE] %*% tryCatch(solve(tww[-1,-1,i]+diag(rep(ridge[i],nc)),
                  tww[-1,1,i]+c(ridge.val,rep(0,nc-1))),
                  error = function(e){
                    ridge[i] <<- ridge[i]+epsilon
                    doridge[i] <<- TRUE
                    return(rep(maxPenalty,nc))
                  })
        }

        while(any(doridge)){
          iloo <- (1:n)[doridge]
          mean.loo[iloo] <- sapply(iloo, ridger)
        }

        cv.console <<- printClear(cv.console)
        ##cv.console <<- printPush(msg = paste("param:", param), console = cv.console)

        stopifnot(all(is.finite(mean.loo)))

        if(!any(mean.loo == maxPenalty)){
          fv <- sum((ydat-mean.loo)^2)/n
          cv.console <<- printPush(msg = paste("fval:", signif(fv, digits = options('digits')$digits)), console = cv.console)
        } else {
          cv.console <<- printPush(msg = "near-singular system encountered, ridging", console = cv.console)
          fv <- maxPenalty
        }

        return(ifelse(is.finite(fv),fv,maxPenalty))

      }

      scoef.looE <- parse(text = paste('npscoef(bws = bws, txdat = xdat, tydat = ydat,',
                            ifelse(miss.z,'','tzdat = zdat,'),
                            'leave.one.out = TRUE, iterate = TRUE,',
                            'maxiter = backfit.maxiter, tol = backfit.tol,',
                            'betas = TRUE)'))
      
      partial.cv.ls <- function(param, partial.index) {
        if (any(param < 0) || ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
          return(maxPenalty)
        
        if (backfit.iterate){
          bws$bw.fitted[,partial.index] <- param
          scoef.loo <- eval(scoef.looE)
          partial.loo <- W[,partial.index]*scoef.loo$beta[,partial.index]
        } else {
          bws$bw <- param
          partial.loo <- W[,partial.index]*
            npksum(txdat = zdat,
                   tydat = partial.orig * W[,partial.index],
                   bws = bws,
                   leave.one.out = TRUE)$ksum/
                     npksum(txdat = zdat,
                            tydat = W[,partial.index]^2,
                            bws = bws,
                            leave.one.out = TRUE)$ksum
        }
        

        fv <- sum((partial.orig - partial.loo)^2)/n
        
        cv.console <<- printClear(cv.console)
        cv.console <<- printPush(msg = paste("fval:",
                                   signif(fv, digits = options('digits')$digits)),
                                 console = cv.console)
        return(ifelse(is.finite(fv),fv,maxPenalty))
      }

      ## Now we implement multistarting

      fval.min <- .Machine$double.xmax
      numimp <- 0
      value.overall <- numeric(nmulti)

      x.scale <- sapply(1:bws$ndim, function(i){
        if (dati$icon[i]){
          if(IQR(zdat[,i]) > 0) {
            return((4/3)^0.2*(ifelse(bws$scaling, 1.0,
                                     min(sd(zdat[,i]), IQR(zdat[,i])/(qnorm(.25,lower.tail=F)*2)) *
                                     n^(-1.0/(2.0*bws$ckerorder+bws$ncon)))))
          } else {
            return((4/3)^0.2*(ifelse(bws$scaling, 1.0,
                                     sd(zdat[,i]) *
                                     n^(-1.0/(2.0*bws$ckerorder+bws$ncon)))))
          }
        }
        
        if (dati$iord[i])
          return(0.5*oMaxL(dati$all.nlev[[i]], kertype = bws$okertype)*
                 ifelse(bws$scaling,n^(2.0/(2.0*bws$ckerorder+bws$ncon)),1.0))
        
        if (dati$iuno[i])
          return(0.5*uMaxL(dati$all.nlev[[i]], kertype = bws$ukertype)*
                 ifelse(bws$scaling,n^(2.0/(2.0*bws$ckerorder+bws$ncon)),1.0))       
      })

      console <- newLineConsole()
      optim.control <- list(abstol = optim.abstol,
                            reltol = optim.reltol,
                            maxit = optim.maxit)

      for (i in 1:nmulti) {

        console <- printPush(msg = paste(sep="", "Multistart ", i, " of ", nmulti, "... "), console)
        cv.console <- newLineConsole(console)
        
        if (i == 1) {
          tbw <- x.scale
          if(all(bws$bw != 0))
            tbw <- bws$bw
        } else {
            tbw <- runif(bws$ndim, min=0.5, max=1.5)*x.scale
        }

        suppressWarnings(optim.return <- optim(tbw,
                                               fn = overall.cv.ls,
                                               control = optim.control))
        attempts <- 0
        while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
          attempts <- attempts + 1
          tbw <- runif(bws$ndim, min=0.5, max=1.5)*x.scale
          optim.control <- lapply(optim.control, '*', 10.0)
          suppressWarnings(optim.return <- optim(tbw,
                                                 fn = overall.cv.ls,
                                                 control = optim.control))

        }

        cv.console <- printClear(cv.console)

        value.overall[i] <- optim.return$value

        if(optim.return$value < fval.min) {
          param <- optim.return$par
          min.overall <- optim.return$value
          fval.min <- min.overall ## Added by jracine Jul 22 2010
          numimp.overall <- numimp + 1
          best.overall <- i
        }

        if(i < nmulti)
          console <- printPop(console)
        else
          console <- printClear(console)
      }

      param.overall <- bws$bw <- param

      if(cv.iterate){
        
        console <- newLineConsole()
        
        n.part <- (ncol(xdat)+1)
        
        bws$bw.fitted <- matrix(data = bws$bw, nrow = length(bws$bw), ncol = n.part)
        ## obtain matrix of alpha.hat | h0 and beta.hat | h0

        scoef <- eval(parse(text = paste('npscoef(bws = bws, txdat = xdat, tydat = ydat,',
                              ifelse(miss.z, '', 'tzdat = zdat,'),
                              'iterate = FALSE, betas = TRUE)')))
        
        resid.full <- ydat - scoef$mean

        
        for(i in 1:cv.num.iterations){
          console <- printPush(msg = paste(sep="", "backfitting iteration ", i, " of ", cv.num.iterations, "... "), console)

          for(j in 1:n.part){
            console <- printPush(msg = paste(sep="", "partial residual ", j, " of ", n.part, "... "), console)
            cv.console <- newLineConsole(console)

            ## estimate partial residuals
            partial.orig <- W[,j] * scoef$beta[,j] + resid.full
            
            ## minimise
            suppressWarnings(optim.return <-
                             optim(tbw, fn = partial.cv.ls,
                                   control = optim.control,
                                   partial.index = j))
            
            cv.console <- printClear(cv.console)
            
            ## grab parameter
            bws$bw.fitted[,j] <- optim.return$par

            if (backfit.iterate){
              ## re-estimate all betas
              scoef <- eval(parse(text = paste('npscoef(bws = bws, txdat = xdat, tydat = ydat,',
                                    ifelse(miss.z, '', 'tzdat = zdat,'),
                                    'iterate = TRUE, maxiter = backfit.maxiter,',
                                    'tol = backfit.tol, betas = TRUE)')))
              resid.full <- ydat - scoef$mean
            } else {
              bws$bw <- bws$bw.fitted[,j]
              ## estimate new beta.hats
              scoef$beta[,j] <-
                npksum(txdat = zdat,
                       tydat = partial.orig * W[,j],
                       bws = bws)$ksum/
                         npksum(txdat = zdat,
                                tydat = W[,j]^2,
                                bws = bws)$ksum
              bws$bw <- param.overall
              ## estimate new full residuals 
              resid.full <- partial.orig - W[,j] * scoef$beta[,j]
            }

            console <- printPop(console)
          }
          if(i < cv.num.iterations)
            console <- printPop(console)
          else
            console <- printClear(console)
        }
        scoef.loo <- eval(parse(text = paste('npscoef(bws = bws, txdat = xdat, tydat = ydat,',
                                  ifelse(miss.z,'', 'tzdat = zdat,'),
                                  'iterate = TRUE, maxiter = backfit.maxiter,',
                                  'tol = backfit.tol, leave.one.out = TRUE)$mean')))
        bws$fval.fitted <- sum((ydat - scoef.loo)^2)/n
      }

      bws$fval = min.overall
      bws$ifval = best.overall
      bws$numimp = numimp.overall
      bws$fval.vector = value.overall
    }
    
    bws$sfactor <- bws$bandwidth <- bws$bw
    nfactor <- nrow^(-2.0/(2.0*bws$ckerorder+bws$ncon))
    dfactor <- EssDee(zdat[, dati$icon, drop = FALSE])*nrow^(-1.0/(2.0*bws$ckerorder+sum(dati$icon)))

    if (bws$scaling) {
      bws$bandwidth[dati$icon] <- bws$bandwidth[dati$icon]*dfactor

      if(bws$nuno > 0)
        bws$bandwidth[dati$iuno] <- bws$bandwidth[dati$iuno]*nfactor

      if(bws$nord > 0)
        bws$bandwidth[dati$iord] <- bws$bandwidth[dati$iord]*nfactor
      
    } else {
      bws$sfactor[dati$icon] <- bws$sfactor[dati$icon]/dfactor

      if(bws$nuno > 0)
        bws$sfactor[dati$iuno] <- bws$sfactor[dati$iuno]/nfactor

      if(bws$nord > 0)
        bws$sfactor[dati$iord] <- bws$sfactor[dati$iord]/nfactor
    }

    ## Restore seed

    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
    
    bws <- scbandwidth(bw = bws$bw,
                       bwmethod = bws$method,
                       bwscaling = bws$scaling,
                       bwtype = bws$type,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       ukertype = bws$ukertype,
                       okertype = bws$okertype,
                       fval = bws$fval,
                       ifval = bws$ifval,
                       numimp = bws$numimp,
                       fval.vector = bws$fval.vector,
                       nobs = bws$nobs,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       zdati = bws$zdati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       znames = bws$znames,
                       sfactor = bws$sfactor,
                       bandwidth = bws$bandwidth,
                       rows.omit = rows.omit,
                       bandwidth.compute = bandwidth.compute,
                       optim.method = optim.method)

    bws
  }

npscoefbw.default <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws,
           nmulti,
           random.seed,
           cv.iterate,
           cv.num.iterations,
           backfit.iterate,
           backfit.maxiter,
           backfit.tol,
           bandwidth.compute = TRUE,
           ## dummy arguments for scbandwidth()
           bwmethod, bwscaling, bwtype,
           ckertype, ckerorder,
           ukertype, okertype,
           optim.method, optim.maxattempts,
           optim.reltol, optim.abstol, optim.maxit,
           ...){


    miss.z <- missing(zdat)
    xdat <- toFrame(xdat)
    
    if(!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector")

    if(!miss.z)
      zdat <- toFrame(zdat)

    ## first grab dummy args for scbandwidth() and perform 'bootstrap'
    ## bandwidth call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
               "ukertype", "okertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    tbw <- eval(parse(text=paste("scbandwidth(bw = bws",
                        ifelse(any.m, ",",""),
                        paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                        ", nobs = dim(xdat)[1],",
                        "xdati = untangle(xdat),",
                        "ydati = untangle(data.frame(ydat)),",
                        "zdati = untangle(zdat),",
                        "xnames = names(xdat),",
                        "ynames = deparse(substitute(ydat)),",
                        "znames = names(zdat),",
                        "bandwidth.compute = bandwidth.compute)")))

    ## next grab dummies for actual bandwidth selection and perform call
    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("zdat", "bandwidth.compute",
               "nmulti",
               "cv.iterate",
               "cv.num.iterations",
               "backfit.iterate",
               "backfit.maxiter",
               "backfit.tol",
               "optim.method", "optim.maxattempts",
               "optim.reltol", "optim.abstol", "optim.maxit")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)


    tbw <- eval(parse(text=paste("npscoefbw.scbandwidth(xdat=xdat, ydat=ydat, bws=tbw",
                        ifelse(any.m, ",",""),
                        paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                        ")")))

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
    
  }

