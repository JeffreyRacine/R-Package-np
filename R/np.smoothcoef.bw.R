npscoefbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    dots <- mc$...

    if (length(dots) == 0L)
      stop("invoked without arguments")

    dot.names <- names(dots)

    if (!is.null(dot.names) && any(dot.names == "formula")) {
      formula.val <- eval(dots[[which(dot.names == "formula")[1L]]], envir = parent.frame())
      return(UseMethod("npscoefbw", formula.val))
    }

    first.val <- eval(dots[[1L]], envir = parent.frame())
    if (inherits(first.val, "formula"))
      return(UseMethod("npscoefbw", first.val))

    if (!is.null(dot.names) && any(dot.names == "bws")) {
      bws.val <- eval(dots[[which(dot.names == "bws")[1L]]], envir = parent.frame())
      return(UseMethod("npscoefbw", bws.val))
    }

    if (!is.null(dot.names) && any(dot.names %in% c("xdat", "ydat", "zdat")))
      return(UseMethod("npscoefbw", NULL))

    UseMethod("npscoefbw", first.val)
  }

npscoefbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    orig.ts <- if (missing(data))
      sapply(eval(attr(terms(formula), "variables"), environment(formula)), inherits, "ts")
    else sapply(eval(attr(terms(formula, data=data), "variables"), data, environment(formula)), inherits, "ts")

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

    mf[[1]] <- as.name("model.frame")

    formula_to_explode <- eval(mf[["formula"]], parent.frame())
    chromoly <- explodePipe(formula_to_explode, env = environment(formula))

    bronze <- sapply(chromoly, paste, collapse = " + ")
    mf[["formula"]] <-
      as.formula(paste(bronze[1]," ~ ",
                       paste(bronze[2:length(bronze)],
                             collapse =" + ")),
                 env = environment(formula))

    mf[["formula"]] <- terms(mf[["formula"]])
    if(all(orig.ts)){
      args <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      attr(mf[["formula"]], "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      attr(mf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
    }
    
    mf <- eval(mf, parent.frame())
    
    ydat <- model.response(mf)
    xdat <- mf[, chromoly[[2]], drop = FALSE]
    if (!(miss.z <- !(length(chromoly) == 3)))
      zdat <- mf[, chromoly[[3]], drop = FALSE]
    
    bw.args <- list(xdat = xdat, ydat = ydat)
    if (!miss.z)
      bw.args$zdat <- zdat
    tbw <- do.call(npscoefbw, c(bw.args, list(...)))

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
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    miss.z <- missing(zdat)

    xdat <- toFrame(xdat)

    if(!miss.z)
      zdat <- toFrame(zdat)

    bws <- double(ifelse(miss.z, ncol(xdat), ncol(zdat)))

    bw.args <- list(xdat = xdat, ydat = ydat, bws = bws)
    if (!miss.z)
      bw.args$zdat <- zdat
    tbw <- do.call(npscoefbw.default, c(bw.args, list(...)))

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
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))
    
    ## Save seed prior to setting

    if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      save.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
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
      stop("'ydat' must be a vector or a factor")

    if (miss.z) {
      bwMatch(xdat, bws$xdati)
    } else {
      bwMatch(zdat, bws$zdati)
    }
    
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
    train.df <- data.frame(xdat, ydat)
    if (!miss.z)
      train.df <- data.frame(train.df, zdat)
    rows.omit <- attr(na.omit(train.df), "na.action")
    
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

    mysd <- EssDee(zdat[, dati$icon, drop = FALSE])
    nconfac <- n^(-1.0/(2.0*bws$ckerorder+bws$ncon))
    ncatfac <- n^(-2.0/(2.0*bws$ckerorder+bws$ncon))

    bws$sdev <- mysd
    bws$nconfac <- nconfac
    bws$ncatfac <- ncatfac

    total.time <-
      system.time({
        if (bandwidth.compute){
          maxPenalty <- sqrt(.Machine$double.xmax)
          overall.cv.ls <- function(param) {
            sbw <- bws
            sbw$bw <- param
            sbw$bandwidth[[1]] <- param
            if (!validateBandwidthTF(sbw) || ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            
            bws$bw <- param

            if(bws$scaling)
              bws$bandwidth[[1]] <- sapply(1:bws$ndim, function(i) { bws$bw[i]*ifelse(bws$icon[i],nconfac*bws$sdev[sum(dati$icon[1:i])], ncatfac) })
            else
              bws$bandwidth[[1]] <- bws$bw
            
            
            tww <- npksum(txdat = zdat, tydat = yW, weights = yW, bws = bws,
                          leave.one.out = TRUE)$ksum

            mean.loo <- rep(maxPenalty,n)
            epsilon <- 1.0/n
            ridge <- double(n)
            doridge <- !logical(n)

            nc <- ncol(tww[-1,-1,1])

            ridger <- function(i) {
              doridge[i] <<- FALSE
              ridge.val <- ridge[i]*tww[-1,1,i][1]/NZD(tww[-1,-1,i][1,1])
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

          scoef.loo.args <- list(
            bws = bws, txdat = xdat, tydat = ydat,
            leave.one.out = TRUE, iterate = TRUE,
            maxiter = backfit.maxiter, tol = backfit.tol,
            betas = TRUE
          )
          if (!miss.z)
            scoef.loo.args$tzdat <- zdat
          
          partial.cv.ls <- function(param, partial.index) {
            sbw <- bws
            sbw$bw <- param
            sbw$bandwidth[[1]] <- param

            if (!validateBandwidthTF(sbw) || ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            
            if (backfit.iterate){
              bws$bw.fitted[,partial.index] <- param
              scoef.loo <- do.call(npscoef, scoef.loo.args)
              partial.loo <- W[,partial.index]*scoef.loo$beta[,partial.index]
            } else {
              bws$bw <- param

              if(bws$scaling)
                bws$bandwidth[[1]] <- sapply(1:bws$ndim, function(i) { bws$bw[i]*ifelse(bws$icon[i],nconfac*bws$sdev[sum(dati$icon[1:i])], ncatfac) })
              else
                bws$bandwidth[[1]] <- bws$bw

              tww <- npksum(txdat=zdat,
                            tydat=cbind(partial.orig * W[,partial.index],W[,partial.index]^2),
                            weights=cbind(partial.orig * W[,partial.index],1),
                            bws=bws,
                            leave.one.out=TRUE)$ksum

              partial.loo <- W[,partial.index]*tww[1,2,]/NZD(tww[2,2,])
              
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
          num.feval.overall <- 0

          x.scale <- sapply(1:bws$ndim, function(i){
            if (dati$icon[i]){
              return(1.059224*(ifelse(bws$scaling, 1.0, mysd[sum(dati$icon[1:i])]*nconfac)))
            }
            
            if (dati$iord[i])
              return(0.5*oMaxL(dati$all.nlev[[i]], kertype = bws$okertype)*
                     ifelse(bws$scaling,ncatfac,1.0))
            
            if (dati$iuno[i])
              return(0.5*uMaxL(dati$all.nlev[[i]], kertype = bws$ukertype)*
                     ifelse(bws$scaling,ncatfac,1.0))       
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
            if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
              num.feval.overall <- num.feval.overall + optim.return$counts[1]
            attempts <- 0
            while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
              attempts <- attempts + 1
              tbw <- runif(bws$ndim, min=0.5, max=1.5)*x.scale
              optim.control <- lapply(optim.control, '*', 10.0)
              suppressWarnings(optim.return <- optim(tbw,
                                                     fn = overall.cv.ls,
                                                     control = optim.control))
              if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                num.feval.overall <- num.feval.overall + optim.return$counts[1]

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

            scoef.args <- list(bws = bws, txdat = xdat, tydat = ydat, iterate = FALSE, betas = TRUE)
            if (!miss.z)
              scoef.args$tzdat <- zdat
            scoef <- do.call(npscoef, scoef.args)
            
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
                if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                  num.feval.overall <- num.feval.overall + optim.return$counts[1]
                
                cv.console <- printClear(cv.console)
                
                ## grab parameter
                bws$bw.fitted[,j] <- optim.return$par

                if (backfit.iterate){
                  ## re-estimate all betas
                  scoef.args <- list(
                    bws = bws, txdat = xdat, tydat = ydat,
                    iterate = TRUE, maxiter = backfit.maxiter,
                    tol = backfit.tol, betas = TRUE
                  )
                  if (!miss.z)
                    scoef.args$tzdat <- zdat
                  scoef <- do.call(npscoef, scoef.args)
                  resid.full <- ydat - scoef$mean
                } else {
                  bws$bw <- bws$bw.fitted[,j]
                  ## estimate new beta.hats

                  if(bws$scaling)
                    bws$bandwidth[[1]] <- sapply(1:bws$ndim, function(i) { bws$bw[i]*ifelse(bws$icon[i],nconfac*bws$sdev[sum(dati$icon[1:i])], ncatfac) })
                  else
                    bws$bandwidth[[1]] <- bws$bw

                  tww <- npksum(txdat=zdat,
                                tydat=cbind(partial.orig * W[,j],W[,j]^2),
                                weights=cbind(partial.orig * W[,j],1),
                                bws=bws)$ksum
                  
                  scoef$beta[,j] <- tww[1,2,]/NZD(tww[2,2,])
                  
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
            scoef.loo.args <- list(
              bws = bws, txdat = xdat, tydat = ydat,
              iterate = TRUE, maxiter = backfit.maxiter,
              tol = backfit.tol, leave.one.out = TRUE
            )
            if (!miss.z)
              scoef.loo.args$tzdat <- zdat
            scoef.loo <- do.call(npscoef, scoef.loo.args)$mean
            bws$fval.fitted <- sum((ydat - scoef.loo)^2)/n
          }

          bws$fval = min.overall
          bws$ifval = best.overall
          bws$num.feval = num.feval.overall
          bws$numimp = numimp.overall
          bws$fval.vector = value.overall
        }
      })[1]
    
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

    if(exists.seed) assign(".Random.seed", save.seed, envir = .GlobalEnv)
    
    bws <- scbandwidth(bw = bws$bw,
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
                       fval = bws$fval,
                       ifval = bws$ifval,
                       num.feval = bws$num.feval,
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
                       optim.method = optim.method,
                       total.time = total.time)

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
           ckertype, ckerorder, ckerbound, ckerlb, ckerub,
           ukertype, okertype,
           optim.method, optim.maxattempts,
           optim.reltol, optim.abstol, optim.maxit,
           ...){
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    if (!missing(bwmethod) && identical(match.arg(bwmethod, c("cv.ls", "manual")), "manual") &&
        missing(bws))
      stop("bwmethod='manual' requires argument 'bws'")


    miss.z <- missing(zdat)
    xdat <- toFrame(xdat)
    
    if(!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector or a factor")

    if(!miss.z)
      zdat <- toFrame(zdat)

    ## first grab dummy args for scbandwidth() and perform 'bootstrap'
    ## bandwidth call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
               "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    sbw.args <- list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      zdati = untangle(zdat),
      xnames = names(xdat),
      ynames = deparse(substitute(ydat)),
      znames = names(zdat),
      bandwidth.compute = bandwidth.compute
    )
    if (any.m) {
      nms <- mc.names[m]
      sbw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    tbw <- do.call(scbandwidth, sbw.args)

    ## next grab dummies for actual bandwidth selection and perform call
    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("zdat", "bandwidth.compute",
               "nmulti",
               "random.seed",
               "cv.iterate",
               "cv.num.iterations",
               "backfit.iterate",
               "backfit.maxiter",
               "backfit.tol",
               "optim.method", "optim.maxattempts",
               "optim.reltol", "optim.abstol", "optim.maxit")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)


    scbw.args <- list(xdat = xdat, ydat = ydat, bws = tbw)
    if (any.m) {
      nms <- mc.names[m]
      scbw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    tbw <- do.call(npscoefbw.scbandwidth, scbw.args)

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
    
  }
