npscoefbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat", "zdat"),
                                     eval_env = parent.frame())
    UseMethod("npscoefbw", target)
  }

npscoefbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = terms(formula),
                        data = environment(formula),
                        eval_env = environment(formula))
    else .np_terms_ts_mask(terms_obj = terms(formula, data = data),
                           data = data,
                           eval_env = environment(formula))

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    formula.call <- .np_bw_formula_from_call(call_obj = call, eval_env = parent.frame())
    if (!is.null(formula.call))
      mf[[2]] <- formula.call

    mf[[1]] <- as.name("model.frame")

    formula.obj <- .np_bw_resolve_formula(formula_obj = formula,
                                        formula_call = formula.call,
                                        eval_env = parent.frame())
    chromoly <- explodePipe(formula.obj, env = environment(formula))

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
    
    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())
    
    ydat <- model.response(mf)
    xdat <- mf[, chromoly[[2]], drop = FALSE]
    miss.z <- !(length(chromoly) == 3)
    if (!miss.z)
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

    n.bw <- if (miss.z) ncol(xdat) else ncol(zdat)
    bws <- double(n.bw)

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
    ## Save seed prior to setting

    seed.state <- .np_seed_enter(random.seed)


    miss.z <- missing(zdat)
    
    xdat <- toFrame(xdat)

    if (!miss.z)
      zdat <- toFrame(zdat)
    
    if (missing(nmulti)){
      nmulti <- min(5,length(bws$bw))
    }
    regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype
    cv.iterate <- npValidateScalarLogical(cv.iterate, "cv.iterate")
    backfit.iterate <- npValidateScalarLogical(backfit.iterate, "backfit.iterate")
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    nmulti <- npValidateNonNegativeInteger(nmulti, "nmulti")
    backfit.maxiter <- npValidatePositiveInteger(backfit.maxiter, "backfit.maxiter")
    backfit.tol <- npValidatePositiveFiniteNumeric(backfit.tol, "backfit.tol")
    optim.maxattempts <- npValidatePositiveInteger(optim.maxattempts, "optim.maxattempts")
    optim.maxit <- npValidatePositiveInteger(optim.maxit, "optim.maxit")
    optim.reltol <- npValidatePositiveFiniteNumeric(optim.reltol, "optim.reltol")
    optim.abstol <- npValidatePositiveFiniteNumeric(optim.abstol, "optim.abstol")
    if (cv.iterate)
      cv.num.iterations <- npValidatePositiveInteger(cv.num.iterations, "cv.num.iterations")
    if (!identical(regtype, "lc") && cv.iterate) {
      warning("cv.iterate currently supports regtype='lc' for npscoefbw; using cv.iterate=FALSE")
      cv.iterate <- FALSE
    }
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    if (!(is.vector(ydat) || is.factor(ydat)))
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

    if (!all(bws$xdati$icon))
      stop("Only continuous 'x' regressors are supported in this version.")

    optim.method <- match.arg(optim.method)
      
    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    train.df <- data.frame(xdat, ydat)
    if (!miss.z)
      train.df <- data.frame(train.df, zdat)
    rows.omit <- attr(na.omit(train.df), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]

    if(!miss.z)
      zdat <- zdat[keep.rows,, drop = FALSE]
    
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

    W <- cbind(1.0, xdat)
    yW <- cbind(ydat, W)
    
    if (miss.z){
      zdat <- xdat
      dati <- bws$xdati
    }
    else
      dati <- bws$zdati
    zdat.df <- if (is.data.frame(zdat)) zdat else as.data.frame(zdat)

    mysd <- EssDee(zdat[, dati$icon, drop = FALSE])
    nconfac <- n^(-1.0/(2.0*bws$ckerorder+bws$ncon))
    ncatfac <- n^(-2.0/(2.0*bws$ckerorder+bws$ncon))

    bws$sdev <- mysd
    bws$nconfac <- nconfac
    bws$ncatfac <- ncatfac
    bw.scale.multiplier <- NULL
    if (bws$scaling) {
      bw.scale.multiplier <- rep(ncatfac, bws$ndim)
      if (any(bws$icon)) {
        icon.cumsum <- cumsum(dati$icon)
        bw.scale.multiplier[bws$icon] <- nconfac * bws$sdev[icon.cumsum[bws$icon]]
      }
    }
    apply_bw_to_scbw <- function(scbw, param) {
      scbw$bw <- param
      if (scbw$scaling)
        scbw$bandwidth[[1]] <- scbw$bw * bw.scale.multiplier
      else
        scbw$bandwidth[[1]] <- scbw$bw
      scbw
    }

    total.time <-
      system.time({
        if (bandwidth.compute){
          maxPenalty <- sqrt(.Machine$double.xmax)
          cv_state <- new.env(parent = emptyenv())
          cv_state$console <- NULL
          overall.cv.ls <- function(param) {
            sbw <- apply_bw_to_scbw(bws, param)
            if (!validateBandwidthTF(sbw) || ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)

            tww <- npksum(txdat = zdat, tydat = yW, weights = yW, bws = sbw,
                          leave.one.out = TRUE)$ksum

            mean.loo <- rep(maxPenalty,n)
            epsilon <- 1.0/n
            ridge <- double(n)
            doridge <- rep.int(TRUE, n)

            nc <- ncol(tww[-1,-1,1])

            while(any(doridge)){
              iloo <- which(doridge)
              for (ii in iloo) {
                doridge[ii] <- FALSE
                ridge.val <- ridge[ii]*tww[-1,1,ii][1]/NZD(tww[-1,-1,ii][1,1])
                beta.ii <- tryCatch(
                  solve(tww[-1,-1,ii] + diag(rep(ridge[ii], nc)),
                        tww[-1,1,ii] + c(ridge.val, rep(0, nc - 1))),
                  error = function(e) e
                )
                if (inherits(beta.ii, "error")) {
                  ridge[ii] <- ridge[ii] + epsilon
                  doridge[ii] <- TRUE
                  beta.ii <- rep(maxPenalty, nc)
                }
                mean.loo[ii] <- W[ii,, drop = FALSE] %*% beta.ii
              }
            }

            cv_state$console <- printClear(cv_state$console)
            ##cv_state$console <- printPush(msg = paste("param:", param), console = cv_state$console)

            stopifnot(all(is.finite(mean.loo)))

            if(!any(mean.loo == maxPenalty)){
              fv <- sum((ydat-mean.loo)^2)/n
              cv_state$console <- printPush(msg = paste("fval:", signif(fv, digits = getOption("digits", 7L))), console = cv_state$console)
            } else {
              cv_state$console <- printPush(msg = "near-singular system encountered, ridging", console = cv_state$console)
              fv <- maxPenalty
            }

            return((if (is.finite(fv)) fv else maxPenalty))

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
            sbw <- apply_bw_to_scbw(bws, param)

            if (!validateBandwidthTF(sbw) || ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            
            if (backfit.iterate){
              bws$bw.fitted[,partial.index] <- param
              scoef.loo <- do.call(npscoef, scoef.loo.args)
              partial.loo <- W[,partial.index]*scoef.loo$beta[,partial.index]
            } else {
              wj <- W[,partial.index]
              if (identical(regtype, "lc")) {
                tww <- npksum(txdat=zdat,
                              tydat=cbind(partial.orig * wj, wj * wj),
                              weights=cbind(partial.orig * wj, 1),
                              bws=sbw,
                              leave.one.out=TRUE)$ksum

                partial.loo <- wj * tww[1,2,]/NZD(tww[2,2,])
              } else {
                kw <- .npscoef_weight_matrix(
                  bws = sbw,
                  tzdat = zdat.df,
                  ezdat = zdat.df,
                  leave.one.out = TRUE
                )
                num <- as.vector(crossprod(kw, partial.orig * wj))
                den <- as.vector(crossprod(kw, wj * wj))
                partial.loo <- wj * num / NZD(den)
              }
            }
            

            fv <- sum((partial.orig - partial.loo)^2)/n
            
            cv_state$console <- printClear(cv_state$console)
            cv_state$console <- printPush(msg = paste("fval:",
                                       signif(fv, digits = getOption("digits", 7L))),
                                     console = cv_state$console)
            return((if (is.finite(fv)) fv else maxPenalty))
          }

          ## Now we implement multistarting

          fval.min <- .Machine$double.xmax
          numimp <- 0
          value.overall <- numeric(nmulti)
          num.feval.overall <- 0

          x.scale <- sapply(seq_len(bws$ndim), function(i){
            if (dati$icon[i]){
              return(1.059224*((if (bws$scaling) 1.0 else mysd[sum(dati$icon[seq_len(i)])]*nconfac)))
            }
            
            if (dati$iord[i])
              return(0.5*oMaxL(dati$all.nlev[[i]], kertype = bws$okertype)*
                     (if (bws$scaling) ncatfac else 1.0))
            
            if (dati$iuno[i])
              return(0.5*uMaxL(dati$all.nlev[[i]], kertype = bws$ukertype)*
                     (if (bws$scaling) ncatfac else 1.0))       
          })

          console <- newLineConsole()
          optim.control <- list(abstol = optim.abstol,
                                reltol = optim.reltol,
                                maxit = optim.maxit)

          for (i in seq_len(nmulti)) {

            console <- printPush(msg = paste(sep="", "Multistart ", i, " of ", nmulti, "... "), console)
            cv_state$console <- newLineConsole(console)
            
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

            cv_state$console <- printClear(cv_state$console)

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

            
            for (i in seq_len(cv.num.iterations)) {
              console <- printPush(msg = paste(sep="", "backfitting iteration ", i, " of ", cv.num.iterations, "... "), console)

              for (j in seq_len(n.part)) {
                console <- printPush(msg = paste(sep="", "partial residual ", j, " of ", n.part, "... "), console)
                cv_state$console <- newLineConsole(console)

                ## estimate partial residuals
                partial.orig <- W[,j] * scoef$beta[,j] + resid.full
                
                ## minimise
                suppressWarnings(optim.return <-
                                 optim(tbw, fn = partial.cv.ls,
                                       control = optim.control,
                                       partial.index = j))
                if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                  num.feval.overall <- num.feval.overall + optim.return$counts[1]
                
                cv_state$console <- printClear(cv_state$console)
                
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

                  bws <- apply_bw_to_scbw(bws, bws$bw)

                  if (identical(regtype, "lc")) {
                    wj <- W[,j]
                    tww <- npksum(txdat=zdat,
                                  tydat=cbind(partial.orig * wj, wj * wj),
                                  weights=cbind(partial.orig * wj, 1),
                                  bws=bws)$ksum
                    scoef$beta[,j] <- tww[1,2,]/NZD(tww[2,2,])
                  } else {
                    wj <- W[,j]
                    kw <- .npscoef_weight_matrix(
                      bws = bws,
                      tzdat = zdat.df,
                      ezdat = zdat.df,
                      leave.one.out = FALSE
                    )
                    num <- as.vector(crossprod(kw, partial.orig * wj))
                    den <- as.vector(crossprod(kw, wj * wj))
                    scoef$beta[,j] <- num / NZD(den)
                  }
                  
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

    .np_seed_exit(seed.state)
    
    bws <- scbandwidth(bw = bws$bw,
                       regtype = regtype,
                       basis = if (is.null(bws$basis)) "glp" else bws$basis,
                       degree = bws$degree,
                       bernstein.basis = bws$bernstein.basis,
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
           backfit.iterate,
           backfit.maxiter,
           backfit.tol,
           bandwidth.compute = TRUE,
           basis,
           bernstein.basis,
           bwmethod,
           bwscaling,
           bwtype,
           ckerbound,
           ckerlb,
           ckerorder,
           ckertype,
           ckerub,
           cv.iterate,
           cv.num.iterations,
           degree,
           nmulti,
           okertype,
           optim.abstol,
           optim.maxattempts,
           optim.maxit,
           optim.method,
           optim.reltol,
           random.seed,
           regtype,
           ukertype,
           ...){
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    if (!missing(bwmethod) && identical(match.arg(bwmethod, c("cv.ls", "manual")), "manual") &&
        missing(bws))
      stop("bwmethod='manual' requires argument 'bws'")


    miss.z <- missing(zdat)
    xdat <- toFrame(xdat)
    
    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector or a factor")

    if(!miss.z)
      zdat <- toFrame(zdat)

    ## first grab dummy args for scbandwidth() and perform 'bootstrap'
    ## bandwidth call

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("regtype", "basis", "degree", "bernstein.basis",
               "bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
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
