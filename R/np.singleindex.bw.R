
npindexbw <-
  function(...){
    args = list(...)
    if (is(args[[1]],"formula"))
      UseMethod("npindexbw",args[[1]])
    else if (!is.null(args$formula))
      UseMethod("npindexbw",args$formula)
    else
      UseMethod("npindexbw",args[[which(names(args)=="bws")[1]]])
  }

npindexbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    orig.class <- if (missing(data))
      sapply(eval(attr(terms(formula), "variables"), environment(formula)),class)
    else sapply(eval(attr(terms(formula), "variables"), data, environment(formula)),class)

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    mf[[1]] <- as.name("model.frame")

    if(all(orig.class == "ts")){
      args <- (as.list(attr(terms(formula), "variables"))[-1])
      formula <- terms(formula)
      attr(formula, "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
      mf[["formula"]] <- formula
    }else if(any(orig.class == "ts")){
      arguments <- (as.list(attr(terms(formula), "variables"))[-1])
      arguments.normal <- arguments[which(orig.class != "ts")]
      arguments.timeseries <- arguments[which(orig.class == "ts")]

      ix <- sort(c(which(orig.class == "ts"),which(orig.class != "ts")),index.return = TRUE)$ix
      formula <- terms(formula)
      attr(formula, "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
      mf[["formula"]] <- formula
    }
    
    mf <- eval(mf, parent.frame())

    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    tbw <- npindexbw(xdat = xdat, ydat = ydat, ...)

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

npindexbw.NULL <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws, ...){

    xdat <- toFrame(xdat)

    bws <- double(ncol(xdat)+1)

    tbw <- npindexbw.default(xdat = xdat,
                             ydat = ydat,
                             bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <- updateBwNameMetadata(nameList =
                                list(ynames = deparse(substitute(ydat))),
                                bws = tbw)

    tbw
  }

npindexbw.default <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws, bandwidth.compute = TRUE,
           nmulti, random.seed, optim.method, optim.maxattempts,
           optim.reltol, optim.abstol, optim.maxit, only.optimize.beta, ...){

    xdat <- toFrame(xdat)

    if(!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (ncol(xdat) < 2) {
      if (coarseclass(xdat[,1]) != "numeric")
        stop("xdat must contain at least one continuous variable")

      warning(paste("xdat has one dimension. Using a single index model to reduce",
                    "dimensionality is unnecessary."))
    }

    if (coarseclass(bws) != "numeric" | length(bws) != ncol(xdat)+1)
      stop(paste("manually specified 'bws' must be a numeric vector of length ncol(xdat)+1.",
                 "See documentation for details."))

    tbw <- sibandwidth(beta = bws[1:ncol(xdat)],
                       h = bws[ncol(xdat)+1], ...,
                       nobs = dim(xdat)[1],
                       xdati = untangle(xdat),
                       ydati = untangle(data.frame(ydat)),
                       xnames = names(xdat),
                       ynames = deparse(substitute(ydat)),
                       bandwidth = bws[ncol(xdat)+1],
                       bandwidth.compute = bandwidth.compute)

    if (tbw$method == "kleinspady" & !setequal(ydat,c(0,1)))
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("nmulti","random.seed", "optim.method", "optim.maxattempts",
               "optim.reltol", "optim.abstol", "optim.maxit", "only.optimize.beta")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (bandwidth.compute)
      tbw <- eval(parse(text=paste("npindexbw.sibandwidth(xdat=xdat, ydat=ydat, bws=tbw",
                          ifelse(any.m, ",",""),
                          paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                          ")")))

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }

npindexbw.sibandwidth <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws, bandwidth.compute = TRUE, nmulti, random.seed = 42,
           optim.method = c("Nelder-Mead", "BFGS", "CG"),
           optim.maxattempts = 10,
           optim.reltol = sqrt(.Machine$double.eps),
           optim.abstol = .Machine$double.eps,
           optim.maxit = 500,
           only.optimize.beta = FALSE,
           ...){

    ## Save seed prior to setting

    if(exists(".Random.seed", .GlobalEnv)) {
      save.seed <- get(".Random.seed", .GlobalEnv)
      exists.seed = TRUE
    } else {
      exists.seed = FALSE
    }

    set.seed(random.seed)

    xdat = toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- min(5,ncol(xdat))
    }

    if (bws$method == "kleinspady" & !setequal(ydat,c(0,1)))
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    if (ncol(xdat) < 2) {
      if (coarseclass(xdat[,1]) != "numeric")
        stop("xdat must contain at least one continuous variable")

      warning(paste("xdat has one dimension. Using a single index model to reduce",
                    "dimensionality is unnecessary."))
    }

    optim.method <- match.arg(optim.method)

    ## catch and destroy NA's
    goodrows = 1:dim(xdat)[1]
    rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows]

    ## convert to numeric
    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat = toMatrix(xdat)

    total.time <-
      system.time({

        if(bandwidth.compute){

          ## Note - there are two methods currently implemented, Ichimura's
          ## least squares approach and Klein and Spady's likelihood approach.

          ##We define ichimura's objective function. Since we normalize beta_1
          ##to equal 1, we only worry about beta_2...beta_k in the index
          ##function for these internals, and when computing the leave-one-out
          ##objective function use c(1,beta). However, we do indeed return
          ##c(1,beta) which can be used in the index.model function above.

          ichimuraMaxPenalty <- 10*mean(ydat^2)

          ichimura <- function(param) {

            ##Define the leave-one-out objective function, sum (y - \hat
            ## G(X\hat\beta))^2. We let beta denote beta_2...beta_k (first k-1
            ## parameters in `param') and then let h denote the kth column.
            if (ncol(xdat) == 1)
              beta = numeric(0)
            else
              beta <- param[1:(ncol(xdat)-1)]

            h <- param[ncol(xdat)]

            ## Next we define the sum of squared leave-one-out residuals

            sum.squares.leave.one.out <- function(xdat,ydat,beta,h) {

              ## Normalize beta_1 = 1 hence multiply X by c(1,beta)

              index <- as.matrix(xdat) %*% c(1,beta)

              ## One call to npksum to avoid repeated computation of the
              ## product kernel (the expensive part)

              W <- as.matrix(data.frame(ydat,1))

              tww <- npksum(txdat=index,
                            tydat=W,
                            weights=W,
                            leave.one.out=TRUE,
                            bandwidth.divide=TRUE,
                            bws=c(h),
                            ckertype = bws$ckertype,
                            ckerorder = bws$ckerorder)$ksum

              fit.loo <- tww[1,2,]/NZD(tww[2,2,])

              t.ret <- mean((ydat-fit.loo)^2)
              return(t.ret)

            }

            ## For the objective function, we require a positive bandwidth, so
            ## return an infinite penalty for negative h

            if(h > 0) {
              return(sum.squares.leave.one.out(xdat,ydat,beta,h))
            } else {
              return(ichimuraMaxPenalty)
            }

          }

          ichimura.nobw <- function(param,h){ return(ichimura(c(param,h))) }

          ## We define ichimura's objective function. Since we normalize beta_1
          ## to equal 1, we only worry about beta_2...beta_k in the index
          ## function for these internals, and when computing the leave-one-out
          ## objective function use c(1,beta). However, we do indeed return
          ## c(1,beta) which can be used in the index.model function above.

          kleinspadyFloor <- sqrt(.Machine$double.eps)

          kleinspady <- function(param) {

            ## Define the leave-one-out objective function, sum (y - \hat
            ## G(X\hat\beta))^2. We let beta denote beta_2...beta_k (first k-1
            ## parameters in `param') and then let h denote the kth column.
            if (ncol(xdat) == 1)
              beta = numeric(0)
            else
              beta <- param[1:(ncol(xdat)-1)]

            h <- param[ncol(xdat)]

            ## Next we define the sum of logs

            sum.log.leave.one.out <- function(xdat,ydat,beta,h) {

              ## Normalize beta_1 = 1 hence multiply X by c(1,beta)

              index <- as.matrix(xdat) %*% c(1,beta)

              ## One call to npksum to avoid repeated computation of the
              ## product kernel (the expensive part)

              W <- as.matrix(data.frame(ydat,1))

              tww <- npksum(txdat=index,
                            tydat=W,
                            weights=W,
                            leave.one.out=TRUE,
                            bandwidth.divide=TRUE,
                            bws=c(h),
                            ckertype = bws$ckertype,
                            ckerorder = bws$ckerorder)$ksum

              ks.loo <- tww[1,2,]/NZD(tww[2,2,])

              ## Avoid taking log of zero (ks.loo = 0 or 1 since we take
              ## the log of ks.loo and the log of 1-ks.loo)

              ks.loo[which(ks.loo < kleinspadyFloor)] <- kleinspadyFloor
              ks.loo[which(ks.loo > 1- kleinspadyFloor)] <- 1-kleinspadyFloor

              ## Maximize the log likelihood, therefore minimize minus...
              t.ret <- - mean(ifelse(ydat==1,log(ks.loo),log(1-ks.loo)))
              return(t.ret)

            }

            ## For the objective function, we require a positive bandwidth, so
            ## return an infinite penalty for negative h

            if(h > 0) {
              return(sum.log.leave.one.out(xdat,ydat,beta,h))
            } else {
              ## No natural counterpart to var of y here, unlike Ichimura above...
              return(sqrt(.Machine$double.xmax))
            }

          }

          kleinspady.nobw <- function(param,h){ return(kleinspady(c(param,h))) }
          ## Now we implement multistarting

          fval.min <- .Machine$double.xmax
          numimp <- 0
          fval.value <- numeric(nmulti)

          console <- newLineConsole()

          if(bws$method == "ichimura"){
            optim.fn <- if(only.optimize.beta) ichimura.nobw else ichimura
            optim.control <- list(abstol=optim.abstol,
                                  reltol=optim.reltol,
                                  maxit=optim.maxit)
          } else if(bws$method == "kleinspady"){
            optim.fn <- if(only.optimize.beta) kleinspady.nobw else  kleinspady
            optim.control <- list(reltol=optim.reltol,maxit=optim.maxit)
          }

          for(i in 1:nmulti) {

            console <- printPush(paste(sep="", "Multistart ", i, " of ", nmulti, "..."), console)
            ##cv.console <- newLineConsole(console)

            n <- nrow(xdat)

            ## We use the nlm command to minimize the objective function using
            ## starting values. Note that since we normalize beta_1=1 here beta
            ## is the k-1 vector containing beta_2...beta_k

            if(i == 1) {

              ## Initial values taken from OLS fit with a constant used for
              ## multistart 1
              ols.fit <- lm(ydat~xdat,x=TRUE)
              fit <- fitted(ols.fit)

              if (ncol(xdat) != 1){
                if (setequal(bws$beta[2:ncol(xdat)],c(0)))
                  beta <- coef(ols.fit)[3:ncol(ols.fit$x)]
                else
                  beta = bws$beta[2:ncol(xdat)]
              } else { beta = numeric(0) }

              if (bws$bw == 0)
                h <- 1.059224*EssDee(fit)*n^(-1/5)
              else
                h <- bws$bw
            } else {
              ## Random initialization used for remaining multistarts

              beta.length <- length(coef(ols.fit)[3:ncol(ols.fit$x)])
              beta <- runif(beta.length,min=0.5,max=1.5)*coef(ols.fit)[3:ncol(ols.fit$x)]
              if (!only.optimize.beta)
                h <- runif(1,min=0.5,max=1.5)*EssDee(fit)*n^(-1/5)
            }

            optim.parm <- if(only.optimize.beta) beta else c(beta,h)

            topt <- parse(text=paste("optim(optim.parm,fn=optim.fn,gr=NULL,method=optim.method,control=optim.control", ifelse(only.optimize.beta, ',h)',')')))
            suppressWarnings(optim.return <- eval(topt))
            attempts <- 0
            while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
              attempts <- attempts + 1
              beta.length <- length(coef(ols.fit)[3:ncol(ols.fit$x)])
              beta <- runif(beta.length,min=0.5,max=1.5)*coef(ols.fit)[3:ncol(ols.fit$x)]
              if(!only.optimize.beta)
                h <- runif(1,min=0.5,max=1.5)*EssDee(fit)*n^(-1/5)

              if(optim.return$convergence == 1){
                if(optim.control$maxit < (2^32/10))
                  optim.control$maxit <- 10*optim.control$maxit
                else
                  stop(paste("optim failed to converge after optim.maxattempts = ", optim.maxattempts, " iterations."))
              }

              if(optim.return$convergence == 10){
                optim.control$reltol <-  10.0*optim.control$reltol
                if(!is.null(optim.control$abstol))
                  optim.control$abstol <-  10.0*optim.control$abstol
              }
              
              suppressWarnings(optim.return <- eval(topt))
            }

            if(optim.return$convergence != 0)
              stop(paste("optim failed to converge after optim.maxattempts = ", optim.maxattempts, " iterations."))

            fval.value[i] <- optim.return$value
            if(optim.return$value < fval.min) {
              param <- if(only.optimize.beta) c(optim.return$par,h) else optim.return$par
              fval.min <- optim.return$value
              numimp <- numimp + 1
              best <- i
            }

            if(i != nmulti)
              console <- printPop(console)
          }
          console <- printClear(console)

          bws$beta <- if(ncol(xdat) == 1) 1.0
          else c(1,param[1:(ncol(xdat)-1)])
          bws$bw <- param[ncol(xdat)]
          bws$fval <- fval.min
          bws$ifval <- best
          bws$numimp <- numimp
          bws$fval.vector <- fval.value
        }
      })[1]
    ## Return a list with beta (we append the restricted value of
    ## beta_1=1), the bandwidth h, the value of the objective function at
    ## its minimum, the number of restarts that resulted in an improved
    ## value of the objective function, the restart that resulted in the
    ## smallest value, and the vector of objective function values.

    ## Restore seed

    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

    bws <- sibandwidth(beta = bws$beta,
                       h = bws$bw,
                       method = bws$method,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       fval = bws$fval,
                       ifval = bws$ifval,
                       numimp = bws$numimp,
                       fval.vector = bws$fval.vector,
                       nobs = bws$nobs,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       bandwidth = bws$bw,
                       rows.omit = rows.omit,
                       bandwidth.compute = bandwidth.compute,
                       optim.method = optim.method,
                       only.optimize.beta = only.optimize.beta,
                       total.time = total.time)

    bws

  }
