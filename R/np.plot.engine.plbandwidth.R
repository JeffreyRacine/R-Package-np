.np_plot_plbandwidth_engine <-
  function(bws,
           xdat,
           ydat,
           zdat,
           data = NULL,
           xq = 0.5,
           zq = 0.5,
           xtrim = 0.0,
           ztrim = 0.0,
           neval = 50,
           common.scale = TRUE,
           perspective = TRUE,
           gradients = FALSE,
           main = NULL,
           type = NULL,
           border = NULL,
           cex.axis = NULL,
           cex.lab = NULL,
           cex.main = NULL,
           cex.sub = NULL,
           col = NULL,
           ylab = NULL,
           xlab = NULL,
           zlab = NULL,
           sub = NULL,
           ylim = NULL,
           xlim = NULL,
           zlim = NULL,
           lty = NULL,
           lwd = NULL,
           theta = 0.0,
           phi = 10.0,
           view = c("rotate","fixed"),
           plot.behavior = c("plot","plot-data","data"),
           plot.errors.method = c("none","bootstrap","asymptotic"),
           plot.errors.boot.method = c("inid", "fixed", "geom"),
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen = NULL,
           plot.errors.boot.num = 399,
           plot.errors.center = c("estimate","bias-corrected"),
           plot.errors.type = c("pmzsd","pointwise","bonferroni","simultaneous","all"),
           plot.errors.alpha = 0.05,
           plot.errors.style = c("band","bar"),
           plot.errors.bar = c("|","I"),
           plot.errors.bar.num = min(neval,25),
           plot.bxp = FALSE,
           plot.bxp.out = TRUE,
           plot.par.mfrow = TRUE,
           ...,
           random.seed){

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    scalar_default <- function(value, default) {
      if (is.null(value)) default else value
    }

    plot.par.mfrow.opt <- getOption("plot.par.mfrow")
    if(!is.null(plot.par.mfrow.opt))
        plot.par.mfrow <- plot.par.mfrow.opt

    if(!missing(gradients))
      stop("gradients not supported with partially linear models. Coefficients may be extracted with coef()")

    miss.xyz = c(missing(xdat), missing(ydat), missing(zdat))
    
    if (any(miss.xyz) && !all(miss.xyz))
      stop("one of, but not both, xdat and ydat was specified")
    else if(all(miss.xyz) && !is.null(bws$formula)){
      tt <- terms(bws)
      m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)

      tmf.xf <- tmf.x <- tmf <- bws$call[c(1,m)]
      tmf.xf[[1]] <- tmf[[1]] <- as.name("model.frame")
      tmf[["formula"]] <- tt
      mf.args <- as.list(tmf)[-1L]
      umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

      bronze <- lapply(bws$chromoly, paste, collapse = " + ")

      tmf.xf[["formula"]] <- as.formula(paste(" ~ ", bronze[[2]]),
                                      env = environment(formula))
      mf.xf.args <- as.list(tmf.xf)[-1L]
      tmf.xf <- do.call(stats::model.frame, mf.xf.args, envir = environment(tt))
      
      ydat <- model.response(tmf)
      xdat <- tmf.xf

      zdat <- tmf[, bws$chromoly[[3]], drop = FALSE]
    } else {
      if(all(miss.xyz) && !is.null(bws$call)){
        xdat <- data.frame(.np_eval_bws_call_arg(bws, "xdat"))
        ydat = .np_eval_bws_call_arg(bws, "ydat")
        zdat <- data.frame(.np_eval_bws_call_arg(bws, "zdat"))
      }
      xdat = toFrame(xdat)
      zdat = toFrame(zdat)
      
      ## catch and destroy NA's
      keep.rows <- rep_len(TRUE, nrow(xdat))
      rows.omit <- attr(na.omit(data.frame(xdat, ydat, zdat)), "na.action")
      if (length(rows.omit) > 0L)
        keep.rows[as.integer(rows.omit)] <- FALSE

      if (!any(keep.rows))
        stop("Training data has no rows without NAs")

      xdat <- xdat[keep.rows,,drop = FALSE]
      ydat <- ydat[keep.rows]
      zdat <- zdat[keep.rows,,drop = FALSE]
    }


    nxcon = sum(bws$xdati$icon)
    nxuno = sum(bws$xdati$iuno)
    nxord = sum(bws$xdati$iord)

    nzcon = sum(bws$zdati$icon)
    nzuno = sum(bws$zdati$iuno)
    nzord = sum(bws$zdati$iord)

    xq = double(bws$xndim)+xq
    zq = double(bws$zndim)+zq
    
    xtrim = double(bws$xndim)+xtrim
    ztrim = double(bws$zndim)+ztrim


    if (missing(plot.errors.method) &
        any(!missing(plot.errors.boot.num), !missing(plot.errors.boot.method),
            !missing(plot.errors.boot.blocklen))){
      warning(paste("plot.errors.method must be set to 'bootstrap' to use bootstrapping.",
                    "\nProceeding without bootstrapping."))
    }
    
    normalized.opts <- .np_plot_normalize_common_options(
      plot.behavior = plot.behavior,
      plot.errors.method = plot.errors.method,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.wild = plot.errors.boot.wild,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      plot.errors.style = plot.errors.style,
      plot.errors.bar = plot.errors.bar,
      xdat = xdat,
      common.scale = common.scale,
      ylim = ylim,
      allow_asymptotic_quantile = FALSE
    )

    plot.behavior <- .npRmpi_plot_behavior_for_rank(normalized.opts$plot.behavior)
    plot.errors.method <- normalized.opts$plot.errors.method
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method,
                                              where = "plot(...)",
                                              allow.direct.bootstrap = TRUE)
    plot.errors.boot.method <- normalized.opts$plot.errors.boot.method
    plot.errors.boot.wild <- normalized.opts$plot.errors.boot.wild
    plot.errors.boot.blocklen <- normalized.opts$plot.errors.boot.blocklen
    plot.errors.center <- normalized.opts$plot.errors.center
    plot.errors.type <- normalized.opts$plot.errors.type
    plot.errors.alpha <- normalized.opts$plot.errors.alpha
    plot.errors.style <- normalized.opts$plot.errors.style
    plot.errors.bar <- normalized.opts$plot.errors.bar
    common.scale <- normalized.opts$common.scale

    if (plot.errors.method == "asymptotic") {
      warning(paste("asymptotic errors are not supported with partially linear regression.\n",
                    "Proceeding without calculating errors"))
      plot.errors.method = "none"
    }

    plot.errors = (plot.errors.method != "none")

    if ((nxcon + nxord == 1) && (nzcon + nzord == 1) && (nxuno + nzuno == 0) &&
        perspective & !gradients & !any(xor(bws$xdati$iord, bws$xdati$inumord)) &
        !any(xor(bws$zdati$iord, bws$zdati$inumord))){

      view = match.arg(view)
      rotate = (view == "rotate")
      
      if (is.ordered(xdat[,1])){
        x1.eval = bws$xdati$all.ulev[[1]]
        x1.neval = length(x1.eval)
      } else {
        x1.neval = neval
        qi = trim.quantiles(xdat[,1], xtrim[1])
        x1.eval = seq(qi[1], qi[2], length.out = x1.neval)
      }

      if (is.ordered(zdat[,1])){
        z1.eval = bws$zdati$all.ulev[[1]]
        z1.neval = length(z1.eval)
      } else {
        z1.neval = neval
        qi = trim.quantiles(zdat[,1], ztrim[1])
        z1.eval = seq(qi[1], qi[2], length.out = z1.neval)
      }

      x.eval <- expand.grid(x1.eval, z1.eval)

      if (bws$xdati$iord[1])
        x1.eval <- (bws$xdati$all.dlev[[1]])[as.integer(x1.eval)]
      
      if (bws$zdati$iord[1])
        z1.eval <- (bws$zdati$all.dlev[[1]])[as.integer(z1.eval)]

      tobj = npplreg(txdat = xdat, tydat = ydat, tzdat = zdat,
        exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,2, drop = FALSE], bws = bws)

      terr = matrix(data = NA, nrow = nrow(x.eval), ncol = 3)
      
      treg = matrix(data = tobj$mean,
        nrow = x1.neval, ncol = z1.neval, byrow = FALSE)

      if (plot.errors.method == "bootstrap"){
        terr <- compute.bootstrap.errors(
          xdat = xdat, ydat = ydat, zdat = zdat,
          exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,2, drop = FALSE],
          gradients = FALSE,
          slice.index = 0,
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.wild = plot.errors.boot.wild,
          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
          plot.errors.boot.num = plot.errors.boot.num,
          plot.errors.center = plot.errors.center,
          plot.errors.type = plot.errors.type,
          plot.errors.alpha = plot.errors.alpha,
          bws = bws)[["boot.err"]]

        pc = (plot.errors.center == "bias-corrected")

        lerr = matrix(data = if(pc) {terr[,3]} else {tobj$mean}
          -terr[,1],
          nrow = x1.neval, ncol = z1.neval, byrow = FALSE)

        herr = matrix(data = if(pc) {terr[,3]} else {tobj$mean}
          +terr[,2],
          nrow = x1.neval, ncol = z1.neval, byrow = FALSE)

      }
      
      if(is.null(zlim)) {
          zlim =
              if (plot.errors)
                  c(min(lerr),max(herr))
              else
                  c(min(tobj$mean),max(tobj$mean))
      }
        
      if (plot.behavior != "plot"){
        r1 = plregression(bws = bws, xcoef = tobj$xcoef,
          xcoefvcov = vcov(tobj),
          xcoeferr = tobj$xcoeferr,
          evalx =  x.eval[,1],
          evalz =  x.eval[,2],
          mean = tobj$mean,
          ntrain = dim(xdat)[1],
          trainiseval = FALSE,
          xtra=c(tobj$RSQ,tobj$MSE,0,0,0,0))

        r1$merr = NA
        r1$bias = NA

        if (plot.errors)
          r1$merr = terr[,1:2]
        
        if (plot.errors.center == "bias-corrected")
          r1$bias = terr[,3] - treg

        if (plot.behavior == "data")
          return ( list(r1 = r1) )

      }

      dtheta = 5.0
      dphi = 10.0

      persp.col = if (plot.errors) FALSE else scalar_default(col, "lightblue")
      
      ##for (j in 0:((50 %/% dphi - 1)*rotate)*dphi+phi){
        for (i in 0:((360 %/% dtheta - 1)*rotate)*dtheta+theta){
          if (plot.errors){
            persp(x1.eval,
                  z1.eval,
                  lerr,
                  zlim = zlim,
                  cex.axis = scalar_default(cex.axis, par()$cex.axis),
                  cex.lab = scalar_default(cex.lab, par()$cex.lab),
                  cex.main = scalar_default(cex.main, par()$cex.main),
                  cex.sub = scalar_default(cex.sub, par()$cex.sub),
                  col = persp.col,
                  border = scalar_default(border, "grey"),
                  ticktype = "detailed",
                  xlab = "",
                  ylab = "",
                  zlab = "",
                  theta = i,
                  phi = phi,
                  lwd = scalar_default(lwd, par()$lwd))
            par(new = TRUE)
          }

          persp(x1.eval,
                z1.eval,
                treg,
                zlim = zlim,
                col = persp.col,
                border = scalar_default(border, "black"),
                ticktype = "detailed",
                cex.axis = scalar_default(cex.axis, par()$cex.axis),
                cex.lab = scalar_default(cex.lab, par()$cex.lab),
                cex.main = scalar_default(cex.main, par()$cex.main),
                cex.sub = scalar_default(cex.sub, par()$cex.sub),
                xlab = scalar_default(xlab, gen.label(names(xdat)[1], "X1")),
                ylab = scalar_default(ylab, gen.label(names(xdat)[2], "Z1")),
                zlab = scalar_default(zlab, gen.label(names(ydat),"Conditional Mean")),
                theta = i,
                phi = phi,
                main = gen.tflabel(!is.null(main), main, paste("[theta= ", i,", phi= ", phi,"]", sep="")))

          if (plot.errors){
            par(new = TRUE)
            persp(x1.eval,
                  z1.eval,
                  herr,
                  zlim = zlim,
                  cex.axis = scalar_default(cex.axis, par()$cex.axis),
                  cex.lab = scalar_default(cex.lab, par()$cex.lab),
                  cex.main = scalar_default(cex.main, par()$cex.main),
                  cex.sub = scalar_default(cex.sub, par()$cex.sub),
                  col = persp.col,
                  border = scalar_default(border, "grey"),
                  ticktype = "detailed",
                  xlab = "",
                  ylab = "",
                  zlab = "",
                  theta = i,
                  phi = phi,
                  lwd = scalar_default(lwd, par()$lwd))
          }

          Sys.sleep(0.5)
        }
      ##}

      if (plot.behavior == "plot-data")
        return ( list(r1 = r1) )

    } else {
##      stop("not yet supported!")
      
      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(bws$xndim + bws$zndim),cex=par()$cex)

      x.ev = xdat[1,,drop = FALSE]
      z.ev = zdat[1,,drop = FALSE]

      for (i in seq_len(bws$xndim))
        x.ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      for (i in seq_len(bws$zndim))
        z.ev[1,i] = uocquantile(zdat[,i], prob=zq[i])


      maxneval = max(c(sapply(xdat,nlevels), sapply(zdat,nlevels), neval))

      ## Preserve original data types (e.g., factors) for evaluation data
      exdat = xdat[rep(1, maxneval), , drop = FALSE]
      ezdat = zdat[rep(1, maxneval), , drop = FALSE]

      for (i in seq_len(bws$xndim))
        exdat[,i] = x.ev[1,i]

      for (i in seq_len(bws$zndim))
        ezdat[,i] = z.ev[1,i]

      if (common.scale){
        data.eval = matrix(data = NA, nrow = maxneval,
          ncol = (bws$xndim + bws$zndim))
        
        data.err = matrix(data = NA, nrow = maxneval,
          ncol = 3*(bws$xndim + bws$zndim))
        data.err.all = vector("list", bws$xndim + bws$zndim)

        allei = as.data.frame(matrix(data = NA, nrow = maxneval,
          ncol = bws$xndim + bws$zndim))

        all.bxp = list()
      }

      all.isFactor = c(sapply(xdat, is.factor), sapply(zdat, is.factor))

      plot.out = list()

      temp.err = matrix(data = NA, nrow = maxneval, ncol = 3)
      temp.mean = replicate(maxneval, NA)

      ## plotting expressions
      
      plot.bootstrap = plot.errors.method == "bootstrap"

      plotOnEstimate = (plot.errors.center == "estimate")

      plot.index = 0
      xOrZ = "x"

      for (i in seq_len(bws$xndim)){
        plot.index = plot.index + 1
        temp.err[,] = NA
        temp.mean[] =  NA
        temp.boot = list()
        temp.all.err <- NULL

        xi.factor = all.isFactor[plot.index]
        
        if (xi.factor){
          ei = levels(xdat[,i])
          ei = factor(ei, levels = ei)
          xi.neval = length(ei)
        } else {
          xi.neval = neval
          qi = trim.quantiles(xdat[,i], xtrim[i])
          ei = seq(qi[1], qi[2], length.out = neval)
        }

        if (xi.neval < maxneval){
          ei[(xi.neval+1):maxneval] = NA
        }


        tobj <- npplreg(txdat = xdat, tydat = ydat, tzdat = zdat,
          exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
          ezdat = ezdat[1:xi.neval,, drop = FALSE],
          bws = bws)

        temp.mean[1:xi.neval] = tobj$mean

        if (plot.errors){
          if (plot.errors.method == "bootstrap"){
            temp.boot <- compute.bootstrap.errors(
                      xdat = xdat,
                      ydat = ydat,
                      zdat = zdat,
                      exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                      ezdat = ezdat[1:xi.neval,, drop = FALSE],
                      gradients = gradients,
                      slice.index = plot.index,
                      plot.errors.boot.method = plot.errors.boot.method,
                      plot.errors.boot.wild = plot.errors.boot.wild,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
                      bws = bws)
            temp.err[1:xi.neval,] <- temp.boot[["boot.err"]]
            temp.all.err <- temp.boot[["boot.all.err"]]
            temp.boot <- temp.boot[["bxp"]]
            if (!plot.bxp.out){
              temp.boot$out <- numeric()
              temp.boot$group <- integer()
            }
          }
        }

        if (common.scale){
          allei[,plot.index] = ei
          data.eval[, plot.index] = temp.mean
          if (plot.errors){
            all.bxp[plot.index] = NA
            all.bxp[[plot.index]] = temp.boot

            data.err[, c(3*plot.index-2,3*plot.index-1,3*plot.index)] = temp.err
            data.err.all[[plot.index]] = temp.all.err
          }
        } else if (plot.behavior != "data") {
          ## plot evaluation
          plot.fun <- if (xi.factor) {
            .np_plot_panel_fun(plot.bootstrap = plot.bootstrap, plot.bxp = plot.bxp)
          } else {
            plot
          }
          plot.args <- list()
          if (xi.factor) {
            if (plot.bootstrap && plot.bxp) plot.args$z <- temp.boot else plot.args$f <- ei
          } else {
            plot.args$x <- ei
          }
          if (!(xi.factor && plot.bootstrap && plot.bxp))
            plot.args$y <- temp.mean
          if (plot.errors)
            plot.args$ylim <- if (plot.errors.type == "all")
              compute.all.error.range(if (plotOnEstimate) temp.mean else temp.err[,3], temp.all.err)
            else
              c(min(na.omit(c(temp.mean - temp.err[,1], temp.err[,3] - temp.err[,1]))),
                max(na.omit(c(temp.mean + temp.err[,2], temp.err[,3] + temp.err[,2]))))
          plot.args$xlab <- gen.label(if (xOrZ == "x") bws$xnames[i] else bws$znames[i],
                                      paste(toupper(xOrZ), i, sep = ""))
          plot.args$ylab <- paste(if (gradients) paste("Gradient Component ", i, " of", sep = "") else "",
                                  gen.label(bws$ynames, "Conditional Mean"))
          if (!xi.factor) {
            plot.args$type <- scalar_default(type, "l")
            plot.args$lty <- scalar_default(lty, par()$lty)
            plot.args$col <- scalar_default(col, par()$col)
            plot.args$lwd <- scalar_default(lwd, par()$lwd)
            plot.args$cex.axis <- scalar_default(cex.axis, par()$cex.axis)
            plot.args$cex.lab <- scalar_default(cex.lab, par()$cex.lab)
            plot.args$cex.main <- scalar_default(cex.main, par()$cex.main)
            plot.args$cex.sub <- scalar_default(cex.sub, par()$cex.sub)
          }
          plot.args$main <- scalar_default(main, "")
          plot.args$sub <- scalar_default(sub, "")
          do.call(plot.fun, plot.args)

          ## error plotting evaluation
          if (plot.errors && !(xi.factor && plot.bootstrap && plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = na.omit(if (plotOnEstimate) temp.mean else temp.err[,3]),
                all.err = temp.all.err,
                xi.factor = xi.factor)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(ei)),
                ely = if (plotOnEstimate) na.omit(temp.mean - temp.err[,1]) else na.omit(temp.err[,3] - temp.err[,1]),
                ehy = if (plotOnEstimate) na.omit(temp.mean + temp.err[,2]) else na.omit(temp.err[,3] + temp.err[,2]),
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = if (xi.factor) 1 else 2
              )
              do.call(draw.errors, draw.args)
            }

          }
        }

        if (plot.behavior != "plot") {
          plot.out[plot.index] = NA
          if (gradients){
          } else {
            plot.out[[plot.index]] =
              plregression(bws = bws, xcoef = tobj$xcoef, xcoefvcov = vcov(tobj),
                           xcoeferr = tobj$xcoeferr,
                           evalx = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                           evalz = ezdat[1:xi.neval,, drop = FALSE],
                           mean = na.omit(temp.mean),
                           ntrain = dim(xdat)[1],
                           trainiseval = FALSE,
                           xtra = c(tobj$RSQ, tobj$MSE, 0, 0, 0, 0))
            plot.out[[plot.index]]$merr = NA
            plot.out[[plot.index]]$bias = NA

            if (plot.errors)
              plot.out[[plot.index]]$merr = temp.err[,1:2]

            if (plot.errors.center == "bias-corrected")
              plot.out[[plot.index]]$bias = temp.err[,3] - temp.mean
            plot.out[[plot.index]]$bxp = temp.boot
          }
        }
      }

      xOrZ = "z"
      for (i in seq_len(bws$zndim)){
        plot.index = plot.index + 1
        temp.err[,] = NA
        temp.mean[] =  NA
        temp.boot = list()
        temp.all.err <- NULL

        xi.factor = all.isFactor[plot.index]
        
        if (xi.factor){
          ei = levels(zdat[,i])
          ei = factor(ei, levels = ei)
          xi.neval = length(ei)
        } else {
          xi.neval = neval
          qi = trim.quantiles(zdat[,i], ztrim[i])
          ei = seq(qi[1], qi[2], length.out = neval)
        }

        if (xi.neval < maxneval){
          ei[(xi.neval+1):maxneval] = NA
        }

        tobj <- npplreg(txdat = xdat, tydat = ydat, tzdat = zdat,
          exdat = exdat[1:xi.neval,, drop = FALSE],
          ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
          bws = bws)

        temp.mean[1:xi.neval] = tobj$mean

        if (plot.errors){
          if (plot.errors.method == "bootstrap"){
            temp.boot <- compute.bootstrap.errors(
                      xdat = xdat,
                      ydat = ydat,
                      zdat = zdat,
                      exdat = exdat[1:xi.neval,, drop = FALSE],
                      ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
                      gradients = gradients,
                      slice.index = plot.index,
                      plot.errors.boot.method = plot.errors.boot.method,
                      plot.errors.boot.wild = plot.errors.boot.wild,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
                      bws = bws)
            temp.err[1:xi.neval,] <- temp.boot[["boot.err"]]
            temp.all.err <- temp.boot[["boot.all.err"]]
            temp.boot <- temp.boot[["bxp"]]
            if (!plot.bxp.out){
              temp.boot$out <- numeric()
              temp.boot$group <- integer()
            }
          }
        }

        if (common.scale){
          allei[,plot.index] = ei
          data.eval[, plot.index] = temp.mean
          if (plot.errors){
            all.bxp[plot.index] = NA
            all.bxp[[plot.index]] = temp.boot

            data.err[, c(3*plot.index-2,3*plot.index-1,3*plot.index)] = temp.err
            data.err.all[[plot.index]] = temp.all.err
          }
        } else if (plot.behavior != "data") {
          ## plot evaluation
          plot.fun <- if (xi.factor) {
            .np_plot_panel_fun(plot.bootstrap = plot.bootstrap, plot.bxp = plot.bxp)
          } else {
            plot
          }
          plot.args <- list()
          if (xi.factor) {
            if (plot.bootstrap && plot.bxp) plot.args$z <- temp.boot else plot.args$f <- ei
          } else {
            plot.args$x <- ei
          }
          if (!(xi.factor && plot.bootstrap && plot.bxp))
            plot.args$y <- temp.mean
          if (plot.errors)
            plot.args$ylim <- if (plot.errors.type == "all")
              compute.all.error.range(if (plotOnEstimate) temp.mean else temp.err[,3], temp.all.err)
            else
              c(min(na.omit(c(temp.mean - temp.err[,1], temp.err[,3] - temp.err[,1]))),
                max(na.omit(c(temp.mean + temp.err[,2], temp.err[,3] + temp.err[,2]))))
          plot.args$xlab <- gen.label(if (xOrZ == "x") bws$xnames[i] else bws$znames[i],
                                      paste(toupper(xOrZ), i, sep = ""))
          plot.args$ylab <- paste(if (gradients) paste("Gradient Component ", i, " of", sep = "") else "",
                                  gen.label(bws$ynames, "Conditional Mean"))
          if (!xi.factor) {
            plot.args$type <- scalar_default(type, "l")
            plot.args$lty <- scalar_default(lty, par()$lty)
            plot.args$col <- scalar_default(col, par()$col)
            plot.args$lwd <- scalar_default(lwd, par()$lwd)
            plot.args$cex.axis <- scalar_default(cex.axis, par()$cex.axis)
            plot.args$cex.lab <- scalar_default(cex.lab, par()$cex.lab)
            plot.args$cex.main <- scalar_default(cex.main, par()$cex.main)
            plot.args$cex.sub <- scalar_default(cex.sub, par()$cex.sub)
          }
          plot.args$main <- scalar_default(main, "")
          plot.args$sub <- scalar_default(sub, "")
          do.call(plot.fun, plot.args)

          ## error plotting evaluation
          if (plot.errors && !(xi.factor && plot.bootstrap && plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = na.omit(if (plotOnEstimate) temp.mean else temp.err[,3]),
                all.err = temp.all.err,
                xi.factor = xi.factor)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(ei)),
                ely = if (plotOnEstimate) na.omit(temp.mean - temp.err[,1]) else na.omit(temp.err[,3] - temp.err[,1]),
                ehy = if (plotOnEstimate) na.omit(temp.mean + temp.err[,2]) else na.omit(temp.err[,3] + temp.err[,2]),
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = if (xi.factor) 1 else 2
              )
              do.call(draw.errors, draw.args)
            }

          }
        }

        
        if (plot.behavior != "plot") {
          plot.out[plot.index] = NA
          if (gradients){
          } else {
            plot.out[[plot.index]] =
              plregression(bws = bws, xcoef = tobj$xcoef,
                           xcoeferr = tobj$xcoeferr,
                           xcoefvcov = vcov(tobj),
                           evalx = exdat[1:xi.neval,, drop = FALSE],
                           evalz = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
                           mean = na.omit(temp.mean),
                           ntrain = dim(zdat)[1],
                           trainiseval = FALSE,
                           xtra = c(tobj$RSQ, tobj$MSE, 0, 0, 0, 0))
            plot.out[[plot.index]]$merr = NA
            plot.out[[plot.index]]$bias = NA

            if (plot.errors)
              plot.out[[plot.index]]$merr = temp.err[,1:2]

            if (plot.errors.center == "bias-corrected")
              plot.out[[plot.index]]$bias = temp.err[,3] - temp.mean
            plot.out[[plot.index]]$bxp = temp.boot
          }
        }
      }

      if (common.scale && (plot.behavior != "data")){
        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in seq_len(bws$xndim + bws$zndim)) {
            if (is.null(data.err.all[[k]])) next
            nkeep.k <- nrow(data.err.all[[k]]$pointwise)
            if (nkeep.k == 0) next
            center.k <- if (plot.errors.center == "estimate")
              na.omit(data.eval[seq_len(nkeep.k), k])
            else
              na.omit(data.err[seq_len(nkeep.k), 3*k])
            range.k <- compute.all.error.range(center.k, data.err.all[[k]])
            y.min <- min(y.min, range.k[1], na.rm = TRUE)
            y.max <- max(y.max, range.k[2], na.rm = TRUE)
          }
          if (!is.finite(y.min) || !is.finite(y.max)) {
            y.min <- min(na.omit(as.double(data.eval)))
            y.max <- max(na.omit(as.double(data.eval)))
          }
        } else {
          jj = seq_len(bws$xndim + bws$zndim)*3
          if (plot.errors.center == "estimate" || !plot.errors) {
            y.max = max(na.omit(as.double(data.eval)) +
              if (plot.errors) na.omit(as.double(data.err[,jj-1]))
              else 0)
            y.min = min(na.omit(as.double(data.eval)) -
              if (plot.errors) na.omit(as.double(data.err[,jj-2]))
              else 0)
          } else if (plot.errors.center == "bias-corrected") {
            y.max = max(na.omit(as.double(data.err[,jj] + data.err[,jj-1])))
            y.min = min(na.omit(as.double(data.err[,jj] - data.err[,jj-2])))
          }
        }

        if(!is.null(ylim)){
          y.min = ylim[1]
          y.max = ylim[2]
        }
        
        xOrZ = "x"
        
        for (plot.index in 1:(bws$xndim + bws$zndim)){
          i = if (plot.index <= bws$xndim) plot.index else plot.index - bws$xndim

          if (plot.index > bws$xndim)
            xOrZ <- "z"
            
          xi.factor = all.isFactor[plot.index]

          ## plot evaluation
          plot.fun <- if (xi.factor) {
            .np_plot_panel_fun(plot.bootstrap = plot.bootstrap, plot.bxp = plot.bxp)
          } else {
            plot
          }
          plot.args <- list()
          if (xi.factor) {
            if (plot.bootstrap && plot.bxp) plot.args$z <- all.bxp[[plot.index]] else plot.args$f <- allei[,plot.index]
          } else {
            plot.args$x <- allei[,plot.index]
          }
          if (!(xi.factor && plot.bootstrap && plot.bxp))
            plot.args$y <- data.eval[,plot.index]
          plot.args$ylim <- c(y.min, y.max)
          plot.args$xlab <- gen.label(if (xOrZ == "x") bws$xnames[i] else bws$znames[i],
                                      paste(toupper(xOrZ), i, sep = ""))
          plot.args$ylab <- paste(if (gradients) paste("Gradient Component ", i, " of", sep = "") else "",
                                  gen.label(bws$ynames, "Conditional Mean"))
          if (!xi.factor) {
            plot.args$type <- scalar_default(type, "l")
            plot.args$lty <- scalar_default(lty, par()$lty)
            plot.args$col <- scalar_default(col, par()$col)
            plot.args$lwd <- scalar_default(lwd, par()$lwd)
            plot.args$cex.axis <- scalar_default(cex.axis, par()$cex.axis)
            plot.args$cex.lab <- scalar_default(cex.lab, par()$cex.lab)
            plot.args$cex.main <- scalar_default(cex.main, par()$cex.main)
            plot.args$cex.sub <- scalar_default(cex.sub, par()$cex.sub)
          }
          plot.args$main <- scalar_default(main, "")
          plot.args$sub <- scalar_default(sub, "")
          do.call(plot.fun, plot.args)

          ## error plotting evaluation
          if (plot.errors && !(xi.factor && plot.bootstrap && plot.bxp)){
            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(allei[,plot.index])),
                center = if (plotOnEstimate)
                  na.omit(data.eval[,plot.index])
                else
                  na.omit(data.err[,3*plot.index]),
                all.err = data.err.all[[plot.index]],
                xi.factor = xi.factor)
            } else {
              if (!xi.factor && !plotOnEstimate)
                lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

              draw.args <- list(
                ex = as.numeric(na.omit(allei[,plot.index])),
                ely = if (plotOnEstimate) na.omit(data.eval[,plot.index] - data.err[,3*plot.index-2]) else na.omit(data.err[,3*plot.index] - data.err[,3*plot.index-2]),
                ehy = if (plotOnEstimate) na.omit(data.eval[,plot.index] + data.err[,3*plot.index-1]) else na.omit(data.err[,3*plot.index] + data.err[,3*plot.index-1]),
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = if (xi.factor) 1 else 2
              )
              do.call(draw.errors, draw.args)
            }
            

          }
        }
      }

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=c(1,1),cex=par()$cex)
      
      if (plot.behavior != "plot"){
        names(plot.out) =
          if (gradients){ }
          else
            paste("plr",1:(bws$xndim+bws$zndim),sep="")
        
        return (plot.out)
      }
    }
  }
