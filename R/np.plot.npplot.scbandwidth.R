npplot.scbandwidth <-
  function(bws,
           xdat,
           ydat,
           zdat = NULL,
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
           plot.errors.boot.num = 399,
           plot.errors.boot.method = c("inid", "fixed", "geom"),
           plot.errors.boot.blocklen = NULL,
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

    if(!is.null(options('plot.par.mfrow')$plot.par.mfrow))
        plot.par.mfrow <- options('plot.par.mfrow')$plot.par.mfrow      

    if(!missing(gradients))
      stop("gradients not supported with smooth coefficient models.")

    miss.xy = c(missing(xdat),missing(ydat))
    miss.z = missing(zdat) & is.null(bws$zdati)
    
    if (any(miss.xy) && !all(miss.xy))
      stop("one of, but not both, xdat and ydat was specified")
    else if(all(miss.xy) && !is.null(bws$formula)){
      tt <- terms(bws)
      m <- match(c("formula", "data", "subset", "na.action"),
                 names(bws$call), nomatch = 0)
      tmf <- bws$call[c(1,m)]
      tmf[[1]] <- as.name("model.frame")
      tmf[["formula"]] <- tt
      umf <- tmf <- eval(tmf, envir = environment(tt))
      
      ydat <- model.response(tmf)
      xdat <- tmf[, bws$chromoly[[2]], drop = FALSE]
      if (!miss.z)
        zdat <- tmf[, bws$chromoly[[3]], drop = FALSE]
    } else {
      if(all(miss.xy) && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["xdat"]], environment(bws$call)))
        ydat <- eval(bws$call[["ydat"]], environment(bws$call))
        if (!miss.z)
          zdat <- data.frame(eval(bws$call[["zdat"]], environment(bws$call)))
      }
          
      ## catch and destroy NA's
      xdat = toFrame(xdat)
      if(!miss.z)
        zdat <- toFrame(zdat)

      goodrows = 1:dim(xdat)[1]
      train.df <- data.frame(xdat, ydat)
      if (!miss.z)
        train.df <- data.frame(train.df, zdat)
      rows.omit <- attr(na.omit(train.df), "na.action")
      
      attr(na.omit(data.frame(xdat,ydat,zdat)), "na.action")
      goodrows[rows.omit] = 0

      if (all(goodrows==0))
        stop("Data has no rows without NAs")

      xdat = xdat[goodrows,,drop = FALSE]

      if(!miss.z)
        zdat <- zdat[goodrows,,drop = FALSE]
      
      ydat = ydat[goodrows]
    }

    ## ydat = as.double(ydat)

    xq = double(ncol(xdat))+xq
    xtrim = double(ncol(xdat))+xtrim

    if (!miss.z){
      zq = double(ncol(zdat))+zq
      ztrim = double(ncol(zdat))+ztrim
    }
    

    if (missing(plot.errors.method) &
        any(!missing(plot.errors.boot.num), !missing(plot.errors.boot.method),
            !missing(plot.errors.boot.blocklen))){
      warning(paste("plot.errors.method must be set to 'bootstrap' to use bootstrapping.",
                    "\nProceeding without bootstrapping."))
    }

    
    normalized.opts <- .npplot_normalize_common_options(
      plot.behavior = plot.behavior,
      plot.errors.method = plot.errors.method,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      plot.errors.style = plot.errors.style,
      plot.errors.bar = plot.errors.bar,
      xdat = xdat,
      common.scale = common.scale,
      ylim = ylim
    )

    plot.behavior <- .npRmpi_plot_behavior_for_rank(normalized.opts$plot.behavior)
    plot.errors.method <- normalized.opts$plot.errors.method
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method,
                                              where = "plot(...)",
                                              allow.direct.bootstrap = TRUE)
    plot.errors.boot.method <- normalized.opts$plot.errors.boot.method
    plot.errors.boot.blocklen <- normalized.opts$plot.errors.boot.blocklen
    plot.errors.center <- normalized.opts$plot.errors.center
    plot.errors.type <- normalized.opts$plot.errors.type
    plot.errors.alpha <- normalized.opts$plot.errors.alpha
    plot.errors.style <- normalized.opts$plot.errors.style
    plot.errors.bar <- normalized.opts$plot.errors.bar
    common.scale <- normalized.opts$common.scale
    plot.errors <- normalized.opts$plot.errors

    if ((sum(c(bws$xdati$icon, bws$xdati$iord, bws$zdati$icon, bws$zdati$iord))== 2) & (sum(c(bws$xdati$iuno, bws$zdati$iuno)) == 0) & perspective & !gradients &
        !any(xor(c(bws$xdati$iord, bws$zdati$iord), c(bws$xdati$inumord, bws$zdati$inumord)))){

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

      if(!miss.z){
        tdat <- zdat[,1]
        ti <- 1
        tdati <-  bws$zdati
        ttrim <- ztrim
        x2.names <- bws$znames
      } else {
        tdat <- xdat[,2]
        ti <- 2
        tdati <- bws$xdati
        ttrim <- xtrim
        x2.names <- bws$xnames
      }
        
      if (is.ordered(tdat)){
        x2.eval = tdati$all.ulev[[ti]]
        x2.neval = length(x2.eval)
      } else {
        x2.neval = neval
        qi = trim.quantiles(tdat, ttrim[ti])
        x2.eval = seq(qi[1], qi[2], length.out = x2.neval)
      }

      x.eval <- expand.grid(x1.eval, x2.eval)
      
      if (is.ordered(xdat[,1]))
        x1.eval <- (bws$xdati$all.dlev[[1]])[as.integer(x1.eval)]
      
      if (is.ordered(tdat))
        x2.eval <- (tdati$all.dlev[[ti]])[as.integer(x2.eval)]

      scoef.args <- list(txdat = xdat, tydat = ydat, bws = bws, iterate = FALSE, errors = plot.errors)
      if (!miss.z)
        scoef.args$tzdat <- zdat
      if (miss.z) {
        scoef.args$exdat <- x.eval
      } else {
        scoef.args$exdat <- x.eval[,1, drop = FALSE]
        scoef.args$ezdat <- x.eval[,2, drop = FALSE]
      }
      tobj <- do.call(npscoef, scoef.args)

      terr = matrix(data = tobj$merr, nrow = dim(x.eval)[1], ncol = 3)
      terr[,3] = NA
      
      treg = matrix(data = tobj$mean,
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      if (plot.errors.method == "bootstrap"){
        boot.args <- list(
          xdat = xdat,
          ydat = ydat,
          gradients = FALSE,
          slice.index = 0,
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
          plot.errors.boot.num = plot.errors.boot.num,
          plot.errors.center = plot.errors.center,
          plot.errors.type = plot.errors.type,
          plot.errors.alpha = plot.errors.alpha,
          bws = bws
        )
        if (!miss.z)
          boot.args$zdat <- zdat
        if (miss.z) {
          boot.args$exdat <- x.eval
        } else {
          boot.args$exdat <- x.eval[,1, drop = FALSE]
          boot.args$ezdat <- x.eval[,1, drop = FALSE]
        }
        terr <- do.call(compute.bootstrap.errors, boot.args)[["boot.err"]]

        pc = (plot.errors.center == "bias-corrected")

        lerr = matrix(data = if(pc) {terr[,3]} else {tobj$mean}
          -terr[,1],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = if(pc) {terr[,3]} else {tobj$mean}
          +terr[,2],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      } else if (plot.errors.method == "asymptotic") {
        lerr = matrix(data = tobj$mean - qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = tobj$mean + qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      }

      if(is.null(zlim)) {
          zlim =
              if (plot.errors)
                  c(min(lerr),max(herr))
              else
                  c(min(tobj$mean),max(tobj$mean))
      }
        
      if (plot.behavior != "plot"){
        r1 <- do.call(smoothcoefficient, list(
          bws = bws,
          eval = if (miss.z) x.eval else list(exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,2, drop = FALSE]),
          mean = as.double(treg),
          merr = terr[,1:2],
          ntrain = dim(xdat)[1]
        ))
        r1$bias = NA

        if (plot.errors.center == "bias-corrected")
          r1$bias = terr[,3] - treg

        if (plot.behavior == "data")
          return ( list(r1 = r1) )

      }

      dtheta = 5.0
      dphi = 10.0

      persp.col = ifelse(plot.errors, FALSE, ifelse(!is.null(col),col,"lightblue"))
      
      ##for (j in 0:((50 %/% dphi - 1)*rotate)*dphi+phi){
        for (i in 0:((360 %/% dtheta - 1)*rotate)*dtheta+theta){
          if (plot.errors){
            persp(x1.eval,
                  x2.eval,
                  lerr,
                  zlim = zlim,
                  cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                  cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                  cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                  cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                  col = persp.col,
                  border = ifelse(!is.null(border),border,"grey"),
                  ticktype = "detailed",
                  xlab = "",
                  ylab = "",
                  zlab = "",
                  theta = i,
                  phi = phi,
                  lwd = ifelse(!is.null(lwd),lwd,par()$lwd))
            par(new = TRUE)
          }

          persp(x1.eval,
                x2.eval,
                treg,
                zlim = zlim,
                col = persp.col,
                border = ifelse(!is.null(border),border,"black"),
                ticktype = "detailed",
                cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                cex.lab =  ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                xlab = ifelse(!is.null(xlab),xlab,gen.label(bws$xnames[1], "X1")),
                ylab = ifelse(!is.null(ylab),ylab,gen.label(x2.names[1], "X2")),
                zlab = ifelse(!is.null(zlab),zlab,gen.label(bws$ynames,"Conditional Mean")),
                theta = i,
                phi = phi,
                main = gen.tflabel(!is.null(main), main, paste("[theta= ", i,", phi= ", phi,"]", sep="")))

          if (plot.errors){
            par(new = TRUE)
            persp(x1.eval,
                  x2.eval,
                  herr,
                  zlim = zlim,
                  cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                  cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                  cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                  cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                  col = persp.col,
                  border = ifelse(!is.null(border),border,"grey"),
                  ticktype = "detailed",
                  xlab = "",
                  ylab = "",
                  zlab = "",
                  theta = i,
                  phi = phi,
                  lwd = ifelse(!is.null(lwd),lwd,par()$lwd))
          }

          Sys.sleep(0.5)
        }
      ##}

      if (plot.behavior == "plot-data")
        return ( list(r1 = r1) )

    } else {

      tot.dim <- (bws$xndim <- length(bws$xdati$icon)) + (bws$zndim <- length(bws$zdati$icon))

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(tot.dim),cex=par()$cex)

      maxneval = max(c(sapply(xdat,nlevels), unlist(sapply(zdat,nlevels)), neval))
      all.isFactor = c(sapply(xdat, is.factor), unlist(sapply(zdat, is.factor)))
      
      x.ev = xdat[1,,drop = FALSE]

      for (i in 1:bws$xndim)
        x.ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      exdat = xdat[rep(1, maxneval), , drop = FALSE]

      for (i in 1:bws$xndim)
        exdat[,i] = x.ev[1,i]

      if (!miss.z){
        z.ev = zdat[1,,drop = FALSE]
        
        for (i in 1:bws$zndim)
          z.ev[1,i] = uocquantile(zdat[,i], prob=zq[i])

        ezdat = zdat[rep(1, maxneval), , drop = FALSE]

        for (i in 1:bws$zndim)
          ezdat[,i] = z.ev[1,i]

      }

      if (common.scale){
        data.eval = matrix(data = NA, nrow = maxneval,
          ncol = tot.dim)
        
        data.err = matrix(data = NA, nrow = maxneval,
          ncol = 3*tot.dim)
        data.err.all = vector("list", tot.dim)

        allei = as.data.frame(matrix(data = NA, nrow = maxneval,
          ncol = tot.dim))

        all.bxp = list()
      }

      plot.out = list()

      temp.err = matrix(data = NA, nrow = maxneval, ncol = 3)
      temp.mean = replicate(maxneval, NA)

      ## plotting expressions
      
      plot.bootstrap = plot.errors.method == "bootstrap"
      
      pfunE = expression(ifelse(xi.factor,
          ifelse(plot.bootstrap & plot.bxp,"bxp","plotFactor"), "plot"))
      pxE = expression(ifelse(common.scale,
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "all.bxp[[plot.index]],",
                        "f = allei[,plot.index],"),
                 "x = allei[,plot.index],"),
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = temp.boot,", "f = ei,"),
                 "x = ei,")))

      pyE = expression(ifelse(xi.factor & plot.bootstrap & plot.bxp, "",
          ifelse(common.scale,"y = data.eval[,plot.index],", "y = temp.mean,")))

      pylimE = ifelse(common.scale, "ylim = c(y.min,y.max),",
        ifelse(plot.errors,
               "ylim = if (plot.errors.type == 'all') compute.all.error.range(if (plotOnEstimate) temp.mean else temp.err[,3], temp.all.err) else c(min(na.omit(c(temp.mean - temp.err[,1], temp.err[,3] - temp.err[,1]))), max(na.omit(c(temp.mean + temp.err[,2], temp.err[,3] + temp.err[,2])))),",
               ""))

      pxlabE = expression(paste("xlab = gen.label(bws$",
          xOrZ, "names[i], paste('", toupper(xOrZ),"', i, sep = '')),",sep=''))

      pylabE = "ylab = paste(ifelse(gradients,
          paste('Gradient Component ', i, ' of', sep=''), ''),
          gen.label(bws$ynames, 'Conditional Mean')),"

      prestE = expression(ifelse(xi.factor,"", "type = ifelse(!is.null(type),type,'l'), lty = ifelse(!is.null(lty),lty,par()$lty), col = ifelse(!is.null(col),col,par()$col), lwd = ifelse(!is.null(lwd),lwd,par()$lwd), cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis), cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab), cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main), cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),"))

      pmainE = "main = ifelse(!is.null(main),main,''), sub = ifelse(!is.null(sub),sub,''),"

      txobj_call <- function(i, ei, xi.neval) {
        tx.args <- list(
          txdat = xdat,
          tydat = ydat,
          exdat = subcol(exdat, ei, i)[1:xi.neval, , drop = FALSE],
          bws = bws,
          errors = plot.errors
        )
        if (!miss.z) {
          tx.args$tzdat <- zdat
          tx.args$ezdat <- ezdat[1:xi.neval, , drop = FALSE]
        }
        do.call(npscoef, tx.args)
      }


      ## error plotting expressions
      plotOnEstimate = (plot.errors.center == "estimate")

      efunE = "draw.errors"
      eexE = expression(ifelse(common.scale, "ex = as.numeric(na.omit(allei[,plot.index])),",
          "ex = as.numeric(na.omit(ei)),"))
      eelyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ely = na.omit(data.eval[,plot.index] - data.err[,3*plot.index-2]),",
                 "ely = na.omit(data.err[,3*plot.index] - data.err[,3*plot.index-2]),"),
          ifelse(plotOnEstimate, "ely = na.omit(temp.mean - temp.err[,1]),",
                 "ely = na.omit(temp.err[,3] - temp.err[,1]),")))
      eehyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ehy = na.omit(data.eval[,plot.index] + data.err[,3*plot.index-1]),",
                 "ehy = na.omit(data.err[,3*plot.index] + data.err[,3*plot.index-1]),"),
          ifelse(plotOnEstimate, "ehy = na.omit(temp.mean + temp.err[,2]),",
                 "ehy = na.omit(temp.err[,3] + temp.err[,2]),")))

      erestE = "plot.errors.style = ifelse(xi.factor,'bar',plot.errors.style),
                plot.errors.bar = ifelse(xi.factor,'I',plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(xi.factor,1,2)"

      plot.index = 0
      xOrZ = "x"

      for (i in 1:bws$xndim){
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

        tobj <- txobj_call(i, ei, xi.neval)

        temp.mean[1:xi.neval] = tobj$mean

        if (plot.errors){
          if (plot.errors.method == "asymptotic")
            temp.err[1:xi.neval,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr
          else if (plot.errors.method == "bootstrap"){
            boot.args <- list(
              xdat = xdat,
              ydat = ydat,
              exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
              gradients = gradients,
              slice.index = plot.index,
              plot.errors.boot.method = plot.errors.boot.method,
              plot.errors.boot.blocklen = plot.errors.boot.blocklen,
              plot.errors.boot.num = plot.errors.boot.num,
              plot.errors.center = plot.errors.center,
              plot.errors.type = plot.errors.type,
              plot.errors.alpha = plot.errors.alpha,
              bws = bws
            )
            if (!miss.z) {
              boot.args$zdat <- zdat
              boot.args$ezdat <- ezdat[1:xi.neval,, drop = FALSE]
            }
            temp.boot.raw <- do.call(compute.bootstrap.errors, boot.args)
            temp.err[1:xi.neval,] <- temp.boot.raw[["boot.err"]]
            temp.all.err <- temp.boot.raw[["boot.all.err"]]
            temp.boot <- temp.boot.raw[["bxp"]]
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
            if (plot.bootstrap && plot.bxp) bxp else plotFactor
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
          plot.args$ylab <- gen.label(bws$ynames, "Conditional Mean")
          if (!xi.factor) {
            plot.args$type <- ifelse(!is.null(type), type, "l")
            plot.args$lty <- ifelse(!is.null(lty), lty, par()$lty)
            plot.args$col <- ifelse(!is.null(col), col, par()$col)
            plot.args$lwd <- ifelse(!is.null(lwd), lwd, par()$lwd)
            plot.args$cex.axis <- ifelse(!is.null(cex.axis), cex.axis, par()$cex.axis)
            plot.args$cex.lab <- ifelse(!is.null(cex.lab), cex.lab, par()$cex.lab)
            plot.args$cex.main <- ifelse(!is.null(cex.main), cex.main, par()$cex.main)
            plot.args$cex.sub <- ifelse(!is.null(cex.sub), cex.sub, par()$cex.sub)
          }
          plot.args$main <- ifelse(!is.null(main), main, "")
          plot.args$sub <- ifelse(!is.null(sub), sub, "")
          do.call(plot.fun, plot.args)

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = as.numeric(na.omit(if (plotOnEstimate) temp.mean else temp.err[,3])),
                all.err = temp.all.err,
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(ei)),
                ely = if (plotOnEstimate) na.omit(temp.mean - temp.err[,1]) else na.omit(temp.err[,3] - temp.err[,1]),
                ehy = if (plotOnEstimate) na.omit(temp.mean + temp.err[,2]) else na.omit(temp.err[,3] + temp.err[,2]),
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2
              )
              do.call(draw.errors, draw.args)
            }

          }
        }

        if (plot.behavior != "plot") {
          plot.out[plot.index] = NA
          if (gradients){
          } else {
            eval.obj <- if (miss.z) {
              subcol(exdat, ei, i)[1:xi.neval,, drop = FALSE]
            } else {
              list(exdat = subcol(exdat, ei, i)[1:xi.neval,, drop = FALSE],
                   ezdat = ezdat[1:xi.neval,, drop = FALSE])
            }
            plot.out[[plot.index]] <-
              smoothcoefficient(bws = bws,
                                eval = eval.obj,
                                mean = na.omit(temp.mean),
                                ntrain = dim(xdat)[1],
                                trainiseval = FALSE,
                                xtra = c(0, 0, 0, 0, 0, 0))
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

      if (!miss.z){
        xOrZ = "z"
        for (i in 1:bws$zndim){
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

          tobj <- npscoef(txdat = xdat, tydat = ydat, tzdat = zdat,
                          exdat = exdat[1:xi.neval,, drop = FALSE],
                          ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
                          bws = bws)

          temp.mean[1:xi.neval] = tobj$mean

          if (plot.errors){
            if (plot.errors.method == "asymptotic")
              temp.err[1:xi.neval,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr
          else if (plot.errors.method == "bootstrap"){
              temp.boot.raw <- compute.bootstrap.errors(
                                                    xdat = xdat,
                                                    ydat = ydat,
                                                    zdat = zdat,
                                                    exdat = exdat[1:xi.neval,, drop = FALSE],
                                                    ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
                                                    gradients = gradients,
                                                    slice.index = plot.index,
                                                    plot.errors.boot.method = plot.errors.boot.method,
                                                    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                                                    plot.errors.boot.num = plot.errors.boot.num,
                                                    plot.errors.center = plot.errors.center,
                                                    plot.errors.type = plot.errors.type,
                                                    plot.errors.alpha = plot.errors.alpha,
                                                    bws = bws)
              temp.err[1:xi.neval,] <- temp.boot.raw[["boot.err"]]
              temp.all.err <- temp.boot.raw[["boot.all.err"]]
              temp.boot <- temp.boot.raw[["bxp"]]
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
              if (plot.bootstrap && plot.bxp) bxp else plotFactor
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
            plot.args$ylab <- gen.label(bws$ynames, "Conditional Mean")
            if (!xi.factor) {
              plot.args$type <- ifelse(!is.null(type), type, "l")
              plot.args$lty <- ifelse(!is.null(lty), lty, par()$lty)
              plot.args$col <- ifelse(!is.null(col), col, par()$col)
              plot.args$lwd <- ifelse(!is.null(lwd), lwd, par()$lwd)
              plot.args$cex.axis <- ifelse(!is.null(cex.axis), cex.axis, par()$cex.axis)
              plot.args$cex.lab <- ifelse(!is.null(cex.lab), cex.lab, par()$cex.lab)
              plot.args$cex.main <- ifelse(!is.null(cex.main), cex.main, par()$cex.main)
              plot.args$cex.sub <- ifelse(!is.null(cex.sub), cex.sub, par()$cex.sub)
            }
            plot.args$main <- ifelse(!is.null(main), main, "")
            plot.args$sub <- ifelse(!is.null(sub), sub, "")
            do.call(plot.fun, plot.args)

            ## error plotting evaluation
            if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
              if (!xi.factor && !plotOnEstimate)
                lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

              if (plot.errors.type == "all") {
                draw.all.error.types(
                  ex = as.numeric(na.omit(ei)),
                  center = as.numeric(na.omit(if (plotOnEstimate) temp.mean else temp.err[,3])),
                  all.err = temp.all.err,
                  plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                  plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                  plot.errors.bar.num = plot.errors.bar.num,
                  lty = 2,
                  add.legend = TRUE)
              } else {
                draw.args <- list(
                  ex = as.numeric(na.omit(ei)),
                  ely = if (plotOnEstimate) na.omit(temp.mean - temp.err[,1]) else na.omit(temp.err[,3] - temp.err[,1]),
                  ehy = if (plotOnEstimate) na.omit(temp.mean + temp.err[,2]) else na.omit(temp.err[,3] + temp.err[,2]),
                  plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                  plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                  plot.errors.bar.num = plot.errors.bar.num,
                  lty = 2
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
                smoothcoefficient(bws = bws, 
                                  eval = list(exdat = exdat[1:xi.neval,, drop = FALSE],
                                    ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE]),
                                  mean = na.omit(temp.mean),
                                  ntrain = dim(zdat)[1],
                                  trainiseval = FALSE)

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
      }
      
      if (common.scale & (plot.behavior != "data")){
        jj = 1:(bws$xndim + bws$zndim)*3
        
        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in 1:(bws$xndim + bws$zndim)) {
            if (is.null(data.err.all[[k]])) next
            nkeep.k <- nrow(data.err.all[[k]]$pointwise)
            center.k <- if (plot.errors.center == "estimate")
              data.eval[1:nkeep.k,k]
            else
              data.err[1:nkeep.k,3*k]
            range.k <- compute.all.error.range(center.k, data.err.all[[k]])
            y.min <- min(y.min, range.k[1], na.rm = TRUE)
            y.max <- max(y.max, range.k[2], na.rm = TRUE)
          }
          if (!is.finite(y.min) || !is.finite(y.max)) {
            if (plot.errors.center == "estimate") {
              y.max = max(na.omit(as.double(data.eval)) + na.omit(as.double(data.err[,jj-1])))
              y.min = min(na.omit(as.double(data.eval)) - na.omit(as.double(data.err[,jj-2])))
            } else {
              y.max = max(na.omit(as.double(data.err[,jj] + data.err[,jj-1])))
              y.min = min(na.omit(as.double(data.err[,jj] - data.err[,jj-2])))
            }
          }
        } else if (plot.errors.center == "estimate" | !plot.errors) {
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

        if(!is.null(ylim)){
          y.min = ylim[1]
          y.max = ylim[2]
        }
        
        xOrZ = "x"
        
        for (plot.index in 1:(bws$xndim + bws$zndim)){
          i = ifelse(plot.index <= bws$xndim, plot.index, plot.index - bws$xndim)

          if (plot.index > bws$xndim)
            xOrZ <- "z"
            
          xi.factor = all.isFactor[plot.index]

          ## plot evaluation
          plot.fun <- if (xi.factor) {
            if (plot.bootstrap && plot.bxp) bxp else plotFactor
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
          plot.args$ylab <- gen.label(bws$ynames, "Conditional Mean")
          if (!xi.factor) {
            plot.args$type <- ifelse(!is.null(type), type, "l")
            plot.args$lty <- ifelse(!is.null(lty), lty, par()$lty)
            plot.args$col <- ifelse(!is.null(col), col, par()$col)
            plot.args$lwd <- ifelse(!is.null(lwd), lwd, par()$lwd)
            plot.args$cex.axis <- ifelse(!is.null(cex.axis), cex.axis, par()$cex.axis)
            plot.args$cex.lab <- ifelse(!is.null(cex.lab), cex.lab, par()$cex.lab)
            plot.args$cex.main <- ifelse(!is.null(cex.main), cex.main, par()$cex.main)
            plot.args$cex.sub <- ifelse(!is.null(cex.sub), cex.sub, par()$cex.sub)
          }
          plot.args$main <- ifelse(!is.null(main), main, "")
          plot.args$sub <- ifelse(!is.null(sub), sub, "")
          do.call(plot.fun, plot.args)

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(allei[,plot.index])),
                center = as.numeric(na.omit(if (plotOnEstimate) data.eval[,plot.index] else data.err[,3*plot.index])),
                all.err = data.err.all[[plot.index]],
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(allei[,plot.index])),
                ely = if (plotOnEstimate) na.omit(data.eval[,plot.index] - data.err[,3*plot.index-2]) else na.omit(data.err[,3*plot.index] - data.err[,3*plot.index-2]),
                ehy = if (plotOnEstimate) na.omit(data.eval[,plot.index] + data.err[,3*plot.index-1]) else na.omit(data.err[,3*plot.index] + data.err[,3*plot.index-1]),
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2
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
            paste("sc",1:(bws$xndim+bws$zndim),sep="")
        
        return (plot.out)
      }
    }

      
  }
