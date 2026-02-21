npplot.bandwidth <-
  function(bws,
           xdat,
           data = NULL,
           xq = 0.5,
           xtrim = 0.0,
           neval = 50,
           common.scale = TRUE,
           perspective = TRUE,
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

    if(!is.null(options('plot.par.mfrow')$plot.par.mfrow))
        plot.par.mfrow <- options('plot.par.mfrow')$plot.par.mfrow
      
    miss.x <- missing(xdat)

    if(miss.x && !is.null(bws$formula)){
      tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf, envir = environment(tt))

      xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]
    } else {
      if(miss.x && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["dat"]], environment(bws$call)))
      }
      xdat = toFrame(xdat)
      xdat = na.omit(xdat)
    }

    xq = double(bws$ndim)+xq
    xtrim = double(bws$ndim)+xtrim

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

    plot.behavior <- normalized.opts$plot.behavior
    plot.errors.method <- normalized.opts$plot.errors.method
    plot.errors.boot.method <- normalized.opts$plot.errors.boot.method
    plot.errors.boot.blocklen <- normalized.opts$plot.errors.boot.blocklen
    plot.errors.center <- normalized.opts$plot.errors.center
    plot.errors.type <- normalized.opts$plot.errors.type
    plot.errors.alpha <- normalized.opts$plot.errors.alpha
    plot.errors.style <- normalized.opts$plot.errors.style
    plot.errors.bar <- normalized.opts$plot.errors.bar
    common.scale <- normalized.opts$common.scale
    plot.errors <- normalized.opts$plot.errors

    if ((bws$ncon + bws$nord == 2) & (bws$nuno == 0) & perspective &
        !any(xor(bws$xdati$iord, bws$xdati$inumord))){
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

      if (is.ordered(xdat[,2])){
        x2.eval = bws$xdati$all.ulev[[2]]
        x2.neval = length(x2.eval)
      } else {
        x2.neval = neval
        qi = trim.quantiles(xdat[,2], xtrim[2])
        x2.eval = seq(qi[1], qi[2], length.out = x2.neval)
      }

      x.eval <- expand.grid(x1.eval, x2.eval)

      ##x.eval = as.data.frame(matrix(data = NA, ncol = 2, nrow = x1.neval*x2.neval))
      ##x.eval[,1] = rep(x1.eval, times = x2.neval)
      ##x.eval[,2] = rep(x2.eval, each = x1.neval)

      if (is.ordered(xdat[,1]))
        x1.eval <- (bws$xdati$all.dlev[[1]])[as.integer(x1.eval)]

      if (is.ordered(xdat[,2]))
        x2.eval <- (bws$xdati$all.dlev[[2]])[as.integer(x2.eval)]

      tobj =  npudens(tdat = xdat, edat = x.eval, bws = bws)

      tdens = matrix(data = tobj$dens,
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      terr = matrix(data = tobj$derr, nrow = nrow(x.eval), ncol = 3)
      terr[,3] = NA
      
      lerr.all <- NULL
      herr.all <- NULL

      if (plot.errors.method == "bootstrap"){
        terr.obj <- compute.bootstrap.errors(xdat = xdat, 
          exdat = x.eval,
          cdf = FALSE,
          slice.index = 0,
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
          plot.errors.boot.num = plot.errors.boot.num,
          plot.errors.center = plot.errors.center,
          plot.errors.type = plot.errors.type,
          plot.errors.alpha = plot.errors.alpha,
          bws = bws)
        terr <- terr.obj[["boot.err"]]
        terr.all <- terr.obj[["boot.all.err"]]

        pc = (plot.errors.center == "bias-corrected")
        center.val <- if (pc) terr[,3] else tobj$dens

        lerr = matrix(data = center.val - terr[,1],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = center.val + terr[,2],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        if (plot.errors.type == "all" && !is.null(terr.all)) {
          lerr.all <- lapply(terr.all, function(te)
            matrix(data = center.val - te[,1], nrow = x1.neval, ncol = x2.neval, byrow = FALSE))
          herr.all <- lapply(terr.all, function(te)
            matrix(data = center.val + te[,2], nrow = x1.neval, ncol = x2.neval, byrow = FALSE))
        }

      } else if (plot.errors.method == "asymptotic") {
        lerr = matrix(data = tobj$dens - qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = tobj$dens + qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      }

      if(is.null(zlim)) {
          zlim =
              if (plot.errors){
                  if (plot.errors.type == "all" && !is.null(lerr.all))
                  c(min(c(unlist(lerr.all), lerr)), max(c(unlist(herr.all), herr)))
                  else
                      c(min(lerr), max(herr))
              } else
                  c(min(tobj$dens),max(tobj$dens))
      }

      if (plot.behavior != "plot"){
        d1 <- npdensity(bws = bws, eval = x.eval, dens = tobj$dens, derr = terr[,1:2], ntrain = nrow(xdat))
        d1$bias = NA

        if (plot.errors.center == "bias-corrected")
          d1$bias = terr[,3] - tobj$dens
        
        if (plot.behavior == "data")
          return ( list(d1 = d1) )
      }
      
      # rows = constant x2
      # cols = constant x1

      dtheta = 5.0
      dphi = 10.0

      persp.col = ifelse(plot.errors, FALSE, ifelse(!is.null(col),col,"lightblue"))
      
      for (i in 0:((360 %/% dtheta - 1)*rotate)*dtheta+theta){
          if (plot.errors){
            if (plot.errors.type == "all" && !is.null(lerr.all)) {
              band.cols <- c(pointwise = "red", simultaneous = "green3", bonferroni = "blue")
              for (bn in c("pointwise", "simultaneous", "bonferroni")) {
                persp(x1.eval,
                      x2.eval,
                      lerr.all[[bn]],
                      zlim = zlim,
                      cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                      cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                      cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                      cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                      col = persp.col,
                      border = band.cols[bn],
                      ticktype = "detailed",
                      xlab = "",
                      ylab = "",
                      zlab = "",
                      theta = i,
                      phi = phi,
                      lwd = ifelse(!is.null(lwd),lwd,par()$lwd))
                par(new = TRUE)
              }
            } else {
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
          }
          
          persp(x1.eval,
                x2.eval,
                tdens,
                zlim = zlim,
                col = persp.col,
                border = ifelse(!is.null(border),border,"black"),
                ticktype = "detailed",
                cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                cex.lab =  ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                xlab = ifelse(!is.null(xlab),xlab,gen.label(names(xdat)[1], "X1")),
                ylab = ifelse(!is.null(ylab),ylab,gen.label(names(xdat)[2], "X2")),
                zlab = ifelse(!is.null(zlab),zlab,"Joint Density"),
                theta = i,
                phi = phi,
                main = gen.tflabel(!is.null(main), main, paste("[theta= ", i,", phi= ", phi,"]", sep="")))

          if (plot.errors){
            par(new = TRUE)
            if (plot.errors.type == "all" && !is.null(herr.all)) {
              band.cols <- c(pointwise = "red", simultaneous = "green3", bonferroni = "blue")
              for (bn in c("pointwise", "simultaneous", "bonferroni")) {
                persp(x1.eval,
                      x2.eval,
                      herr.all[[bn]],
                      zlim = zlim,
                      cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                      cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                      cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                      cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                      col = persp.col,
                      border = band.cols[bn],
                      ticktype = "detailed",
                      xlab = "",
                      ylab = "",
                      zlab = "",
                      theta = i,
                      phi = phi,
                      lwd = ifelse(!is.null(lwd),lwd,par()$lwd))
                if (bn != "bonferroni") par(new = TRUE)
              }
              legend("topleft",
                     legend = c("Pointwise","Simultaneous","Bonferroni"),
                     lty = 1, col = c("red","green3","blue"), lwd = 2, bty = "n")
            } else {
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
          }

          Sys.sleep(0.5)
      }
      

      if (plot.behavior == "plot-data")
        return ( list(d1 = d1) )        

    } else {

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(bws$ndim),cex=par()$cex)

      ev = xdat[1,,drop = FALSE]

      for (i in 1:bws$ndim)
        ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      maxneval = max(c(sapply(xdat,nlevels),neval))

      exdat = xdat[rep(1, maxneval), , drop = FALSE]

      for (i in 1:bws$ndim)
        exdat[,i] = ev[1,i]
      
      if (common.scale){
        data.eval = matrix(data = NA, nrow = maxneval, ncol = bws$ndim)
        data.err = matrix(data = NA, nrow = maxneval, ncol = 3*bws$ndim)
        data.err.all = vector("list", bws$ndim)
        allei = as.data.frame(matrix(data = NA, nrow = maxneval, ncol = bws$ndim))
        all.bxp = list()
      }

      plot.out = list()

      temp.err = matrix(data = NA, nrow = maxneval, ncol = 3)
      temp.dens = replicate(maxneval, NA)

      ## plotting expressions
      
      plot.bootstrap = plot.errors.method == "bootstrap"

      pfunE = expression(ifelse(xi.factor,
          ifelse(plot.bootstrap & plot.bxp,"bxp","plotFactor"), "plot"))

      pxE = expression(ifelse(common.scale,
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = all.bxp[[i]],", "f = allei[,i],"),
                 "x = allei[,i],"),
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = temp.boot,", "f = ei,"), "x = ei,")))

      pyE = expression(ifelse(xi.factor & plot.bootstrap & plot.bxp, "",
          ifelse(common.scale,"y = data.eval[,i],", "y = temp.dens,")))

      pylimE = ifelse(common.scale, "ylim = c(y.min,y.max),",
        ifelse(plot.errors,
               "ylim = if (plot.errors.type == 'all') compute.all.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.all.err) else compute.default.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.err),",
               ""))

      pxlabE = "xlab = ifelse(!is.null(xlab),xlab,gen.label(bws$xnames[i], paste('X', i, sep = ''))),"

      pylabE = paste("ylab = ", "ifelse(!is.null(ylab),ylab,'Density')", ",")

      prestE = expression(ifelse(xi.factor,"", "type = ifelse(!is.null(type),type,'l'), lty = ifelse(!is.null(lty),lty,par()$lty), col = ifelse(!is.null(col),col,par()$col), lwd = ifelse(!is.null(lwd),lwd,par()$lwd), cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis), cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab), cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main), cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),"))

      pmainE = "main = ifelse(!is.null(main),main,''), sub = ifelse(!is.null(sub),sub,''),"

      ## error plotting expressions
      plotOnEstimate = (plot.errors.center == "estimate")

      efunE = "draw.errors"
      eexE = ifelse(common.scale, "ex = as.numeric(na.omit(allei[,i])),",
          "ex = as.numeric(na.omit(ei)),")
      eelyE = ifelse(common.scale,
        ifelse(plotOnEstimate, "ely = na.omit(data.eval[,i] - data.err[,3*i-2]),",
               "ely = na.omit(data.err[,3*i] - data.err[,3*i-2]),"),
        ifelse(plotOnEstimate, "ely = na.omit(temp.dens - temp.err[,1]),",
               "ely = na.omit(temp.err[,3] - temp.err[,1]),"))
      eehyE = ifelse(common.scale,
        ifelse(plotOnEstimate, "ehy = na.omit(data.eval[,i] + data.err[,3*i-1]),",
               "ehy = na.omit(data.err[,3*i] + data.err[,3*i-1]),"),
        ifelse(plotOnEstimate, "ehy = na.omit(temp.dens + temp.err[,2]),",
               "ehy = na.omit(temp.err[,3] + temp.err[,2]),"))

      erestE = "plot.errors.style = ifelse(xi.factor,'bar',plot.errors.style),
                plot.errors.bar = ifelse(xi.factor,'I',plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(!is.null(lty),lty,ifelse(xi.factor,1,2))"

      ## density / distribution expressions

      for (i in 1:bws$ndim){
        temp.err[,] = NA
        temp.dens[] =  NA
        temp.boot = list()
        temp.all.err <- NULL

        xi.factor = is.factor(xdat[,i])

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

        eval.slice <- subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE]
        tobj <- npudens(tdat = xdat, edat = eval.slice, bws = bws)
        temp.dens[1:xi.neval] <- tobj$dens

        if (plot.errors){
          if (plot.errors.method == "asymptotic")
            temp.err[1:xi.neval,1:2] = replicate(2,qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr)
          else if (plot.errors.method == "bootstrap"){
            temp.boot.raw <- compute.bootstrap.errors(
                      xdat = xdat,
                      exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                      cdf = FALSE,
                      slice.index = i,
                      plot.errors.boot.method = plot.errors.boot.method,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
                      bws = bws)
            temp.err[1:xi.neval,] = temp.boot.raw[["boot.err"]]
            temp.all.err <- temp.boot.raw[["boot.all.err"]]
            temp.boot <- temp.boot.raw[["bxp"]]
            if (!plot.bxp.out){
              temp.boot$out <- numeric()
              temp.boot$group <- integer()
            }
          }
        }
        
        if (common.scale){
          allei[,i] = ei
          data.eval[,i] = temp.dens
          if (plot.errors){
            all.bxp[i] = NA
            all.bxp[[i]] = temp.boot

            data.err[,c(3*i-2,3*i-1,3*i)] = temp.err
            data.err.all[[i]] = temp.all.err
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
            plot.args$y <- temp.dens
          if (plot.errors)
            plot.args$ylim <- if (plot.errors.type == "all")
              compute.all.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.all.err)
            else
              compute.default.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.err)
          plot.args$xlab <- ifelse(!is.null(xlab), xlab, gen.label(bws$xnames[i], paste("X", i, sep = "")))
          plot.args$ylab <- ifelse(!is.null(ylab), ylab, "Density")
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
                center = as.numeric(na.omit(if (plotOnEstimate) temp.dens else temp.err[,3])),
                all.err = temp.all.err,
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(ei)),
                ely = if (plotOnEstimate) na.omit(temp.dens - temp.err[,1]) else na.omit(temp.err[,3] - temp.err[,1]),
                ehy = if (plotOnEstimate) na.omit(temp.dens + temp.err[,2]) else na.omit(temp.err[,3] + temp.err[,2]),
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(!is.null(lty), lty, ifelse(xi.factor, 1, 2))
              )
              do.call(draw.errors, draw.args)
            }
          }
        }

        if (plot.behavior != "plot") {
          plot.out[i] = NA
          plot.out[[i]] <- npdensity(
            bws = bws,
            eval = eval.slice,
            dens = na.omit(temp.dens),
            derr = na.omit(cbind(-temp.err[,1], temp.err[,2])),
            ntrain = bws$nobs
          )
          plot.out[[i]]$bias = na.omit(temp.dens - temp.err[,3])
          plot.out[[i]]$bxp = temp.boot
        }
      }
      
      if (common.scale & (plot.behavior != "data")){
        jj = 1:bws$ndim*3

        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in 1:bws$ndim) {
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
            else 0
            )
          y.min = min(na.omit(as.double(data.eval)) -
            if (plot.errors) na.omit(as.double(data.err[,jj-2]))
            else 0
            )
        } else if (plot.errors.center == "bias-corrected") {
          y.max = max(na.omit(as.double(data.err[,jj] + data.err[,jj-1])))
          y.min = min(na.omit(as.double(data.err[,jj] - data.err[,jj-2])))
        }

        if(!is.null(ylim)){
          y.min = ylim[1]
          y.max = ylim[2]
        }
        
        for (i in 1:bws$ndim){
          xi.factor = is.factor(xdat[,i])

          ## plot evaluation
          plot.fun <- if (xi.factor) {
            .np_plot_panel_fun(plot.bootstrap = plot.bootstrap, plot.bxp = plot.bxp)
          } else {
            plot
          }
          plot.args <- list()
          if (xi.factor) {
            if (plot.bootstrap && plot.bxp) plot.args$z <- all.bxp[[i]] else plot.args$f <- allei[,i]
          } else {
            plot.args$x <- allei[,i]
          }
          if (!(xi.factor && plot.bootstrap && plot.bxp))
            plot.args$y <- data.eval[,i]
          plot.args$ylim <- c(y.min, y.max)
          plot.args$xlab <- ifelse(!is.null(xlab), xlab, gen.label(bws$xnames[i], paste("X", i, sep = "")))
          plot.args$ylab <- ifelse(!is.null(ylab), ylab, "Density")
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
                ex = as.numeric(na.omit(allei[,i])),
                center = as.numeric(na.omit(if (plotOnEstimate) data.eval[,i] else data.err[,3*i])),
                all.err = data.err.all[[i]],
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(allei[,i])),
                ely = if (plotOnEstimate) na.omit(data.eval[,i] - data.err[,3*i-2]) else na.omit(data.err[,3*i] - data.err[,3*i-2]),
                ehy = if (plotOnEstimate) na.omit(data.eval[,i] + data.err[,3*i-1]) else na.omit(data.err[,3*i] + data.err[,3*i-1]),
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(!is.null(lty), lty, ifelse(xi.factor, 1, 2))
              )
              do.call(draw.errors, draw.args)
            }
          }
        }
      }

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=c(1,1),cex=par()$cex)

      if (plot.behavior != "plot"){
        names(plot.out) = paste("d",1:bws$ndim,sep="")
        return (plot.out)
      }
    }
  }
