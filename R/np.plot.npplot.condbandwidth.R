npplot.condbandwidth <-
  function(bws,
           xdat,
           ydat,
           data = NULL,
           xq = 0.5,
           yq = 0.5,
           xtrim = 0.0,
           ytrim = 0.0,
           neval = 50,
           quantreg = FALSE,
           gradients = FALSE,
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
           tau = 0.5,
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
      
    cdf <- TRUE
    miss.xy = c(missing(xdat),missing(ydat))
    
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

      ydat <- tmf[, bws$variableNames[["response"]], drop = FALSE]
      xdat <- tmf[, bws$variableNames[["terms"]], drop = FALSE]
    } else {
      if(all(miss.xy) && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["xdat"]], environment(bws$call)))
        ydat <- data.frame(eval(bws$call[["ydat"]], environment(bws$call)))
      }

      ## catch and destroy NA's
      xdat = toFrame(xdat)
      ydat = toFrame(ydat)
      
      goodrows = 1:dim(xdat)[1]
      rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
      goodrows[rows.omit] = 0

      if (all(goodrows==0))
        stop("Data has no rows without NAs")

      xdat = xdat[goodrows,,drop = FALSE]
      ydat = ydat[goodrows,,drop = FALSE]

    }

    if (quantreg & dim(ydat)[2] != 1)
      stop("'ydat' must have one column for quantile regression")
    
    xq = double(bws$xndim)+xq
    yq = double(bws$yndim)+yq
    
    xtrim = double(bws$xndim)+xtrim
    ytrim = double(bws$yndim)+ytrim

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

    if (plot.errors.method == "asymptotic" && quantreg && gradients){
      warning(paste("no asymptotic errors available for quantile regression gradients.",
                    "\nOne must instead use bootstrapping."))
      plot.errors.method = "none"
    }

    plot.errors = (plot.errors.method != "none")

    if ((bws$xncon + bws$xnord + bws$yncon + bws$ynord - quantreg == 2) &
        (bws$xnuno + bws$ynuno == 0) & perspective & !gradients &
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

      ## if we are doing quantile regression then we are dealing
      ## with 2 x variables ...

      if (quantreg){
        tx2 <- xdat[,2]
        txi <- 2
        txdati <- bws$xdati
        txtrim <- xtrim
      }
      else{
        tx2 <- ydat[,1]
        txi <- 1
        txdati <- bws$ydati
        txtrim <- ytrim
      }
      
      if (txdati$iord[txi]){
        x2.eval = txdati$all.ulev[[txi]]
        x2.neval = length(x2.eval)
      } else {
        x2.neval = neval
        qi = trim.quantiles(tx2, txtrim[txi])
        x2.eval = seq(qi[1], qi[2], length.out = x2.neval)
      }

      x.eval <- expand.grid(x1.eval, x2.eval)

      if (bws$xdati$iord[1])
        x1.eval <- (bws$xdati$all.dlev[[1]])[as.integer(x1.eval)]
      
      if (txdati$iord[txi])
        x2.eval <- (txdati$all.dlev[[txi]])[as.integer(x2.eval)]


      tboo =
        if(quantreg) "quant"
        else if (cdf) "dist"
        else "dens"

      method.fun <- switch(tboo,
                           "quant" = npqreg,
                           "dist" = npcdist,
                           "dens" = npcdens)
      margs <- list(txdat = xdat, tydat = ydat, bws = bws)
      if (quantreg) {
        margs$exdat <- x.eval
        margs$tau <- tau
      } else {
        margs$exdat <- x.eval[,1, drop = FALSE]
        margs$eydat <- x.eval[,2, drop = FALSE]
      }
      tobj <- do.call(method.fun, margs)

      tcomp = parse(text=paste("tobj$",
                      switch(tboo,
                             "quant" = "quantile",
                             "dist" = "condist",
                             "dens" = "condens"), sep=""))

      tcerr = parse(text=paste(ifelse(quantreg, "tobj$quanterr", "tobj$conderr")))

      tex = parse(text=paste(ifelse(quantreg, "x.eval", "x.eval[,1]")))
      tey = parse(text=paste(ifelse(quantreg, "NA", "x.eval[,2]")))

      tdens = matrix(data = eval(tcomp),
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      terr = matrix(data = eval(tcerr), nrow = length(eval(tcomp)), ncol = 3)
      terr[,3] = NA
      lerr.all <- NULL
      herr.all <- NULL
      
      if (plot.errors.method == "bootstrap"){
        terr.obj <- compute.bootstrap.errors(xdat = xdat, ydat = ydat,
          exdat = eval(tex), eydat = eval(tey),
          cdf = cdf,
          quantreg = quantreg,
          tau = tau,
          gradients = FALSE,
          gradient.index = 0,
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
        center.val <- if(pc) terr[,3] else eval(tcomp)

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
        lerr = matrix(data = eval(tcomp) - qnorm(plot.errors.alpha/2, lower.tail = FALSE)*eval(tcerr),
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = eval(tcomp) + qnorm(plot.errors.alpha/2, lower.tail = FALSE)*eval(tcerr),
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      }

      if(is.null(zlim)) {
        zlim =
          if (plot.errors){
            if (plot.errors.type == "all" && !is.null(lerr.all))
              c(min(c(unlist(lerr.all), lerr)), max(c(unlist(herr.all), herr)))
            else
              c(min(lerr),max(herr))
          } else
            c(min(eval(tcomp)),max(eval(tcomp)))
      }

      ## I am sorry it had to come to this ...
      tret = parse(text=paste(
                     switch(tboo,
                            "quant" = "qregression",
                            "dist" = "condistribution",
                            "dens" = "condensity"),
                     "(bws = bws, xeval = eval(tex),",
                     ifelse(quantreg, "tau = tau, quantile = eval(tcomp), quanterr = terr[,1:2]",
                            paste("yeval = eval(tey),", ifelse(cdf, "condist = ", "condens = "),
                                  "eval(tcomp), conderr = terr[,1:2]")),
                     ", ntrain = dim(xdat)[1])", sep=""))

      if (plot.behavior != "plot"){
        cd1 = eval(tret)
        cd1$bias = NA

        if (plot.errors.center == "bias-corrected")
          cd1$bias = terr[,3] - eval(tcomp)
        
        if (plot.behavior == "data")
          return ( list(cd1 = cd1) )
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
                xlab = ifelse(!is.null(xlab),xlab,gen.label(names(xdat)[1], "X")),
                ylab = ifelse(!is.null(ylab),ylab,gen.label(names(ydat)[1], "Y")),
                zlab = ifelse(!is.null(zlab),zlab,paste("Conditional", ifelse(cdf,"Distribution", "Density"))),
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
    } else {

      dsf = ifelse(gradients,bws$xndim,1)
      tot.dim = bws$xndim + bws$yndim - quantreg

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(dsf*tot.dim),cex=par()$cex)

      x.ev = xdat[1,,drop = FALSE]
      y.ev = ydat[1,,drop = FALSE]

      for (i in 1:bws$xndim)
        x.ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      for (i in 1:bws$yndim)
        y.ev[1,i] = uocquantile(ydat[,i], prob=yq[i])


      maxneval = max(c(sapply(xdat,nlevels), sapply(ydat,nlevels), neval))

      exdat = xdat[rep(1, maxneval), , drop = FALSE]
      eydat = as.data.frame(matrix(data = 0, nrow = maxneval, ncol = bws$yndim))

      for (i in 1:bws$xndim)
        exdat[,i] = x.ev[1,i]

      for (i in 1:bws$yndim)
        eydat[,i] = y.ev[1,i]

      if (common.scale){
        data.eval = matrix(data = NA, nrow = maxneval,
          ncol = tot.dim*dsf)
        
        data.err = matrix(data = NA, nrow = maxneval,
          ncol = 3*tot.dim*dsf)
        data.err.all = vector("list", tot.dim*dsf)

        allei = as.data.frame(matrix(data = NA, nrow = maxneval,
          ncol = tot.dim))

        all.bxp = list()
      }

      all.isFactor = c(sapply(xdat, is.factor), sapply(ydat, is.factor))

      plot.out = list()

      temp.err = matrix(data = NA, nrow = maxneval, ncol = 3)
      temp.dens = replicate(maxneval, NA)

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
          ifelse(common.scale,"y = data.eval[,plot.index],", "y = temp.dens,")))

      pylimE = ifelse(common.scale, "ylim = c(y.min,y.max),",
        ifelse(plot.errors, "ylim = if (plot.errors.type == 'all') compute.all.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.all.err) else c(min(na.omit(c(temp.dens - temp.err[,1], temp.err[,3] - temp.err[,1]))), max(na.omit(c(temp.dens + temp.err[,2], temp.err[,3] + temp.err[,2])))),", ""))

      pxlabE = expression(paste("xlab = ifelse(!is.null(xlab),xlab,gen.label(bws$",
          xOrY, "names[i], paste('", toupper(xOrY),"', i, sep = ''))),",sep=''))

      tylabE = ifelse(quantreg, paste(tau, 'quantile'),
        paste('Conditional', ifelse(cdf,'Distribution', 'Density')))

      pylabE = paste("ylab =", "paste(", ifelse(gradients,"'GC',j,'of',",''), "tylabE),")

      prestE = expression(ifelse(xi.factor,"", "type = ifelse(!is.null(type),type,'l'), lty = ifelse(!is.null(lty),lty,par()$lty), col = ifelse(!is.null(col),col,par()$col), lwd = ifelse(!is.null(lwd),lwd,par()$lwd), cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis), cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab), cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main), cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),"))

      pmainE = "main = ifelse(!is.null(main),main,''), sub = ifelse(!is.null(sub),sub,''),"

      ## error plotting expressions
      plotOnEstimate = (plot.errors.center == "estimate")

      efunE = "draw.errors"
      eexE = expression(ifelse(common.scale, "ex = as.numeric(na.omit(allei[,plot.index])),",
          "ex = as.numeric(na.omit(ei)),"))
      eelyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ely = na.omit(data.eval[,plot.index] - data.err[,3*plot.index-2]),",
                 "ely = na.omit(data.err[,3*plot.index] - data.err[,3*plot.index-2]),"),
          ifelse(plotOnEstimate, "ely = na.omit(temp.dens - temp.err[,1]),",
                 "ely = na.omit(temp.err[,3] - temp.err[,1]),")))
      eehyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ehy = na.omit(data.eval[,plot.index] + data.err[,3*plot.index-1]),",
                 "ehy = na.omit(data.err[,3*plot.index] + data.err[,3*plot.index-1]),"),
          ifelse(plotOnEstimate, "ehy = na.omit(temp.dens + temp.err[,2]),",
                 "ehy = na.omit(temp.err[,3] + temp.err[,2]),")))

      erestE = "plot.errors.style = ifelse(xi.factor,'bar',plot.errors.style),
                plot.errors.bar = ifelse(xi.factor,'I',plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(xi.factor,1,2)"


      plot.index = 0
      xOrY = "x"

      for (i in 1:bws$xndim){
        plot.index = plot.index + 1
        temp.err[,] = NA
        temp.dens[] =  NA
        temp.boot = list()

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

        method.fun <- if (cdf) npcdist else if (quantreg) npqreg else npcdens
        margs <- list(
          txdat = xdat,
          tydat = ydat,
          exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
          gradients = gradients,
          bws = bws
        )
        if (quantreg)
          margs$tau <- tau
        else
          margs$eydat <- eydat[1:xi.neval,, drop = FALSE]
        tobj <- do.call(method.fun, margs)

        
        ## if there are gradients then we need to repeat the process for each component

        tevalexpr = parse(text=paste("tobj$",ifelse(gradients,
                            ifelse(quantreg, "quantgrad[,j]","congrad[,j]"),
                            ifelse(cdf, "condist", ifelse(quantreg, "quantile",
                                                          "condens"))), sep=""))
        terrexpr = parse(text=paste("tobj$",ifelse(gradients,
                           "congerr[,j]", ifelse(quantreg,"quanterr",
                                                 "conderr")), sep=""))

        if (gradients & quantreg)
          terrexpr = parse(text="NA")

        if (plot.behavior != "plot"){
          plot.out[plot.index] = NA
          plot.out[[plot.index]] = tobj
        }

        for (j in 1:dsf){
          temp.boot = list()
          temp.all.err <- NULL
          temp.dens[1:xi.neval] = eval(tevalexpr) 
          
          if (plot.errors){
            if (plot.errors.method == "asymptotic")
              temp.err[1:xi.neval,1:2] = replicate(2,qnorm(plot.errors.alpha/2, lower.tail = FALSE)*eval(terrexpr))
            else if (plot.errors.method == "bootstrap"){
              temp.boot <- compute.bootstrap.errors(
                        xdat = xdat,
                        ydat = ydat,
                        exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                        eydat = eydat[1:xi.neval,, drop = FALSE],
                        cdf = cdf,
                        quantreg = quantreg,
                        tau = tau,
                        gradients = gradients,
                        gradient.index = j,
                        slice.index = plot.index,
                        plot.errors.boot.method = plot.errors.boot.method,
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
            data.eval[,(plot.index-1)*dsf+j] = temp.dens
            if (plot.errors){
              all.bxp[plot.index] = NA
              all.bxp[[plot.index]] = temp.boot

              data.err[,seq(3*((plot.index-1)*dsf+j)-2,length=3)] = temp.err
              data.err.all[[(plot.index-1)*dsf+j]] = temp.all.err
            }
          } else if (plot.behavior != "data") {
            ## plot evaluation
            eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                         eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                         eval(pmainE), ")")))

            ## error plotting evaluation
            if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
              if (plot.errors.type == "all") {
                draw.all.error.types(
                  ex = as.numeric(na.omit(ei)),
                  center = na.omit(if (plotOnEstimate) temp.dens else temp.err[,3]),
                  all.err = temp.all.err,
                  xi.factor = xi.factor)
              } else {
                if (!xi.factor && !plotOnEstimate)
                  lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

                eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                             eval(eehyE), eval(erestE), ")")))
              }

            }
          }
        
          if (plot.behavior != "plot" & plot.errors) {
            err.name <- ifelse(gradients,
                               paste("gc", j, "err", sep = ""),
                               ifelse(quantreg, "quanterr", "conderr"))
            bias.name <- ifelse(gradients, paste("gc", j, "bias", sep = ""), "bias")
            plot.out[[plot.index]][[err.name]] <- na.omit(cbind(-temp.err[,1], temp.err[,2]))
            plot.out[[plot.index]][[bias.name]] <- na.omit(temp.dens - temp.err[,3])
            plot.out[[plot.index]]$bxp <- temp.boot
          }
        }
      }

      if (!quantreg){
        xOrY = "y"
        for (i in 1:bws$yndim){
          plot.index = plot.index + 1
          temp.err[,] = NA
          temp.dens[] =  NA
          temp.boot = list()

          xi.factor = all.isFactor[plot.index]
          
          if (xi.factor){
            ei = bws$ydati$all.ulev[[i]]
            xi.neval = length(ei)
          } else {
            xi.neval = neval
            qi = trim.quantiles(ydat[,i], ytrim[i])
            ei = seq(qi[1], qi[2], length.out = neval)
          }

          if (xi.neval < maxneval){
            ei[(xi.neval+1):maxneval] = NA
          }

          method.fun <- if (cdf) npcdist else if (quantreg) npqreg else npcdens
          margs <- list(
            txdat = xdat,
            tydat = ydat,
            eydat = subcol(eydat,ei,i)[1:xi.neval,, drop = FALSE],
            gradients = gradients,
            bws = bws
          )
          if (quantreg)
            margs$tau <- tau
          else
            margs$exdat <- exdat[1:xi.neval,, drop = FALSE]
          tobj <- do.call(method.fun, margs)

          
          ## if there are gradients then we need to repeat the process for each component

          tevalexpr = parse(text=paste("tobj$",ifelse(gradients,
                              ifelse(quantreg, "quantgrad[,j]","congrad[,j]"),
                              ifelse(cdf, "condist", ifelse(quantreg, "quantile",
                                                            "condens"))), sep=""))
          terrexpr = parse(text=paste("tobj$",ifelse(gradients,
                             "congerr[,j]", ifelse(quantreg,"quanterr",
                                                   "conderr")), sep=""))

          if (gradients & quantreg)
            terrexpr = parse(text="NA")

          if (plot.behavior != "plot"){
            plot.out[plot.index] = NA
            plot.out[[plot.index]] = tobj
          }

          for (j in 1:dsf){
            temp.boot = list()
            temp.all.err <- NULL
            temp.dens[1:xi.neval] = eval(tevalexpr) 
            
            if (plot.errors){
              if (plot.errors.method == "asymptotic")
                temp.err[1:xi.neval,1:2] = replicate(2,qnorm(plot.errors.alpha/2, lower.tail = FALSE)*eval(terrexpr))
              else if (plot.errors.method == "bootstrap"){
                temp.boot <- compute.bootstrap.errors(
                          xdat = xdat,
                          ydat = ydat,
                          exdat = exdat[1:xi.neval,, drop = FALSE],
                          eydat = subcol(eydat,ei,i)[1:xi.neval,, drop = FALSE],
                          cdf = cdf,
                          quantreg = quantreg,
                          tau = tau,
                          gradients = gradients,
                          gradient.index = j,
                          slice.index = plot.index,
                          plot.errors.boot.method = plot.errors.boot.method,
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
              data.eval[,(plot.index-1)*dsf+j] = temp.dens
              if (plot.errors){
                all.bxp[plot.index] = NA
                all.bxp[[plot.index]] = temp.boot

                data.err[,seq(3*((plot.index-1)*dsf+j)-2,length=3)] = temp.err
                data.err.all[[(plot.index-1)*dsf+j]] = temp.all.err
              }
            } else if (plot.behavior != "data") {
              ## plot evaluation
              eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                           eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                           eval(pmainE), ")")))

              ## error plotting evaluation
              if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
                if (plot.errors.type == "all") {
                  draw.all.error.types(
                    ex = as.numeric(na.omit(ei)),
                    center = na.omit(if (plotOnEstimate) temp.dens else temp.err[,3]),
                    all.err = temp.all.err,
                    xi.factor = xi.factor)
                } else {
                  if (!xi.factor && !plotOnEstimate)
                    lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

                  eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                               eval(eehyE), eval(erestE), ")")))
                }

              }
            }
              
            if (plot.behavior != "plot" & plot.errors) {
              err.name <- ifelse(gradients,
                                 paste("gc", j, "err", sep = ""),
                                 ifelse(quantreg, "quanterr", "conderr"))
              bias.name <- ifelse(gradients, paste("gc", j, "bias", sep = ""), "bias")
              plot.out[[plot.index]][[err.name]] <- na.omit(cbind(-temp.err[,1], temp.err[,2]))
              plot.out[[plot.index]][[bias.name]] <- na.omit(temp.dens - temp.err[,3])
              plot.out[[plot.index]]$bxp <- temp.boot
            }
          }
        }
      }
      
      if (common.scale & (plot.behavior != "data")){
        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in 1:(tot.dim*dsf)) {
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
          jj = 1:(dsf*tot.dim)*3
          if (plot.errors.center == "estimate" | !plot.errors) {
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
        
        xOrY = "x"
        
        for (plot.index in 1:tot.dim){
          i = ifelse(plot.index <= bws$xndim, plot.index, plot.index - bws$xndim)

          if (plot.index > bws$xndim)
            xOrY <- "y"
            
          xi.factor = all.isFactor[plot.index]

          for (j in 1:dsf){
            ## plot evaluation
            eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                         eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                         eval(pmainE), ")")))

            ## error plotting evaluation
            if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
              if (plot.errors.type == "all") {
                idx <- (plot.index-1)*dsf+j
                draw.all.error.types(
                  ex = as.numeric(na.omit(allei[,plot.index])),
                  center = if (plotOnEstimate)
                    na.omit(data.eval[,idx])
                  else
                    na.omit(data.err[,3*idx]),
                  all.err = data.err.all[[idx]],
                  xi.factor = xi.factor)
              } else {
                if (!xi.factor && !plotOnEstimate)
                  lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

                eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                             eval(eehyE), eval(erestE), ")")))
              }
              
            }
          }
        }
      }

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=c(1,1),cex=par()$cex)

      if (plot.behavior != "plot"){
        names(plot.out) = paste("cd", 1:tot.dim, sep="")
        return (plot.out)
      }
    }
  }
