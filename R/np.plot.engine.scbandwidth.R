.np_plot_scbandwidth_engine <-
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
           coef = FALSE,
           coef.index = 1L,
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
           plot.errors.boot.num = 1999,
           plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
           plot.errors.boot.nonfixed = c("exact", "frozen"),
           plot.errors.boot.wild = c("rademacher", "mammen"),
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

    oldpar <- .np_plot_capture_par(c("mfrow", "cex"))
    on.exit(.np_plot_restore_par(oldpar), add = TRUE)

    scalar_default <- function(value, default) {
      if (is.null(value)) default else value
    }

    opt.plot.par.mfrow <- getOption("plot.par.mfrow")
    if (!is.null(opt.plot.par.mfrow))
        plot.par.mfrow <- opt.plot.par.mfrow

    if(!missing(gradients))
      stop("gradients not supported with smooth coefficient models.")
    coef <- isTRUE(coef)
    coef.index <- as.integer(coef.index)[1L]
    if (!is.finite(coef.index) || is.na(coef.index) || coef.index < 1L)
      coef.index <- 1L
    extract_scoef_value <- function(obj) {
      if (!coef)
        return(as.double(obj$mean))
      beta.plot <- obj$beta
      if (is.null(beta.plot))
        stop("coef=TRUE requires npscoef(..., betas=TRUE) output.")
      if (is.matrix(beta.plot))
        return(as.double(beta.plot[, min(coef.index, ncol(beta.plot))]))
      as.double(beta.plot)
    }

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
      mf.args <- as.list(tmf)[-1L]
      umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))
      
      ydat <- model.response(tmf)
      xdat <- tmf[, bws$chromoly[[2]], drop = FALSE]
      if (!miss.z)
        zdat <- tmf[, bws$chromoly[[3]], drop = FALSE]
    } else {
      if(all(miss.xy) && !is.null(bws$call)){
        xdat <- data.frame(.np_eval_bws_call_arg(bws, "xdat"))
        ydat <- .np_eval_bws_call_arg(bws, "ydat")
        if (!miss.z)
          zdat <- data.frame(.np_eval_bws_call_arg(bws, "zdat"))
      }
          
      ## catch and destroy NA's
      xdat = toFrame(xdat)
      if(!miss.z)
        zdat <- toFrame(zdat)

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

      if(!miss.z)
        zdat <- zdat[keep.rows,,drop = FALSE]
      
      ydat <- ydat[keep.rows]
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
            !missing(plot.errors.boot.nonfixed),
            !missing(plot.errors.boot.blocklen))){
      stop(
        "plot.errors.method must be set to 'bootstrap' when bootstrap error arguments are supplied",
        call. = FALSE
      )
    }
    normalized.opts <- .np_plot_normalize_common_options(
      plot.behavior = plot.behavior,
      plot.errors.method = plot.errors.method,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
      plot.errors.boot.wild = plot.errors.boot.wild,
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
    plot.errors.boot.nonfixed <- normalized.opts$plot.errors.boot.nonfixed
    plot.errors.boot.wild <- normalized.opts$plot.errors.boot.wild
    plot.errors.boot.blocklen <- normalized.opts$plot.errors.boot.blocklen
    plot.errors.center <- normalized.opts$plot.errors.center
    plot.errors.type <- normalized.opts$plot.errors.type
    plot.errors.alpha <- normalized.opts$plot.errors.alpha
    plot.errors.style <- normalized.opts$plot.errors.style
    plot.errors.bar <- normalized.opts$plot.errors.bar
    common.scale <- normalized.opts$common.scale
    plot.errors <- normalized.opts$plot.errors
    if (plot.errors.method == "bootstrap" &&
        identical(plot.errors.boot.nonfixed, "frozen") &&
        !identical(bws$type, "fixed")) {
      stop(
        "plot.errors.boot.nonfixed='frozen' is currently supported only for nonfixed unconditional/conditional density and distribution bootstrap routes",
        call. = FALSE
      )
    }
    if (coef && plot.errors.method != "none") {
      .np_warning("coef=TRUE currently disables plot errors for smooth coefficient plots.")
      plot.errors.method <- "none"
      plot.errors <- FALSE
    }

    if ((sum(c(bws$xdati$icon, bws$xdati$iord, bws$zdati$icon, bws$zdati$iord))== 2) && (sum(c(bws$xdati$iuno, bws$zdati$iuno)) == 0) && perspective && !gradients &&
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

      scoef.args <- list(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        iterate = FALSE,
        errors = (plot.errors && identical(plot.errors.method, "asymptotic")),
        betas = coef
      )
      if (!miss.z)
        scoef.args$tzdat <- zdat
      if (miss.z) {
        scoef.args$exdat <- x.eval
      } else {
        scoef.args$exdat <- x.eval[,1, drop = FALSE]
        scoef.args$ezdat <- x.eval[,2, drop = FALSE]
      }
      tobj <- do.call(.np_scoef_fit_internal, scoef.args)

      terr = matrix(data = tobj$merr, nrow = dim(x.eval)[1], ncol = 3)
      terr[,3] = NA
      
      treg = matrix(data = extract_scoef_value(tobj),
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      if (plot.errors.method == "bootstrap"){
        boot.args <- list(
          xdat = xdat,
          ydat = ydat,
          gradients = FALSE,
          slice.index = 0,
          progress.target = "surf 1/1",
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.wild = plot.errors.boot.wild,
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
          boot.args$ezdat <- x.eval[,2, drop = FALSE]
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
        terr[,1:2] <- .np_plot_asymptotic_error_from_se(
          se = tobj$merr,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type,
          m = nrow(x.eval)
        )$err
        lerr = matrix(data = tobj$mean - terr[,1],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = tobj$mean + terr[,2],
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

      persp.col = if (plot.errors) FALSE else scalar_default(col, "lightblue")
      
      ##for (j in 0:((50 %/% dphi - 1)*rotate)*dphi+phi){
        for (i in 0:((360 %/% dtheta - 1)*rotate)*dtheta+theta){
          persp(x1.eval,
                x2.eval,
                treg,
                zlim = zlim,
                col = persp.col,
                border = scalar_default(border, "black"),
                ticktype = "detailed",
                cex.axis = scalar_default(cex.axis, par()$cex.axis),
                cex.lab = scalar_default(cex.lab, par()$cex.lab),
                cex.main = scalar_default(cex.main, par()$cex.main),
                cex.sub = scalar_default(cex.sub, par()$cex.sub),
                xlab = scalar_default(xlab, gen.label(bws$xnames[1], "X1")),
                ylab = scalar_default(ylab, gen.label(x2.names[1], "X2")),
                zlab = scalar_default(zlab, gen.label(bws$ynames,"Conditional Mean")),
                theta = i,
                phi = phi,
                main = gen.tflabel(!is.null(main), main, paste("[theta= ", i,", phi= ", phi,"]", sep="")))

          if (plot.errors){
            par(new = TRUE)
            persp(x1.eval,
                  x2.eval,
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
            persp(x1.eval,
                  x2.eval,
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

      tot.dim <- (bws$xndim <- length(bws$xdati$icon)) + (bws$zndim <- length(bws$zdati$icon))

      plot.layout <- .np_plot_layout_begin(
        plot.behavior = plot.behavior,
        plot.par.mfrow = plot.par.mfrow,
        mfrow = n2mfrow(tot.dim)
      )

      maxneval = max(c(sapply(xdat,nlevels), unlist(sapply(zdat,nlevels)), neval))
      all.isFactor = c(vapply(xdat, is.factor, logical(1)), vapply(zdat, is.factor, logical(1)))
      
      x.ev = xdat[1,,drop = FALSE]

      for (i in seq_len(bws$xndim))
        x.ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      exdat = xdat[rep(1, maxneval), , drop = FALSE]

      for (i in seq_len(bws$xndim))
        exdat[,i] = x.ev[1,i]

      if (!miss.z){
        z.ev = zdat[1,,drop = FALSE]
        
        for (i in seq_len(bws$zndim))
          z.ev[1,i] = uocquantile(zdat[,i], prob=zq[i])

        ezdat = zdat[rep(1, maxneval), , drop = FALSE]

        for (i in seq_len(bws$zndim))
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

      plotOnEstimate = (plot.errors.center == "estimate")

      txobj_call <- function(ex.slice, ez.slice = NULL) {
        if (coef || identical(plot.errors.method, "asymptotic")) {
          tx.args <- list(
            txdat = xdat,
            tydat = ydat,
            exdat = ex.slice,
            bws = bws,
            errors = (plot.errors && identical(plot.errors.method, "asymptotic")),
            betas = coef
          )
          if (!miss.z) {
            tx.args$tzdat <- zdat
            tx.args$ezdat <- ez.slice
          }
          return(do.call(.np_scoef_fit_internal, tx.args))
        }

        hat.args <- list(
          bws = bws,
          txdat = xdat,
          exdat = ex.slice,
          y = ydat,
          output = "apply"
        )
        if (!miss.z) {
          hat.args$tzdat <- zdat
          hat.args$ezdat <- ez.slice
        }

        list(
          mean = as.vector(do.call(npscoefhat, hat.args)),
          merr = rep_len(NA_real_, nrow(ex.slice))
        )
      }

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

        ex.slice <- subcol(exdat, ei, i)[seq_len(xi.neval), , drop = FALSE]
        ez.slice <- if (!miss.z) ezdat[seq_len(xi.neval), , drop = FALSE] else NULL
        tobj <- txobj_call(ex.slice = ex.slice, ez.slice = ez.slice)

        temp.mean[seq_len(xi.neval)] = extract_scoef_value(tobj)

        if (plot.errors){
          if (plot.errors.method == "asymptotic") {
            asym.obj <- .np_plot_asymptotic_error_from_se(
              se = tobj$merr,
              alpha = plot.errors.alpha,
              band.type = plot.errors.type,
              m = xi.neval
            )
            temp.err[seq_len(xi.neval),1:2] <- asym.obj$err
            temp.all.err <- asym.obj$all.err
          } else if (plot.errors.method == "bootstrap"){
            boot.args <- list(
              xdat = xdat,
              ydat = ydat,
              exdat = ex.slice,
              gradients = gradients,
              slice.index = plot.index,
              progress.target = .np_plot_scoef_bootstrap_target_label(
                bws = bws,
                slice.index = plot.index
              ),
              plot.errors.boot.method = plot.errors.boot.method,
              plot.errors.boot.wild = plot.errors.boot.wild,
              plot.errors.boot.blocklen = plot.errors.boot.blocklen,
              plot.errors.boot.num = plot.errors.boot.num,
              plot.errors.center = plot.errors.center,
              plot.errors.type = plot.errors.type,
              plot.errors.alpha = plot.errors.alpha,
              bws = bws
            )
            if (!miss.z) {
              boot.args$zdat <- zdat
              boot.args$ezdat <- ez.slice
            }
            temp.boot.raw <- do.call(compute.bootstrap.errors, boot.args)
            temp.err[seq_len(xi.neval),] <- temp.boot.raw[["boot.err"]]
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
          plot.layout <- .np_plot_layout_activate(plot.layout)
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
          plot.args$ylab <- if (coef) paste("Coefficient", min(coef.index, bws$xndim)) else gen.label(bws$ynames, "Conditional Mean")
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
                center = as.numeric(na.omit(if (plotOnEstimate) temp.mean else temp.err[,3])),
                all.err = temp.all.err,
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(ei)),
                ely = if (plotOnEstimate) na.omit(temp.mean - temp.err[,1]) else na.omit(temp.err[,3] - temp.err[,1]),
                ehy = if (plotOnEstimate) na.omit(temp.mean + temp.err[,2]) else na.omit(temp.err[,3] + temp.err[,2]),
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
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
              subcol(exdat, ei, i)[seq_len(xi.neval),, drop = FALSE]
            } else {
              list(exdat = subcol(exdat, ei, i)[seq_len(xi.neval),, drop = FALSE],
                   ezdat = ezdat[seq_len(xi.neval),, drop = FALSE])
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

          ex.slice <- exdat[seq_len(xi.neval),, drop = FALSE]
          ez.slice <- subcol(ezdat,ei,i)[seq_len(xi.neval),, drop = FALSE]
          tobj <- txobj_call(ex.slice = ex.slice, ez.slice = ez.slice)

          temp.mean[seq_len(xi.neval)] = extract_scoef_value(tobj)

          if (plot.errors){
            if (plot.errors.method == "asymptotic") {
              asym.obj <- .np_plot_asymptotic_error_from_se(
                se = tobj$merr,
                alpha = plot.errors.alpha,
                band.type = plot.errors.type,
                m = xi.neval
              )
              temp.err[seq_len(xi.neval),1:2] <- asym.obj$err
              temp.all.err <- asym.obj$all.err
            } else if (plot.errors.method == "bootstrap"){
              temp.boot.raw <- compute.bootstrap.errors(
                                                    xdat = xdat,
                                                    ydat = ydat,
                                                    zdat = zdat,
                                                    exdat = ex.slice,
                                                    ezdat = ez.slice,
                                                    gradients = gradients,
                                                    slice.index = plot.index,
                                                    progress.target = .np_plot_scoef_bootstrap_target_label(
                                                      bws = bws,
                                                      slice.index = plot.index
                                                    ),
                                                    plot.errors.boot.method = plot.errors.boot.method,
                                                    plot.errors.boot.wild = plot.errors.boot.wild,
                                                    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                                                    plot.errors.boot.num = plot.errors.boot.num,
                                                    plot.errors.center = plot.errors.center,
                                                    plot.errors.type = plot.errors.type,
                                                    plot.errors.alpha = plot.errors.alpha,
                                                    bws = bws)
              temp.err[seq_len(xi.neval),] <- temp.boot.raw[["boot.err"]]
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
            plot.layout <- .np_plot_layout_activate(plot.layout)
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
            plot.args$ylab <- if (coef) paste("Coefficient", min(coef.index, bws$xndim)) else gen.label(bws$ynames, "Conditional Mean")
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
                  center = as.numeric(na.omit(if (plotOnEstimate) temp.mean else temp.err[,3])),
                  all.err = temp.all.err,
                  plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                  plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                  plot.errors.bar.num = plot.errors.bar.num,
                  lty = 2,
                  add.legend = TRUE)
              } else {
                draw.args <- list(
                  ex = as.numeric(na.omit(ei)),
                  ely = if (plotOnEstimate) na.omit(temp.mean - temp.err[,1]) else na.omit(temp.err[,3] - temp.err[,1]),
                  ehy = if (plotOnEstimate) na.omit(temp.mean + temp.err[,2]) else na.omit(temp.err[,3] + temp.err[,2]),
                  plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                  plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
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
                                  eval = list(exdat = exdat[seq_len(xi.neval),, drop = FALSE],
                                    ezdat = subcol(ezdat,ei,i)[seq_len(xi.neval),, drop = FALSE]),
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
      
      if (common.scale && (plot.behavior != "data")){
        jj = seq_len(bws$xndim + bws$zndim)*3
        
        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in seq_len(bws$xndim + bws$zndim)) {
            if (is.null(data.err.all[[k]])) next
            nkeep.k <- nrow(data.err.all[[k]]$pointwise)
            center.k <- if (plot.errors.center == "estimate")
              data.eval[seq_len(nkeep.k),k]
            else
              data.err[seq_len(nkeep.k),3*k]
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
        } else if (plot.errors.center == "estimate" || !plot.errors) {
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
        
        for (plot.index in seq_len(bws$xndim + bws$zndim)){
          i = if (plot.index <= bws$xndim) plot.index else plot.index - bws$xndim

          if (plot.index > bws$xndim)
            xOrZ <- "z"
            
          xi.factor = all.isFactor[plot.index]
          plot.layout <- .np_plot_layout_activate(plot.layout)

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
          plot.args$ylab <- if (coef) paste("Coefficient", min(coef.index, bws$xndim)) else gen.label(bws$ynames, "Conditional Mean")
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
                ex = as.numeric(na.omit(allei[,plot.index])),
                center = as.numeric(na.omit(if (plotOnEstimate) data.eval[,plot.index] else data.err[,3*plot.index])),
                all.err = data.err.all[[plot.index]],
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(allei[,plot.index])),
                ely = if (plotOnEstimate) na.omit(data.eval[,plot.index] - data.err[,3*plot.index-2]) else na.omit(data.err[,3*plot.index] - data.err[,3*plot.index-2]),
                ehy = if (plotOnEstimate) na.omit(data.eval[,plot.index] + data.err[,3*plot.index-1]) else na.omit(data.err[,3*plot.index] + data.err[,3*plot.index-1]),
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
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
            paste("sc", seq_len(bws$xndim + bws$zndim), sep = "")
        
        return (plot.out)
      }
    }

      
  }
