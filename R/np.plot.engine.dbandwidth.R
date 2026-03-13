.np_plot_dbandwidth_engine <-
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
           plot.errors.boot.nonfixed = c("exact", "frozen"),
           plot.errors.boot.blocklen = NULL,
           plot.errors.boot.num = 1999,
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

    plot.par.mfrow.opt <- getOption("plot.par.mfrow")
    if(!is.null(plot.par.mfrow.opt))
        plot.par.mfrow <- plot.par.mfrow.opt

    miss.x <- missing(xdat)

    if(miss.x && !is.null(bws$formula)){
      tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

      xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]
    } else {
      if(miss.x && !is.null(bws$call)){
        xdat <- data.frame(.np_eval_bws_call_arg(bws, "dat"))
      }
      xdat = toFrame(xdat)
      xdat = na.omit(xdat)
    }

    xq = double(bws$ndim)+xq
    xtrim = double(bws$ndim)+xtrim

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
    plot.errors.boot.method <- normalized.opts$plot.errors.boot.method
    plot.errors.boot.nonfixed <- normalized.opts$plot.errors.boot.nonfixed
    plot.errors.boot.blocklen <- normalized.opts$plot.errors.boot.blocklen
    plot.errors.center <- normalized.opts$plot.errors.center
    plot.errors.type <- normalized.opts$plot.errors.type
    plot.errors.alpha <- normalized.opts$plot.errors.alpha
    plot.errors.style <- normalized.opts$plot.errors.style
    plot.errors.bar <- normalized.opts$plot.errors.bar
    common.scale <- normalized.opts$common.scale
    plot.errors <- normalized.opts$plot.errors

    if ((bws$ncon + bws$nord == 2) && (bws$nuno == 0) && perspective &&
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

      tobj <- .np_plot_unconditional_eval(
        xdat = xdat,
        exdat = x.eval,
        bws = bws,
        cdf = TRUE,
        need.asymptotic = identical(plot.errors.method, "asymptotic")
      )

      tdens = matrix(data = tobj$dist,
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      terr = matrix(data = tobj$derr, nrow = nrow(x.eval), ncol = 3)
      terr[,3] = NA
      
      lerr.all <- NULL
      herr.all <- NULL

      if (plot.errors.method == "bootstrap"){
        terr.obj <- compute.bootstrap.errors(xdat = xdat, 
          exdat = x.eval,
          slice.index = 0,
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
          plot.errors.boot.num = plot.errors.boot.num,
          plot.errors.center = plot.errors.center,
          plot.errors.type = plot.errors.type,
          plot.errors.alpha = plot.errors.alpha,
          bws = bws)
        terr <- terr.obj[["boot.err"]]
        terr.all <- terr.obj[["boot.all.err"]]

        pc = (plot.errors.center == "bias-corrected")
        center.val <- if (pc) terr[,3] else tobj$dist

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
        terr.obj <- .np_plot_asymptotic_error_from_se(
          se = tobj$derr,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type,
          m = nrow(x.eval)
        )
        terr[,1:2] <- terr.obj$err
        terr.all <- terr.obj$all.err
        center.val <- tobj$dist
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
      }

      if(is.null(zlim)) {
          zlim =
              if (plot.errors){
                  if (plot.errors.type == "all" && !is.null(lerr.all))
                  c(min(c(unlist(lerr.all), lerr)), max(c(unlist(herr.all), herr)))
                  else
                      c(min(lerr),max(herr))
              } else
                  c(min(tobj$dist),max(tobj$dist))
      }

      if (plot.behavior != "plot"){
        d1 <- npdistribution(bws = bws, eval = x.eval,
                             dist = tobj$dist, derr = terr[,1:2], ntrain = nrow(xdat))
        d1$bias = NA

        if (plot.errors.center == "bias-corrected")
          d1$bias = terr[,3] - tobj$dist
        
        if (plot.behavior == "data")
          return ( list(d1 = d1) )
      }
      
      # rows = constant x2
      # cols = constant x1

      dtheta = 5.0
      dphi = 10.0
      persp.col = if (plot.errors) FALSE else scalar_default(col, "lightblue")
      
      for (i in 0:((360 %/% dtheta - 1)*rotate)*dtheta+theta){
          persp(x1.eval,
                x2.eval,
                tdens,
                zlim = zlim,
                col = persp.col,
                border = scalar_default(border, "black"),
                ticktype = "detailed",
                cex.axis = scalar_default(cex.axis, par()$cex.axis),
                cex.lab = scalar_default(cex.lab, par()$cex.lab),
                cex.main = scalar_default(cex.main, par()$cex.main),
                cex.sub = scalar_default(cex.sub, par()$cex.sub),
                xlab = scalar_default(xlab, gen.label(names(xdat)[1], "X1")),
                ylab = scalar_default(ylab, gen.label(names(xdat)[2], "X2")),
                zlab = scalar_default(zlab, "Joint Distribution"),
                theta = i,
                phi = phi,
                main = gen.tflabel(!is.null(main), main, paste("[theta= ", i,", phi= ", phi,"]", sep="")))

          if (plot.errors){
            par(new = TRUE)
            if (plot.errors.type == "all" && !is.null(lerr.all)) {
              band.cols <- c(pointwise = "red", simultaneous = "green3", bonferroni = "blue")
              for (bn in c("pointwise", "simultaneous", "bonferroni")) {
                persp(x1.eval,
                      x2.eval,
                      lerr.all[[bn]],
                      zlim = zlim,
                      cex.axis = scalar_default(cex.axis, par()$cex.axis),
                      cex.lab = scalar_default(cex.lab, par()$cex.lab),
                      cex.main = scalar_default(cex.main, par()$cex.main),
                      cex.sub = scalar_default(cex.sub, par()$cex.sub),
                      col = persp.col,
                      border = band.cols[bn],
                      ticktype = "detailed",
                      xlab = "",
                      ylab = "",
                      zlab = "",
                      theta = i,
                      phi = phi,
                      lwd = scalar_default(lwd, par()$lwd))
                par(new = TRUE)
              }
            } else {
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
            }

            par(new = TRUE)
            if (plot.errors.type == "all" && !is.null(herr.all)) {
              band.cols <- c(pointwise = "red", simultaneous = "green3", bonferroni = "blue")
              for (bn in c("pointwise", "simultaneous", "bonferroni")) {
                persp(x1.eval,
                      x2.eval,
                      herr.all[[bn]],
                      zlim = zlim,
                      cex.axis = scalar_default(cex.axis, par()$cex.axis),
                      cex.lab = scalar_default(cex.lab, par()$cex.lab),
                      cex.main = scalar_default(cex.main, par()$cex.main),
                      cex.sub = scalar_default(cex.sub, par()$cex.sub),
                      col = persp.col,
                      border = band.cols[bn],
                      ticktype = "detailed",
                      xlab = "",
                      ylab = "",
                      zlab = "",
                      theta = i,
                      phi = phi,
                      lwd = scalar_default(lwd, par()$lwd))
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
          }

          Sys.sleep(0.5)
      }


      if (plot.behavior == "plot-data")
        return ( list(d1 = d1) )        

    } else {

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(bws$ndim),cex=par()$cex)

      ev = xdat[1,,drop = FALSE]

      for (i in seq_len(bws$ndim))
        ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      maxneval = max(c(sapply(xdat,nlevels),neval))

      exdat = xdat[rep(1, maxneval), , drop = FALSE]

      for (i in seq_len(bws$ndim))
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

      plotOnEstimate = (plot.errors.center == "estimate")

      ## density / distribution expressions

      for (i in seq_len(bws$ndim)){
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

        eval.slice <- subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE]
        tobj <- .np_plot_unconditional_eval(
          xdat = xdat,
          exdat = eval.slice,
          bws = bws,
          cdf = TRUE,
          need.asymptotic = identical(plot.errors.method, "asymptotic")
        )
        temp.dens[seq_len(xi.neval)] <- tobj$dist

        if (plot.errors){
          if (plot.errors.method == "asymptotic") {
            asym.obj <- .np_plot_asymptotic_error_from_se(
              se = tobj$derr,
              alpha = plot.errors.alpha,
              band.type = plot.errors.type,
              m = xi.neval
            )
            temp.err[seq_len(xi.neval),1:2] <- asym.obj$err
            temp.all.err <- asym.obj$all.err
          } else if (plot.errors.method == "bootstrap"){
            temp.boot.raw <- compute.bootstrap.errors(
                      xdat = xdat,
                      exdat = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
                      slice.index = i,
                      plot.errors.boot.method = plot.errors.boot.method,
                      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
                      bws = bws)
            temp.err[seq_len(xi.neval),] = temp.boot.raw[["boot.err"]]
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
              c(min(na.omit(c(temp.dens - temp.err[,1], temp.err[,3] - temp.err[,1]))),
                max(na.omit(c(temp.dens + temp.err[,2], temp.err[,3] + temp.err[,2]))))
          plot.args$xlab <- scalar_default(xlab, gen.label(bws$xnames[i], paste("X", i, sep = "")))
          plot.args$ylab <- scalar_default(ylab, "Distribution")
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
                center = as.numeric(na.omit(if (plotOnEstimate) temp.dens else temp.err[,3])),
                all.err = temp.all.err,
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(ei)),
                ely = if (plotOnEstimate) na.omit(temp.dens - temp.err[,1]) else na.omit(temp.err[,3] - temp.err[,1]),
                ehy = if (plotOnEstimate) na.omit(temp.dens + temp.err[,2]) else na.omit(temp.err[,3] + temp.err[,2]),
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
          plot.out[i] = NA
          plot.out[[i]] <- npdistribution(
            bws = bws,
            eval = eval.slice,
            dist = na.omit(temp.dens),
            derr = na.omit(cbind(-temp.err[,1], temp.err[,2])),
            ntrain = bws$nobs
          )
          plot.out[[i]]$bias = na.omit(temp.dens - temp.err[,3])
          plot.out[[i]]$bxp = temp.boot
        }
      }
      
      if (common.scale && (plot.behavior != "data")){
        jj = seq_len(bws$ndim)*3

        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in seq_len(bws$ndim)) {
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
        
        for (i in seq_len(bws$ndim)){
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
          plot.args$xlab <- scalar_default(xlab, gen.label(bws$xnames[i], paste("X", i, sep = "")))
          plot.args$ylab <- scalar_default(ylab, "Distribution")
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
                ex = as.numeric(na.omit(allei[,i])),
                center = as.numeric(na.omit(if (plotOnEstimate) data.eval[,i] else data.err[,3*i])),
                all.err = data.err.all[[i]],
                plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              draw.args <- list(
                ex = as.numeric(na.omit(allei[,i])),
                ely = if (plotOnEstimate) na.omit(data.eval[,i] - data.err[,3*i-2]) else na.omit(data.err[,3*i] - data.err[,3*i-2]),
                ehy = if (plotOnEstimate) na.omit(data.eval[,i] + data.err[,3*i-1]) else na.omit(data.err[,3*i] + data.err[,3*i-1]),
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
        names(plot.out) = paste("d",seq_len(bws$ndim),sep="")
        return (plot.out)
      }
    }
  }
