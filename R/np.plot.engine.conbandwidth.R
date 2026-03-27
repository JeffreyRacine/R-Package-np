.np_plot_conbandwidth_engine <-
  function(bws,
           xdat,
           ydat,
           data = NULL,
           xq = 0.5,
           yq = 0.5,
           xtrim = 0.0,
           ytrim = 0.0,
           neval = 50,
           gradients = FALSE,
           common.scale = TRUE,
           perspective = TRUE,
           renderer = c("base", "rgl"),
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
           phi = 20.0,
           tau = 0.5,
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
           proper = FALSE,
           proper.method = c("project"),
           proper.control = list(),
           plot.rug = FALSE,
           ...,
           random.seed){

    engine.ctx <- .np_plot_engine_begin(plot.par.mfrow = plot.par.mfrow)
    on.exit(.np_plot_restore_par(engine.ctx$oldpar), add = TRUE)
    plot.par.mfrow <- engine.ctx$plot.par.mfrow
    scalar_default <- .np_plot_scalar_default

    dots <- list(...)
    plot.user.args <- .np_plot_user_args(dots, "plot")
    bxp.user.args <- .np_plot_user_args(dots, "bxp")
    rgl.persp3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.persp3d", "rgl.persp3d.")
    rgl.view3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.view3d", "rgl.view3d.")
    rgl.par3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.par3d", "rgl.par3d.")
    rgl.grid3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.grid3d", "rgl.grid3d.")
    rgl.widget.user.args <- .np_plot_collect_rgl_args(dots, "rgl.widget", "rgl.widget.")
    rgl.legend3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.legend3d", "rgl.legend3d.")
    rgl.surface3d.user.args <- .np_plot_extract_prefixed_args(dots, "rgl.surface3d.")
    bxp.args <- bxp.user.args
    if (!is.null(col)) bxp.args$col <- col
    if (!is.null(lty)) bxp.args$lty <- lty
    if (!is.null(lwd)) bxp.args$lwd <- lwd
    if (!is.null(border)) bxp.args$border <- border

    cdf <- FALSE
    quantreg <- FALSE
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
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

      ydat <- tmf[, bws$variableNames[["response"]], drop = FALSE]
      xdat <- tmf[, bws$variableNames[["terms"]], drop = FALSE]
    } else {
      if(all(miss.xy) && !is.null(bws$call)){
        xdat <- data.frame(.np_eval_bws_call_arg(bws, "xdat"))
        ydat <- data.frame(.np_eval_bws_call_arg(bws, "ydat"))
      }

      ## catch and destroy NA's
      xdat = toFrame(xdat)
      ydat = toFrame(ydat)
      
      keep.rows <- rep_len(TRUE, nrow(xdat))
      rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
      if (length(rows.omit) > 0L)
        keep.rows[as.integer(rows.omit)] <- FALSE

      if (!any(keep.rows))
        stop("Data has no rows without NAs")

      xdat <- xdat[keep.rows,,drop = FALSE]
      ydat <- ydat[keep.rows,,drop = FALSE]

    }

    if (quantreg && dim(ydat)[2] != 1)
      stop("'ydat' must have one column for quantile regression")
    
    xq = double(bws$xndim)+xq
    yq = double(bws$yndim)+yq
    
    xtrim = double(bws$xndim)+xtrim
    ytrim = double(bws$yndim)+ytrim

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

    if (plot.errors.method == "asymptotic" && quantreg && gradients){
      stop(
        "asymptotic errors are unsupported for quantile regression gradients; use bootstrap errors",
        call. = FALSE
      )
    }

    plot.errors = (plot.errors.method != "none")
    proper.args <- .np_condens_validate_proper_args(
      proper = proper,
      proper.method = proper.method,
      proper.control = proper.control
    )
    surface.supported <- isTRUE((bws$xncon + bws$xnord + bws$yncon + bws$ynord - quantreg == 2) &&
                                (bws$xnuno + bws$ynuno == 0) &&
                                !any(xor(bws$xdati$iord, bws$xdati$inumord)))
    renderer <- .np_plot_validate_renderer_request(
      renderer = renderer,
      route = "plot.conbandwidth()",
      perspective = perspective,
      supported.route = surface.supported,
      view = as.character(view)[1L],
      gradients = gradients,
      plot.errors.method = plot.errors.method,
      plot.behavior = plot.behavior,
      allow.plot.errors = TRUE
    )
    plot.rug <- .np_plot_validate_rug_request(
      plot.rug = plot.rug,
      route = "plot.conbandwidth()",
      supported.route = if (identical(renderer, "rgl")) {
        isTRUE(surface.supported && perspective)
      } else {
        TRUE
      },
      renderer = renderer,
      reason = if (identical(renderer, "rgl")) {
        "supported rgl surface routes"
      } else {
        "supported base plot routes"
      }
    )

    if (surface.supported &
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

      if (quantreg) {
        tobj <- .np_plot_quantile_eval(
          txdat = xdat,
          tydat = ydat,
          exdat = x.eval,
          tau = tau,
          bws = bws
        )
      } else {
        tobj <- .np_plot_conditional_eval(
          bws = bws,
          xdat = xdat,
          ydat = ydat,
          exdat = x.eval[,1, drop = FALSE],
          eydat = x.eval[,2, drop = FALSE],
          cdf = FALSE,
          gradients = FALSE,
          proper = isTRUE(proper.args$proper.requested),
          proper.method = proper.args$proper.method,
          proper.control = proper.args$proper.control
        )
      }
      tcomp <- switch(tboo,
                      "quant" = tobj$quantile,
                      "dist" = tobj$condist,
                      "dens" = tobj$condens)
      tcerr <- if (quantreg) tobj$quanterr else tobj$conderr
      tex <- if (quantreg) x.eval else x.eval[,1]
      tey <- if (quantreg) NA else x.eval[,2]

      tdens = matrix(data = tcomp,
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      terr = matrix(data = tcerr, nrow = length(tcomp), ncol = 3)
      terr[,3] = NA
      lerr.all <- NULL
      herr.all <- NULL
      
      if (plot.errors.method == "bootstrap"){
        terr.obj <- compute.bootstrap.errors(xdat = xdat, ydat = ydat,
          exdat = tex, eydat = tey,
          cdf = cdf,
          quantreg = quantreg,
          tau = tau,
          gradients = FALSE,
          gradient.index = 0,
          slice.index = 0,
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
          plot.errors.boot.num = plot.errors.boot.num,
          plot.errors.center = plot.errors.center,
          plot.errors.type = plot.errors.type,
          plot.errors.alpha = plot.errors.alpha,
          progress.target = NULL,
          proper = isTRUE(proper.args$proper.requested),
          proper.method = proper.args$proper.method,
          proper.control = proper.args$proper.control,
          bws = bws)
        terr <- terr.obj[["boot.err"]]
        terr.all <- terr.obj[["boot.all.err"]]

        pc = (plot.errors.center == "bias-corrected")
        center.val <- if(pc) terr[,3] else tcomp

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
          se = tcerr,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type,
          m = nrow(tex)
        )
        terr[,1:2] <- terr.obj$err
        terr.all <- terr.obj$all.err
        center.val <- tcomp
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
                c(min(tcomp),max(tcomp))
      }

      if (plot.behavior != "plot"){
        ret.fun <- switch(tboo,
                          "quant" = qregression,
                          "dist" = condistribution,
                          "dens" = condensity)
        ret.args <- list(bws = bws, xeval = tex, ntrain = dim(xdat)[1])
        if (quantreg) {
          ret.args$tau <- tau
          ret.args$quantile <- tcomp
          ret.args$quanterr <- terr[,1:2]
        } else {
          ret.args$yeval <- tey
          if (cdf) ret.args$condist <- tcomp else ret.args$condens <- tcomp
          ret.args$conderr <- terr[,1:2]
          ret.args$proper.requested <- tobj$proper.requested
          ret.args$proper.applied <- tobj$proper.applied
          ret.args$proper.method <- tobj$proper.method
          ret.args$condens.raw <- tobj$condens.raw
          ret.args$proper.info <- tobj$proper.info
        }
        cd1 <- do.call(ret.fun, ret.args)
        cd1$bias = NA

        if (plot.errors.center == "bias-corrected")
          cd1$bias = terr[,3] - tcomp
        
        if (plot.behavior == "data")
          return ( list(cd1 = cd1) )
      }


      # rows = constant x2
      # cols = constant x1

      xlab.val <- scalar_default(xlab, gen.label(names(xdat)[1], "X"))
      ylab.val <- scalar_default(ylab, gen.label(names(ydat)[1], "Y"))
      zlab.val <- scalar_default(zlab, paste("Conditional", if (cdf) "Distribution" else "Density"))

      if (identical(renderer, "rgl")) {
        rgl.first.render <- .np_plot_first_render_state()
        on.exit(.np_plot_activity_end(rgl.first.render$activity), add = TRUE)
        rgl.view <- .np_plot_rgl_view_angles(theta = theta, phi = phi)
        main.val <- if (!is.null(main)) main else NULL
        .np_plot_first_render_begin(rgl.first.render)
        rgl.out <- .np_plot_render_surface_rgl(
          x = x1.eval,
          y = x2.eval,
          z = tdens,
          zlim = zlim,
          col = col,
          border = scalar_default(border, "black"),
          xlab = xlab.val,
          ylab = ylab.val,
          zlab = zlab.val,
          theta = rgl.view$theta,
          phi = rgl.view$phi,
          main = main.val,
          par3d.args = rgl.par3d.user.args,
          view3d.args = rgl.view3d.user.args,
          persp3d.args = rgl.persp3d.user.args,
          grid3d.args = rgl.grid3d.user.args,
          widget.args = rgl.widget.user.args,
          draw.extras = function() {
            if (plot.rug) {
              .np_plot_draw_floor_rug_rgl(
                x1 = xdat[,1],
                x2 = if (quantreg) xdat[,2] else ydat[,1],
                zlim = zlim
              )
            }
            if (plot.errors) {
              .np_plot_error_surfaces_rgl(
                x = x1.eval,
                y = x2.eval,
                plot.errors.type = plot.errors.type,
                lerr = lerr,
                herr = herr,
                lerr.all = lerr.all,
                herr.all = herr.all,
                surface3d.args = rgl.surface3d.user.args,
                legend3d.args = rgl.legend3d.user.args
              )
            }
          }
        )
        .np_plot_first_render_end(rgl.first.render)
        return(.np_plot_rgl_finalize(
          rgl.out = rgl.out,
          plot.behavior = plot.behavior,
          plot.data = list(cd1 = cd1)
        ))
      }

      rotate.defaults <- .np_plot_rotate_defaults()
      dtheta = rotate.defaults$dtheta
      persp.col = grDevices::adjustcolor(
        .np_plot_persp_surface_colors(z = tdens, col = col),
        alpha.f = 0.5
      )
      first.render.activity <- NULL
      first.render.pending <- TRUE
      on.exit(.np_plot_activity_end(first.render.activity), add = TRUE)
      frame.theta <- (0:((360 %/% dtheta - 1L) * rotate)) * dtheta + theta
      rotation.progress <- .np_plot_rotation_progress_begin(length(frame.theta))
      on.exit(.np_plot_rotation_progress_end(rotation.progress), add = TRUE)
      
      for (frame.idx in seq_along(frame.theta)){
          i <- frame.theta[[frame.idx]]
          if (isTRUE(first.render.pending))
            first.render.activity <- .np_plot_activity_begin("Rendering plot surface")
          persp.mat <- persp(x1.eval,
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
                             lwd = 0.8 * scalar_default(lwd, par()$lwd),
                             xlab = xlab.val,
                             ylab = ylab.val,
                             zlab = zlab.val,
                             theta = i,
                             phi = phi,
                             main = scalar_default(main, ""))
          .np_plot_draw_box_grid_persp(
            xlim = range(x1.eval, finite = TRUE),
            ylim = range(x2.eval, finite = TRUE),
            zlim = zlim,
            persp.mat = persp.mat
          )
          if (plot.rug) {
            .np_plot_draw_floor_rug_persp(
              x1 = xdat[,1],
              x2 = if (quantreg) xdat[,2] else ydat[,1],
              zlim = zlim,
              persp.mat = persp.mat
            )
          }

          if (isTRUE(first.render.pending)) {
            .np_plot_activity_end(first.render.activity)
            first.render.activity <- NULL
            first.render.pending <- FALSE
          }

          if (plot.errors){
            .np_plot_draw_error_wireframes_persp(
              x = x1.eval,
              y = x2.eval,
              persp.mat = persp.mat,
              plot.errors.type = plot.errors.type,
              lerr = lerr,
              herr = herr,
              lerr.all = lerr.all,
              herr.all = herr.all,
              border = scalar_default(border, "grey"),
              lwd = scalar_default(lwd, par()$lwd)
            )
            if (plot.errors.type == "all" && !is.null(lerr.all) && !is.null(herr.all)) {
              band.cols <- .np_plot_all_band_colors()
              legend("topright",
                     legend = c("Pointwise","Simultaneous","Bonferroni"),
                     lty = 1,
                     col = unname(band.cols[c("pointwise", "simultaneous", "bonferroni")]),
                     lwd = 2.15 * scalar_default(lwd, par()$lwd),
                     bty = "n")
            }
          }

          rotation.progress <- .np_plot_rotation_progress_tick(rotation.progress, done = frame.idx)
          Sys.sleep(if (isTRUE(rotate)) rotate.defaults$sleep else 0.5)
      }

      if (plot.behavior == "plot-data")
        return ( list(cd1 = cd1) )
    } else {

      dsf = if (gradients) bws$xndim else 1
      tot.dim = bws$xndim + bws$yndim - quantreg

      plot.layout <- .np_plot_layout_begin(
        plot.behavior = plot.behavior,
        plot.par.mfrow = plot.par.mfrow,
        mfrow = n2mfrow(dsf * tot.dim)
      )

      x.ev = xdat[1,,drop = FALSE]
      y.ev = ydat[1,,drop = FALSE]

      for (i in seq_len(bws$xndim))
        x.ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      for (i in seq_len(bws$yndim))
        y.ev[1,i] = uocquantile(ydat[,i], prob=yq[i])


      maxneval = max(c(sapply(xdat,nlevels), sapply(ydat,nlevels), neval))

      exdat = xdat[rep(1, maxneval), , drop = FALSE]
      eydat = as.data.frame(matrix(data = 0, nrow = maxneval, ncol = bws$yndim))

      for (i in seq_len(bws$xndim))
        exdat[,i] = x.ev[1,i]

      for (i in seq_len(bws$yndim))
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

      all.isFactor = c(vapply(xdat, is.factor, logical(1)), vapply(ydat, is.factor, logical(1)))

      plot.out = list()

      temp.err = matrix(data = NA, nrow = maxneval, ncol = 3)
      temp.dens = replicate(maxneval, NA)
      ## plotting controls
      plot.bootstrap <- plot.errors.method == "bootstrap"
      tylabE <- if (quantreg) paste(tau, "quantile") else paste("Conditional", if (cdf) "Distribution" else "Density")
      plotOnEstimate <- (plot.errors.center == "estimate")

      plot.index = 0
      xOrY = "x"

      for (i in seq_len(bws$xndim)){
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

        if (quantreg) {
          tobj <- .np_plot_quantile_eval(
            txdat = xdat,
            tydat = ydat,
            exdat = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
            gradients = gradients,
            tau = tau,
            bws = bws
          )
        } else {
          tobj <- .np_plot_conditional_eval(
            bws = bws,
            xdat = xdat,
            ydat = ydat,
            exdat = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
            eydat = eydat[seq_len(xi.neval),, drop = FALSE],
            cdf = cdf,
            gradients = gradients
          )
        }
        if (!quantreg && isTRUE(proper.args$proper.requested)) {
          tobj$proper.requested <- TRUE
          tobj$proper.applied <- FALSE
          tobj$proper.method <- proper.args$proper.method
          tobj$proper.info <- .np_condens_make_reason_info(
            reason = "x_slices_not_repeated",
            supported = FALSE
          )
        }

        
        ## if there are gradients then we need to repeat the process for each component

        eval.extract <- function(obj, jj){
          if (gradients) {
            if (quantreg) obj$quantgrad[,jj] else obj$congrad[,jj]
          } else if (quantreg) {
            obj$quantile
          } else if (cdf) {
            obj$condist
          } else {
            obj$condens
          }
        }
        err.extract <- function(obj, jj){
          if (gradients) {
            if (quantreg) rep(NA_real_, length(eval.extract(obj, jj))) else obj$congerr[,jj]
          } else if (quantreg) {
            obj$quanterr
          } else {
            obj$conderr
          }
        }

        if (plot.behavior != "plot"){
          plot.out[plot.index] = NA
          plot.out[[plot.index]] = tobj
        }

        for (j in seq_len(dsf)){
          temp.boot = list()
          temp.all.err <- NULL
          temp.dens[seq_len(xi.neval)] <- eval.extract(tobj, j)
          
          if (plot.errors){
            if (plot.errors.method == "asymptotic") {
              terr.j <- err.extract(tobj, j)
              asym.obj <- .np_plot_asymptotic_error_from_se(
                  se = terr.j,
                  alpha = plot.errors.alpha,
                  band.type = plot.errors.type,
                  m = xi.neval
                )
                temp.err[seq_len(xi.neval),1:2] <- asym.obj$err
                temp.all.err <- asym.obj$all.err
            }
            else if (plot.errors.method == "bootstrap"){
              temp.boot <- compute.bootstrap.errors(
                        xdat = xdat,
                        ydat = ydat,
                        exdat = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
                        eydat = eydat[seq_len(xi.neval),, drop = FALSE],
                        cdf = cdf,
                        quantreg = quantreg,
                        tau = tau,
                        gradients = gradients,
                        gradient.index = j,
                        slice.index = plot.index,
                        plot.errors.boot.method = plot.errors.boot.method,
                        plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
                        plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                        plot.errors.boot.num = plot.errors.boot.num,
                        plot.errors.center = plot.errors.center,
                        plot.errors.type = plot.errors.type,
                        plot.errors.alpha = plot.errors.alpha,
                        progress.target = .np_plot_conditional_bootstrap_target_label(
                          bws = bws,
                          slice.index = plot.index,
                          gradients = gradients,
                          gradient.index = if (gradients) i else 0L
                        ),
                        bws = bws)
              temp.err[seq_len(xi.neval),] <- temp.boot[["boot.err"]]
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
              plot.args$y <- temp.dens
            if (plot.errors)
              plot.args$ylim <- if (plot.errors.type == "all")
                compute.all.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.all.err)
              else
                c(min(na.omit(c(temp.dens - temp.err[,1], temp.err[,3] - temp.err[,1]))),
                  max(na.omit(c(temp.dens + temp.err[,2], temp.err[,3] + temp.err[,2]))))
            plot.args$xlab <- scalar_default(xlab, gen.label(if (xOrY == "x") bws$xnames[i] else bws$ynames[i], paste(toupper(xOrY), i, sep = "")))
            plot.args$ylab <- scalar_default(ylab, if (gradients) paste("GC", j, "of", tylabE) else tylabE)
            if (!xi.factor) {
              plot.args$type <- scalar_default(type, "l")
              plot.args$lty <- scalar_default(lty, par()$lty)
              plot.args$col <- scalar_default(col, par()$col)
              plot.args$lwd <- scalar_default(lwd, par()$lwd)
              plot.args$cex.axis <- scalar_default(cex.axis, par()$cex.axis)
              plot.args$cex.lab <- scalar_default(cex.lab, par()$cex.lab)
              plot.args$cex.main <- scalar_default(cex.main, par()$cex.main)
              plot.args$cex.sub <- scalar_default(cex.sub, par()$cex.sub)
            } else {
              if (!is.null(col)) plot.args$col <- col
              if (!is.null(lty)) plot.args$lty <- lty
              if (!is.null(lwd)) plot.args$lwd <- lwd
            }
            plot.args$main <- scalar_default(main, "")
            plot.args$sub <- scalar_default(sub, "")
            plot.args <- .np_plot_merge_user_args(
              plot.args,
              if (xi.factor && plot.bootstrap && plot.bxp) bxp.args else plot.user.args
            )
            do.call(plot.fun, plot.args)
            if (plot.rug && !xi.factor)
              .np_plot_draw_rug_1d(if (xOrY == "x") xdat[,i] else ydat[,i])

            ## error plotting evaluation
            if (plot.errors && !(xi.factor && plot.bootstrap && plot.bxp)){
              if (plot.errors.type == "all") {
                draw.all.error.types(
                  ex = as.numeric(na.omit(ei)),
                  center = na.omit(if (plotOnEstimate) temp.dens else temp.err[,3]),
                  all.err = temp.all.err,
                  xi.factor = xi.factor)
              } else {
                if (!xi.factor && !plotOnEstimate)
                  lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)
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
        
          if (plot.behavior != "plot" && plot.errors) {
            err.name <- if (gradients) paste("gc", j, "err", sep = "") else if (quantreg) "quanterr" else "conderr"
            bias.name <- if (gradients) paste("gc", j, "bias", sep = "") else "bias"
            plot.out[[plot.index]][[err.name]] <- na.omit(cbind(-temp.err[,1], temp.err[,2]))
            plot.out[[plot.index]][[bias.name]] <- na.omit(temp.dens - temp.err[,3])
            plot.out[[plot.index]]$bxp <- temp.boot
          }
        }
      }

      if (!quantreg){
        xOrY = "y"
        for (i in seq_len(bws$yndim)){
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

          if (quantreg) {
            tobj <- .np_plot_quantile_eval(
              txdat = xdat,
              tydat = ydat,
              eydat = subcol(eydat,ei,i)[seq_len(xi.neval),, drop = FALSE],
              gradients = gradients,
              tau = tau,
              bws = bws
            )
          } else {
            tobj <- .np_plot_conditional_eval(
              bws = bws,
              xdat = xdat,
              ydat = ydat,
              exdat = exdat[seq_len(xi.neval),, drop = FALSE],
              eydat = subcol(eydat,ei,i)[seq_len(xi.neval),, drop = FALSE],
              cdf = cdf,
              gradients = gradients,
              proper = isTRUE(proper.args$proper.requested),
              proper.method = proper.args$proper.method,
              proper.control = proper.args$proper.control
            )
          }

          
          ## if there are gradients then we need to repeat the process for each component

          eval.extract <- function(obj, jj){
            if (gradients) {
              if (quantreg) obj$quantgrad[,jj] else obj$congrad[,jj]
            } else if (quantreg) {
              obj$quantile
            } else if (cdf) {
              obj$condist
            } else {
              obj$condens
            }
          }
          err.extract <- function(obj, jj){
            if (gradients) {
              if (quantreg) rep(NA_real_, length(eval.extract(obj, jj))) else obj$congerr[,jj]
            } else if (quantreg) {
              obj$quanterr
            } else {
              obj$conderr
            }
          }

          if (plot.behavior != "plot"){
            plot.out[plot.index] = NA
            plot.out[[plot.index]] = tobj
          }

          for (j in seq_len(dsf)){
            temp.boot = list()
            temp.all.err <- NULL
            temp.dens[seq_len(xi.neval)] <- eval.extract(tobj, j)
            
            if (plot.errors){
              if (plot.errors.method == "asymptotic") {
                terr.j <- err.extract(tobj, j)
                asym.obj <- .np_plot_asymptotic_error_from_se(
                  se = terr.j,
                  alpha = plot.errors.alpha,
                  band.type = plot.errors.type,
                  m = xi.neval
                )
                temp.err[seq_len(xi.neval),1:2] <- asym.obj$err
                temp.all.err <- asym.obj$all.err
              }
              else if (plot.errors.method == "bootstrap"){
                temp.boot <- compute.bootstrap.errors(
                          xdat = xdat,
                          ydat = ydat,
                          exdat = exdat[seq_len(xi.neval),, drop = FALSE],
                          eydat = subcol(eydat,ei,i)[seq_len(xi.neval),, drop = FALSE],
                          cdf = cdf,
                          quantreg = quantreg,
                          tau = tau,
                          gradients = gradients,
                          gradient.index = j,
                          slice.index = plot.index,
                          plot.errors.boot.method = plot.errors.boot.method,
                          plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
                          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                          plot.errors.boot.num = plot.errors.boot.num,
                          plot.errors.center = plot.errors.center,
                          plot.errors.type = plot.errors.type,
                          plot.errors.alpha = plot.errors.alpha,
                          progress.target = .np_plot_conditional_bootstrap_target_label(
                            bws = bws,
                            slice.index = plot.index,
                            gradients = gradients,
                            gradient.index = if (gradients) i else 0L
                          ),
                          proper = isTRUE(proper.args$proper.requested),
                          proper.method = proper.args$proper.method,
                          proper.control = proper.args$proper.control,
                          bws = bws)
                temp.err[seq_len(xi.neval),] <- temp.boot[["boot.err"]]
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
                plot.args$y <- temp.dens
              if (plot.errors)
                plot.args$ylim <- if (plot.errors.type == "all")
                  compute.all.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.all.err)
                else
                  c(min(na.omit(c(temp.dens - temp.err[,1], temp.err[,3] - temp.err[,1]))),
                    max(na.omit(c(temp.dens + temp.err[,2], temp.err[,3] + temp.err[,2]))))
              plot.args$xlab <- scalar_default(xlab, gen.label(if (xOrY == "x") bws$xnames[i] else bws$ynames[i], paste(toupper(xOrY), i, sep = "")))
              plot.args$ylab <- scalar_default(ylab, if (gradients) paste("GC", j, "of", tylabE) else tylabE)
              if (!xi.factor) {
                plot.args$type <- scalar_default(type, "l")
                plot.args$lty <- scalar_default(lty, par()$lty)
                plot.args$col <- scalar_default(col, par()$col)
                plot.args$lwd <- scalar_default(lwd, par()$lwd)
                plot.args$cex.axis <- scalar_default(cex.axis, par()$cex.axis)
                plot.args$cex.lab <- scalar_default(cex.lab, par()$cex.lab)
                plot.args$cex.main <- scalar_default(cex.main, par()$cex.main)
                plot.args$cex.sub <- scalar_default(cex.sub, par()$cex.sub)
              } else {
                if (!is.null(col)) plot.args$col <- col
                if (!is.null(lty)) plot.args$lty <- lty
                if (!is.null(lwd)) plot.args$lwd <- lwd
              }
              plot.args$main <- scalar_default(main, "")
              plot.args$sub <- scalar_default(sub, "")
              plot.args <- .np_plot_merge_user_args(
                plot.args,
                if (xi.factor && plot.bootstrap && plot.bxp) bxp.args else plot.user.args
              )
              do.call(plot.fun, plot.args)
              if (plot.rug && !xi.factor)
                .np_plot_draw_rug_1d(if (xOrY == "x") xdat[,i] else ydat[,i])

              ## error plotting evaluation
              if (plot.errors && !(xi.factor && plot.bootstrap && plot.bxp)){
                if (plot.errors.type == "all") {
                  draw.all.error.types(
                    ex = as.numeric(na.omit(ei)),
                    center = na.omit(if (plotOnEstimate) temp.dens else temp.err[,3]),
                    all.err = temp.all.err,
                    xi.factor = xi.factor)
                } else {
                  if (!xi.factor && !plotOnEstimate)
                    lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)
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
              
            if (plot.behavior != "plot" && plot.errors) {
              err.name <- if (gradients) paste("gc", j, "err", sep = "") else if (quantreg) "quanterr" else "conderr"
              bias.name <- if (gradients) paste("gc", j, "bias", sep = "") else "bias"
              plot.out[[plot.index]][[err.name]] <- na.omit(cbind(-temp.err[,1], temp.err[,2]))
              plot.out[[plot.index]][[bias.name]] <- na.omit(temp.dens - temp.err[,3])
              plot.out[[plot.index]]$bxp <- temp.boot
            }
          }
        }
      }
      
      if (common.scale && (plot.behavior != "data")){
        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in seq_len(tot.dim*dsf)) {
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
          jj = seq_len(dsf*tot.dim)*3
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
        
        xOrY = "x"
        
        for (plot.index in seq_len(tot.dim)){
          i = if (plot.index <= bws$xndim) plot.index else plot.index - bws$xndim

          if (plot.index > bws$xndim)
            xOrY <- "y"
            
          xi.factor = all.isFactor[plot.index]

          for (j in seq_len(dsf)){
            plot.layout <- .np_plot_layout_activate(plot.layout)
            ## plot evaluation
            idx <- (plot.index-1)*dsf+j
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
              plot.args$y <- data.eval[,idx]
            plot.args$ylim <- c(y.min, y.max)
            plot.args$xlab <- scalar_default(xlab, gen.label(if (xOrY == "x") bws$xnames[i] else bws$ynames[i], paste(toupper(xOrY), i, sep = "")))
            plot.args$ylab <- scalar_default(ylab, if (gradients) paste("GC", j, "of", tylabE) else tylabE)
            if (!xi.factor) {
              plot.args$type <- scalar_default(type, "l")
              plot.args$lty <- scalar_default(lty, par()$lty)
              plot.args$col <- scalar_default(col, par()$col)
              plot.args$lwd <- scalar_default(lwd, par()$lwd)
              plot.args$cex.axis <- scalar_default(cex.axis, par()$cex.axis)
              plot.args$cex.lab <- scalar_default(cex.lab, par()$cex.lab)
              plot.args$cex.main <- scalar_default(cex.main, par()$cex.main)
              plot.args$cex.sub <- scalar_default(cex.sub, par()$cex.sub)
            } else {
              if (!is.null(col)) plot.args$col <- col
              if (!is.null(lty)) plot.args$lty <- lty
              if (!is.null(lwd)) plot.args$lwd <- lwd
            }
            plot.args$main <- scalar_default(main, "")
            plot.args$sub <- scalar_default(sub, "")
            plot.args <- .np_plot_merge_user_args(
              plot.args,
              if (xi.factor && plot.bootstrap && plot.bxp) bxp.args else plot.user.args
            )
            do.call(plot.fun, plot.args)
            if (plot.rug && !xi.factor)
              .np_plot_draw_rug_1d(if (xOrY == "x") xdat[,i] else ydat[,i])

            ## error plotting evaluation
            if (plot.errors && !(xi.factor && plot.bootstrap && plot.bxp)){
              if (plot.errors.type == "all") {
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
                draw.args <- list(
                  ex = as.numeric(na.omit(allei[,plot.index])),
                  ely = if (plotOnEstimate) na.omit(data.eval[,idx] - data.err[,3*idx-2]) else na.omit(data.err[,3*idx] - data.err[,3*idx-2]),
                  ehy = if (plotOnEstimate) na.omit(data.eval[,idx] + data.err[,3*idx-1]) else na.omit(data.err[,3*idx] + data.err[,3*idx-1]),
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
      }

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=c(1,1),cex=par()$cex)

      if (plot.behavior != "plot"){
        names(plot.out) = paste("cd", seq_len(tot.dim), sep="")
        return (plot.out)
      }
    }
  }
