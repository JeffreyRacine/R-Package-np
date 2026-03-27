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
           renderer = c("base", "rgl"),
           gradients = FALSE,
           coef = FALSE,
           main = NULL,
           type = NULL,
           border = NULL,
           cex = NULL,
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
           view = c("rotate","fixed"),
           plot.behavior = c("plot","plot-data","data"),
           plot.errors.method = c("none","bootstrap","asymptotic"),
           plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
           plot.errors.boot.nonfixed = c("exact", "frozen"),
           plot.errors.boot.wild = c("rademacher", "mammen"),
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
           plot.data.overlay = TRUE,
           plot.rug = FALSE,
           ...,
           random.seed){

    engine.ctx <- .np_plot_engine_begin(plot.par.mfrow = plot.par.mfrow)
    on.exit(.np_plot_restore_par(engine.ctx$oldpar), add = TRUE)
    plot.par.mfrow <- engine.ctx$plot.par.mfrow
    scalar_default <- .np_plot_scalar_default

    dots <- list(...)
    plot.user.args <- .np_plot_user_args(dots, "plot")
    points.user.args <- .np_plot_user_args(dots, "points")
    persp.user.args <- .np_plot_user_args(dots, "persp")
    bxp.user.args <- .np_plot_user_args(dots, "bxp")
    rgl.persp3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.persp3d", "rgl.persp3d.")
    rgl.view3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.view3d", "rgl.view3d.")
    rgl.par3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.par3d", "rgl.par3d.")
    rgl.grid3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.grid3d", "rgl.grid3d.")
    rgl.widget.user.args <- .np_plot_collect_rgl_args(dots, "rgl.widget", "rgl.widget.")
    rgl.legend3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.legend3d", "rgl.legend3d.")
    rgl.points3d.user.args <- .np_plot_extract_prefixed_args(dots, "rgl.points3d.")
    rgl.surface3d.user.args <- .np_plot_extract_prefixed_args(dots, "rgl.surface3d.")
    if (!is.null(cex)) {
      if (is.null(plot.user.args$cex)) plot.user.args$cex <- cex
      if (is.null(points.user.args$cex)) points.user.args$cex <- cex
      if (is.null(bxp.user.args$cex)) bxp.user.args$cex <- cex
    }
    overlay.col <- points.user.args$col
    overlay.points.args <- points.user.args

    bxp.args <- bxp.user.args
    if (!is.null(col)) bxp.args$col <- col
    if (!is.null(lty)) bxp.args$lty <- lty
    if (!is.null(lwd)) bxp.args$lwd <- lwd
    if (!is.null(border)) bxp.args$border <- border

    if(!missing(gradients))
      stop("gradients not supported with partially linear models. Coefficients may be extracted with coef()")
    coef <- isTRUE(coef)

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
                                      env = environment(tt))
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
      ylim = ylim,
      allow_asymptotic_quantile = FALSE
    )

    plot.behavior <- .npRmpi_plot_behavior_for_rank(normalized.opts$plot.behavior)
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

    plot.errors = (plot.errors.method != "none")
    first.render <- .np_plot_first_render_state()
    on.exit(.np_plot_activity_end(first.render$activity), add = TRUE)
    if (identical(.np_plot_match_renderer(renderer), "rgl") && coef) {
      .np_plot_validate_renderer_request(
        renderer = renderer,
        route = "npplreg/npplregbw",
        perspective = perspective,
        supported.route = TRUE,
        view = if (missing(view)) "fixed" else view,
        gradients = FALSE,
        coef = coef,
        plot.errors.method = plot.errors.method,
        plot.data.overlay = FALSE,
        plot.behavior = plot.behavior
      )
    }
    if (coef) {
      fit.coef <- if (identical(plot.errors.method, "asymptotic")) {
        npplreg(
          bws = bws,
          txdat = xdat,
          tydat = ydat,
          tzdat = zdat
        )
      } else {
        .np_plot_plreg_local_fit(
          bws = bws,
          xdat = xdat,
          ydat = ydat,
          zdat = zdat
        )
      }
      cf <- as.double(fit.coef$xcoef)
      se <- as.double(fit.coef$xcoeferr)
      cf.names <- if (!is.null(names(fit.coef$xcoef))) names(fit.coef$xcoef) else bws$xnames
      if (is.null(cf.names) || length(cf.names) != length(cf))
        cf.names <- paste0("x", seq_along(cf))

      if (plot.behavior != "data") {
        bp <- barplot(cf,
                      names.arg = cf.names,
                      ylab = scalar_default(ylab, "Linear Coefficient"),
                      xlab = scalar_default(xlab, "Regressor"),
                      main = scalar_default(main, ""),
                      col = scalar_default(col, "gray70"),
                      border = scalar_default(border, par("fg")))
        if (plot.errors && length(se) == length(cf) && all(is.finite(se))) {
          arrows(bp, cf - se, bp, cf + se, angle = 90, code = 3, length = 0.05, lwd = 1)
        }
      }

      if (plot.behavior == "plot")
        return(invisible(NULL))

      return(list(coefficients = cf, coefficient.stderr = se, fit = fit.coef))
    }

    overlay.ok <- .np_plot_overlay_enabled(
      plot.data.overlay = plot.data.overlay,
      plot.behavior = plot.behavior,
      gradients = FALSE,
      coef = coef,
      plot.data.overlay.missing = missing(plot.data.overlay)
    )

    surface.supported <-
      (nxcon + nxord == 1) &&
      (nzcon + nzord == 1) &&
      (nxuno + nzuno == 0) &&
      !any(xor(bws$xdati$iord, bws$xdati$inumord)) &&
      !any(xor(bws$zdati$iord, bws$zdati$inumord))

    renderer <- .np_plot_validate_renderer_request(
      renderer = renderer,
      route = "npplreg/npplregbw",
      perspective = perspective,
      supported.route = surface.supported,
      view = if (missing(view)) "fixed" else view,
      gradients = FALSE,
      coef = coef,
      plot.errors.method = plot.errors.method,
      plot.data.overlay = overlay.ok,
      plot.behavior = plot.behavior,
      allow.plot.errors = TRUE,
      allow.plot.data.overlay = TRUE
    )
    plot.rug <- .np_plot_validate_rug_request(
      plot.rug = plot.rug,
      route = "plot.plbandwidth()",
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

    if (surface.supported &&
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

      tobj <- if (identical(plot.errors.method, "asymptotic")) {
        .np_plot_plreg_asymptotic_fit(
          bws = bws,
          xdat = xdat,
          ydat = ydat,
          zdat = zdat,
          exdat = x.eval[,1, drop = FALSE],
          ezdat = x.eval[,2, drop = FALSE]
        )
      } else {
        .np_plot_plreg_local_fit(
          bws = bws,
          xdat = xdat,
          ydat = ydat,
          zdat = zdat,
          exdat = x.eval[,1, drop = FALSE],
          ezdat = x.eval[,2, drop = FALSE]
        )
      }

      terr = matrix(data = NA, nrow = nrow(x.eval), ncol = 3)
      lerr.all <- NULL
      herr.all <- NULL
      
      treg = matrix(data = tobj$mean,
        nrow = x1.neval, ncol = z1.neval, byrow = FALSE)

      if (plot.errors.method == "asymptotic") {
        terr.obj <- .np_plot_asymptotic_error_from_se(
          se = tobj$merr,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type,
          m = nrow(x.eval)
        )
        terr[,1:2] <- terr.obj$err
        terr.all <- terr.obj$all.err
        lerr = matrix(data = tobj$mean - terr[,1],
          nrow = x1.neval, ncol = z1.neval, byrow = FALSE)
        herr = matrix(data = tobj$mean + terr[,2],
          nrow = x1.neval, ncol = z1.neval, byrow = FALSE)
        if (plot.errors.type == "all" && !is.null(terr.all)) {
          lerr.all <- lapply(terr.all, function(te)
            matrix(data = tobj$mean - te[,1], nrow = x1.neval, ncol = z1.neval, byrow = FALSE))
          herr.all <- lapply(terr.all, function(te)
            matrix(data = tobj$mean + te[,2], nrow = x1.neval, ncol = z1.neval, byrow = FALSE))
        }

      } else if (plot.errors.method == "bootstrap"){
        terr.obj <- compute.bootstrap.errors(
          xdat = xdat, ydat = ydat, zdat = zdat,
          exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,2, drop = FALSE],
          gradients = FALSE,
          slice.index = 0,
          progress.target = NULL,
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
          plot.errors.boot.wild = plot.errors.boot.wild,
          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
          plot.errors.boot.num = plot.errors.boot.num,
          plot.errors.center = plot.errors.center,
          plot.errors.type = plot.errors.type,
          plot.errors.alpha = plot.errors.alpha,
          bws = bws)
        terr <- terr.obj[["boot.err"]]
        terr.all <- terr.obj[["boot.all.err"]]

        pc = (plot.errors.center == "bias-corrected")
        center.val <- if(pc) terr[,3] else tobj$mean

        lerr = matrix(data = center.val - terr[,1],
          nrow = x1.neval, ncol = z1.neval, byrow = FALSE)

        herr = matrix(data = center.val + terr[,2],
          nrow = x1.neval, ncol = z1.neval, byrow = FALSE)
        if (plot.errors.type == "all" && !is.null(terr.all)) {
          lerr.all <- lapply(terr.all, function(te)
            matrix(data = center.val - te[,1], nrow = x1.neval, ncol = z1.neval, byrow = FALSE))
          herr.all <- lapply(terr.all, function(te)
            matrix(data = center.val + te[,2], nrow = x1.neval, ncol = z1.neval, byrow = FALSE))
        }

      }
      
      if(is.null(zlim)) {
          zlim =
              if (plot.errors && plot.errors.type == "all" && !is.null(lerr.all))
                  c(min(c(unlist(lerr.all), lerr)), max(c(unlist(herr.all), herr)))
              else if (plot.errors)
                  c(min(lerr),max(herr))
              else
                  c(min(tobj$mean),max(tobj$mean))
      }
      if (overlay.ok)
        zlim <- .np_plot_overlay_range(zlim, ydat)
        
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

      xlab.val <- scalar_default(xlab, gen.label(names(xdat)[1], "X1"))
      ylab.val <- scalar_default(ylab, gen.label(names(zdat)[1], "Z1"))
      zlab.val <- scalar_default(zlab, gen.label(names(ydat),"Conditional Mean"))

      if (identical(renderer, "rgl")) {
        rgl.view <- .np_plot_rgl_view_angles(theta = theta, phi = phi)
        main.val <- if (!is.null(main)) main else NULL
        overlay.x1 <- xdat[,1]
        if (is.factor(overlay.x1) || is.ordered(overlay.x1))
          overlay.x1 <- (bws$xdati$all.dlev[[1]])[as.integer(overlay.x1)]
        overlay.x2 <- zdat[,1]
        if (is.factor(overlay.x2) || is.ordered(overlay.x2))
          overlay.x2 <- (bws$zdati$all.dlev[[1]])[as.integer(overlay.x2)]
        .np_plot_first_render_begin(first.render)
        rgl.out <- .np_plot_render_surface_rgl(
          x = x1.eval,
          y = z1.eval,
          z = treg,
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
                x1 = overlay.x1,
                x2 = overlay.x2,
                zlim = zlim
              )
            }
            if (plot.errors) {
              .np_plot_error_surfaces_rgl(
                x = x1.eval,
                y = z1.eval,
                plot.errors.type = plot.errors.type,
                lerr = lerr,
                herr = herr,
                lerr.all = lerr.all,
                herr.all = herr.all,
                surface3d.args = rgl.surface3d.user.args,
                legend3d.args = rgl.legend3d.user.args
              )
            }
            if (overlay.ok) {
              .np_plot_overlay_points_rgl(
                x1 = overlay.x1,
                x2 = overlay.x2,
                y = ydat,
                points3d.args = rgl.points3d.user.args
              )
            }
          }
        )
        .np_plot_first_render_end(first.render)
        return(.np_plot_rgl_finalize(
          rgl.out = rgl.out,
          plot.behavior = plot.behavior,
          plot.data = list(r1 = r1)
        ))
      }

      rotate.defaults <- .np_plot_rotate_defaults()
      dtheta = rotate.defaults$dtheta

      persp.col = grDevices::adjustcolor(
        .np_plot_persp_surface_colors(z = treg, col = col),
        alpha.f = 0.5
      )
      overlay.x1 <- xdat[,1]
      if (is.factor(overlay.x1) || is.ordered(overlay.x1))
        overlay.x1 <- (bws$xdati$all.dlev[[1]])[as.integer(overlay.x1)]
      overlay.x2 <- zdat[,1]
      if (is.factor(overlay.x2) || is.ordered(overlay.x2))
        overlay.x2 <- (bws$zdati$all.dlev[[1]])[as.integer(overlay.x2)]
      
      ##for (j in 0:((50 %/% dphi - 1)*rotate)*dphi+phi){
        frame.theta <- (0:((360 %/% dtheta - 1L) * rotate)) * dtheta + theta
        rotation.progress <- .np_plot_rotation_progress_begin(length(frame.theta))
        on.exit(.np_plot_rotation_progress_end(rotation.progress), add = TRUE)
        for (frame.idx in seq_along(frame.theta)){
          i <- frame.theta[[frame.idx]]
          .np_plot_first_render_begin(first.render)
          persp.args <- list(x = x1.eval,
                             y = z1.eval,
                             z = treg,
                             zlim = zlim,
                             col = persp.col,
                             border = scalar_default(border, "black"),
                             ticktype = "detailed",
                             cex.axis = scalar_default(cex.axis, par()$cex.axis),
                             cex.lab = scalar_default(cex.lab, par()$cex.lab),
                             cex.main = scalar_default(cex.main, par()$cex.main),
                             cex.sub = scalar_default(cex.sub, par()$cex.sub),
                             lwd = 0.8 * scalar_default(lwd, par()$lwd),
                             xlab = scalar_default(xlab, gen.label(names(xdat)[1], "X1")),
                             ylab = scalar_default(ylab, gen.label(names(xdat)[2], "Z1")),
                             zlab = scalar_default(zlab, gen.label(names(ydat),"Conditional Mean")),
                             theta = i,
                             phi = phi,
                             main = scalar_default(main, ""))
          persp.args <- .np_plot_merge_user_args(persp.args, persp.user.args)
          persp.mat <- do.call(persp, persp.args)
          .np_plot_first_render_end(first.render)
          .np_plot_draw_box_grid_persp(
            xlim = range(x1.eval, finite = TRUE),
            ylim = range(z1.eval, finite = TRUE),
            zlim = zlim,
            persp.mat = persp.mat
          )
          if (plot.rug) {
            .np_plot_draw_floor_rug_persp(
              x1 = overlay.x1,
              x2 = overlay.x2,
              zlim = zlim,
              persp.mat = persp.mat
            )
          }

          if (plot.errors){
            .np_plot_draw_error_wireframes_persp(
              x = x1.eval,
              y = z1.eval,
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
          if (overlay.ok)
            do.call(.np_plot_overlay_points_persp,
                    c(list(x1 = overlay.x1,
                           x2 = overlay.x2,
                           y = ydat,
                           persp.mat = persp.mat),
                      overlay.points.args))

          rotation.progress <- .np_plot_rotation_progress_tick(rotation.progress, done = frame.idx)
          Sys.sleep(if (isTRUE(rotate)) rotate.defaults$sleep else 0.5)
        }
      ##}

      if (plot.behavior == "plot-data")
        return ( list(r1 = r1) )

    } else {
      plot.layout <- .np_plot_layout_begin(
        plot.behavior = plot.behavior,
        plot.par.mfrow = plot.par.mfrow,
        mfrow = n2mfrow(bws$xndim + bws$zndim)
      )

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

      all.isFactor = c(vapply(xdat, is.factor, logical(1)), vapply(zdat, is.factor, logical(1)))

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


        tobj <- if (identical(plot.errors.method, "asymptotic")) {
          .np_plot_plreg_asymptotic_fit(
            bws = bws,
            xdat = xdat,
            ydat = ydat,
            zdat = zdat,
            exdat = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
            ezdat = ezdat[seq_len(xi.neval),, drop = FALSE]
          )
        } else {
          .np_plot_plreg_local_fit(
            bws = bws,
            xdat = xdat,
            ydat = ydat,
            zdat = zdat,
            exdat = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
            ezdat = ezdat[seq_len(xi.neval),, drop = FALSE]
          )
        }

        temp.mean[seq_len(xi.neval)] = tobj$mean

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
            temp.boot <- compute.bootstrap.errors(
                      xdat = xdat,
                      ydat = ydat,
                      zdat = zdat,
                      exdat = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
                      ezdat = ezdat[seq_len(xi.neval),, drop = FALSE],
                      gradients = gradients,
                      slice.index = plot.index,
                      progress.target = .np_plot_scoef_bootstrap_target_label(
                        bws = bws,
                        slice.index = plot.index
                      ),
                      plot.errors.boot.method = plot.errors.boot.method,
                      t0.override = as.vector(tobj$mean),
                      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
                      plot.errors.boot.wild = plot.errors.boot.wild,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
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
          panel.ylim <- NULL
          if (plot.errors)
            panel.ylim <- if (plot.errors.type == "all")
              compute.all.error.range(if (plotOnEstimate) temp.mean else temp.err[,3], temp.all.err)
            else
              c(min(na.omit(c(temp.mean - temp.err[,1], temp.err[,3] - temp.err[,1]))),
                max(na.omit(c(temp.mean + temp.err[,2], temp.err[,3] + temp.err[,2]))))
          if (overlay.ok)
            panel.ylim <- .np_plot_overlay_range(panel.ylim, ydat)
          if (!is.null(panel.ylim))
            plot.args$ylim <- panel.ylim
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
          if (overlay.ok && !xi.factor) {
            type.val <- plot.args$type
            plot.args$type <- "n"
            .np_plot_first_render_begin(first.render)
            do.call(plot.fun, plot.args)
            .np_plot_first_render_end(first.render)
            overlay.x <- if (xOrZ == "x") xdat[,i] else zdat[,i]
            do.call(.np_plot_overlay_points_1d,
                    c(list(x = overlay.x, y = ydat),
                      overlay.points.args))
            if (!identical(type.val, "n")) {
              ok.line <- is.finite(ei) & is.finite(temp.mean)
              line.args <- list(x = ei[ok.line],
                                y = temp.mean[ok.line],
                                type = type.val,
                                lty = plot.args$lty,
                                lwd = plot.args$lwd,
                                col = plot.args$col)
              do.call(lines, line.args)
            }
          } else if (overlay.ok && xi.factor) {
            axis.labels <- levels(ei)
            axis.at <- seq_along(axis.labels)
            add.axis <- is.null(plot.user.args$xaxt)
            base.args <- list(x = as.numeric(ei),
                              y = temp.mean,
                              type = "n",
                              xlab = plot.args$xlab,
                              ylab = plot.args$ylab,
                              main = plot.args$main,
                              sub = plot.args$sub)
            if (!is.null(plot.args$ylim))
              base.args$ylim <- plot.args$ylim
            axis.lim <- if (xOrZ == "x") xlim else zlim
            if (!is.null(axis.lim)) {
              base.args$xlim <- axis.lim
            } else {
              base.args$xlim <- c(0.5, length(axis.labels) + 0.5)
            }
            if (add.axis)
              base.args$xaxt <- "n"
            base.args <- .np_plot_merge_user_args(base.args, plot.user.args)
            .np_plot_first_render_begin(first.render)
            do.call(graphics::plot.default, base.args)
            .np_plot_first_render_end(first.render)
            if (add.axis)
              axis(1, at = axis.at, labels = axis.labels)
            overlay.x <- if (xOrZ == "x") xdat[,i] else zdat[,i]
            do.call(.np_plot_overlay_points_factor,
                    c(list(x = overlay.x, y = ydat),
                      overlay.points.args))
            if (plot.bootstrap && plot.bxp) {
              do.call(bxp, c(list(z = temp.boot, add = TRUE), bxp.args))
            } else {
              l.f <- rep(ei, each = 3)
              l.f[3 * seq_along(ei)] <- NA
              l.y <- unlist(lapply(temp.mean, function(p) c(0, p, NA)))
              lines(x = l.f, y = l.y, lty = 2)
              point.args <- list(x = ei, y = temp.mean)
              if (!is.null(col)) point.args$col <- col
              if (!is.null(points.user.args$pch)) point.args$pch <- points.user.args$pch
              if (!is.null(points.user.args$cex)) point.args$cex <- points.user.args$cex
              if (!is.null(points.user.args$bg)) point.args$bg <- points.user.args$bg
              do.call(points, point.args)
            }
          } else {
            .np_plot_first_render_begin(first.render)
            do.call(plot.fun, plot.args)
            .np_plot_first_render_end(first.render)
            if (plot.rug)
              .np_plot_draw_rug_1d(if (xOrZ == "x") xdat[,i] else zdat[,i])
          }

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
                           evalx = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
                           evalz = ezdat[seq_len(xi.neval),, drop = FALSE],
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

        tobj <- if (identical(plot.errors.method, "asymptotic")) {
          .np_plot_plreg_asymptotic_fit(
            bws = bws,
            xdat = xdat,
            ydat = ydat,
            zdat = zdat,
            exdat = exdat[seq_len(xi.neval),, drop = FALSE],
            ezdat = subcol(ezdat,ei,i)[seq_len(xi.neval),, drop = FALSE]
          )
        } else {
          .np_plot_plreg_local_fit(
            bws = bws,
            xdat = xdat,
            ydat = ydat,
            zdat = zdat,
            exdat = exdat[seq_len(xi.neval),, drop = FALSE],
            ezdat = subcol(ezdat,ei,i)[seq_len(xi.neval),, drop = FALSE]
          )
        }

        temp.mean[seq_len(xi.neval)] = tobj$mean

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
            temp.boot <- compute.bootstrap.errors(
                      xdat = xdat,
                      ydat = ydat,
                      zdat = zdat,
                      exdat = exdat[seq_len(xi.neval),, drop = FALSE],
                      ezdat = subcol(ezdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
                      gradients = gradients,
                      slice.index = plot.index,
                      progress.target = .np_plot_scoef_bootstrap_target_label(
                        bws = bws,
                        slice.index = plot.index
                      ),
                      plot.errors.boot.method = plot.errors.boot.method,
                      t0.override = as.vector(tobj$mean),
                      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
                      plot.errors.boot.wild = plot.errors.boot.wild,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
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
          .np_plot_first_render_begin(first.render)
          do.call(plot.fun, plot.args)
          .np_plot_first_render_end(first.render)
          if (plot.rug && !xi.factor)
            .np_plot_draw_rug_1d(if (xOrZ == "x") xdat[,i] else zdat[,i])

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
                           evalx = exdat[seq_len(xi.neval),, drop = FALSE],
                           evalz = subcol(ezdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
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

        if (overlay.ok) {
          y.range <- .np_plot_overlay_range(c(y.min, y.max), ydat)
          y.min <- y.range[1L]
          y.max <- y.range[2L]
        }

        if(!is.null(ylim)){
          y.min = ylim[1]
          y.max = ylim[2]
        }
        
        xOrZ = "x"
        
        for (plot.index in seq_len(bws$xndim + bws$zndim)){
          plot.layout <- .np_plot_layout_activate(plot.layout)
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
          if (overlay.ok && !xi.factor) {
            type.val <- plot.args$type
            plot.args$type <- "n"
            .np_plot_first_render_begin(first.render)
            do.call(plot.fun, plot.args)
            .np_plot_first_render_end(first.render)
            overlay.x <- if (xOrZ == "x") xdat[,i] else zdat[,i]
            do.call(.np_plot_overlay_points_1d,
                    c(list(x = overlay.x, y = ydat),
                      overlay.points.args))
            if (!identical(type.val, "n")) {
              ok.line <- is.finite(allei[,plot.index]) & is.finite(data.eval[,plot.index])
              line.args <- list(x = allei[ok.line, plot.index],
                                y = data.eval[ok.line, plot.index],
                                type = type.val,
                                lty = plot.args$lty,
                                lwd = plot.args$lwd,
                                col = plot.args$col)
              do.call(lines, line.args)
            }
          } else if (overlay.ok && xi.factor) {
            axis.labels <- levels(allei[,plot.index])
            axis.at <- seq_along(axis.labels)
            add.axis <- is.null(plot.user.args$xaxt)
            base.args <- list(x = as.numeric(allei[,plot.index]),
                              y = data.eval[,plot.index],
                              type = "n",
                              xlab = plot.args$xlab,
                              ylab = plot.args$ylab,
                              main = plot.args$main,
                              sub = plot.args$sub,
                              ylim = plot.args$ylim)
            axis.lim <- if (xOrZ == "x") xlim else zlim
            if (!is.null(axis.lim)) {
              base.args$xlim <- axis.lim
            } else {
              base.args$xlim <- c(0.5, length(axis.labels) + 0.5)
            }
            if (add.axis)
              base.args$xaxt <- "n"
            base.args <- .np_plot_merge_user_args(base.args, plot.user.args)
            .np_plot_first_render_begin(first.render)
            do.call(graphics::plot.default, base.args)
            .np_plot_first_render_end(first.render)
            if (add.axis)
              axis(1, at = axis.at, labels = axis.labels)
            overlay.x <- if (xOrZ == "x") xdat[,i] else zdat[,i]
            do.call(.np_plot_overlay_points_factor,
                    c(list(x = overlay.x, y = ydat),
                      overlay.points.args))
            if (plot.bootstrap && plot.bxp) {
              do.call(bxp, c(list(z = all.bxp[[plot.index]], add = TRUE), bxp.args))
            } else {
              l.f <- rep(allei[,plot.index], each = 3)
              l.f[3 * seq_along(allei[,plot.index])] <- NA
              l.y <- unlist(lapply(data.eval[,plot.index], function(p) c(0, p, NA)))
              lines(x = l.f, y = l.y, lty = 2)
              point.args <- list(x = allei[,plot.index], y = data.eval[,plot.index])
              if (!is.null(col)) point.args$col <- col
              if (!is.null(points.user.args$pch)) point.args$pch <- points.user.args$pch
              if (!is.null(points.user.args$cex)) point.args$cex <- points.user.args$cex
              if (!is.null(points.user.args$bg)) point.args$bg <- points.user.args$bg
              do.call(points, point.args)
            }
          } else {
            .np_plot_first_render_begin(first.render)
            do.call(plot.fun, plot.args)
            .np_plot_first_render_end(first.render)
            if (plot.rug)
              .np_plot_draw_rug_1d(if (xOrZ == "x") xdat[,i] else zdat[,i])
          }

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
            paste("plr", seq_len(bws$xndim + bws$zndim), sep = "")
        
        return (plot.out)
      }
    }
  }
