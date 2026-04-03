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
           renderer = c("base", "rgl"),
           gradients = FALSE,
           coef = FALSE,
           coef.index = 1L,
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
    plot.errors <- normalized.opts$plot.errors
    first.render <- .np_plot_first_render_state()
    on.exit(.np_plot_activity_end(first.render$activity), add = TRUE)
    if (coef && plot.errors.method != "none") {
      .np_warning("coef=TRUE currently disables plot errors for smooth coefficient plots.")
      plot.errors.method <- "none"
      plot.errors <- FALSE
    }
    overlay.ok <- .np_plot_overlay_enabled(
      plot.data.overlay = plot.data.overlay,
      plot.behavior = plot.behavior,
      gradients = FALSE,
      coef = coef,
      plot.data.overlay.missing = missing(plot.data.overlay)
    )

    surface.supported <- isTRUE(
      (sum(c(bws$xdati$icon, bws$xdati$iord, bws$zdati$icon, bws$zdati$iord)) == 2) &&
      (sum(c(bws$xdati$iuno, bws$zdati$iuno)) == 0) &&
      !any(xor(c(bws$xdati$iord, bws$zdati$iord),
               c(bws$xdati$inumord, bws$zdati$inumord)))
    )
    renderer <- .np_plot_validate_renderer_request(
      renderer = renderer,
      route = "npscoef/npscoefbw",
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
      route = "plot.scbandwidth()",
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

    if (surface.supported && perspective && !gradients &&
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
      colnames(x.eval) <- c(
        colnames(xdat)[1L],
        if (!miss.z) colnames(zdat)[1L] else colnames(xdat)[2L]
      )
      
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
      tobj <- .np_plot_with_local_compiled_eval(do.call(.np_scoef_fit_internal, scoef.args))

      terr = matrix(data = tobj$merr, nrow = dim(x.eval)[1], ncol = 3)
      terr[,3] = NA
      lerr.all <- NULL
      herr.all <- NULL
      
      treg = matrix(data = extract_scoef_value(tobj),
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      if (plot.errors.method == "bootstrap"){
        boot.args <- list(
          xdat = xdat,
          ydat = ydat,
          gradients = FALSE,
          slice.index = 0,
          progress.target = NULL,
          plot.errors.boot.method = plot.errors.boot.method,
          t0.override = as.vector(tobj$mean),
          plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
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
        terr.obj <- do.call(compute.bootstrap.errors, boot.args)
        terr <- terr.obj[["boot.err"]]
        terr.all <- terr.obj[["boot.all.err"]]

        pc = (plot.errors.center == "bias-corrected")
        center.val <- if (pc) terr[,3] else tobj$mean

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
          se = tobj$merr,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type,
          m = nrow(x.eval)
        )
        terr[,1:2] <- terr.obj$err
        terr.all <- terr.obj$all.err
        lerr = matrix(data = tobj$mean - terr[,1],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = tobj$mean + terr[,2],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)
        if (plot.errors.type == "all" && !is.null(terr.all)) {
          lerr.all <- lapply(terr.all, function(te)
            matrix(data = tobj$mean - te[,1], nrow = x1.neval, ncol = x2.neval, byrow = FALSE))
          herr.all <- lapply(terr.all, function(te)
            matrix(data = tobj$mean + te[,2], nrow = x1.neval, ncol = x2.neval, byrow = FALSE))
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

      xlab.val <- scalar_default(xlab, gen.label(bws$xnames[1], "X1"))
      ylab.val <- scalar_default(ylab, gen.label(x2.names[1], "X2"))
      zlab.val <- scalar_default(zlab, gen.label(bws$ynames,"Conditional Mean"))

      if (identical(renderer, "rgl")) {
        rgl.view <- .np_plot_rgl_view_angles(theta = theta, phi = phi)
        main.val <- if (!is.null(main)) main else NULL
        overlay.x1 <- xdat[,1]
        if (is.factor(overlay.x1) || is.ordered(overlay.x1))
          overlay.x1 <- (bws$xdati$all.dlev[[1]])[as.integer(overlay.x1)]
        if (miss.z) {
          overlay.x2 <- xdat[,2]
          if (is.factor(overlay.x2) || is.ordered(overlay.x2))
            overlay.x2 <- (bws$xdati$all.dlev[[2]])[as.integer(overlay.x2)]
        } else {
          overlay.x2 <- zdat[,1]
          if (is.factor(overlay.x2) || is.ordered(overlay.x2))
            overlay.x2 <- (bws$zdati$all.dlev[[1]])[as.integer(overlay.x2)]
        }
        .np_plot_first_render_begin(first.render)
        rgl.out <- .np_plot_render_surface_rgl(
          x = x1.eval,
          y = x2.eval,
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
      if (miss.z) {
        overlay.x2 <- xdat[,2]
        if (is.factor(overlay.x2) || is.ordered(overlay.x2))
          overlay.x2 <- (bws$xdati$all.dlev[[2]])[as.integer(overlay.x2)]
      } else {
        overlay.x2 <- zdat[,1]
        if (is.factor(overlay.x2) || is.ordered(overlay.x2))
          overlay.x2 <- (bws$zdati$all.dlev[[1]])[as.integer(overlay.x2)]
      }
      
      ##for (j in 0:((50 %/% dphi - 1)*rotate)*dphi+phi){
        frame.theta <- (0:((360 %/% dtheta - 1L) * rotate)) * dtheta + theta
        rotation.progress <- .np_plot_rotation_progress_begin(length(frame.theta))
        on.exit(.np_plot_rotation_progress_end(rotation.progress), add = TRUE)
        for (frame.idx in seq_along(frame.theta)){
          i <- frame.theta[[frame.idx]]
          .np_plot_first_render_begin(first.render)
          persp.args <- list(x = x1.eval,
                             y = x2.eval,
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
                             xlab = scalar_default(xlab, gen.label(bws$xnames[1], "X1")),
                             ylab = scalar_default(ylab, gen.label(x2.names[1], "X2")),
                             zlab = scalar_default(zlab, gen.label(bws$ynames,"Conditional Mean")),
                             theta = i,
                             phi = phi,
                             main = scalar_default(main, ""))
          persp.args <- .np_plot_merge_user_args(persp.args, persp.user.args)
          persp.mat <- do.call(persp, persp.args)
          .np_plot_first_render_end(first.render)
          .np_plot_draw_box_grid_persp(
            xlim = range(x1.eval, finite = TRUE),
            ylim = range(x2.eval, finite = TRUE),
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
          return(.np_plot_with_local_compiled_eval(do.call(.np_scoef_fit_internal, tx.args)))
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
              t0.override = as.vector(tobj$mean),
              plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
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
                                                    t0.override = as.vector(tobj$mean),
                                                    plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
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
              if (plot.rug && !xi.factor)
                .np_plot_draw_rug_1d(if (xOrZ == "x") xdat[,i] else zdat[,i])
            }

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
            if (add.axis)
              base.args$xaxt <- "n"
            axis.lim <- if (xOrZ == "x") xlim else zlim
            if (!is.null(axis.lim)) {
              base.args$xlim <- axis.lim
            } else {
              base.args$xlim <- c(0.5, length(axis.labels) + 0.5)
            }
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
            if (plot.rug && !xi.factor)
              .np_plot_draw_rug_1d(if (xOrZ == "x") xdat[,i] else zdat[,i])
          }

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
