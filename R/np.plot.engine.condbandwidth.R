.np_plot_condbandwidth_engine <-
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
           gradient.order = 1L,
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
           legend = TRUE,
           view = c("rotate","fixed"),
           plot.behavior = c("plot","plot-data","data"),
           plot.errors.method = c("none","bootstrap","asymptotic"),
           plot.errors.boot.method = c("inid", "fixed", "geom"),
           plot.errors.boot.nonfixed = c("exact", "frozen"),
           plot.errors.boot.blocklen = NULL,
           plot.errors.boot.num = 1999,
           plot.errors.center = c("estimate", "bias-corrected"),
           plot.errors.type = c("pmzsd","pointwise","bonferroni","simultaneous","all"),
           plot.errors.alpha = 0.05,
           plot.errors.style = c("band","bar"),
           plot.errors.bar = c("|","I"),
           plot.errors.bar.num = min(neval,25),
           plot.bxp = FALSE,
           plot.bxp.out = TRUE,
           plot.par.mfrow = TRUE,
           plot.data.overlay = FALSE,
           proper = FALSE,
           proper.method = c("isotonic"),
           proper.control = list(),
           plot.rug = FALSE,
           ...,
           random.seed){

    plot.gradient.available <- rep.int(TRUE, bws$xndim)
    if (isTRUE(gradients)) {
      reg.spec <- npConditionalRegEngineSpec(bws, where = "plot.condbandwidth()")
      if (identical(reg.spec$reg.engine, "lp") && bws$xncon > 0L) {
        gradient.order <- npConditionalGradientOrder(
          bws = bws,
          reg.engine = reg.spec$reg.engine,
          gradient.order = gradient.order,
          where = "plot.condbandwidth()"
        )
        available <- npGlpGradientAvailability(
          regtype.engine = reg.spec$reg.engine,
          degree.engine = reg.spec$degree.engine,
          gradient.order = gradient.order,
          ncon = bws$xncon
        )
        plot.gradient.available[which(bws$ixcon)] <- available
        if (!any(plot.gradient.available)) {
          npStopGlpGradientNoneAvailable(
            where = "plot.condbandwidth()",
            action = "plot",
            degree.engine = reg.spec$degree.engine,
            gradient.order = gradient.order,
            available = available,
            con.names = bws$xnames[bws$ixcon]
          )
        }
        if (any(!available)) {
          npWarnGlpGradientPartialAvailability(
            where = "plot.condbandwidth()",
            degree.engine = reg.spec$degree.engine,
            gradient.order = gradient.order,
            available = available,
            con.names = bws$xnames[bws$ixcon]
          )
        }
      }
    }

    engine.ctx <- .np_plot_engine_begin(plot.par.mfrow = plot.par.mfrow)
    on.exit(.np_plot_restore_par(engine.ctx$oldpar), add = TRUE)
    plot.par.mfrow <- engine.ctx$plot.par.mfrow
    scalar_default <- .np_plot_scalar_default

    dots <- list(...)
    plot.legend <- legend
    plot.user.args <- .np_plot_user_args(dots, "plot")
    bxp.user.args <- .np_plot_user_args(dots, "bxp")
    rgl.persp3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.persp3d", "rgl.persp3d.")
    rgl.view3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.view3d", "rgl.view3d.")
    rgl.par3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.par3d", "rgl.par3d.")
    rgl.grid3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.grid3d", "rgl.grid3d.")
    rgl.widget.user.args <- .np_plot_collect_rgl_args(dots, "rgl.widget", "rgl.widget.")
    rgl.legend3d.user.args <- .np_plot_collect_rgl_args(dots, "rgl.legend3d", "rgl.legend3d.")
    rgl.legend3d.user.args <- .np_plot_merge_rgl_legend_control(rgl.legend3d.user.args, plot.legend)
    rgl.surface3d.user.args <- .np_plot_extract_prefixed_args(dots, "rgl.surface3d.")
    bxp.args <- bxp.user.args
    if (!is.null(col)) bxp.args$col <- col
    if (!is.null(lty)) bxp.args$lty <- lty
    if (!is.null(lwd)) bxp.args$lwd <- lwd
    if (!is.null(border)) bxp.args$border <- border

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
        "errors must be set to \"bootstrap\" when bootstrap controls are supplied",
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

    if (quantreg) {
      tau <- .npqreg_validate_tau(tau)
    }
    plot.data.overlay <- .np_plot_match_flag(plot.data.overlay, "plot.data.overlay")
    multi.tau <- isTRUE(quantreg) && length(tau) > 1L
    if (multi.tau) {
      common.scale <- FALSE
      perspective <- FALSE
    }
    if (isTRUE(quantreg) && isTRUE(plot.data.overlay) && !isTRUE(gradients))
      common.scale <- FALSE

    plot.errors = (plot.errors.method != "none")
    proper.args <- .np_condist_validate_proper_args(
      proper = proper,
      proper.method = proper.method,
      proper.control = proper.control
    )
    surface.supported <- isTRUE((bws$xncon + bws$xnord + bws$yncon + bws$ynord - quantreg == 2) &&
                                (bws$xnuno + bws$ynuno == 0) &&
                                !any(xor(bws$xdati$iord, bws$xdati$inumord)) &&
                                !isTRUE(quantreg))
    renderer <- .np_plot_validate_renderer_request(
      renderer = renderer,
      route = "plot.condbandwidth()",
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
      route = "plot.condbandwidth()",
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
          cdf = cdf,
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

        pc = (.np_plot_center_is_bias_corrected(plot.errors.center))
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
          if (cdf) {
            ret.args$proper.requested <- tobj$proper.requested
            ret.args$proper.applied <- tobj$proper.applied
            ret.args$proper.method <- tobj$proper.method
            ret.args$condist.raw <- tobj$condist.raw
            ret.args$proper.info <- tobj$proper.info
          }
        }
        cd1 <- do.call(ret.fun, ret.args)
        cd1$bias = NA
        cd1$bias.corrected = NA

        if (.np_plot_center_is_bias_corrected(plot.errors.center))
          cd1 <- .np_plot_add_bias_fields(cd1, tcomp, terr[,3])
        
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
          border = scalar_default(border, .np_plot_color("surface_border")),
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
                             border = scalar_default(border, .np_plot_color("surface_border")),
                             ticktype = "detailed",
                             cex.axis = scalar_default(cex.axis, par()$cex.axis),
                             cex.lab = scalar_default(cex.lab, par()$cex.lab),
                             cex.main = scalar_default(cex.main, par()$cex.main),
                             cex.sub = scalar_default(cex.sub, par()$cex.sub),
                             lwd = .np_plot_lwd("surface_border", scalar_default(lwd, par()$lwd)),
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
              border = scalar_default(border, .np_plot_color("context_border")),
              lwd = scalar_default(lwd, par()$lwd)
            )
            if (plot.errors.type == "all" && !is.null(lerr.all) && !is.null(herr.all)) {
              .np_plot_draw_all_band_legend(
                legend = plot.legend,
                x = "topright",
                lty = .np_plot_lty("solid"),
                lwd = .np_plot_lwd("band_all_surface", scalar_default(lwd, par()$lwd))
              )
            }
          }

          rotation.progress <- .np_plot_rotation_progress_tick(rotation.progress, done = frame.idx)
          Sys.sleep(if (isTRUE(rotate)) rotate.defaults$sleep else 0.5)
      }

      if (plot.behavior == "plot-data")
        return ( list(cd1 = cd1) )
    } else {

      quantreg.gradient.diagonal <- isTRUE(quantreg) && isTRUE(gradients)
      dsf <- if (gradients && !quantreg.gradient.diagonal) bws$xndim else 1L
      tot.dim = bws$xndim + bws$yndim - quantreg

      gradient_component_index <- function(slice.index, component.index) {
        if (quantreg.gradient.diagonal) slice.index else component.index
      }

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
      tylabE <- if (quantreg) {
        if (multi.tau) "Conditional quantile" else paste(tau, "quantile")
      } else paste("Conditional", if (cdf) "Distribution" else "Density")
      plotOnEstimate <- (plot.errors.center == "estimate")

      conditional_gradient_axis_label <- function(component) {
        .np_plot_conditional_gradient_axis_label(
          target = tylabE,
          predictor = bws$xnames[component],
          component = component,
          continuous = bws$ixcon,
          gradient.order = gradient.order
        )
      }

      quantile_matrix_extract <- function(obj, jj) {
        if (gradients) {
          if (multi.tau) {
            g <- obj$quantgrad
            matrix(
              g[, jj, , drop = FALSE],
              nrow = dim(g)[1L],
              dimnames = list(NULL, dimnames(g)[[3L]])
            )
          } else {
            obj$quantgrad[, jj]
          }
        } else if (multi.tau) {
          obj$quantile
        } else {
          obj$quantile
        }
      }

      quantile_error_extract <- function(obj, jj) {
        if (gradients) {
          if (multi.tau) {
            gerr <- obj$quantgerr
            matrix(
              gerr[, jj, , drop = FALSE],
              nrow = dim(gerr)[1L],
              dimnames = list(NULL, dimnames(gerr)[[3L]])
            )
          } else {
            obj$quantgerr[, jj]
          }
        } else if (multi.tau) {
          obj$quanterr
        } else {
          obj$quanterr
        }
      }

      multi_tau_plotout_err <- function(err) {
        if (is.null(err))
          return(NULL)
        out <- array(
          NA_real_,
          dim = c(dim(err)[1L], 2L, dim(err)[3L]),
          dimnames = list(NULL, c("lower", "upper"), dimnames(err)[[3L]])
        )
        out[, 1L, ] <- -err[, 1L, ]
        out[, 2L, ] <- err[, 2L, ]
        out
      }

      plot_multi_tau <- function(ei, value, xi.factor, xlab.value, ylab.value,
                                 err = NULL, all.err = NULL) {
        tau.labels <- colnames(value)
        if (is.null(tau.labels))
          tau.labels <- .npqreg_tau_labels(tau)
        .np_plot_quantile_overlay_1d(
          ei = ei,
          value = value,
          xi.factor = xi.factor,
          xlab.value = xlab.value,
          ylab.value = ylab.value,
          err = err,
          all.err = all.err,
          tau.labels = tau.labels,
          tau = tau,
          overlay.x = if (xOrY == "x") xdat[, i] else ydat[, i],
          overlay.y = ydat[, 1L],
          plot.data.overlay = isTRUE(quantreg) && isTRUE(plot.data.overlay) && !isTRUE(gradients),
          gradients = gradients,
          plotOnEstimate = plotOnEstimate,
          plot.errors.type = plot.errors.type,
          plot.errors.style = plot.errors.style,
          plot.errors.bar = plot.errors.bar,
          plot.errors.bar.num = plot.errors.bar.num,
          plot.user.args = plot.user.args,
          legend = plot.legend,
          col = col,
          lty = lty,
          lwd = lwd,
          type = type,
          xlim = xlim,
          ylim = ylim,
          main = main,
          sub = sub,
          cex.axis = cex.axis,
          cex.lab = cex.lab,
          cex.main = cex.main,
          cex.sub = cex.sub
        )
      }

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
            gradient.order = gradient.order,
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
            gradients = gradients,
            gradient.order = gradient.order
          )
        }
        if (cdf && isTRUE(proper.args$proper.requested)) {
          tobj$proper.requested <- TRUE
          tobj$proper.applied <- FALSE
          tobj$proper.method <- proper.args$proper.method
          tobj$proper.info <- .np_condist_make_reason_info(
            reason = "x_slices_not_repeated",
            supported = FALSE
          )
        }

        
        ## if there are gradients then we need to repeat the process for each component

        eval.extract <- function(obj, jj){
          if (gradients) {
            if (quantreg) quantile_matrix_extract(obj, jj) else obj$congrad[,jj]
          } else if (quantreg) {
            quantile_matrix_extract(obj, jj)
          } else if (cdf) {
            obj$condist
          } else {
            obj$condens
          }
        }
        err.extract <- function(obj, jj){
          if (gradients) {
            if (quantreg) quantile_error_extract(obj, jj) else obj$congerr[,jj]
          } else if (quantreg) {
            quantile_error_extract(obj, jj)
          } else {
            obj$conderr
          }
        }

        if (plot.behavior != "plot"){
          plot.out[plot.index] = NA
          plot.out[[plot.index]] = tobj
        }

        for (j in seq_len(dsf)){
          grad.j <- gradient_component_index(i, j)
          component.available <- !isTRUE(gradients) || plot.gradient.available[grad.j]
          temp.boot = list()
          temp.all.err <- NULL
          temp.err <- matrix(data = NA, nrow = maxneval, ncol = 3)
          temp.err.arr <- NULL
          temp.eval <- eval.extract(tobj, grad.j)
          multi.eval <- is.matrix(temp.eval)
          if (multi.eval) {
            temp.dens.mat <- matrix(NA_real_, nrow = maxneval,
                                    ncol = ncol(temp.eval),
                                    dimnames = list(NULL, colnames(temp.eval)))
            temp.dens.mat[seq_len(xi.neval), ] <- temp.eval
            temp.dens[seq_len(xi.neval)] <- temp.eval[, 1L]
          } else {
            temp.dens[seq_len(xi.neval)] <- temp.eval
          }
          
          if (plot.errors && component.available){
            if (plot.errors.method == "asymptotic") {
              terr.j <- err.extract(tobj, grad.j)
              if (multi.eval) {
                ntau.j <- ncol(terr.j)
                temp.err.arr <- array(
                  NA_real_,
                  dim = c(maxneval, 3L, ntau.j),
                  dimnames = list(NULL, c("lower", "upper", "center"), colnames(terr.j))
                )
                temp.all.err <- vector("list", ntau.j)
                names(temp.all.err) <- colnames(terr.j)
                for (kk in seq_len(ntau.j)) {
                  asym.obj <- .np_plot_asymptotic_error_from_se(
                    se = terr.j[, kk],
                    alpha = plot.errors.alpha,
                    band.type = plot.errors.type,
                    m = xi.neval
                  )
                  temp.err.arr[seq_len(xi.neval), 1:2, kk] <- asym.obj$err
                  temp.all.err[[kk]] <- asym.obj$all.err
                }
              } else {
                asym.obj <- .np_plot_asymptotic_error_from_se(
                  se = terr.j,
                  alpha = plot.errors.alpha,
                  band.type = plot.errors.type,
                  m = xi.neval
                )
                temp.err[seq_len(xi.neval),1:2] <- asym.obj$err
                temp.all.err <- asym.obj$all.err
              }
            }
            else if (plot.errors.method == "bootstrap"){
              temp.boot.raw <- compute.bootstrap.errors(
                        xdat = xdat,
                        ydat = ydat,
                        exdat = subcol(exdat,ei,i)[seq_len(xi.neval),, drop = FALSE],
                        eydat = eydat[seq_len(xi.neval),, drop = FALSE],
                        cdf = cdf,
                        quantreg = quantreg,
                        tau = tau,
                        gradients = gradients,
                        gradient.order = gradient.order,
                        gradient.index = grad.j,
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
                          gradient.index = grad.j
                        ),
                        bws = bws)
              if (multi.eval) {
                temp.err.arr <- array(
                  NA_real_,
                  dim = c(maxneval, dim(temp.boot.raw[["boot.err"]])[2:3]),
                  dimnames = dimnames(temp.boot.raw[["boot.err"]])
                )
                temp.err.arr[seq_len(xi.neval), , ] <- temp.boot.raw[["boot.err"]]
                temp.all.err <- temp.boot.raw[["boot.all.err"]]
                temp.boot <- temp.boot.raw[["bxp"]]
              } else {
                temp.err[seq_len(xi.neval),] <- temp.boot.raw[["boot.err"]]
                temp.all.err <- temp.boot.raw[["boot.all.err"]]
                temp.boot <- temp.boot.raw[["bxp"]]
              }
              if (!multi.eval && !plot.bxp.out){
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
              data.err.all[(plot.index-1)*dsf+j] = list(temp.all.err)
            }
          } else if (plot.behavior == "data" && multi.eval && plot.errors) {
            err.name <- if (gradients) paste("gc", grad.j, "err", sep = "") else "quanterr"
            bias.name <- if (gradients) paste("gc", grad.j, "bias", sep = "") else "bias"
            multi.err <- temp.err.arr[seq_len(xi.neval), , , drop = FALSE]
            plot.out[[plot.index]][[err.name]] <- multi_tau_plotout_err(multi.err)
            plot.out[[plot.index]][[bias.name]] <-
              temp.dens.mat[seq_len(xi.neval), , drop = FALSE] -
              multi.err[, 3L, , drop = TRUE]
            plot.out[[plot.index]]$bxp <- temp.boot
          } else if (plot.behavior != "data") {
            plot.layout <- .np_plot_layout_activate(plot.layout)
            if (multi.eval) {
              multi.err <- if (plot.errors) temp.err.arr[seq_len(xi.neval), , , drop = FALSE] else NULL
              plot_multi_tau(
                ei = ei,
                value = temp.dens.mat[seq_len(xi.neval), , drop = FALSE],
                xi.factor = xi.factor,
                xlab.value = scalar_default(xlab, gen.label(if (xOrY == "x") bws$xnames[i] else bws$ynames[i], paste(toupper(xOrY), i, sep = ""))),
                ylab.value = scalar_default(ylab, if (gradients) conditional_gradient_axis_label(grad.j) else tylabE),
                err = multi.err,
                all.err = temp.all.err
              )
              if (plot.rug && !xi.factor)
                .np_plot_draw_rug_1d(if (xOrY == "x") xdat[,i] else ydat[,i])
              if (plot.behavior != "plot" && plot.errors) {
                err.name <- if (gradients) paste("gc", grad.j, "err", sep = "") else "quanterr"
                bias.name <- if (gradients) paste("gc", grad.j, "bias", sep = "") else "bias"
                plot.out[[plot.index]][[err.name]] <- multi_tau_plotout_err(multi.err)
                plot.out[[plot.index]][[bias.name]] <-
                  temp.dens.mat[seq_len(xi.neval), , drop = FALSE] -
                  multi.err[, 3L, , drop = TRUE]
                plot.out[[plot.index]]$bxp <- temp.boot
              }
              next
            }
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
              plot.args$y <- .np_plot_geometry_values(temp.dens)
            panel.ylim <- if (plot.errors) {
              .np_plot_panel_error_range(
                estimate = temp.dens,
                err = temp.err,
                all.err = temp.all.err,
                plot.errors.type = plot.errors.type,
                plotOnEstimate = plotOnEstimate
              )
            } else {
              range(temp.dens, finite = TRUE)
            }
            if (!all(is.finite(panel.ylim)))
              panel.ylim <- c(-1, 1)
            panel.ylim <- .np_plot_resolve_requested_ylim(panel.ylim, ylim)
            plot.args$ylim <- panel.ylim
            overlay.ok <- isTRUE(quantreg) && isTRUE(plot.data.overlay) && !isTRUE(gradients)
            if (overlay.ok && is.null(ylim))
              plot.args$ylim <- .np_plot_overlay_range(if (is.null(plot.args$ylim)) range(temp.dens, finite = TRUE) else plot.args$ylim,
                                                       ydat[, 1L])
            plot.args$xlab <- scalar_default(xlab, gen.label(if (xOrY == "x") bws$xnames[i] else bws$ynames[i], paste(toupper(xOrY), i, sep = "")))
            plot.args$ylab <- scalar_default(ylab, if (gradients) conditional_gradient_axis_label(grad.j) else tylabE)
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
            if (overlay.ok) {
              overlay.x <- if (xOrY == "x") xdat[, i] else ydat[, i]
              overlay.y <- ydat[, 1L]
              if (xi.factor) {
                data.args <- plot.args
                data.args$f <- NULL
                data.args$z <- NULL
                data.args$x <- overlay.x
                data.args$y <- overlay.y
                data.args$type <- NULL
                data.args$lty <- NULL
                data.args$lwd <- NULL
                data.args$col <- NULL
                do.call(graphics::plot, data.args)
              } else {
                data.args <- plot.args
                data.args$x <- overlay.x
                data.args$y <- overlay.y
                data.args$type <- "p"
                data.args$lty <- NULL
                data.args$lwd <- NULL
                data.args$pch <- scalar_default(data.args$pch, 20)
                data.args$cex <- scalar_default(data.args$cex, 0.5)
                data.args$col <- .np_plot_color("data_overlay")
                do.call(graphics::plot, data.args)
              }
            } else {
              do.call(plot.fun, plot.args)
            }
            if (plot.rug && !xi.factor)
              .np_plot_draw_rug_1d(if (xOrY == "x") xdat[,i] else ydat[,i])

            ## error plotting evaluation
            if (plot.errors && component.available &&
                !(xi.factor && plot.bootstrap && plot.bxp)){
              if (plot.errors.type == "all") {
                draw.all.error.types(
                  ex = as.numeric(ei),
                  center = .np_plot_error_center(temp.dens, temp.err, plotOnEstimate),
                  all.err = temp.all.err,
                  xi.factor = xi.factor,
                  legend = plot.legend)
                .np_plot_draw_bias_center_1d(
                  x = ei,
                  center = .np_plot_error_center(temp.dens, temp.err, plotOnEstimate),
                  xi.factor = xi.factor,
                  plotOnEstimate = plotOnEstimate,
                  legend = plot.legend,
                  estimate.col = plot.args$col,
                  estimate.lty = plot.args$lty,
                  estimate.lwd = plot.args$lwd,
                  draw.legend = FALSE)
              } else {
                .np_plot_draw_bias_center_1d(
                  x = ei,
                  center = .np_plot_error_center(temp.dens, temp.err, plotOnEstimate),
                  xi.factor = xi.factor,
                  plotOnEstimate = plotOnEstimate,
                  legend = plot.legend,
                  estimate.col = plot.args$col,
                  estimate.lty = plot.args$lty,
                  estimate.lwd = plot.args$lwd)
                draw.vec <- .np_plot_error_draw_vectors(ei, temp.dens, temp.err, plotOnEstimate)
                draw.args <- list(
                  ex = draw.vec$x,
                  ely = draw.vec$lower,
                  ehy = draw.vec$upper,
                  plot.errors.style = if (xi.factor) "bar" else plot.errors.style,
                  plot.errors.bar = if (xi.factor) "I" else plot.errors.bar,
                  plot.errors.bar.num = plot.errors.bar.num,
                  lty = if (xi.factor) 1 else 2
                )
                do.call(draw.errors, draw.args)
              }

            }
            if (overlay.ok) {
              line.x <- if (xi.factor) seq_len(length(ei)) else ei
              lines(line.x, temp.dens, type = scalar_default(type, "l"),
                    col = scalar_default(col, par()$col),
                    lty = scalar_default(lty, par()$lty),
                    lwd = scalar_default(lwd, par()$lwd))
            }
          }
        
          if (plot.behavior != "plot" && plot.errors && !multi.eval) {
            err.name <- if (gradients) paste("gc", grad.j, "err", sep = "") else if (quantreg) "quanterr" else "conderr"
            bias.name <- if (gradients) paste("gc", grad.j, "bias", sep = "") else "bias"
            corrected.name <- if (gradients) paste("gc", grad.j, "bias.corrected", sep = "") else "bias.corrected"
            plot.out[[plot.index]][[err.name]] <- if (component.available) {
              na.omit(cbind(-temp.err[, 1], temp.err[, 2]))
            } else {
              cbind(-temp.err[seq_len(xi.neval), 1],
                    temp.err[seq_len(xi.neval), 2])
            }
            plot.out[[plot.index]] <- .np_plot_add_named_bias_fields(
              plot.out[[plot.index]],
              estimate = temp.dens,
              bias.corrected = temp.err[,3],
              bias.name = bias.name,
              corrected.name = corrected.name)
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
              gradient.order = gradient.order,
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
            component.available <- !isTRUE(gradients) || plot.gradient.available[j]
            temp.err[,] <- NA_real_
            temp.boot = list()
            temp.all.err <- NULL
            temp.dens[seq_len(xi.neval)] <- eval.extract(tobj, j)
            
            if (plot.errors && component.available){
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
                          gradient.order = gradient.order,
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
                            gradient.index = j
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
                data.err.all[(plot.index-1)*dsf+j] = list(temp.all.err)
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
              panel.ylim <- if (plot.errors) {
                .np_plot_panel_error_range(
                  estimate = temp.dens,
                  err = temp.err,
                  all.err = temp.all.err,
                  plot.errors.type = plot.errors.type,
                  plotOnEstimate = plotOnEstimate
                )
              } else {
                range(temp.dens, finite = TRUE)
              }
              if (!all(is.finite(panel.ylim)))
                panel.ylim <- c(-1, 1)
              panel.ylim <- .np_plot_resolve_requested_ylim(panel.ylim, ylim)
              plot.args$ylim <- panel.ylim
              plot.args$xlab <- scalar_default(xlab, gen.label(if (xOrY == "x") bws$xnames[i] else bws$ynames[i], paste(toupper(xOrY), i, sep = "")))
              plot.args$ylab <- scalar_default(ylab, if (gradients) conditional_gradient_axis_label(j) else tylabE)
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
              if (plot.errors && component.available &&
                  !(xi.factor && plot.bootstrap && plot.bxp)){
                if (plot.errors.type == "all") {
                  draw.all.error.types(
                    ex = as.numeric(ei),
                    center = .np_plot_error_center(temp.dens, temp.err, plotOnEstimate),
                    all.err = temp.all.err,
                    xi.factor = xi.factor,
                  legend = plot.legend)
                  .np_plot_draw_bias_center_1d(
                    x = ei,
                    center = .np_plot_error_center(temp.dens, temp.err, plotOnEstimate),
                    xi.factor = xi.factor,
                    plotOnEstimate = plotOnEstimate,
                    legend = plot.legend,
                    estimate.col = plot.args$col,
                    estimate.lty = plot.args$lty,
                    estimate.lwd = plot.args$lwd,
                    draw.legend = FALSE)
                } else {
                  .np_plot_draw_bias_center_1d(
                    x = ei,
                    center = .np_plot_error_center(temp.dens, temp.err, plotOnEstimate),
                    xi.factor = xi.factor,
                    plotOnEstimate = plotOnEstimate,
                    legend = plot.legend,
                    estimate.col = plot.args$col,
                    estimate.lty = plot.args$lty,
                    estimate.lwd = plot.args$lwd)
                  draw.vec <- .np_plot_error_draw_vectors(ei, temp.dens, temp.err, plotOnEstimate)
                  draw.args <- list(
                    ex = draw.vec$x,
                    ely = draw.vec$lower,
                    ehy = draw.vec$upper,
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
              corrected.name <- if (gradients) paste("gc", j, "bias.corrected", sep = "") else "bias.corrected"
              plot.out[[plot.index]][[err.name]] <- if (component.available) {
                na.omit(cbind(-temp.err[, 1], temp.err[, 2]))
              } else {
                cbind(-temp.err[seq_len(xi.neval), 1],
                      temp.err[seq_len(xi.neval), 2])
              }
              plot.out[[plot.index]] <- .np_plot_add_named_bias_fields(
                plot.out[[plot.index]],
                estimate = temp.dens,
                bias.corrected = temp.err[,3],
                bias.name = bias.name,
                corrected.name = corrected.name)
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
            err.block <- data.err[, seq(3*k - 2L, length.out = 3L), drop = FALSE]
            range.k <- .np_plot_panel_error_range(
              estimate = data.eval[, k],
              err = err.block,
              all.err = data.err.all[[k]],
              plot.errors.type = plot.errors.type,
              plotOnEstimate = plotOnEstimate
            )
            y.min <- min(y.min, range.k[1], na.rm = TRUE)
            y.max <- max(y.max, range.k[2], na.rm = TRUE)
          }
          if (!is.finite(y.min) || !is.finite(y.max)) {
            y.range <- .np_plot_finite_range(data.eval)
            y.min <- y.range[1L]
            y.max <- y.range[2L]
          }
        } else {
          y.min <- Inf
          y.max <- -Inf
          for (k in seq_len(tot.dim*dsf)) {
            err.block <- if (plot.errors)
              data.err[, seq(3*k - 2L, length.out = 3L), drop = FALSE]
            else
              NULL
            range.k <- .np_plot_panel_error_range(
              estimate = data.eval[, k],
              err = err.block,
              plot.errors.type = plot.errors.type,
              plotOnEstimate = plotOnEstimate
            )
            y.min <- min(y.min, range.k[1L], na.rm = TRUE)
            y.max <- max(y.max, range.k[2L], na.rm = TRUE)
          }
          if (!is.finite(y.min) || !is.finite(y.max)) {
            y.range <- .np_plot_finite_range(data.eval)
            y.min <- y.range[1L]
            y.max <- y.range[2L]
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
            grad.j <- if (plot.index <= bws$xndim) gradient_component_index(i, j) else j
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
              plot.args$y <- .np_plot_geometry_values(data.eval[,idx])
            plot.args$ylim <- c(y.min, y.max)
            plot.args$xlab <- scalar_default(xlab, gen.label(if (xOrY == "x") bws$xnames[i] else bws$ynames[i], paste(toupper(xOrY), i, sep = "")))
            plot.args$ylab <- scalar_default(ylab, if (gradients) conditional_gradient_axis_label(grad.j) else tylabE)
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
              err.block <- data.err[, seq(3*idx - 2L, length.out = 3L), drop = FALSE]
              if (plot.errors.type == "all") {
                draw.all.error.types(
                  ex = as.numeric(allei[,plot.index]),
                  center = .np_plot_error_center(data.eval[,idx], err.block, plotOnEstimate),
                  all.err = data.err.all[[idx]],
                  xi.factor = xi.factor,
                  legend = plot.legend)
                .np_plot_draw_bias_center_1d(
                  x = allei[,plot.index],
                  center = .np_plot_error_center(data.eval[,idx], err.block, plotOnEstimate),
                  xi.factor = xi.factor,
                  plotOnEstimate = plotOnEstimate,
                  legend = plot.legend,
                  estimate.col = plot.args$col,
                  estimate.lty = plot.args$lty,
                  estimate.lwd = plot.args$lwd,
                  draw.legend = FALSE)
              } else {
                .np_plot_draw_bias_center_1d(
                  x = allei[,plot.index],
                  center = .np_plot_error_center(data.eval[,idx], err.block, plotOnEstimate),
                  xi.factor = xi.factor,
                  plotOnEstimate = plotOnEstimate,
                  legend = plot.legend,
                  estimate.col = plot.args$col,
                  estimate.lty = plot.args$lty,
                  estimate.lwd = plot.args$lwd)
                draw.vec <- .np_plot_error_draw_vectors(allei[,plot.index], data.eval[,idx], err.block, plotOnEstimate)
                draw.args <- list(
                  ex = draw.vec$x,
                  ely = draw.vec$lower,
                  ehy = draw.vec$upper,
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
