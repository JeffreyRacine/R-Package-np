.np_plot_sibandwidth_engine <-
  function(bws,
           xdat,
           ydat,
           data = NULL,
           common.scale = TRUE,
           gradients = FALSE,
           main = NULL,
           type = NULL,
           cex.axis = NULL,
           cex.lab = NULL,
           cex.main = NULL,
           cex.sub = NULL,
           col = NULL,
           ylab = NULL,
           xlab = NULL,
           sub = NULL,
           ylim = NULL,
           xlim = NULL,
           lty = NULL,
           lwd = NULL,
           plot.behavior = c("plot","plot-data","data"),
           neval = 50,
           plot.errors.method = c("none","bootstrap","asymptotic"),
           plot.errors.boot.num = 1999,
           plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
           plot.errors.boot.nonfixed = c("exact", "frozen"),
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen = NULL,
           plot.errors.center = c("estimate", "bias-corrected"),
           plot.errors.type = c("pmzsd","pointwise","bonferroni","simultaneous","all"),
           plot.errors.alpha = 0.05,
           plot.errors.style = c("band","bar"),
           plot.errors.bar = c("|","I"),
           plot.errors.bar.num = NULL,
           plot.par.mfrow = TRUE,
           plot.data.overlay = TRUE,
           plot.rug = FALSE,
           ...,
           random.seed){

    dots <- list(...)
    .np_singleindex_reject_higher_gradient_order(
      dots,
      where = "plot.sibandwidth"
    )
    engine.ctx <- .np_plot_engine_begin(plot.par.mfrow = plot.par.mfrow)
    on.exit(.np_plot_restore_par(engine.ctx$oldpar), add = TRUE)
    plot.par.mfrow <- engine.ctx$plot.par.mfrow

    miss.xy = c(missing(xdat),missing(ydat))
    xy <- .np_plot_resolve_xydat(bws = bws, xdat = xdat, ydat = ydat, miss.xy = miss.xy)
    xdat <- xy$xdat
    ydat <- xy$ydat


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

    plot.errors = (plot.errors.method != "none")
    first.render <- .np_plot_first_render_state()
    on.exit(.np_plot_activity_end(first.render$activity), add = TRUE)
    overlay.ok <- .np_plot_overlay_enabled(
      plot.data.overlay = plot.data.overlay,
      plot.behavior = plot.behavior,
      gradients = gradients,
      plot.data.overlay.missing = missing(plot.data.overlay)
    )
    plot.rug <- .np_plot_validate_rug_request(
      plot.rug = plot.rug,
      route = "plot.sibandwidth()",
      supported.route = TRUE,
      renderer = "base",
      reason = "supported base plot routes"
    )

    plot.layout <- .np_plot_layout_begin(
      plot.behavior = plot.behavior,
      plot.par.mfrow = plot.par.mfrow,
      mfrow = if (gradients) n2mfrow(bws$ndim) else c(1, 1)
    )

    plot.out = list()

    scalar_default <- function(value, default) {
      if (is.null(value)) default else value
    }

    plot.legend <- if (!is.null(dots$legend)) dots$legend else TRUE
    plot.user.args <- .np_plot_user_args(dots, "plot")
    points.user.args <- .np_plot_user_args(dots, "points")
    bxp.user.args <- .np_plot_user_args(dots, "bxp")
    if (!is.null(dots$cex) && is.null(points.user.args$cex))
      points.user.args$cex <- dots$cex
    overlay.points.args <- points.user.args
    bxp.args <- bxp.user.args
    if (!is.null(col)) bxp.args$col <- col
    if (!is.null(lty)) bxp.args$lty <- lty
    if (!is.null(lwd)) bxp.args$lwd <- lwd

    make_singleindex_payload <- function(index,
                                         mean,
                                         grad = NA,
                                         gradients.flag = FALSE,
                                         trainiseval = FALSE) {
      do.call(singleindex, list(
        bws = bws,
        index = index,
        mean = mean,
        grad = grad,
        ntrain = nrow(xdat),
        trainiseval = trainiseval,
        gradients = gradients.flag
      ))
    }

    xdat <- toFrame(xdat)
    neval <- .np_plot_validate_neval(neval, where = "plot.sibandwidth")
    eval.info <- .np_plot_singleindex_eval_grid(
      bws = bws,
      xdat = xdat,
      neval = neval,
      trim = 0.0,
      where = "plot.sibandwidth"
    )
    maxneval <- nrow(eval.info$idx.eval)

    if (is.null(plot.errors.bar.num))
      plot.errors.bar.num = min(maxneval, 25L)

    if (identical(plot.errors.method, "asymptotic")) {
      tobj <- .np_plot_singleindex_asymptotic_eval(
        bws = bws,
        txdat = xdat,
        tydat = ydat,
        gradients = gradients,
        index.eval = eval.info$index.eval
      )
    } else {
      tobj <- .np_plot_singleindex_local_eval(
        bws = bws,
        idx.train = eval.info$idx.train,
        idx.eval = eval.info$idx.eval,
        ydat = ydat,
        gradients = gradients
      )
    }
    
    temp.err = matrix(data = NA, nrow = maxneval, ncol = 3)
    temp.mean = rep(NA_real_, maxneval)
    temp.all.err <- NULL
    

    temp.mean[] = if(gradients) tobj$grad.index else tobj$mean

    if (plot.errors){
      if (plot.errors.method == "asymptotic") {
        asym.obj <- .np_plot_asymptotic_error_from_se(
          se = if (gradients) tobj$gerr.index else tobj$merr,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type,
          m = maxneval
        )
        temp.err[,1:2] <- asym.obj$err
        temp.all.err <- asym.obj$all.err
      } else if (plot.errors.method == "bootstrap") {
        temp.boot.raw <- compute.bootstrap.errors(
                  xdat = xdat, ydat = ydat,
                  gradients = gradients,
                  plot.errors.boot.method = plot.errors.boot.method,
                  plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
                  plot.errors.boot.wild = plot.errors.boot.wild,
                  plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                  plot.errors.boot.num = plot.errors.boot.num,
                  plot.errors.center = plot.errors.center,
                  plot.errors.type = plot.errors.type,
                  plot.errors.alpha = plot.errors.alpha,
                  idx.eval = eval.info$idx.eval,
                  bws = bws)
        temp.err[,] <- temp.boot.raw[["boot.err"]]
        temp.all.err <- temp.boot.raw[["boot.all.err"]]
      }
    }

    i.sort = sort(tobj$index, index.return=TRUE)$ix
    train.index <- as.numeric(as.matrix(xdat) %*% as.numeric(bws$beta))

    finite_range <- function(...) {
      vals <- unlist(list(...), use.names = FALSE)
      vals <- vals[is.finite(vals)]
      if (!length(vals))
        return(c(NA_real_, NA_real_))
      range(vals)
    }

    singleindex_level_range <- function() {
      if (plot.errors && plot.errors.type == "all") {
        center <- if (plot.errors.center == "estimate") temp.mean else temp.err[, 3L]
        yr <- compute.all.error.range(center, temp.all.err)
        vals <- c(yr, temp.mean)
        if (.np_plot_center_is_bias_corrected(plot.errors.center))
          vals <- c(vals, temp.err[, 3L])
        out <- finite_range(vals)
      } else if (plot.errors) {
        vals <- c(temp.mean,
                  temp.mean + temp.err[, 2L],
                  temp.mean - temp.err[, 1L])
        if (.np_plot_center_is_bias_corrected(plot.errors.center)) {
          vals <- c(vals,
                    temp.err[, 3L],
                    temp.err[, 3L] + temp.err[, 2L],
                    temp.err[, 3L] - temp.err[, 1L])
        }
        out <- finite_range(vals)
      } else {
        out <- finite_range(temp.mean)
      }
      if (overlay.ok) {
        overlay.range <- .np_plot_overlay_range(out, ydat)
        out <- overlay.range
      }
      out
    }

    singleindex_gradient_component_range <- function(beta.i, idx = seq_along(temp.mean)) {
      estimate.i <- beta.i * temp.mean[idx]
      if (plot.errors && plot.errors.type == "all" && !is.null(temp.all.err)) {
        center.base <- if (plot.errors.center == "estimate") temp.mean[idx] else temp.err[idx, 3L]
        center.i <- beta.i * center.base
        scaled.all.err <- lapply(temp.all.err, function(err) {
          if (is.null(err)) return(NULL)
          abs(beta.i) * err[idx,,drop = FALSE]
        })
        yr <- compute.all.error.range(center.i, scaled.all.err)
        vals <- c(estimate.i, yr)
        if (.np_plot_center_is_bias_corrected(plot.errors.center))
          vals <- c(vals, center.i)
        finite_range(vals)
      } else if (plot.errors) {
        center.base <- if (plot.errors.center == "estimate") temp.mean[idx] else temp.err[idx, 3L]
        lo.i <- beta.i * (center.base - temp.err[idx, 1L])
        hi.i <- beta.i * (center.base + temp.err[idx, 2L])
        vals <- c(estimate.i, pmin(lo.i, hi.i), pmax(lo.i, hi.i))
        if (.np_plot_center_is_bias_corrected(plot.errors.center))
          vals <- c(vals, beta.i * temp.err[idx, 3L])
        finite_range(vals)
      } else {
        finite_range(estimate.i)
      }
    }

    singleindex_gradient_common_range <- function() {
      ranges <- lapply(seq_len(ncol(xdat)), function(j) {
        singleindex_gradient_component_range(bws$beta[j])
      })
      finite_range(ranges)
    }

    if (!gradients){
      if(!is.null(ylim)){
        ymin = ylim[1]
        ymax = ylim[2]
      } else {
        yr <- singleindex_level_range()
        ymin <- yr[1L]
        ymax <- yr[2L]
      }


      if (plot.behavior != "data"){      
        plot.layout <- .np_plot_layout_activate(plot.layout)
        if (plot.errors){
          plot.args <- list(x = tobj$index[i.sort],
                            y = temp.mean[i.sort],
                            ylim = if (!is.null(ylim)) ylim else c(ymin, ymax),
                            xlim = xlim,
                            cex.axis = scalar_default(cex.axis, par()$cex.axis),
                            cex.lab =  scalar_default(cex.lab, par()$cex.lab),
                            cex.main = scalar_default(cex.main, par()$cex.main),
                            cex.sub = scalar_default(cex.sub, par()$cex.sub),
                            xlab = scalar_default(xlab, "index"),
                            ylab = scalar_default(ylab, gen.label(bws$ynames, "Conditional Mean")),
                            type = scalar_default(type, "l"),
                            lty = scalar_default(lty, par()$lty),
                            col = scalar_default(col, par()$col),
                            main = main,
                            sub = sub,
                            lwd = scalar_default(lwd, par()$lwd))
          plot.args <- .np_plot_merge_user_args(plot.args, plot.user.args)
          .np_plot_first_render_begin(first.render)
          do.call(plot, plot.args)
          .np_plot_first_render_end(first.render)
          if (overlay.ok)
            do.call(.np_plot_overlay_points_1d,
                    c(list(x = train.index, y = ydat), overlay.points.args))
          if (plot.rug)
            .np_plot_draw_rug_1d(train.index)
          if (plot.errors.type == "all") {
            sorted.all.err <- lapply(temp.all.err, function(err) {
              if (is.null(err)) return(NULL)
              err[i.sort,,drop = FALSE]
            })
            draw.all.error.types(
              ex = na.omit(tobj$index[i.sort]),
              center = na.omit(if (plot.errors.center == "estimate") temp.mean[i.sort] else temp.err[i.sort,3]),
              all.err = sorted.all.err,
              plot.errors.style = plot.errors.style,
              plot.errors.bar = plot.errors.bar,
              plot.errors.bar.num = plot.errors.bar.num,
              lty = .np_plot_lty("interval"),
              legend = plot.legend)
          } else if (plot.errors.center == "estimate") {
            draw.errors(ex = na.omit(tobj$index[i.sort]),
                        ely = na.omit(temp.mean[i.sort] - temp.err[i.sort,1]),
                        ehy = na.omit(temp.mean[i.sort] + temp.err[i.sort,2]),
                        plot.errors.style = plot.errors.style,
                        plot.errors.bar = plot.errors.bar,
                        plot.errors.bar.num = plot.errors.bar.num,
                        lty = .np_plot_lty("interval"))
          } else if (.np_plot_center_is_bias_corrected(plot.errors.center)) {
            lines(na.omit(tobj$index[i.sort]), na.omit(temp.err[i.sort,3]), lty = .np_plot_lty("center"))
            draw.errors(ex = na.omit(tobj$index[i.sort]),
                        ely = na.omit(temp.err[i.sort,3] - temp.err[i.sort,1]),
                        ehy = na.omit(temp.err[i.sort,3] + temp.err[i.sort,2]),
                        plot.errors.style  = plot.errors.style,
                        plot.errors.bar = plot.errors.bar,
                        plot.errors.bar.num = plot.errors.bar.num,
                        lty = .np_plot_lty("interval"))
            .np_plot_draw_bias_center_legend(
              legend = plot.legend,
              estimate.col = scalar_default(col, par()$col),
              estimate.lty = scalar_default(lty, par()$lty),
              estimate.lwd = scalar_default(lwd, par()$lwd),
              center.col = scalar_default(col, par()$col),
              center.lwd = scalar_default(lwd, par()$lwd)
            )
          }
        } else {
          plot.args <- list(x = tobj$index[i.sort],
                            y = temp.mean[i.sort],
                            cex.axis = scalar_default(cex.axis, par()$cex.axis),
                            cex.lab =  scalar_default(cex.lab, par()$cex.lab),
                            cex.main = scalar_default(cex.main, par()$cex.main),
                            cex.sub = scalar_default(cex.sub, par()$cex.sub),
                            xlab = scalar_default(xlab, "Index"),
                            ylab = scalar_default(ylab, gen.label(bws$ynames, "Conditional Mean")),
                            type = scalar_default(type, "l"),
                            lty = scalar_default(lty, par()$lty),
                            col = scalar_default(col, par()$col),
                            main = main,
                            sub = sub,
                            xlim = xlim,
                            ylim = if (!is.null(ylim)) ylim else if (overlay.ok) c(ymin, ymax) else ylim,
                            lwd = scalar_default(lwd, par()$lwd))
          plot.args <- .np_plot_merge_user_args(plot.args, plot.user.args)
          .np_plot_first_render_begin(first.render)
          do.call(plot, plot.args)
          .np_plot_first_render_end(first.render)
          if (overlay.ok)
            do.call(.np_plot_overlay_points_1d,
                    c(list(x = train.index, y = ydat), overlay.points.args))
          if (plot.rug)
            .np_plot_draw_rug_1d(train.index)
        }
      }

      if (plot.behavior != "plot") {
        plot.out[1] = NA
        plot.out[[1]] <- make_singleindex_payload(
          index = tobj$index,
          mean = tobj$mean,
          gradients.flag = FALSE,
          trainiseval = eval.info$trainiseval
        )
        if (plot.errors) {
          plot.out[[1]]$merr <- cbind(-temp.err[,1], temp.err[,2])
          if (.np_plot_center_is_bias_corrected(plot.errors.center)) {
            plot.out[[1]] <- .np_plot_add_bias_fields(
              object = plot.out[[1]],
              estimate = temp.mean,
              bias.corrected = temp.err[,3]
            )
          }
        }
      }

    } else {

        if (is.null(ylim)) {
          common.ylim <- singleindex_gradient_common_range()
        } else {
          common.ylim <- ylim
        }


        if (plot.behavior != "plot"){
          plot.out[1] = NA
          plot.out[[1]] <- make_singleindex_payload(
            index = tobj$index[i.sort],
            mean = tobj$mean[i.sort],
            grad = matrix(data = 0, nrow = maxneval, ncol = ncol(xdat)),
            gradients.flag = TRUE,
            trainiseval = eval.info$trainiseval
          )
          plot.out[[1]]$glerr = matrix(data=0,nrow = maxneval, ncol = ncol(xdat))
          plot.out[[1]]$gherr = matrix(data=0,nrow = maxneval, ncol = ncol(xdat))
          if (.np_plot_center_is_bias_corrected(plot.errors.center)) {
            plot.out[[1]]$gbias = matrix(data=0,nrow = maxneval, ncol = ncol(xdat))
            plot.out[[1]]$gradient.bias.corrected = matrix(data=0,nrow = maxneval, ncol = ncol(xdat))
          }
          
        }


        for (i in seq_len(ncol(xdat))) {
          if (plot.behavior != "data"){
            plot.layout <- .np_plot_layout_activate(plot.layout)

            if (is.null(ylim)) {
              if (!common.scale) {
                  panel.ylim <- singleindex_gradient_component_range(bws$beta[i], idx = i.sort)
              } else {
                panel.ylim <- common.ylim
              }
            } else {
              panel.ylim <- common.ylim
            }
            
            plot.args <- list(x = tobj$index[i.sort],
                              y = temp.mean[i.sort] * bws$beta[i],
                              ylim = panel.ylim,
                              cex.axis = scalar_default(cex.axis, par()$cex.axis),
                              cex.lab =  scalar_default(cex.lab, par()$cex.lab),
                              cex.main = scalar_default(cex.main, par()$cex.main),
                              cex.sub = scalar_default(cex.sub, par()$cex.sub),
                              xlab = scalar_default(xlab, "index"),
                              ylab = paste("Gradient Component", i, "of", gen.label(bws$ynames, "Conditional Mean")),
                              lty = scalar_default(lty, par()$lty),
                              col = scalar_default(col, par()$col),
                              type = scalar_default(type, "l"),
                              main = main,
                              sub = sub,
                              lwd = scalar_default(lwd, par()$lwd))
            plot.args <- .np_plot_merge_user_args(plot.args, plot.user.args)
            .np_plot_first_render_begin(first.render)
            do.call(plot, plot.args)
            .np_plot_first_render_end(first.render)
            if (plot.rug)
              .np_plot_draw_rug_1d(train.index)
            
            if (plot.errors){
              if (plot.errors.type == "all") {
                scaled.all.err <- lapply(temp.all.err, function(err) {
                  if (is.null(err)) return(NULL)
                  abs(bws$beta[i]) * err[i.sort,,drop = FALSE]
                })
                draw.all.error.types(
                  ex = na.omit(tobj$index[i.sort]),
                  center = na.omit(bws$beta[i] * if (plot.errors.center == "estimate") temp.mean[i.sort] else temp.err[i.sort,3]),
                  all.err = scaled.all.err,
                  plot.errors.style = plot.errors.style,
                  plot.errors.bar = plot.errors.bar,
                  plot.errors.bar.num = plot.errors.bar.num,
                  lty = .np_plot_lty("interval"),
                  legend = plot.legend)
              } else if (plot.errors.center == "estimate") {
                lo.i <- bws$beta[i] * (temp.mean[i.sort] - temp.err[i.sort,1])
                hi.i <- bws$beta[i] * (temp.mean[i.sort] + temp.err[i.sort,2])
                draw.errors(ex = na.omit(tobj$index[i.sort]),
                            ely = na.omit(pmin(lo.i, hi.i)),
                            ehy = na.omit(pmax(lo.i, hi.i)),
                            plot.errors.style = plot.errors.style,
                            plot.errors.bar = plot.errors.bar,
                            plot.errors.bar.num = plot.errors.bar.num,
                            lty = .np_plot_lty("interval"))
              } else if (.np_plot_center_is_bias_corrected(plot.errors.center)) {
                lo.i <- bws$beta[i] * (temp.err[i.sort,3] - temp.err[i.sort,1])
                hi.i <- bws$beta[i] * (temp.err[i.sort,3] + temp.err[i.sort,2])
                lines(na.omit(tobj$index[i.sort]), na.omit(bws$beta[i] * temp.err[i.sort,3]), lty = .np_plot_lty("center"))
                draw.errors(ex = na.omit(tobj$index[i.sort]),
                            ely = na.omit(pmin(lo.i, hi.i)),
                            ehy = na.omit(pmax(lo.i, hi.i)),
                            plot.errors.style  = plot.errors.style,
                            plot.errors.bar = plot.errors.bar,
                            plot.errors.bar.num = plot.errors.bar.num,
                            lty = .np_plot_lty("interval"))
                .np_plot_draw_bias_center_legend(
                  legend = plot.legend,
                  estimate.col = scalar_default(col, par()$col),
                  estimate.lty = scalar_default(lty, par()$lty),
                  estimate.lwd = scalar_default(lwd, par()$lwd),
                  center.col = scalar_default(col, par()$col),
                  center.lwd = scalar_default(lwd, par()$lwd)
                )
              }
            }
          }

          if (plot.behavior != "plot"){
            grad.i <- bws$beta[i]*temp.mean[i.sort]
            plot.out[[1]]$grad[,i] = grad.i
            center.out.i <- if (plot.errors.center == "estimate")
              temp.mean[i.sort] else temp.err[i.sort,3]
            lo.i <- bws$beta[i] * (center.out.i - temp.err[i.sort,1])
            hi.i <- bws$beta[i] * (center.out.i + temp.err[i.sort,2])
            plot.out[[1]]$glerr[,i] = pmin(lo.i, hi.i)
            plot.out[[1]]$gherr[,i] = pmax(lo.i, hi.i)
            if (.np_plot_center_is_bias_corrected(plot.errors.center)) {
              gradient.bias.corrected.i <- bws$beta[i]*temp.err[i.sort,3]
              plot.out[[1]]$gbias[,i] = grad.i - gradient.bias.corrected.i
              plot.out[[1]]$gradient.bias.corrected[,i] = gradient.bias.corrected.i
            }
          }

        }
      }
    
    if (plot.behavior != "data" && plot.par.mfrow)
      par(mfrow=c(1,1),cex=par()$cex)
    
    if (plot.behavior != "plot"){
      names(plot.out) = paste(if (gradients) "si.grad" else "si", seq_along(plot.out), sep = "")
      
      return (plot.out)
    }

    
  }
