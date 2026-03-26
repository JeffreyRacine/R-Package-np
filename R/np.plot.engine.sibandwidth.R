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
           plot.errors.center = c("estimate","bias-corrected"),
           plot.errors.type = c("pmzsd","pointwise","bonferroni","simultaneous","all"),
           plot.errors.alpha = 0.05,
           plot.errors.style = c("band","bar"),
           plot.errors.bar = c("|","I"),
           plot.errors.bar.num = NULL,
           plot.par.mfrow = TRUE,
           plot.rug = FALSE,
           ...,
           random.seed){

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

    scalar_default <- .np_plot_scalar_default

    dots <- list(...)
    plot.user.args <- .np_plot_user_args(dots, "plot")
    bxp.user.args <- .np_plot_user_args(dots, "bxp")
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

    if (!gradients){
      if(!is.null(ylim)){
        ymin = ylim[1]
        ymax = ylim[2]
      } else {
        if (plot.errors && plot.errors.type == "all") {
          yr <- compute.all.error.range(if (plot.errors.center == "estimate") temp.mean else temp.err[,3],
                                        temp.all.err)
          ymin <- yr[1]
          ymax <- yr[2]
        } else {
          if (plot.errors){
            ymin <- min(na.omit(c(temp.mean - temp.err[,1],
                                  temp.err[,3] - temp.err[,1])))
            ymax <- max(na.omit(c(temp.mean + temp.err[,2],
                                  temp.err[,3] + temp.err[,2])))
          } else {
            ymin <- min(c(temp.mean, temp.err[,3]))
            ymax <- max(c(temp.mean, temp.err[,3]))
          }
        }
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
          if (plot.rug)
            .np_plot_draw_rug_1d(tobj$index)
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
              lty = 2)
          } else if (plot.errors.center == "estimate") {
            draw.errors(ex = na.omit(tobj$index[i.sort]),
                        ely = na.omit(temp.mean[i.sort] - temp.err[i.sort,1]),
                        ehy = na.omit(temp.mean[i.sort] + temp.err[i.sort,2]),
                        plot.errors.style = plot.errors.style,
                        plot.errors.bar = plot.errors.bar,
                        plot.errors.bar.num = plot.errors.bar.num,
                        lty = 2)
          } else if (plot.errors.center == "bias-corrected") {
            lines(na.omit(tobj$index[i.sort]), na.omit(temp.err[i.sort,3]), lty = 3)
            draw.errors(ex = na.omit(tobj$index[i.sort]),
                        ely = na.omit(temp.err[i.sort,3] - temp.err[i.sort,1]),
                        ehy = na.omit(temp.err[i.sort,3] + temp.err[i.sort,2]),
                        plot.errors.style  = plot.errors.style,
                        plot.errors.bar = plot.errors.bar,
                        plot.errors.bar.num = plot.errors.bar.num,
                        lty = 2)
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
                            ylim = ylim,
                            lwd = scalar_default(lwd, par()$lwd))
          plot.args <- .np_plot_merge_user_args(plot.args, plot.user.args)
          .np_plot_first_render_begin(first.render)
          do.call(plot, plot.args)
          .np_plot_first_render_end(first.render)
          if (plot.rug)
            .np_plot_draw_rug_1d(tobj$index)
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
          if (plot.errors.center == "bias-corrected") {
            plot.out[[1]]$bias <- temp.err[,3] - temp.mean
          }
        }
      }

    } else {

        bmax = max(bws$beta)
        bmin = min(bws$beta)

        if (is.null(ylim)) {
          if (plot.errors && plot.errors.type == "all"){
            yr <- compute.all.error.range(if (plot.errors.center == "estimate") temp.mean else temp.err[,3],
                                          temp.all.err)
            ymax = yr[2]
            ymin = yr[1]
          } else if (plot.errors){
            ymax = max(temp.mean + temp.err[,2], na.rm = TRUE)
            ymin = min(temp.mean - temp.err[,1], na.rm = TRUE)
          } else {
            ymax = max(temp.mean, na.rm = TRUE)
            ymin = min(temp.mean, na.rm = TRUE)
          }
          common.ylim <- c(min(bmin*ymax,bmax*ymin),max(bmax*ymax,bmin*ymin))
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
          plot.out[[1]]$gbias = matrix(data=0,nrow = maxneval, ncol = ncol(xdat))
          
        }


        for (i in seq_len(ncol(xdat))) {
          if (plot.behavior != "data"){
            plot.layout <- .np_plot_layout_activate(plot.layout)

            if (is.null(ylim)) {
              if (!common.scale) {
                center.base.i <- if (plot.errors.center == "estimate")
                  temp.mean[i.sort] else temp.err[i.sort,3]
                center.i <- bws$beta[i] * center.base.i
                if (plot.errors && plot.errors.type == "all" && !is.null(temp.all.err)) {
                  sorted.all.err <- lapply(temp.all.err, function(err) {
                    if (is.null(err)) return(NULL)
                    abs(bws$beta[i]) * err[i.sort,,drop = FALSE]
                  })
                  panel.range <- compute.all.error.range(center.i, sorted.all.err)
                } else if (plot.errors) {
                  lo.i <- bws$beta[i] * (center.base.i - temp.err[i.sort,1])
                  hi.i <- bws$beta[i] * (center.base.i + temp.err[i.sort,2])
                  panel.range <- c(min(pmin(lo.i, hi.i), na.rm = TRUE),
                                   max(pmax(lo.i, hi.i), na.rm = TRUE))
                } else {
                  panel.range <- c(min(center.i, na.rm = TRUE), max(center.i, na.rm = TRUE))
                }
                panel.ylim <- panel.range
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
              .np_plot_draw_rug_1d(tobj$index)
            
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
                  lty = 2)
              } else if (plot.errors.center == "estimate") {
                lo.i <- bws$beta[i] * (temp.mean[i.sort] - temp.err[i.sort,1])
                hi.i <- bws$beta[i] * (temp.mean[i.sort] + temp.err[i.sort,2])
                draw.errors(ex = na.omit(tobj$index[i.sort]),
                            ely = na.omit(pmin(lo.i, hi.i)),
                            ehy = na.omit(pmax(lo.i, hi.i)),
                            plot.errors.style = plot.errors.style,
                            plot.errors.bar = plot.errors.bar,
                            plot.errors.bar.num = plot.errors.bar.num,
                            lty = 2)
              } else if (plot.errors.center == "bias-corrected") {
                lo.i <- bws$beta[i] * (temp.err[i.sort,3] - temp.err[i.sort,1])
                hi.i <- bws$beta[i] * (temp.err[i.sort,3] + temp.err[i.sort,2])
                lines(na.omit(tobj$index[i.sort]), na.omit(bws$beta[i] * temp.err[i.sort,3]), lty = 3)
                draw.errors(ex = na.omit(tobj$index[i.sort]),
                            ely = na.omit(pmin(lo.i, hi.i)),
                            ehy = na.omit(pmax(lo.i, hi.i)),
                            plot.errors.style  = plot.errors.style,
                            plot.errors.bar = plot.errors.bar,
                            plot.errors.bar.num = plot.errors.bar.num,
                            lty = 2)
              }
            }
          }

          if (plot.behavior != "plot"){
            plot.out[[1]]$grad[,i] = bws$beta[i]*temp.mean[i.sort]
            center.out.i <- if (plot.errors.center == "estimate")
              temp.mean[i.sort] else temp.err[i.sort,3]
            lo.i <- bws$beta[i] * (center.out.i - temp.err[i.sort,1])
            hi.i <- bws$beta[i] * (center.out.i + temp.err[i.sort,2])
            plot.out[[1]]$glerr[,i] = pmin(lo.i, hi.i)
            plot.out[[1]]$gherr[,i] = pmax(lo.i, hi.i)
            plot.out[[1]]$gbias[,i] = bws$beta[i]*temp.err[i.sort,3]
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
