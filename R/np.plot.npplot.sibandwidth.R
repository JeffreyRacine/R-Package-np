npplot.sibandwidth <-
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
           plot.errors.method = c("none","bootstrap","asymptotic"),
           plot.errors.boot.num = 399,
           plot.errors.boot.method = c("inid", "fixed", "geom"),
           plot.errors.boot.blocklen = NULL,
           plot.errors.center = c("estimate","bias-corrected"),
           plot.errors.type = c("pmzsd","pointwise","bonferroni","simultaneous","all"),
           plot.errors.alpha = 0.05,
           plot.errors.style = c("band","bar"),
           plot.errors.bar = c("|","I"),
           plot.errors.bar.num = NULL,
           plot.par.mfrow = TRUE,
           ...,
           random.seed){

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    miss.xy = c(missing(xdat),missing(ydat))
    xy <- .npplot_resolve_xydat(bws = bws, xdat = xdat, ydat = ydat, miss.xy = miss.xy)
    xdat <- xy$xdat
    ydat <- xy$ydat


    if (is.null(plot.errors.bar.num))
      plot.errors.bar.num = min(length(ydat),25)

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
      ylim = ylim,
      allow_asymptotic_quantile = FALSE
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

    if (plot.errors.method == "asymptotic") {
      warning(paste("asymptotic errors are not supported with single index regression.\n",
                    "Proceeding without calculating errors"))
      plot.errors.method = "none"
    }

    plot.errors = (plot.errors.method != "none")


    if (plot.behavior != "data" && plot.par.mfrow)
      par(mfrow=if(gradients) n2mfrow(bws$ndim) else c(1,1),cex=par()$cex)

    plot.out = list()

    neval = maxneval = length(ydat)
    
    tobj = npindex(txdat = xdat, tydat = ydat,
      bws = bws, gradients = gradients)
    
    temp.err = matrix(data = NA, nrow = maxneval, ncol = 3)
    temp.mean = replicate(maxneval, NA)
    temp.all.err <- NULL
    

    temp.mean[] = if(gradients) tobj$grad[,1] else tobj$mean

    if (plot.errors){
      if (plot.errors.method == "bootstrap") {
        temp.boot.raw <- compute.bootstrap.errors(
                  xdat = xdat, ydat = ydat,
                  gradients = gradients,
                  plot.errors.boot.method = plot.errors.boot.method,
                  plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                  plot.errors.boot.num = plot.errors.boot.num,
                  plot.errors.center = plot.errors.center,
                  plot.errors.type = plot.errors.type,
                  plot.errors.alpha = plot.errors.alpha,
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
        if (plot.errors){
          plot(tobj$index[i.sort], temp.mean[i.sort],
               ylim = if (!is.null(ylim)) ylim else c(ymin,ymax),
               xlim = xlim,
               cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
               cex.lab =  ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
               cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
               cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
               xlab = ifelse(!is.null(xlab),xlab,"index"),
               ylab = ifelse(!is.null(ylab),ylab,gen.label(bws$ynames, 'Conditional Mean')),
               type = ifelse(!is.null(type),type,'l'),
               lty = ifelse(!is.null(lty),lty,par()$lty),
               col = ifelse(!is.null(col),col,par()$col),
               main = main,
               sub = sub,
               lwd = ifelse(!is.null(lwd),lwd,par()$lwd))
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
          plot(tobj$index[i.sort], temp.mean[i.sort],
               cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
               cex.lab =  ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
               cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
               cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
               xlab = ifelse(!is.null(xlab),xlab,"Index"),
               ylab = ifelse(!is.null(ylab),ylab,gen.label(bws$ynames, 'Conditional Mean')),
               type = ifelse(!is.null(type),type,'l'),
               lty = ifelse(!is.null(lty),lty,par()$lty),
               col = ifelse(!is.null(col),col,par()$col),
               main = main,
               sub = sub,
               xlim = xlim,
               ylim = ylim,
               lwd = ifelse(!is.null(lwd),lwd,par()$lwd))
        }
      }

      if (plot.behavior != "plot") {
        plot.out[1] = NA
        plot.out[[1]] = tobj
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
          plot.out[[1]] = tobj

          plot.out[[1]]$index = tobj$index[i.sort]
          plot.out[[1]]$mean = tobj$mean[i.sort]
          plot.out[[1]]$grad = matrix(data=0,nrow = nrow(xdat), ncol = ncol(xdat))
          plot.out[[1]]$glerr = matrix(data=0,nrow = nrow(xdat), ncol = ncol(xdat))
          plot.out[[1]]$gherr = matrix(data=0,nrow = nrow(xdat), ncol = ncol(xdat))
          plot.out[[1]]$gbias = matrix(data=0,nrow = nrow(xdat), ncol = ncol(xdat))
          
        }


        for (i in 1:ncol(xdat)){
          if (plot.behavior != "data"){

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
            
            plot(tobj$index[i.sort], temp.mean[i.sort]*bws$beta[i],
                 ylim = panel.ylim,
                 cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                 cex.lab =  ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                 cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                 cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                 xlab = ifelse(!is.null(xlab),xlab,"index"),
                 ylab = paste("Gradient Component",i, "of", gen.label(bws$ynames, 'Conditional Mean')),
                 lty = ifelse(!is.null(lty),lty,par()$lty),
                 col = ifelse(!is.null(col),col,par()$col),
                 type = ifelse(!is.null(type),type,'l'),
                 main = main,
                 sub = sub,
                 lwd = ifelse(!is.null(lwd),lwd,par()$lwd))
            
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
      names(plot.out) = paste(ifelse(gradients, "si.grad", "si"),1:length(plot.out),sep="")
      
      return (plot.out)
    }

    
  }
