## the idea is that you have done bandwidth selection
## you just need to supply training data and the bandwidth
## this tool will help you visualize the result

gen.label = function(label, altlabel){
  paste(ifelse(is.null(label), altlabel, label))
}

gen.tflabel = function(condition, tlabel, flabel){
  paste(ifelse(condition,tlabel,flabel))
}

draw.error.bands = function(ex, ely, ehy, lty = 2, col = par("col")){
  lines(ex,ely,lty=lty,col=col)
  lines(ex,ehy,lty=lty,col=col)
}

draw.error.bars = function(ex, ely, ehy, hbar = TRUE, hbarscale = 0.3, lty = 2, col = par("col")){
  yy = double(3*length(ex))
  jj = 1:length(ex)*3

  yy[jj-2] = ely
  yy[jj-1] = ehy
  yy[jj] = NA
  
  xx = double(3*length(ex))
  xx[jj-2] = ex
  xx[jj-1] = ex
  xx[jj] = NA

  lines(xx,yy,lty=lty,col=col)

  if (hbar){
    ## hbars look silly if they are too wide in relation to their height
    ## this only matters in the limit of few points, since that is when
    ## hbardist may get relatively large

    golden = (1+sqrt(5))/2
    hbardist = abs(max(ex) - min(ex))/length(ex)*hbarscale

    yg = abs(yy[jj-2]-yy[jj-1])/golden
    htest = (hbardist >= yg)
    
    xx[jj-2] = ex - ifelse(htest, yg/2, hbardist/2) 
    xx[jj-1] = ex + ifelse(htest, yg/2, hbardist/2) 
    
    ty = yy[jj-1]
    yy[jj-1] = yy[jj-2]

    lines(xx,yy,col=col)

    yy[jj-2] = ty
    yy[jj-1] = ty

    lines(xx,yy,col=col)
  }
}

draw.errors =
  function(ex, ely, ehy,
           plot.errors.style,
           plot.errors.bar,
           plot.errors.bar.num,
           lty,
           col = par("col")){
    if (plot.errors.style == "bar"){
      ei = seq(1,length(ex),length.out = min(length(ex),plot.errors.bar.num))
      draw.error.bars(ex = ex[ei],
                      ely = ely[ei],
                      ehy = ehy[ei],
                      hbar = (plot.errors.bar == "I"),
                      lty = lty,
                      col = col)
    } else if (plot.errors.style == "band") {
      draw.error.bands(ex = ex,
                       ely = ely,
                       ehy = ehy,
                       lty = lty,
                       col = col)
    }
  }

draw.all.error.types <- function(ex, center, all.err,
                                 plot.errors.style = "band",
                                 plot.errors.bar = "|",
                                 plot.errors.bar.num = min(length(ex), 25),
                                 lty = 2, add.legend = TRUE, legend.loc = "topleft",
                                 xi.factor = FALSE){
  if (is.null(all.err)) return(invisible(NULL))

  if (xi.factor) {
    plot.errors.style <- "bar"
    plot.errors.bar <- "I"
  }

  draw_one <- function(err, col) {
    if (is.null(err)) return(invisible(NULL))
    lower <- center - err[,1]
    upper <- center + err[,2]
    good <- complete.cases(ex, lower, upper)
    if (!any(good)) return(invisible(NULL))
    draw.errors(ex = ex[good], ely = lower[good], ehy = upper[good],
                plot.errors.style = plot.errors.style,
                plot.errors.bar = plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = lty, col = col)
  }

  draw_one(all.err$pointwise, "red")
  draw_one(all.err$simultaneous, "green3")
  draw_one(all.err$bonferroni, "blue")

  if (add.legend) {
    legend(legend.loc,
           legend = c("Pointwise","Simultaneous","Bonferroni"),
           lty = 2, col = c("red","green3","blue"), lwd = 2, bty = "n")
  }
}

plotFactor <- function(f, y, ...){
  dot.expr <- as.list(substitute(list(...)))[-1L]
  dot.names <- names(dot.expr)
  has.user.lty <- !is.null(dot.names) && any(dot.names == "lty")

  if (has.user.lty) {
    plot.call <- as.call(c(list(quote(plot), x = quote(f), y = quote(y)), dot.expr))
    eval(plot.call, parent.frame())
  } else {
    plot(x = f, y = y, lty = "blank", ...)
  }

  l.f = rep(f, each=3)
  l.f[3*(1:length(f))] = NA

  l.y = unlist(lapply(y, function (p) { c(0,p,NA) }))

  lines(x = l.f, y = l.y, lty = 2)
  points(x = f, y = y)
}

.npRmpi_guard_bootstrap_plot_autodispatch <- function(plot.errors.method,
                                                      where = "plot()") {
  if (.npRmpi_autodispatch_active() && identical(plot.errors.method, "bootstrap")) {
    stop(sprintf("%s with plot.errors.method='bootstrap' does not currently support npRmpi.autodispatch; use mpi.bcast.cmd(%s, caller.execute=TRUE)", where, where))
  }
}

## Rank-based simultaneous confidence set helper, vendored from
## MCPAN::SCSrank (MCPAN 1.1-21, GPL-2; Schaarschmidt, Gerhard, Sill).
np.plot.SCSrank <- function(x, conf.level = 0.95, alternative = "two.sided", ...) {
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))

  DataMatrix <- x
  N <- nrow(DataMatrix)
  k <- round(conf.level * N, 0)
  RankDat <- apply(DataMatrix, 2, rank)

  switch(alternative,
    "two.sided" = {
      W1 <- apply(RankDat, 1, max)
      W2 <- N + 1 - apply(RankDat, 1, min)

      Wmat <- cbind(W1, W2)
      w <- apply(Wmat, 1, max)
      tstar <- round(sort(w)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(sortx[N + 1 - tstar], sortx[tstar])
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    },
    "less" = {
      W1 <- apply(RankDat, 1, max)
      tstar <- round(sort(W1)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(-Inf, sortx[tstar])
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    },
    "greater" = {
      W2 <- N + 1 - apply(RankDat, 1, min)
      tstar <- round(sort(W2)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(sortx[N + 1 - tstar], Inf)
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    }
  )

  colnames(SCS) <- c("lower", "upper")

  attr(SCS, which = "k") <- k
  attr(SCS, which = "N") <- N
  OUT <- list(conf.int = SCS, conf.level = conf.level, alternative = alternative)
  return(OUT)
}

compute.bootstrap.quantile.bounds <- function(boot.t, alpha, band.type) {
  neval <- ncol(boot.t)

  if (band.type == "pointwise") {
    probs <- c(alpha / 2.0, 1.0 - alpha / 2.0)
    return(t(apply(boot.t, 2, quantile, probs = probs)))
  }

  if (band.type == "bonferroni") {
    probs <- c(alpha / (2.0 * neval), 1.0 - alpha / (2.0 * neval))
    return(t(apply(boot.t, 2, quantile, probs = probs)))
  }

  if (band.type == "simultaneous") {
    return(np.plot.SCSrank(boot.t, conf.level = 1.0 - alpha)$conf.int)
  }

  if (band.type == "all") {
    return(list(
      pointwise = compute.bootstrap.quantile.bounds(boot.t, alpha, "pointwise"),
      bonferroni = compute.bootstrap.quantile.bounds(boot.t, alpha, "bonferroni"),
      simultaneous = compute.bootstrap.quantile.bounds(boot.t, alpha, "simultaneous")
    ))
  }

  stop("'band.type' must be one of pointwise, bonferroni, simultaneous, all")
}

compute.all.error.range <- function(center, all.err) {
  if (is.null(all.err)) {
    return(c(NA_real_, NA_real_))
  }
  lower <- c(center - all.err$pointwise[,1],
             center - all.err$simultaneous[,1],
             center - all.err$bonferroni[,1])
  upper <- c(center + all.err$pointwise[,2],
             center + all.err$simultaneous[,2],
             center + all.err$bonferroni[,2])
  c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
}

compute.default.error.range <- function(center, err) {
  lower <- c(center - err[,1], err[,3] - err[,1])
  upper <- c(center + err[,2], err[,3] + err[,2])
  c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
}


compute.bootstrap.errors = function(...,bws){
  UseMethod("compute.bootstrap.errors",bws)
}

compute.bootstrap.errors.rbandwidth =
  function(xdat, ydat,
           exdat,
           gradients,
           gradient.order,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .npRmpi_guard_bootstrap_plot_autodispatch("bootstrap", where = "compute.bootstrap.errors(...)")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"

    strf = ifelse(is.inid, "function(data,indices){", "function(tsb){")
    strtx = ifelse(is.inid, "txdat = xdat[indices,],",
      "txdat = tsb[,1:(ncol(tsb)-1),drop=FALSE],")
    strty = ifelse(is.inid, "tydat = ydat[indices],",
      "tydat = tsb[,ncol(tsb)],")
    
    boofun = eval(parse(text=paste(strf, "suppressWarnings(npreg(", strtx, strty,
                          "exdat = exdat, bws = bws,",
                          "gradients = gradients,",
                          "gradient.order = gradient.order,",
                          "warn.glp.gradient = FALSE))$",
                          ifelse(gradients, "grad[,slice.index]", "mean"), "}", sep="")))

    if (is.inid){
      boot.out = boot(data = data.frame(xdat,ydat), statistic = boofun,
        R = plot.errors.boot.num)
    } else {
      boot.out = tsboot(tseries = data.frame(xdat,ydat), statistic = boofun,
        R = plot.errors.boot.num,
        l = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method)
    }

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord | bws$xdati$iuno)[slice.index]){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in 1:length(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
  }

compute.bootstrap.errors.scbandwidth =
  function(xdat, ydat, zdat,
           exdat, ezdat,
           gradients,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .npRmpi_guard_bootstrap_plot_autodispatch("bootstrap", where = "compute.bootstrap.errors(...)")
    miss.z <- missing(zdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"

    xi <- 1:ncol(xdat)
    yi <- ncol(xdat)+1
    if (!miss.z)
      zi <- yi+1:ncol(zdat)

    strf = ifelse(is.inid, "function(data,indices){", "function(tsb){")
    strtx = ifelse(is.inid, "txdat = xdat[indices,, drop = FALSE],",
      "txdat = tsb[,xi,drop=FALSE],")
    strty = ifelse(is.inid, "tydat = ydat[indices],",
      "tydat = tsb[,yi],")
    strtz <-
      ifelse(miss.z, '',
             ifelse(is.inid, 'tzdat = zdat[indices,, drop = FALSE],',
                    'tzdat = tsb[,zi, drop = FALSE],'))
    
    boofun = eval(parse(text=paste(strf, "npscoef(", strtx, strty, strtz,
                          "exdat = exdat,", ifelse(miss.z,"", "ezdat = ezdat,"),
                          "bws = bws, iterate = FALSE)$",
                          "mean", "}", sep="")))


    boot.out <-
      eval(parse(text = paste(ifelse(is.inid, 'boot(data = ',
                   'tsboot(tseries ='), 'data.frame(xdat,ydat',
                   ifelse(miss.z,'', ',zdat'),
                   '), statistic = boofun, R = plot.errors.boot.num',
                   ifelse(is.inid,'','l = plot.errors.boot.blocklen, sim = plot.errors.boot.method'),')')))

    all.bp <- list()

    if ((slice.index > 0) && (((slice.index <= ncol(xdat)) && (bws$xdati$iord | bws$xdati$iuno)[slice.index]) ||
                              ((slice.index > ncol(xdat)) && (bws$zdati$iord | bws$zdati$iuno)[slice.index-ncol(xdat)]))) {
      boot.frame <- as.data.frame(boot.out$t)

      if(slice.index <= ncol(xdat))
          u.lev <- bws$xdati$all.ulev[[slice.index]]
      else
          u.lev <- bws$zdati$all.ulev[[slice.index-ncol(xdat)]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in 1:length(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))

      if(slice.index <= ncol(xdat))
          all.bp$names <- bws$xdati$all.lev[[slice.index]]
      else
          all.bp$names <- bws$zdati$all.lev[[slice.index-ncol(xdat)]]
      rm(boot.frame)
    }
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
  }

compute.bootstrap.errors.plbandwidth =
  function(xdat, ydat, zdat,
           exdat, ezdat,
           gradients,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .npRmpi_guard_bootstrap_plot_autodispatch("bootstrap", where = "compute.bootstrap.errors(...)")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"

    strf = ifelse(is.inid, "function(data,indices){", "function(tsb){")
    strtx = ifelse(is.inid, "txdat = xdat[indices,],",
      "txdat = tsb[,1:ncol(xdat),drop=FALSE],")
    strty = ifelse(is.inid, "tydat = ydat[indices],",
      "tydat = tsb[,ncol(xdat)+1],")
    strtz = ifelse(is.inid, "tzdat = zdat[indices,],",
      "tzdat = tsb[,(ncol(xdat)+2):ncol(tsb), drop=FALSE],")


    boofun = eval(parse(text=paste(strf, "npplreg(", strtx, strty, strtz,
                          "exdat = exdat, ezdat = ezdat, bws = bws)$mean}", sep="")))

    if (is.inid){
      boot.out = boot(data = data.frame(xdat,ydat,zdat), statistic = boofun,
        R = plot.errors.boot.num)
    } else {
      boot.out = tsboot(tseries = data.frame(xdat,ydat,zdat), statistic = boofun,
        R = plot.errors.boot.num,
        l = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method)
    }

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$zdati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord | tdati$iuno)[ti]){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in 1:length(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
  }

compute.bootstrap.errors.bandwidth =
  function(xdat, 
           exdat,
           cdf,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .npRmpi_guard_bootstrap_plot_autodispatch("bootstrap", where = "compute.bootstrap.errors(...)")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"

    strf = ifelse(is.inid, "function(data,indices){", "function(tsb){")
    strt = ifelse(is.inid, "tdat = xdat[indices,],",
      "tdat = tsb,")

    
    boofun = eval(parse(text=paste(strf,
                          ifelse(cdf, "npudist(", "npudens("),
                          strt,
                          "edat = exdat, bws = bws)$",
                          ifelse(cdf, "dist", "dens"),
                          "}", sep="")))

    if (is.inid) {
      boot.out = boot(data = data.frame(xdat), statistic = boofun,
        R = plot.errors.boot.num)
    } else {
      boot.out = tsboot(tseries = data.frame(xdat), statistic = boofun,
        R = plot.errors.boot.num,
        l = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method)
    }

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord | bws$xdati$iuno)[slice.index]){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in 1:length(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
  }

compute.bootstrap.errors.dbandwidth =
  function(xdat, 
           exdat,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .npRmpi_guard_bootstrap_plot_autodispatch("bootstrap", where = "compute.bootstrap.errors(...)")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"

    strf = ifelse(is.inid, "function(data,indices){", "function(tsb){")
    strt = ifelse(is.inid, "tdat = xdat[indices,],",
      "tdat = tsb,")

    
    boofun = eval(parse(text=paste(strf, "npudist(",
                          strt, "edat = exdat, bws = bws)$dist}", sep="")))

    if (is.inid) {
      boot.out = boot(data = data.frame(xdat), statistic = boofun,
        R = plot.errors.boot.num)
    } else {
      boot.out = tsboot(tseries = data.frame(xdat), statistic = boofun,
        R = plot.errors.boot.num,
        l = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method)
    }

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord | bws$xdati$iuno)[slice.index]){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in 1:length(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
  }

compute.bootstrap.errors.conbandwidth =
  function(xdat, ydat,
           exdat, eydat,
           cdf,
           quantreg,
           tau,
           gradients,
           gradient.index,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .npRmpi_guard_bootstrap_plot_autodispatch("bootstrap", where = "compute.bootstrap.errors(...)")
    exdat = toFrame(exdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    tboo =
      if(quantreg) "quant"
      else if (cdf) "dist"
      else "dens"

    is.inid = plot.errors.boot.method=="inid"

    strf = ifelse(is.inid, "function(data,indices){", "function(tsb){")
    strtx = ifelse(is.inid, "txdat = xdat[indices,],",
      "txdat = tsb[,1:ncol(xdat),drop=FALSE],")
    strty = ifelse(is.inid, "tydat = ydat[indices,],",
      "tydat = tsb[,(ncol(xdat)+1):ncol(tsb), drop=FALSE],")
    
    
    boofun = eval(parse(text=paste(strf, 
                          switch(tboo,
                                 "quant" = "npqreg(",
                                 "dist" = "npcdist(",
                                 "dens" = "npcdens("),
                          strtx, strty,
                          "exdat = exdat,",
                          ifelse(quantreg, "tau = tau", "eydat = eydat"),
                          ", bws = bws, gradients = gradients)$",
                          switch(tboo,
                                 "quant" = ifelse(gradients, "yqgrad[,gradient.index]", "quantile"),
                                 "dist" = ifelse(gradients, "congrad[,gradient.index]", "condist"),
                                 "dens" = ifelse(gradients, "congrad[,gradient.index]", "condens")),
                          "}", sep="")))
    if (is.inid){
      boot.out = boot(data = data.frame(xdat,ydat), statistic = boofun,
        R = plot.errors.boot.num)
    } else {
      boot.out = tsboot(tseries = data.frame(xdat,ydat), statistic = boofun,
        R = plot.errors.boot.num,
        l = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method)
    }

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$ydati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord | tdati$iuno)[ti]){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in 1:length(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
  }

compute.bootstrap.errors.condbandwidth =
  function(xdat, ydat,
           exdat, eydat,
           cdf,
           quantreg,
           tau,
           gradients,
           gradient.index,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .npRmpi_guard_bootstrap_plot_autodispatch("bootstrap", where = "compute.bootstrap.errors(...)")
    exdat = toFrame(exdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    tboo =
      if(quantreg) "quant"
      else if (cdf) "dist"
      else "dens"

    is.inid = plot.errors.boot.method=="inid"

    strf = ifelse(is.inid, "function(data,indices){", "function(tsb){")
    strtx = ifelse(is.inid, "txdat = xdat[indices,],",
      "txdat = tsb[,1:ncol(xdat),drop=FALSE],")
    strty = ifelse(is.inid, "tydat = ydat[indices,],",
      "tydat = tsb[,(ncol(xdat)+1):ncol(tsb), drop=FALSE],")
    
    
    boofun = eval(parse(text=paste(strf, 
                          switch(tboo,
                                 "quant" = "npqreg(",
                                 "dist" = "npcdist(",
                                 "dens" = "npcdens("),
                          strtx, strty,
                          "exdat = exdat,",
                          ifelse(quantreg, "tau = tau", "eydat = eydat"),
                          ", bws = bws, gradients = gradients)$",
                          switch(tboo,
                                 "quant" = ifelse(gradients, "yqgrad[,gradient.index]", "quantile"),
                                 "dist" = ifelse(gradients, "congrad[,gradient.index]", "condist"),
                                 "dens" = ifelse(gradients, "congrad[,gradient.index]", "condens")),
                          "}", sep="")))
    if (is.inid){
      boot.out = boot(data = data.frame(xdat,ydat), statistic = boofun,
        R = plot.errors.boot.num)
    } else {
      boot.out = tsboot(tseries = data.frame(xdat,ydat), statistic = boofun,
        R = plot.errors.boot.num,
        l = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method)
    }

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$ydati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord | tdati$iuno)[ti]){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in 1:length(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = all.bp, boot.all.err = boot.all.err)
  }

compute.bootstrap.errors.sibandwidth =
  function(xdat, ydat,
           gradients,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .npRmpi_guard_bootstrap_plot_autodispatch("bootstrap", where = "compute.bootstrap.errors(...)")

    boot.err = matrix(data = NA, nrow = nrow(xdat), ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"

    strf = ifelse(is.inid, "function(data,indices){", "function(tsb){")
    strtx = ifelse(is.inid, "txdat = xdat[indices,],",
      "txdat = tsb[,1:(ncol(tsb)-1),drop=FALSE],")
    strty = ifelse(is.inid, "tydat = ydat[indices],",
      "tydat = tsb[,ncol(tsb)],")
    
    ## beta[1] is always 1.0, so use first column of gradients matrix ... 
    
    boofun = eval(parse(text=paste(strf, "npindex(", strtx, strty,
                          "exdat = xdat, bws = bws,",
                          "gradients = gradients)$",
                          ifelse(gradients, "grad[,1]", "mean"), "}", sep="")))

    if (is.inid){
      boot.out = boot(data = data.frame(xdat,ydat), statistic = boofun,
        R = plot.errors.boot.num)
    } else {
      boot.out = tsboot(tseries = data.frame(xdat,ydat), statistic = boofun,
        R = plot.errors.boot.num,
        l = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method)
    }
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise")
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    list(boot.err = boot.err, bxp = list(), boot.all.err = boot.all.err)
  }


uocquantile <- function(x, prob) {
  if(any(prob < 0 | prob > 1)) stop("'prob' outside [0,1]")
  if(any(is.na(x) | is.nan(x))) stop("missing values and NaN's not allowed")
  if (is.ordered(x)){
    x <- droplevels(x)
    tq = unclass(table(x))
    tq = tq / sum(tq)
    tq[length(tq)] <- 1.0
    bscape <- levels(x)
    tq <- sapply(1:length(tq), function(y){ sum(tq[1:y]) })
    j <- sapply(prob, function(p){ which(tq >= p)[1] })
    bscape[j]
  } else if (is.factor(x)) {
    ## just returns mode
    x <- droplevels(x)
    tq = unclass(table(x))
    j = which(tq == max(tq))[1]
    levels(x)[j]
  } else {
    quantile(x, probs = prob)
  }
}


trim.quantiles <- function(dat, trim){
  if (sign(trim) == sign(-1)){
    trim <- abs(trim)
    tq <- quantile(dat, probs = c(0.0, 0.0+trim, 1.0-trim,1.0))
    tq <- c(2.0*tq[1]-tq[2], 2.0*tq[4]-tq[3])
  }
  else {
    tq <- quantile(dat, probs = c(0.0+trim, 1.0-trim))
  }
  tq
}



npplot <- function(bws = stop("'bws' has not been set"), ..., random.seed = 42){
  .npRmpi_require_active_slave_pool(where = "npplot()")
  ## Save seed prior to setting
  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  ## Restore seed
  on.exit(if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv))

  UseMethod("npplot",bws)
}

npplot.rbandwidth <-
  function(bws,
           xdat,
           ydat,
           data = NULL,
           xq = 0.5,
           xtrim = 0.0,
           neval = 50,
           common.scale = TRUE,
           perspective = TRUE,
           gradients = FALSE,
           gradient.order = 1L,
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
           plot.errors.boot.num = 399,
           plot.errors.boot.method = c("inid", "fixed", "geom"),
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

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    if(!is.null(options('plot.par.mfrow')$plot.par.mfrow))
        plot.par.mfrow <- options('plot.par.mfrow')$plot.par.mfrow
      
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

      ydat <- model.response(tmf)
      xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]
    } else {
      if(all(miss.xy) && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["xdat"]], environment(bws$call)))
        ydat = eval(bws$call[["ydat"]], environment(bws$call))
      }
          
      ## catch and destroy NA's
      xdat = toFrame(xdat)
      goodrows = 1:dim(xdat)[1]
      rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
      goodrows[rows.omit] = 0

      if (all(goodrows==0))
        stop("Data has no rows without NAs")

      xdat = xdat[goodrows,,drop = FALSE]
      ydat = ydat[goodrows]
    }

    ## ydat = as.double(ydat)

    xq = double(bws$ndim)+xq
    xtrim = double(bws$ndim)+xtrim

    if (missing(plot.errors.method) &
        any(!missing(plot.errors.boot.num), !missing(plot.errors.boot.method),
            !missing(plot.errors.boot.blocklen))){
      warning(paste("plot.errors.method must be set to 'bootstrap' to use bootstrapping.",
                    "\nProceeding without bootstrapping."))
    }

    
    plot.behavior = match.arg(plot.behavior)
    plot.errors.method = match.arg(plot.errors.method)
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method, where = "plot(...)")
    plot.errors.boot.method = match.arg(plot.errors.boot.method)
    plot.errors.center = match.arg(plot.errors.center)
    plot.errors.type = match.arg(plot.errors.type)

    if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
        is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
      stop("the tail probability plot.errors.alpha must lie in (0,0.5)")
    plot.errors.style = match.arg(plot.errors.style)
    plot.errors.bar = match.arg(plot.errors.bar)

    common.scale = common.scale | (!is.null(ylim))

    if (plot.errors.method == "none" && plot.errors.type == "all") {
      warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
      plot.errors.method <- "bootstrap"
    }

    if (plot.errors.method == "asymptotic") {
      if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")){
        warning("bootstrap quantile bands cannot be calculated with asymptotics, calculating pmzsd errors")
        plot.errors.type = "pmzsd"
      }

      if (plot.errors.center == "bias-corrected") {
        warning("no bias corrections can be calculated with asymptotics, centering on estimate")
        plot.errors.center = "estimate"
      }
    }

    if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
        is.null(plot.errors.boot.blocklen))
      plot.errors.boot.blocklen = b.star(xdat,round=TRUE)[1,1]    

    plot.errors = (plot.errors.method != "none")
    plot.gradient.order.label <- rep.int(1L, bws$ndim)
    if (gradients && identical(bws$regtype, "glp")) {
      go <- npValidateGlpGradientOrder(regtype = bws$regtype,
                                       gradient.order = gradient.order,
                                       ncon = bws$ncon)
      if (length(go))
        plot.gradient.order.label[which(bws$icon)] <- go
      if (any(go > bws$glp.degree))
        warning("some requested glp derivatives exceed polynomial degree; plotting NA for those components")
      if (any(go > 1L))
        warning("higher-order glp derivatives are not yet available at C level; plotting NA for requested orders > 1")
    }

    if ((bws$ncon + bws$nord == 2) & (bws$nuno == 0) & perspective & !gradients &
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
      
      if (is.ordered(xdat[,1]))
        x1.eval <- (bws$xdati$all.dlev[[1]])[as.integer(x1.eval)]
      
      if (is.ordered(xdat[,2]))
        x2.eval <- (bws$xdati$all.dlev[[2]])[as.integer(x2.eval)]

      tobj = npreg(txdat = xdat, tydat = ydat,
        exdat = x.eval, bws = bws,
        gradient.order = gradient.order,
        warn.glp.gradient = FALSE)

      terr = matrix(data = tobj$merr, nrow = dim(x.eval)[1], ncol = 3)
      terr[,3] = NA
      lerr.all <- NULL
      herr.all <- NULL
      
      treg = matrix(data = tobj$mean,
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      if (plot.errors.method == "bootstrap"){
        terr.obj <- compute.bootstrap.errors(xdat = xdat, ydat = ydat,
          exdat = x.eval,
          gradients = FALSE,
          gradient.order = gradient.order,
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
        center.val <- if(pc) terr[,3] else tobj$mean

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
        lerr = matrix(data = tobj$mean - qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = tobj$mean + qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      }

      if(is.null(zlim)) {
          zlim =
              if (plot.errors){
                  if (plot.errors.type == "all" && !is.null(lerr.all))
                    c(min(c(unlist(lerr.all), lerr)), max(c(unlist(herr.all), herr)))
                  else
                    c(min(lerr),max(herr))
              }
              else
                  c(min(tobj$mean),max(tobj$mean))
      }
        
      if (plot.behavior != "plot"){
        r1 = npregression(bws = bws,
          eval = x.eval,
          mean = as.double(treg),
          merr = terr[,1:2],
          ntrain = dim(xdat)[1])
          r1$bias = NA

        if (plot.errors.center == "bias-corrected")
          r1$bias = terr[,3] - treg

        if (plot.behavior == "data")
          return ( list(r1 = r1) )

      }

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
                treg,
                zlim = zlim,
                col = persp.col,
                border = ifelse(!is.null(border),border,"black"),
                ticktype = "detailed",
                cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                cex.lab =  ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                xlab = ifelse(!is.null(xlab),xlab,gen.label(bws$xnames[1], "X1")),
                ylab = ifelse(!is.null(ylab),ylab,gen.label(bws$xnames[2], "X2")),
                zlab = ifelse(!is.null(zlab),zlab,gen.label(bws$ynames,"Conditional Mean")),
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

      if (plot.behavior == "plot-data")
        return ( list(r1 = r1) )

    } else {

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(bws$ndim),cex=par()$cex)

      ev = xdat[1,,drop = FALSE]

      for (i in 1:bws$ndim)
        ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      maxneval = max(c(sapply(xdat,nlevels),neval))

      exdat = xdat[rep(1, maxneval), , drop = FALSE]

      for (i in 1:bws$ndim)
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
      temp.mean = replicate(maxneval, NA)

      ## plotting expressions
      
      plot.bootstrap = plot.errors.method == "bootstrap"
      
      pfunE = expression(ifelse(xi.factor,
          ifelse(plot.bootstrap & plot.bxp,"bxp","plotFactor"), "plot"))
      pxE = expression(ifelse(common.scale,
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = all.bxp[[i]],", "f = allei[,i],"),
                 "x = allei[,i],"),
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = temp.boot,", "f = ei,"), "x = ei,")))

      pyE = expression(ifelse(xi.factor & plot.bootstrap & plot.bxp, "",
          ifelse(common.scale,"y = data.eval[,i],", "y = temp.mean,")))

      pylimE = ifelse(common.scale, "ylim = c(y.min,y.max),",
        ifelse(plot.errors,
               "ylim = if (plot.errors.type == 'all') compute.all.error.range(if (plotOnEstimate) temp.mean else temp.err[,3], temp.all.err) else c(min(na.omit(c(temp.mean - temp.err[,1], temp.err[,3] - temp.err[,1]))), max(na.omit(c(temp.mean + temp.err[,2], temp.err[,3] + temp.err[,2])))),",
               ""))

      pxlabE = "xlab = ifelse(!is.null(xlab),xlab,gen.label(bws$xnames[i], paste('X', i, sep = ''))),"
      pylabE = "ylab = ifelse(!is.null(ylab),ylab,paste(ifelse(gradients,
          paste('Derivative order ', plot.gradient.order.label[i], ' component ', i, ' of', sep=''), ''),
          gen.label(bws$ynames, 'Conditional Mean'))),"

      prestE = expression(ifelse(xi.factor,"", "type = ifelse(!is.null(type),type,'l'), lty = ifelse(!is.null(lty),lty,par()$lty), col = ifelse(!is.null(col),col,par()$col), lwd = ifelse(!is.null(lwd),lwd,par()$lwd), cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis), cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab), cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main), cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),"))

      pmainE = "main = ifelse(!is.null(main),main,''), sub = ifelse(!is.null(sub),sub,''),"

      ## error plotting expressions
      plotOnEstimate = (plot.errors.center == "estimate")

      efunE = "draw.errors"
      eexE = expression(ifelse(common.scale, "ex = as.numeric(na.omit(allei[,i])),",
          "ex = as.numeric(na.omit(ei)),"))
      eelyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ely = na.omit(data.eval[,i] - data.err[,3*i-2]),",
                 "ely = na.omit(data.err[,3*i] - data.err[,3*i-2]),"),
          ifelse(plotOnEstimate, "ely = na.omit(temp.mean - temp.err[,1]),",
                 "ely = na.omit(temp.err[,3] - temp.err[,1]),")))
      eehyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ehy = na.omit(data.eval[,i] + data.err[,3*i-1]),",
                 "ehy = na.omit(data.err[,3*i] + data.err[,3*i-1]),"),
          ifelse(plotOnEstimate, "ehy = na.omit(temp.mean + temp.err[,2]),",
                 "ehy = na.omit(temp.err[,3] + temp.err[,2]),")))

      erestE = "plot.errors.style = ifelse(xi.factor,'bar',plot.errors.style),
                plot.errors.bar = ifelse(xi.factor,'I',plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(xi.factor,1,2)"


      for (i in 1:bws$ndim){
        temp.err[,] = NA
        temp.mean[] =  NA
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
        
        tr <- if (gradients && identical(bws$regtype, "glp")) {
          suppressWarnings(npreg(txdat = xdat, tydat = ydat,
            exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE], bws = bws,
            gradients = gradients,
            gradient.order = gradient.order,
            warn.glp.gradient = FALSE))
        } else {
          npreg(txdat = xdat, tydat = ydat,
            exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE], bws = bws,
            gradients = gradients,
            gradient.order = gradient.order,
            warn.glp.gradient = FALSE)
        }

        temp.mean[1:xi.neval] = if(gradients) tr$grad[,i] else tr$mean

        if (plot.errors){
          if (plot.errors.method == "asymptotic")
            temp.err[1:xi.neval,1:2] = replicate(2,qnorm(plot.errors.alpha/2, lower.tail = FALSE)*(if(gradients) tr$gerr[,i] else tr$merr))
          else if (plot.errors.method == "bootstrap"){
            temp.boot.raw <- compute.bootstrap.errors(
                      xdat = xdat, ydat = ydat,
                      exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                      gradients = gradients,
                      gradient.order = gradient.order,
                      slice.index = i,
                      plot.errors.boot.method = plot.errors.boot.method,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
                      bws = bws)
            temp.err[1:xi.neval,] = temp.boot.raw[["boot.err"]]
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
          data.eval[,i] = temp.mean
          if (plot.errors){
            all.bxp[i] = NA
            all.bxp[[i]] = temp.boot

            data.err[,c(3*i-2,3*i-1,3*i)] = temp.err
            data.err.all[[i]] = temp.all.err
          }
        } else if (plot.behavior != "data") {
          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = as.numeric(na.omit(if (plotOnEstimate) temp.mean else temp.err[,3])),
                all.err = temp.all.err,
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
            }
          }
                            
        }
        
        if (plot.behavior != "plot") {
          plot.out[i] = NA
          if (gradients){
            plot.out[[i]] =
              npregression(bws = bws,
                           eval = as.data.frame(subcol(exdat,ei,i)[1:xi.neval,]),
                           mean = tr$mean,
                           merr = tr$merr,
                           grad = na.omit(temp.mean),
                           gerr = na.omit(cbind(-temp.err[,1],
                             temp.err[,2])),
                           ntrain = dim(xdat)[1])
            plot.out[[i]]$gbias = na.omit(temp.mean - temp.err[,3])
          } else {
            plot.out[[i]] =
              npregression(bws = bws,
                           eval = as.data.frame(subcol(exdat,ei,i)[1:xi.neval,]),
                           mean = na.omit(temp.mean),
                           merr = na.omit(cbind(-temp.err[,1],
                             temp.err[,2])),
                           ntrain = dim(xdat)[1])
            plot.out[[i]]$bias = na.omit(temp.mean - temp.err[,3])
          }
          plot.out[[i]]$bxp = temp.boot
        }
      }
      
      if (common.scale & (plot.behavior != "data")){
        jj = 1:bws$ndim*3

        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in 1:bws$ndim) {
            if (is.null(data.err.all[[k]])) next
            nkeep.k <- nrow(data.err.all[[k]]$pointwise)
            center.k <- if (plot.errors.center == "estimate")
              data.eval[1:nkeep.k,k]
            else
              data.err[1:nkeep.k,3*k]
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
        } else if (plot.errors.center == "estimate" | !plot.errors) {
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
        
        for (i in 1:bws$ndim){
          xi.factor = is.factor(xdat[,i])

          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(allei[,i]), na.omit(data.err[,3*i]), lty = 3)
            
            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(allei[,i])),
                center = as.numeric(na.omit(if (plotOnEstimate) data.eval[,i] else data.err[,3*i])),
                all.err = data.err.all[[i]],
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
            }
          }
        }
      }

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=c(1,1),cex=par()$cex)
      
      if (plot.behavior != "plot"){
        names(plot.out) =
          if (gradients)
            paste("rg",1:bws$ndim,sep="")
          else
            paste("r",1:bws$ndim,sep="")
        
        return (plot.out)
      }

    }
  }

npplot.scbandwidth <-
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
           plot.errors.boot.num = 399,
           plot.errors.boot.method = c("inid", "fixed", "geom"),
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

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    if(!is.null(options('plot.par.mfrow')$plot.par.mfrow))
        plot.par.mfrow <- options('plot.par.mfrow')$plot.par.mfrow      

    if(!missing(gradients))
      stop("gradients not supported with smooth coefficient models.")

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
      umf <- tmf <- eval(tmf, envir = environment(tt))
      
      ydat <- model.response(tmf)
      xdat <- tmf[, bws$chromoly[[2]], drop = FALSE]
      if (!miss.z)
        zdat <- tmf[, bws$chromoly[[3]], drop = FALSE]
    } else {
      if(all(miss.xy) && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["xdat"]], environment(bws$call)))
        ydat <- eval(bws$call[["ydat"]], environment(bws$call))
        if (!miss.z)
          zdat <- data.frame(eval(bws$call[["zdat"]], environment(bws$call)))
      }
          
      ## catch and destroy NA's
      xdat = toFrame(xdat)
      if(!miss.z)
        zdat <- toFrame(zdat)

      goodrows = 1:dim(xdat)[1]
      rows.omit =
        eval(parse(text = paste("attr(na.omit(data.frame(xdat, ydat",
                     ifelse(miss.z,'',',zdat'),')), "na.action")')))
      
      attr(na.omit(data.frame(xdat,ydat,zdat)), "na.action")
      goodrows[rows.omit] = 0

      if (all(goodrows==0))
        stop("Data has no rows without NAs")

      xdat = xdat[goodrows,,drop = FALSE]

      if(!miss.z)
        zdat <- zdat[goodrows,,drop = FALSE]
      
      ydat = ydat[goodrows]
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
            !missing(plot.errors.boot.blocklen))){
      warning(paste("plot.errors.method must be set to 'bootstrap' to use bootstrapping.",
                    "\nProceeding without bootstrapping."))
    }

    
    plot.behavior = match.arg(plot.behavior)
    plot.errors.method = match.arg(plot.errors.method)
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method, where = "plot(...)")
    plot.errors.boot.method = match.arg(plot.errors.boot.method)
    plot.errors.center = match.arg(plot.errors.center)
    plot.errors.type = match.arg(plot.errors.type)

    if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
        is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
      stop("the tail probability plot.errors.alpha must lie in (0,0.5)")
    plot.errors.style = match.arg(plot.errors.style)
    plot.errors.bar = match.arg(plot.errors.bar)

    common.scale = common.scale | (!is.null(ylim))

    if (plot.errors.method == "none" && plot.errors.type == "all") {
      warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
      plot.errors.method <- "bootstrap"
    }

    if (plot.errors.method == "asymptotic") {
      if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")){
        warning("bootstrap quantile bands cannot be calculated with asymptotics, calculating pmzsd errors")
        plot.errors.type = "pmzsd"
      }

      if (plot.errors.center == "bias-corrected") {
        warning("no bias corrections can be calculated with asymptotics, centering on estimate")
        plot.errors.center = "estimate"
      }
    }

    if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
        is.null(plot.errors.boot.blocklen))
      plot.errors.boot.blocklen = b.star(xdat,round=TRUE)[1,1]

    plot.errors = (plot.errors.method != "none")

    if ((sum(c(bws$xdati$icon, bws$xdati$iord, bws$zdati$icon, bws$zdati$iord))== 2) & (sum(c(bws$xdati$iuno, bws$zdati$iuno)) == 0) & perspective & !gradients &
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

      tobj <- eval(parse(text = paste('npscoef(txdat = xdat, tydat = ydat,',
                           ifelse(miss.z,'','tzdat = zdat,'),
                           ifelse(miss.z,'exdat = x.eval,','exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,2, drop = FALSE],'),
                           'bws = bws, iterate = FALSE, errors = plot.errors)')))

      terr = matrix(data = tobj$merr, nrow = dim(x.eval)[1], ncol = 3)
      terr[,3] = NA
      
      treg = matrix(data = tobj$mean,
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      if (plot.errors.method == "bootstrap"){
        terr <-
          eval(parse(text = paste(
                       'compute.bootstrap.errors(xdat = xdat, ydat = ydat,',
                       ifelse(miss.z,'','zdat = zdat,'),
                       ifelse(miss.z, 'exdat = x.eval,',
                              'exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,1, drop = FALSE],'),
                       ' gradients = FALSE, slice.index = 0,',
                       'plot.errors.boot.method = plot.errors.boot.method,',
                       'plot.errors.boot.blocklen = plot.errors.boot.blocklen,',
                       'plot.errors.boot.num = plot.errors.boot.num,',
                       'plot.errors.center = plot.errors.center,',
                       'plot.errors.type = plot.errors.type,',
                       'plot.errors.alpha = plot.errors.alpha,',
                       'bws = bws)[["boot.err"]]')))

        pc = (plot.errors.center == "bias-corrected")

        lerr = matrix(data = if(pc) {terr[,3]} else {tobj$mean}
          -terr[,1],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = if(pc) {terr[,3]} else {tobj$mean}
          +terr[,2],
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      } else if (plot.errors.method == "asymptotic") {
        lerr = matrix(data = tobj$mean - qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = tobj$mean + qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr,
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
        r1 <-
          eval(parse(text = paste("smoothcoefficient(bws = bws,",
                       ifelse(miss.z, "eval = x.eval,", "eval = list(exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,2, drop = FALSE])"),
                       "mean = as.double(treg),",
                       "merr = terr[,1:2],",
                       "ntrain = dim(xdat)[1])")))
        r1$bias = NA

        if (plot.errors.center == "bias-corrected")
          r1$bias = terr[,3] - treg

        if (plot.behavior == "data")
          return ( list(r1 = r1) )

      }

      dtheta = 5.0
      dphi = 10.0

      persp.col = ifelse(plot.errors, FALSE, ifelse(!is.null(col),col,"lightblue"))
      
      ##for (j in 0:((50 %/% dphi - 1)*rotate)*dphi+phi){
        for (i in 0:((360 %/% dtheta - 1)*rotate)*dtheta+theta){
          if (plot.errors){
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

          persp(x1.eval,
                x2.eval,
                treg,
                zlim = zlim,
                col = persp.col,
                border = ifelse(!is.null(border),border,"black"),
                ticktype = "detailed",
                cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                cex.lab =  ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                xlab = ifelse(!is.null(xlab),xlab,gen.label(bws$xnames[1], "X1")),
                ylab = ifelse(!is.null(ylab),ylab,gen.label(x2.names[1], "X2")),
                zlab = ifelse(!is.null(zlab),zlab,gen.label(bws$ynames,"Conditional Mean")),
                theta = i,
                phi = phi,
                main = gen.tflabel(!is.null(main), main, paste("[theta= ", i,", phi= ", phi,"]", sep="")))

          if (plot.errors){
            par(new = TRUE)
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

          Sys.sleep(0.5)
        }
      ##}

      if (plot.behavior == "plot-data")
        return ( list(r1 = r1) )

    } else {

      tot.dim <- (bws$xndim <- length(bws$xdati$icon)) + (bws$zndim <- length(bws$zdati$icon))

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(tot.dim),cex=par()$cex)

      maxneval = max(c(sapply(xdat,nlevels), unlist(sapply(zdat,nlevels)), neval))
      all.isFactor = c(sapply(xdat, is.factor), unlist(sapply(zdat, is.factor)))
      
      x.ev = xdat[1,,drop = FALSE]

      for (i in 1:bws$xndim)
        x.ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      exdat = xdat[rep(1, maxneval), , drop = FALSE]

      for (i in 1:bws$xndim)
        exdat[,i] = x.ev[1,i]

      if (!miss.z){
        z.ev = zdat[1,,drop = FALSE]
        
        for (i in 1:bws$zndim)
          z.ev[1,i] = uocquantile(zdat[,i], prob=zq[i])

        ezdat = zdat[rep(1, maxneval), , drop = FALSE]

        for (i in 1:bws$zndim)
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
          ifelse(common.scale,"y = data.eval[,plot.index],", "y = temp.mean,")))

      pylimE = ifelse(common.scale, "ylim = c(y.min,y.max),",
        ifelse(plot.errors,
               "ylim = if (plot.errors.type == 'all') compute.all.error.range(if (plotOnEstimate) temp.mean else temp.err[,3], temp.all.err) else c(min(na.omit(c(temp.mean - temp.err[,1], temp.err[,3] - temp.err[,1]))), max(na.omit(c(temp.mean + temp.err[,2], temp.err[,3] + temp.err[,2])))),",
               ""))

      pxlabE = expression(paste("xlab = gen.label(bws$",
          xOrZ, "names[i], paste('", toupper(xOrZ),"', i, sep = '')),",sep=''))

      pylabE = "ylab = paste(ifelse(gradients,
          paste('Gradient Component ', i, ' of', sep=''), ''),
          gen.label(bws$ynames, 'Conditional Mean')),"

      prestE = expression(ifelse(xi.factor,"", "type = ifelse(!is.null(type),type,'l'), lty = ifelse(!is.null(lty),lty,par()$lty), col = ifelse(!is.null(col),col,par()$col), lwd = ifelse(!is.null(lwd),lwd,par()$lwd), cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis), cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab), cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main), cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),"))

      pmainE = "main = ifelse(!is.null(main),main,''), sub = ifelse(!is.null(sub),sub,''),"

      txobjE <-
        parse(text = paste("npscoef(txdat = xdat, tydat = ydat,",
                ifelse(miss.z,"","tzdat = zdat,"),
                "exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],",
                ifelse(miss.z,"","ezdat = ezdat[1:xi.neval,, drop = FALSE],"),
                "bws = bws, errors = plot.errors)"))


      ## error plotting expressions
      plotOnEstimate = (plot.errors.center == "estimate")

      efunE = "draw.errors"
      eexE = expression(ifelse(common.scale, "ex = as.numeric(na.omit(allei[,plot.index])),",
          "ex = as.numeric(na.omit(ei)),"))
      eelyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ely = na.omit(data.eval[,plot.index] - data.err[,3*plot.index-2]),",
                 "ely = na.omit(data.err[,3*plot.index] - data.err[,3*plot.index-2]),"),
          ifelse(plotOnEstimate, "ely = na.omit(temp.mean - temp.err[,1]),",
                 "ely = na.omit(temp.err[,3] - temp.err[,1]),")))
      eehyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ehy = na.omit(data.eval[,plot.index] + data.err[,3*plot.index-1]),",
                 "ehy = na.omit(data.err[,3*plot.index] + data.err[,3*plot.index-1]),"),
          ifelse(plotOnEstimate, "ehy = na.omit(temp.mean + temp.err[,2]),",
                 "ehy = na.omit(temp.err[,3] + temp.err[,2]),")))

      erestE = "plot.errors.style = ifelse(xi.factor,'bar',plot.errors.style),
                plot.errors.bar = ifelse(xi.factor,'I',plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(xi.factor,1,2)"

      plot.index = 0
      xOrZ = "x"

      for (i in 1:bws$xndim){
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

        tobj <- eval(txobjE)

        temp.mean[1:xi.neval] = tobj$mean

        if (plot.errors){
          if (plot.errors.method == "asymptotic")
            temp.err[1:xi.neval,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr
          else if (plot.errors.method == "bootstrap"){
            temp.boot.raw <- eval(parse(text = paste("compute.bootstrap.errors(",
                                      "xdat = xdat, ydat = ydat,",
                                      ifelse(miss.z, "", "zdat = zdat,"),
                                      "exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],",
                                      ifelse(miss.z,"","ezdat = ezdat[1:xi.neval,, drop = FALSE],"),
                                      "gradients = gradients,",
                                      "slice.index = plot.index,",
                                      "plot.errors.boot.method = plot.errors.boot.method,",
                                      "plot.errors.boot.blocklen = plot.errors.boot.blocklen,",
                                      "plot.errors.boot.num = plot.errors.boot.num,",
                                      "plot.errors.center = plot.errors.center,",
                                      "plot.errors.type = plot.errors.type,",
                                      "plot.errors.alpha = plot.errors.alpha,",
                                      "bws = bws)")))
            temp.err[1:xi.neval,] <- temp.boot.raw[["boot.err"]]
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
          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = as.numeric(na.omit(if (plotOnEstimate) temp.mean else temp.err[,3])),
                all.err = temp.all.err,
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
            }

          }
        }

        if (plot.behavior != "plot") {
          plot.out[plot.index] = NA
          if (gradients){
          } else {
            plot.out[[plot.index]] =
              eval(parse(text = paste("smoothcoefficient(bws = bws,",
                           "eval = ", ifelse(miss.z, "subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],",
                                             "list(exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE], ezdat = ezdat[1:xi.neval,, drop = FALSE]),"),
                           "mean = na.omit(temp.mean),",
                           "ntrain = dim(xdat)[1],",
                           "trainiseval = FALSE,",
                           "xtra = c(0, 0, 0, 0, 0, 0))")))
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
        for (i in 1:bws$zndim){
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

          tobj <- npscoef(txdat = xdat, tydat = ydat, tzdat = zdat,
                          exdat = exdat[1:xi.neval,, drop = FALSE],
                          ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
                          bws = bws)

          temp.mean[1:xi.neval] = tobj$mean

          if (plot.errors){
            if (plot.errors.method == "asymptotic")
              temp.err[1:xi.neval,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$merr
          else if (plot.errors.method == "bootstrap"){
              temp.boot.raw <- compute.bootstrap.errors(
                                                    xdat = xdat,
                                                    ydat = ydat,
                                                    zdat = zdat,
                                                    exdat = exdat[1:xi.neval,, drop = FALSE],
                                                    ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
                                                    gradients = gradients,
                                                    slice.index = plot.index,
                                                    plot.errors.boot.method = plot.errors.boot.method,
                                                    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                                                    plot.errors.boot.num = plot.errors.boot.num,
                                                    plot.errors.center = plot.errors.center,
                                                    plot.errors.type = plot.errors.type,
                                                    plot.errors.alpha = plot.errors.alpha,
                                                    bws = bws)
              temp.err[1:xi.neval,] <- temp.boot.raw[["boot.err"]]
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
            ## plot evaluation
            eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                         eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                         eval(pmainE), ")")))

            ## error plotting evaluation
            if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
              if (!xi.factor && !plotOnEstimate)
                lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

              if (plot.errors.type == "all") {
                draw.all.error.types(
                  ex = as.numeric(na.omit(ei)),
                  center = as.numeric(na.omit(if (plotOnEstimate) temp.mean else temp.err[,3])),
                  all.err = temp.all.err,
                  plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                  plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                  plot.errors.bar.num = plot.errors.bar.num,
                  lty = 2,
                  add.legend = TRUE)
              } else {
                eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                             eval(eehyE), eval(erestE), ")")))
              }

            }
          }

          
          if (plot.behavior != "plot") {
            plot.out[plot.index] = NA
            if (gradients){
            } else {
              plot.out[[plot.index]] =
                smoothcoefficient(bws = bws, 
                                  eval = list(exdat = exdat[1:xi.neval,, drop = FALSE],
                                    ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE]),
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
      
      if (common.scale & (plot.behavior != "data")){
        jj = 1:(bws$xndim + bws$zndim)*3
        
        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in 1:(bws$xndim + bws$zndim)) {
            if (is.null(data.err.all[[k]])) next
            nkeep.k <- nrow(data.err.all[[k]]$pointwise)
            center.k <- if (plot.errors.center == "estimate")
              data.eval[1:nkeep.k,k]
            else
              data.err[1:nkeep.k,3*k]
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
        } else if (plot.errors.center == "estimate" | !plot.errors) {
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
        
        for (plot.index in 1:(bws$xndim + bws$zndim)){
          i = ifelse(plot.index <= bws$xndim, plot.index, plot.index - bws$xndim)

          if (plot.index > bws$xndim)
            xOrZ <- "z"
            
          xi.factor = all.isFactor[plot.index]

          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(allei[,plot.index])),
                center = as.numeric(na.omit(if (plotOnEstimate) data.eval[,plot.index] else data.err[,3*plot.index])),
                all.err = data.err.all[[plot.index]],
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
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
            paste("sc",1:(bws$xndim+bws$zndim),sep="")
        
        return (plot.out)
      }
    }

      
  }

npplot.plbandwidth <-
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
           gradients = FALSE,
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
      
    if(!missing(gradients))
      stop("gradients not supported with partially linear models. Coefficients may be extracted with coef()")

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
      umf <- tmf <- eval(tmf, envir = environment(tt))

      bronze <- lapply(bws$chromoly, paste, collapse = " + ")

      tmf.xf[["formula"]] <- as.formula(paste(" ~ ", bronze[[2]]),
                                      env = environment(formula))
      tmf.xf <- eval(tmf.xf,parent.frame())
      
      ydat <- model.response(tmf)
      xdat <- tmf.xf

      zdat <- tmf[, bws$chromoly[[3]], drop = FALSE]
    } else {
      if(all(miss.xyz) && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["xdat"]], environment(bws$call)))
        ydat = eval(bws$call[["ydat"]], environment(bws$call))
        zdat <- data.frame(eval(bws$call[["zdat"]], environment(bws$call)))
      }
      xdat = toFrame(xdat)
      zdat = toFrame(zdat)
      
      ## catch and destroy NA's
      goodrows = 1:dim(xdat)[1]
      rows.omit = attr(na.omit(data.frame(xdat,ydat,zdat)), "na.action")
      goodrows[rows.omit] = 0

      if (all(goodrows==0))
        stop("Training data has no rows without NAs")

      xdat = xdat[goodrows,,drop = FALSE]
      ydat = ydat[goodrows]
      zdat = zdat[goodrows,,drop = FALSE]
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
            !missing(plot.errors.boot.blocklen))){
      warning(paste("plot.errors.method must be set to 'bootstrap' to use bootstrapping.",
                    "\nProceeding without bootstrapping."))
    }
    
    plot.behavior = match.arg(plot.behavior)
    plot.errors.method = match.arg(plot.errors.method)
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method, where = "plot(...)")
    plot.errors.boot.method = match.arg(plot.errors.boot.method)
    plot.errors.center = match.arg(plot.errors.center)
    plot.errors.type = match.arg(plot.errors.type)

    if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
        is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
      stop("the tail probability plot.errors.alpha must lie in (0,0.5)")
    plot.errors.style = match.arg(plot.errors.style)
    plot.errors.bar = match.arg(plot.errors.bar)

    common.scale = common.scale | (!is.null(ylim))

    if (plot.errors.method == "none" && plot.errors.type == "all") {
      warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
      plot.errors.method <- "bootstrap"
    }

    if (plot.errors.method == "asymptotic") {
      warning(paste("asymptotic errors are not supported with partially linear regression.\n",
                    "Proceeding without calculating errors"))
      plot.errors.method = "none"
    }

    if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
        is.null(plot.errors.boot.blocklen))
      plot.errors.boot.blocklen = b.star(xdat,round=TRUE)[1,1]


    plot.errors = (plot.errors.method != "none")

    if ((nxcon + nxord == 1) & (nzcon + nzord == 1) & (nxuno + nzuno == 0) &
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

      tobj = npplreg(txdat = xdat, tydat = ydat, tzdat = zdat,
        exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,2, drop = FALSE], bws = bws)

      terr = matrix(data = NA, nrow = nrow(x.eval), ncol = 3)
      
      treg = matrix(data = tobj$mean,
        nrow = x1.neval, ncol = z1.neval, byrow = FALSE)

      if (plot.errors.method == "bootstrap"){
        terr <- compute.bootstrap.errors(
          xdat = xdat, ydat = ydat, zdat = zdat,
          exdat = x.eval[,1, drop = FALSE], ezdat = x.eval[,2, drop = FALSE],
          gradients = FALSE,
          slice.index = 0,
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
          plot.errors.boot.num = plot.errors.boot.num,
          plot.errors.center = plot.errors.center,
          plot.errors.type = plot.errors.type,
          plot.errors.alpha = plot.errors.alpha,
          bws = bws)[["boot.err"]]

        pc = (plot.errors.center == "bias-corrected")

        lerr = matrix(data = if(pc) {terr[,3]} else {tobj$mean}
          -terr[,1],
          nrow = x1.neval, ncol = z1.neval, byrow = FALSE)

        herr = matrix(data = if(pc) {terr[,3]} else {tobj$mean}
          +terr[,2],
          nrow = x1.neval, ncol = z1.neval, byrow = FALSE)

      }
      
      if(is.null(zlim)) {
          zlim =
              if (plot.errors)
                  c(min(lerr),max(herr))
              else
                  c(min(tobj$mean),max(tobj$mean))
      }
        
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

      dtheta = 5.0
      dphi = 10.0

      persp.col = ifelse(plot.errors, FALSE, ifelse(!is.null(col),col,"lightblue"))
      
      ##for (j in 0:((50 %/% dphi - 1)*rotate)*dphi+phi){
        for (i in 0:((360 %/% dtheta - 1)*rotate)*dtheta+theta){
          if (plot.errors){
            persp(x1.eval,
                  z1.eval,
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

          persp(x1.eval,
                z1.eval,
                treg,
                zlim = zlim,
                col = persp.col,
                border = ifelse(!is.null(border),border,"black"),
                ticktype = "detailed",
                cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis),
                cex.lab =  ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab),
                cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main),
                cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),
                xlab = ifelse(!is.null(xlab),xlab,gen.label(names(xdat)[1], "X1")),
                ylab = ifelse(!is.null(ylab),ylab,gen.label(names(xdat)[2], "Z1")),
                zlab = ifelse(!is.null(zlab),zlab,gen.label(names(ydat),"Conditional Mean")),
                theta = i,
                phi = phi,
                main = gen.tflabel(!is.null(main), main, paste("[theta= ", i,", phi= ", phi,"]", sep="")))

          if (plot.errors){
            par(new = TRUE)
            persp(x1.eval,
                  z1.eval,
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

          Sys.sleep(0.5)
        }
      ##}

      if (plot.behavior == "plot-data")
        return ( list(r1 = r1) )

    } else {
##      stop("not yet supported!")
      
      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(bws$xndim + bws$zndim),cex=par()$cex)

      x.ev = xdat[1,,drop = FALSE]
      z.ev = zdat[1,,drop = FALSE]

      for (i in 1:bws$xndim)
        x.ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      for (i in 1:bws$zndim)
        z.ev[1,i] = uocquantile(zdat[,i], prob=zq[i])


      maxneval = max(c(sapply(xdat,nlevels), sapply(zdat,nlevels), neval))

      ## Preserve original data types (e.g., factors) for evaluation data
      exdat = xdat[rep(1, maxneval), , drop = FALSE]
      ezdat = zdat[rep(1, maxneval), , drop = FALSE]

      for (i in 1:bws$xndim)
        exdat[,i] = x.ev[1,i]

      for (i in 1:bws$zndim)
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

      all.isFactor = c(sapply(xdat, is.factor), sapply(zdat, is.factor))

      plot.out = list()

      temp.err = matrix(data = NA, nrow = maxneval, ncol = 3)
      temp.mean = replicate(maxneval, NA)

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
          ifelse(common.scale,"y = data.eval[,plot.index],", "y = temp.mean,")))

      pylimE = ifelse(common.scale, "ylim = c(y.min,y.max),",
        ifelse(plot.errors, "ylim = if (plot.errors.type == 'all') compute.all.error.range(if (plotOnEstimate) temp.mean else temp.err[,3], temp.all.err) else c(min(na.omit(c(temp.mean - temp.err[,1], temp.err[,3] - temp.err[,1]))), max(na.omit(c(temp.mean + temp.err[,2], temp.err[,3] + temp.err[,2])))),", ""))

      pxlabE = expression(paste("xlab = gen.label(bws$",
          xOrZ, "names[i], paste('", toupper(xOrZ),"', i, sep = '')),",sep=''))

      pylabE = "ylab = paste(ifelse(gradients,
          paste('Gradient Component ', i, ' of', sep=''), ''),
          gen.label(bws$ynames, 'Conditional Mean')),"

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
          ifelse(plotOnEstimate, "ely = na.omit(temp.mean - temp.err[,1]),",
                 "ely = na.omit(temp.err[,3] - temp.err[,1]),")))
      eehyE = expression(ifelse(common.scale,
          ifelse(plotOnEstimate, "ehy = na.omit(data.eval[,plot.index] + data.err[,3*plot.index-1]),",
                 "ehy = na.omit(data.err[,3*plot.index] + data.err[,3*plot.index-1]),"),
          ifelse(plotOnEstimate, "ehy = na.omit(temp.mean + temp.err[,2]),",
                 "ehy = na.omit(temp.err[,3] + temp.err[,2]),")))

      erestE = "plot.errors.style = ifelse(xi.factor,'bar',plot.errors.style),
                plot.errors.bar = ifelse(xi.factor,'I',plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(xi.factor,1,2)"

      plot.index = 0
      xOrZ = "x"

      for (i in 1:bws$xndim){
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


        tobj <- npplreg(txdat = xdat, tydat = ydat, tzdat = zdat,
          exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
          ezdat = ezdat[1:xi.neval,, drop = FALSE],
          bws = bws)

        temp.mean[1:xi.neval] = tobj$mean

        if (plot.errors){
          if (plot.errors.method == "bootstrap"){
            temp.boot <- compute.bootstrap.errors(
                      xdat = xdat,
                      ydat = ydat,
                      zdat = zdat,
                      exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                      ezdat = ezdat[1:xi.neval,, drop = FALSE],
                      gradients = gradients,
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
          data.eval[, plot.index] = temp.mean
          if (plot.errors){
            all.bxp[plot.index] = NA
            all.bxp[[plot.index]] = temp.boot

            data.err[, c(3*plot.index-2,3*plot.index-1,3*plot.index)] = temp.err
            data.err.all[[plot.index]] = temp.all.err
          }
        } else if (plot.behavior != "data") {
          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = na.omit(if (plotOnEstimate) temp.mean else temp.err[,3]),
                all.err = temp.all.err,
                xi.factor = xi.factor)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
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
                           evalx = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                           evalz = ezdat[1:xi.neval,, drop = FALSE],
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
      for (i in 1:bws$zndim){
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

        tobj <- npplreg(txdat = xdat, tydat = ydat, tzdat = zdat,
          exdat = exdat[1:xi.neval,, drop = FALSE],
          ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
          bws = bws)

        temp.mean[1:xi.neval] = tobj$mean

        if (plot.errors){
          if (plot.errors.method == "bootstrap"){
            temp.boot <- compute.bootstrap.errors(
                      xdat = xdat,
                      ydat = ydat,
                      zdat = zdat,
                      exdat = exdat[1:xi.neval,, drop = FALSE],
                      ezdat = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
                      gradients = gradients,
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
          data.eval[, plot.index] = temp.mean
          if (plot.errors){
            all.bxp[plot.index] = NA
            all.bxp[[plot.index]] = temp.boot

            data.err[, c(3*plot.index-2,3*plot.index-1,3*plot.index)] = temp.err
            data.err.all[[plot.index]] = temp.all.err
          }
        } else if (plot.behavior != "data") {
          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = na.omit(if (plotOnEstimate) temp.mean else temp.err[,3]),
                all.err = temp.all.err,
                xi.factor = xi.factor)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
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
                           evalx = exdat[1:xi.neval,, drop = FALSE],
                           evalz = subcol(ezdat,ei,i)[1:xi.neval,, drop = FALSE],
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

      if (common.scale & (plot.behavior != "data")){
        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in 1:(bws$xndim + bws$zndim)) {
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
          jj = 1:(bws$xndim + bws$zndim)*3
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
        
        xOrZ = "x"
        
        for (plot.index in 1:(bws$xndim + bws$zndim)){
          i = ifelse(plot.index <= bws$xndim, plot.index, plot.index - bws$xndim)

          if (plot.index > bws$xndim)
            xOrZ <- "z"
            
          xi.factor = all.isFactor[plot.index]

          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
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

              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
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
            paste("plr",1:(bws$xndim+bws$zndim),sep="")
        
        return (plot.out)
      }
    }
  }


npplot.bandwidth <-
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
      
    miss.x <- missing(xdat)

    if(miss.x && !is.null(bws$formula)){
      tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf, envir = environment(tt))

      xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]
    } else {
      if(miss.x && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["dat"]], environment(bws$call)))
      }
      xdat = toFrame(xdat)
      xdat = na.omit(xdat)
    }

    xq = double(bws$ndim)+xq
    xtrim = double(bws$ndim)+xtrim

    if (missing(plot.errors.method) &
        any(!missing(plot.errors.boot.num), !missing(plot.errors.boot.method),
            !missing(plot.errors.boot.blocklen))){
      warning(paste("plot.errors.method must be set to 'bootstrap' to use bootstrapping.",
                    "\nProceeding without bootstrapping."))
    }

    plot.behavior = match.arg(plot.behavior)
    plot.errors.method = match.arg(plot.errors.method)
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method, where = "plot(...)")
    plot.errors.boot.method = match.arg(plot.errors.boot.method)
    plot.errors.center = match.arg(plot.errors.center)
    plot.errors.type = match.arg(plot.errors.type)

    if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
        is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
      stop("the tail probability plot.errors.alpha must lie in (0,0.5)")
    plot.errors.style = match.arg(plot.errors.style)
    plot.errors.bar = match.arg(plot.errors.bar)

    common.scale = common.scale | (!is.null(ylim))

    if (plot.errors.method == "none" && plot.errors.type == "all") {
      warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
      plot.errors.method <- "bootstrap"
    }

    if (plot.errors.method == "asymptotic") {
      if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")){
        warning("bootstrap quantile bands cannot be calculated with asymptotics, calculating pmzsd errors")
        plot.errors.type = "pmzsd"
      }

      if (plot.errors.center == "bias-corrected") {
        warning("no bias corrections can be calculated with asymptotics, centering on estimate")
        plot.errors.center = "estimate"
      }
    }

    if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
        is.null(plot.errors.boot.blocklen))
      plot.errors.boot.blocklen = b.star(xdat,round=TRUE)[1,1]


    plot.errors = (plot.errors.method != "none")

    if ((bws$ncon + bws$nord == 2) & (bws$nuno == 0) & perspective &
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

      tobj =  npudens(tdat = xdat, edat = x.eval, bws = bws)

      tcomp = parse(text="tobj$dens")

      tdens = matrix(data = eval(tcomp),
        nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      terr = matrix(data = tobj$derr, nrow = nrow(x.eval), ncol = 3)
      terr[,3] = NA
      
      lerr.all <- NULL
      herr.all <- NULL

      if (plot.errors.method == "bootstrap"){
        terr.obj <- compute.bootstrap.errors(xdat = xdat, 
          exdat = x.eval,
          cdf = FALSE,
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
        center.val <- if (pc) terr[,3] else eval(tcomp)

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
        lerr = matrix(data = eval(tcomp) - qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = eval(tcomp) + qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

      }

      if(is.null(zlim)) {
          zlim =
              if (plot.errors){
                  if (plot.errors.type == "all" && !is.null(lerr.all))
                  c(min(c(unlist(lerr.all), lerr)), max(c(unlist(herr.all), herr)))
                  else
                      c(min(lerr), max(herr))
              } else
                  c(min(eval(tcomp)),max(eval(tcomp)))
      }
          
      tret = parse(text=paste(
                     "npdensity",
                     "(bws = bws, eval = x.eval,",
                     "dens",
                     " = eval(tcomp), derr = terr[,1:2], ntrain = nrow(xdat))", sep=""))

      if (plot.behavior != "plot"){
        d1 = eval(tret)
        d1$bias = NA

        if (plot.errors.center == "bias-corrected")
          d1$bias = terr[,3] - eval(tcomp)
        
        if (plot.behavior == "data")
          return ( list(d1 = d1) )
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
                xlab = ifelse(!is.null(xlab),xlab,gen.label(names(xdat)[1], "X1")),
                ylab = ifelse(!is.null(ylab),ylab,gen.label(names(xdat)[2], "X2")),
                zlab = ifelse(!is.null(zlab),zlab,"Joint Density"),
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
      

      if (plot.behavior == "plot-data")
        return ( list(d1 = d1) )        

    } else {

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(bws$ndim),cex=par()$cex)

      ev = xdat[1,,drop = FALSE]

      for (i in 1:bws$ndim)
        ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      maxneval = max(c(sapply(xdat,nlevels),neval))

      exdat = xdat[rep(1, maxneval), , drop = FALSE]

      for (i in 1:bws$ndim)
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

      pfunE = expression(ifelse(xi.factor,
          ifelse(plot.bootstrap & plot.bxp,"bxp","plotFactor"), "plot"))

      pxE = expression(ifelse(common.scale,
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = all.bxp[[i]],", "f = allei[,i],"),
                 "x = allei[,i],"),
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = temp.boot,", "f = ei,"), "x = ei,")))

      pyE = expression(ifelse(xi.factor & plot.bootstrap & plot.bxp, "",
          ifelse(common.scale,"y = data.eval[,i],", "y = temp.dens,")))

      pylimE = ifelse(common.scale, "ylim = c(y.min,y.max),",
        ifelse(plot.errors,
               "ylim = if (plot.errors.type == 'all') compute.all.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.all.err) else compute.default.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.err),",
               ""))

      pxlabE = "xlab = ifelse(!is.null(xlab),xlab,gen.label(bws$xnames[i], paste('X', i, sep = ''))),"

      pylabE = paste("ylab = ", "ifelse(!is.null(ylab),ylab,'Density')", ",")

      prestE = expression(ifelse(xi.factor,"", "type = ifelse(!is.null(type),type,'l'), lty = ifelse(!is.null(lty),lty,par()$lty), col = ifelse(!is.null(col),col,par()$col), lwd = ifelse(!is.null(lwd),lwd,par()$lwd), cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis), cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab), cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main), cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),"))

      pmainE = "main = ifelse(!is.null(main),main,''), sub = ifelse(!is.null(sub),sub,''),"

      ## error plotting expressions
      plotOnEstimate = (plot.errors.center == "estimate")

      efunE = "draw.errors"
      eexE = ifelse(common.scale, "ex = as.numeric(na.omit(allei[,i])),",
          "ex = as.numeric(na.omit(ei)),")
      eelyE = ifelse(common.scale,
        ifelse(plotOnEstimate, "ely = na.omit(data.eval[,i] - data.err[,3*i-2]),",
               "ely = na.omit(data.err[,3*i] - data.err[,3*i-2]),"),
        ifelse(plotOnEstimate, "ely = na.omit(temp.dens - temp.err[,1]),",
               "ely = na.omit(temp.err[,3] - temp.err[,1]),"))
      eehyE = ifelse(common.scale,
        ifelse(plotOnEstimate, "ehy = na.omit(data.eval[,i] + data.err[,3*i-1]),",
               "ehy = na.omit(data.err[,3*i] + data.err[,3*i-1]),"),
        ifelse(plotOnEstimate, "ehy = na.omit(temp.dens + temp.err[,2]),",
               "ehy = na.omit(temp.err[,3] + temp.err[,2]),"))

      erestE = "plot.errors.style = ifelse(xi.factor,'bar',plot.errors.style),
                plot.errors.bar = ifelse(xi.factor,'I',plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(!is.null(lty),lty,ifelse(xi.factor,1,2))"

      ## density / distribution expressions

      devalE = parse(text=paste("npudens",
                       "(tdat = xdat, edat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE], bws = bws)",
                       sep=""))
      
      dcompE = parse(text="tobj$dens")


      doutE = parse(text=paste("npdensity",
                      "(bws = bws, eval = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],",
                      "dens",
                      "= na.omit(temp.dens), derr = na.omit(cbind(-temp.err[,1], temp.err[,2])), ntrain = bws$nobs)",
                      sep=""))

      for (i in 1:bws$ndim){
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

        tobj = eval(devalE)

        temp.dens[1:xi.neval] = eval(dcompE)

        if (plot.errors){
          if (plot.errors.method == "asymptotic")
            temp.err[1:xi.neval,1:2] = replicate(2,qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr)
          else if (plot.errors.method == "bootstrap"){
            temp.boot.raw <- compute.bootstrap.errors(
                      xdat = xdat,
                      exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                      cdf = FALSE,
                      slice.index = i,
                      plot.errors.boot.method = plot.errors.boot.method,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
                      bws = bws)
            temp.err[1:xi.neval,] = temp.boot.raw[["boot.err"]]
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
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = as.numeric(na.omit(if (plotOnEstimate) temp.dens else temp.err[,3])),
                all.err = temp.all.err,
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
            }
          }
        }

        if (plot.behavior != "plot") {
          plot.out[i] = NA
          plot.out[[i]] = eval(doutE)
          plot.out[[i]]$bias = na.omit(temp.dens - temp.err[,3])
          plot.out[[i]]$bxp = temp.boot
        }
      }
      
      if (common.scale & (plot.behavior != "data")){
        jj = 1:bws$ndim*3

        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in 1:bws$ndim) {
            if (is.null(data.err.all[[k]])) next
            nkeep.k <- nrow(data.err.all[[k]]$pointwise)
            center.k <- if (plot.errors.center == "estimate")
              data.eval[1:nkeep.k,k]
            else
              data.err[1:nkeep.k,3*k]
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
        } else if (plot.errors.center == "estimate" | !plot.errors) {
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
        
        for (i in 1:bws$ndim){
          xi.factor = is.factor(xdat[,i])

          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(allei[,i])),
                center = as.numeric(na.omit(if (plotOnEstimate) data.eval[,i] else data.err[,3*i])),
                all.err = data.err.all[[i]],
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
            }
          }
        }
      }

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=c(1,1),cex=par()$cex)

      if (plot.behavior != "plot"){
        names(plot.out) = paste("d",1:bws$ndim,sep="")
        return (plot.out)
      }
    }
  }

npplot.dbandwidth <-
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

    miss.x <- missing(xdat)

    if(miss.x && !is.null(bws$formula)){
      tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf, envir = environment(tt))

      xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]
    } else {
      if(miss.x && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["dat"]], environment(bws$call)))
      }
      xdat = toFrame(xdat)
      xdat = na.omit(xdat)
    }

    xq = double(bws$ndim)+xq
    xtrim = double(bws$ndim)+xtrim

    if (missing(plot.errors.method) &
        any(!missing(plot.errors.boot.num), !missing(plot.errors.boot.method),
            !missing(plot.errors.boot.blocklen))){
      warning(paste("plot.errors.method must be set to 'bootstrap' to use bootstrapping.",
                    "\nProceeding without bootstrapping."))
    }

    plot.behavior = match.arg(plot.behavior)
    plot.errors.method = match.arg(plot.errors.method)
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method, where = "plot(...)")
    plot.errors.boot.method = match.arg(plot.errors.boot.method)
    plot.errors.center = match.arg(plot.errors.center)
    plot.errors.type = match.arg(plot.errors.type)

    if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
        is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
      stop("the tail probability plot.errors.alpha must lie in (0,0.5)")
    plot.errors.style = match.arg(plot.errors.style)
    plot.errors.bar = match.arg(plot.errors.bar)

    common.scale = common.scale | (!is.null(ylim))

    if (plot.errors.method == "none" && plot.errors.type == "all") {
      warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
      plot.errors.method <- "bootstrap"
    }

    if (plot.errors.method == "asymptotic") {
      if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")){
        warning("bootstrap quantile bands cannot be calculated with asymptotics, calculating pmzsd errors")
        plot.errors.type = "pmzsd"
      }

      if (plot.errors.center == "bias-corrected") {
        warning("no bias corrections can be calculated with asymptotics, centering on estimate")
        plot.errors.center = "estimate"
      }
    }

    if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
        is.null(plot.errors.boot.blocklen))
      plot.errors.boot.blocklen = b.star(xdat,round=TRUE)[1,1]


    plot.errors = (plot.errors.method != "none")

    if ((bws$ncon + bws$nord == 2) & (bws$nuno == 0) & perspective &
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

      tobj =  npudist(tdat = xdat, edat = x.eval, bws = bws)

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
        lerr = matrix(data = tobj$dist - qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr,
          nrow = x1.neval, ncol = x2.neval, byrow = FALSE)

        herr = matrix(data = tobj$dist + qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr,
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
                xlab = ifelse(!is.null(xlab),xlab,gen.label(names(xdat)[1], "X1")),
                ylab = ifelse(!is.null(ylab),ylab,gen.label(names(xdat)[2], "X2")),
                zlab = ifelse(!is.null(zlab),zlab,"Joint Distribution"),
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


      if (plot.behavior == "plot-data")
        return ( list(d1 = d1) )        

    } else {

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=n2mfrow(bws$ndim),cex=par()$cex)

      ev = xdat[1,,drop = FALSE]

      for (i in 1:bws$ndim)
        ev[1,i] = uocquantile(xdat[,i], prob=xq[i])

      maxneval = max(c(sapply(xdat,nlevels),neval))

      exdat = xdat[rep(1, maxneval), , drop = FALSE]

      for (i in 1:bws$ndim)
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

      pfunE = expression(ifelse(xi.factor,
          ifelse(plot.bootstrap & plot.bxp,"bxp","plotFactor"), "plot"))

      pxE = expression(ifelse(common.scale,
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = all.bxp[[i]],", "f = allei[,i],"),
                 "x = allei[,i],"),
          ifelse(xi.factor,
                 ifelse(plot.bootstrap & plot.bxp, "z = temp.boot,", "f = ei,"), "x = ei,")))

      pyE = expression(ifelse(xi.factor & plot.bootstrap & plot.bxp, "",
          ifelse(common.scale,"y = data.eval[,i],", "y = temp.dens,")))

      pylimE = ifelse(common.scale, "ylim = c(y.min,y.max),",
        ifelse(plot.errors,
               "ylim = if (plot.errors.type == 'all') compute.all.error.range(if (plotOnEstimate) temp.dens else temp.err[,3], temp.all.err) else c(min(na.omit(c(temp.dens - temp.err[,1], temp.err[,3] - temp.err[,1]))), max(na.omit(c(temp.dens + temp.err[,2], temp.err[,3] + temp.err[,2])))),",
               ""))

      pxlabE = "xlab = ifelse(!is.null(xlab),xlab,gen.label(bws$xnames[i], paste('X', i, sep = ''))),"

      pylabE = "ylab = ifelse(!is.null(ylab),ylab,'Distribution'),"

      prestE = expression(ifelse(xi.factor,"", "type = ifelse(!is.null(type),type,'l'), lty = ifelse(!is.null(lty),lty,par()$lty), col = ifelse(!is.null(col),col,par()$col), lwd = ifelse(!is.null(lwd),lwd,par()$lwd), cex.axis = ifelse(!is.null(cex.axis),cex.axis,par()$cex.axis), cex.lab = ifelse(!is.null(cex.lab),cex.lab,par()$cex.lab), cex.main = ifelse(!is.null(cex.main),cex.main,par()$cex.main), cex.sub = ifelse(!is.null(cex.sub),cex.sub,par()$cex.sub),"))

      pmainE = "main = ifelse(!is.null(main),main,''), sub = ifelse(!is.null(sub),sub,''),"

      ## error plotting expressions
      plotOnEstimate = (plot.errors.center == "estimate")

      efunE = "draw.errors"
      eexE = ifelse(common.scale, "ex = as.numeric(na.omit(allei[,i])),",
          "ex = as.numeric(na.omit(ei)),")
      eelyE = ifelse(common.scale,
        ifelse(plotOnEstimate, "ely = na.omit(data.eval[,i] - data.err[,3*i-2]),",
               "ely = na.omit(data.err[,3*i] - data.err[,3*i-2]),"),
        ifelse(plotOnEstimate, "ely = na.omit(temp.dens - temp.err[,1]),",
               "ely = na.omit(temp.err[,3] - temp.err[,1]),"))
      eehyE = ifelse(common.scale,
        ifelse(plotOnEstimate, "ehy = na.omit(data.eval[,i] + data.err[,3*i-1]),",
               "ehy = na.omit(data.err[,3*i] + data.err[,3*i-1]),"),
        ifelse(plotOnEstimate, "ehy = na.omit(temp.dens + temp.err[,2]),",
               "ehy = na.omit(temp.err[,3] + temp.err[,2]),"))

      erestE = "plot.errors.style = ifelse(xi.factor,'bar',plot.errors.style),
                plot.errors.bar = ifelse(xi.factor,'I',plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = ifelse(xi.factor,1,2)"

      ## density / distribution expressions

      devalE = parse(text="npudist(tdat = xdat, edat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE], bws = bws)")
      
      dcompE = parse(text="tobj$dist")

      doutE = parse(text="npdistribution(bws = bws, eval = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE], dist = na.omit(temp.dens), derr = na.omit(cbind(-temp.err[,1], temp.err[,2])), ntrain = bws$nobs)")

      for (i in 1:bws$ndim){
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

        tobj = eval(devalE)

        temp.dens[1:xi.neval] = eval(dcompE)

        if (plot.errors){
          if (plot.errors.method == "asymptotic")
            temp.err[1:xi.neval,1:2] = replicate(2,qnorm(plot.errors.alpha/2, lower.tail = FALSE)*tobj$derr)
          else if (plot.errors.method == "bootstrap"){
            temp.boot.raw <- compute.bootstrap.errors(
                      xdat = xdat,
                      exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],
                      slice.index = i,
                      plot.errors.boot.method = plot.errors.boot.method,
                      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
                      plot.errors.boot.num = plot.errors.boot.num,
                      plot.errors.center = plot.errors.center,
                      plot.errors.type = plot.errors.type,
                      plot.errors.alpha = plot.errors.alpha,
                      bws = bws)
            temp.err[1:xi.neval,] = temp.boot.raw[["boot.err"]]
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
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(ei)),
                center = as.numeric(na.omit(if (plotOnEstimate) temp.dens else temp.err[,3])),
                all.err = temp.all.err,
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
            }
          }
        }

        if (plot.behavior != "plot") {
          plot.out[i] = NA
          plot.out[[i]] = eval(doutE)
          plot.out[[i]]$bias = na.omit(temp.dens - temp.err[,3])
          plot.out[[i]]$bxp = temp.boot
        }
      }
      
      if (common.scale & (plot.behavior != "data")){
        jj = 1:bws$ndim*3

        if (plot.errors && plot.errors.type == "all") {
          y.min <- Inf
          y.max <- -Inf
          for (k in 1:bws$ndim) {
            if (is.null(data.err.all[[k]])) next
            nkeep.k <- nrow(data.err.all[[k]]$pointwise)
            center.k <- if (plot.errors.center == "estimate")
              data.eval[1:nkeep.k,k]
            else
              data.err[1:nkeep.k,3*k]
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
        } else if (plot.errors.center == "estimate" | !plot.errors) {
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
        
        for (i in 1:bws$ndim){
          xi.factor = is.factor(xdat[,i])

          ## plot evaluation
          eval(parse(text = paste(eval(pfunE), "(", eval(pxE), eval(pyE),
                       eval(pylimE), eval(pxlabE), eval(pylabE), eval(prestE),
                       eval(pmainE), ")")))

          ## error plotting evaluation
          if (plot.errors && !(xi.factor & plot.bootstrap & plot.bxp)){
            if (!xi.factor && !plotOnEstimate)
              lines(na.omit(ei), na.omit(temp.err[,3]), lty = 3)

            if (plot.errors.type == "all") {
              draw.all.error.types(
                ex = as.numeric(na.omit(allei[,i])),
                center = as.numeric(na.omit(if (plotOnEstimate) data.eval[,i] else data.err[,3*i])),
                all.err = data.err.all[[i]],
                plot.errors.style = ifelse(xi.factor, "bar", plot.errors.style),
                plot.errors.bar = ifelse(xi.factor, "I", plot.errors.bar),
                plot.errors.bar.num = plot.errors.bar.num,
                lty = 2,
                add.legend = TRUE)
            } else {
              eval(parse(text = paste(eval(efunE), "(", eval(eexE), eval(eelyE),
                           eval(eehyE), eval(erestE), ")")))
            }
          }
        }
      }

      if (plot.behavior != "data" && plot.par.mfrow)
        par(mfrow=c(1,1),cex=par()$cex)

      if (plot.behavior != "plot"){
        names(plot.out) = paste("d",1:bws$ndim,sep="")
        return (plot.out)
      }
    }
  }


npplot.conbandwidth <-
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
    
    plot.behavior = match.arg(plot.behavior)
    plot.errors.method = match.arg(plot.errors.method)
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method, where = "plot(...)")
    plot.errors.boot.method = match.arg(plot.errors.boot.method)
    plot.errors.center = match.arg(plot.errors.center)
    plot.errors.type = match.arg(plot.errors.type)

    if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
        is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
      stop("the tail probability plot.errors.alpha must lie in (0,0.5)")
    plot.errors.style = match.arg(plot.errors.style)
    plot.errors.bar = match.arg(plot.errors.bar)

    common.scale = common.scale | (!is.null(ylim))

    if (plot.errors.method == "none" && plot.errors.type == "all") {
      warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
      plot.errors.method <- "bootstrap"
    }

    if (plot.errors.method == "asymptotic") {
      if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")){
        warning("bootstrap quantile bands cannot be calculated with asymptotics, calculating pmzsd errors")
        plot.errors.type = "pmzsd"
      }

      if (plot.errors.center == "bias-corrected") {
        warning("no bias corrections can be calculated with asymptotics, centering on estimate")
        plot.errors.center = "estimate"
      }

      if (quantreg & gradients){
        warning(paste("no asymptotic errors available for quantile regression gradients.",
                      "\nOne must instead use bootstrapping."))
        plot.errors.method = "none"
      }
    }

    if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
        is.null(plot.errors.boot.blocklen))
      plot.errors.boot.blocklen = b.star(xdat,round=TRUE)[1,1]


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

      tobj = eval(parse(text = paste(
                          switch(tboo,
                                 "quant" = "npqreg",
                                 "dist" = "npcdist",
                                 "dens" = "npcdens"),
                          "(txdat = xdat, tydat = ydat, exdat =",
                          ifelse(quantreg, "x.eval, tau = tau",
                                 "x.eval[,1, drop = FALSE], eydat = x.eval[,2, drop = FALSE]"),
                          ", bws = bws)", sep="")))

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

      pylabE = paste("ylab =", "ifelse(!is.null(ylab),ylab,paste(", ifelse(gradients,"'GC',j,'of',",''), "tylabE)),")

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

        tobj = eval(parse(text=paste(ifelse(cdf, "npcdist",
                          ifelse(quantreg, "npqreg", "npcdens")),
                          "(txdat = xdat, tydat = ydat,",
                          "exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],",
                          ifelse(quantreg, "tau = tau,",
                                 "eydat = eydat[1:xi.neval,, drop = FALSE],"),
                          "gradients = gradients, bws = bws)",sep="")))

        
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
            eval(parse(text=paste("plot.out[[plot.index]]$",ifelse(gradients,
                         paste("gc",j,"err",sep=""),
                         ifelse(quantreg, "quanterr", "conderr")),
                         "= na.omit(cbind(-temp.err[,1], temp.err[,2]))", sep="")))
            eval(parse(text=paste("plot.out[[plot.index]]$",
                         ifelse(gradients, paste("gc",j,"bias",sep=""), "bias"),
                         "= na.omit(temp.dens - temp.err[,3])", sep="")))
            plot.out[[plot.index]]$bxp = temp.boot
            
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

          tobj = eval(parse(text=paste(ifelse(cdf, "npcdist",
                              ifelse(quantreg, "npqreg", "npcdens")),
                              "(txdat = xdat, tydat = ydat,",
                              ifelse(quantreg, "tau = tau,",
                                     "exdat = exdat[1:xi.neval,, drop = FALSE],"),
                              "eydat = subcol(eydat,ei,i)[1:xi.neval,, drop = FALSE],",
                              "gradients = gradients, bws = bws)",sep="")))

          
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
              eval(parse(text=paste("plot.out[[plot.index]]$",ifelse(gradients,
                           paste("gc",j,"err",sep=""),
                           ifelse(quantreg, "quanterr", "conderr")),
                           "= na.omit(cbind(-temp.err[,1], temp.err[,2]))", sep="")))
              eval(parse(text=paste("plot.out[[plot.index]]$",
                           ifelse(gradients, paste("gc",j,"bias",sep=""), "bias"),
                           "= na.omit(temp.dens - temp.err[,3])", sep="")))
              plot.out[[plot.index]]$bxp = temp.boot
                
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
    
    plot.behavior = match.arg(plot.behavior)
    plot.errors.method = match.arg(plot.errors.method)
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method, where = "plot(...)")
    plot.errors.boot.method = match.arg(plot.errors.boot.method)
    plot.errors.center = match.arg(plot.errors.center)
    plot.errors.type = match.arg(plot.errors.type)

    if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
        is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
      stop("the tail probability plot.errors.alpha must lie in (0,0.5)")
    plot.errors.style = match.arg(plot.errors.style)
    plot.errors.bar = match.arg(plot.errors.bar)

    common.scale = common.scale | (!is.null(ylim))

    if (plot.errors.method == "none" && plot.errors.type == "all") {
      warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
      plot.errors.method <- "bootstrap"
    }

    if (plot.errors.method == "asymptotic") {
      if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")){
        warning("bootstrap quantile bands cannot be calculated with asymptotics, calculating pmzsd errors")
        plot.errors.type = "pmzsd"
      }

      if (plot.errors.center == "bias-corrected") {
        warning("no bias corrections can be calculated with asymptotics, centering on estimate")
        plot.errors.center = "estimate"
      }

      if (quantreg & gradients){
        warning(paste("no asymptotic errors available for quantile regression gradients.",
                      "\nOne must instead use bootstrapping."))
        plot.errors.method = "none"
      }
    }

    if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
        is.null(plot.errors.boot.blocklen))
      plot.errors.boot.blocklen = b.star(xdat,round=TRUE)[1,1]


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

      tobj = eval(parse(text = paste(
                          switch(tboo,
                                 "quant" = "npqreg",
                                 "dist" = "npcdist",
                                 "dens" = "npcdens"),
                          "(txdat = xdat, tydat = ydat, exdat =",
                          ifelse(quantreg, "x.eval, tau = tau",
                                 "x.eval[,1, drop = FALSE], eydat = x.eval[,2, drop = FALSE]"),
                          ", bws = bws)", sep="")))

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

        tobj = eval(parse(text=paste(ifelse(cdf, "npcdist",
                          ifelse(quantreg, "npqreg", "npcdens")),
                          "(txdat = xdat, tydat = ydat,",
                          "exdat = subcol(exdat,ei,i)[1:xi.neval,, drop = FALSE],",
                          ifelse(quantreg, "tau = tau,",
                                 "eydat = eydat[1:xi.neval,, drop = FALSE],"),
                          "gradients = gradients, bws = bws)",sep="")))

        
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
            eval(parse(text=paste("plot.out[[plot.index]]$",ifelse(gradients,
                         paste("gc",j,"err",sep=""),
                         ifelse(quantreg, "quanterr", "conderr")),
                         "= na.omit(cbind(-temp.err[,1], temp.err[,2]))", sep="")))
            eval(parse(text=paste("plot.out[[plot.index]]$",
                         ifelse(gradients, paste("gc",j,"bias",sep=""), "bias"),
                         "= na.omit(temp.dens - temp.err[,3])", sep="")))
            plot.out[[plot.index]]$bxp = temp.boot
            
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

          tobj = eval(parse(text=paste(ifelse(cdf, "npcdist",
                              ifelse(quantreg, "npqreg", "npcdens")),
                              "(txdat = xdat, tydat = ydat,",
                              ifelse(quantreg, "tau = tau,",
                                     "exdat = exdat[1:xi.neval,, drop = FALSE],"),
                              "eydat = subcol(eydat,ei,i)[1:xi.neval,, drop = FALSE],",
                              "gradients = gradients, bws = bws)",sep="")))

          
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
              eval(parse(text=paste("plot.out[[plot.index]]$",ifelse(gradients,
                           paste("gc",j,"err",sep=""),
                           ifelse(quantreg, "quanterr", "conderr")),
                           "= na.omit(cbind(-temp.err[,1], temp.err[,2]))", sep="")))
              eval(parse(text=paste("plot.out[[plot.index]]$",
                           ifelse(gradients, paste("gc",j,"bias",sep=""), "bias"),
                           "= na.omit(temp.dens - temp.err[,3])", sep="")))
              plot.out[[plot.index]]$bxp = temp.boot
                
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

      ydat <- model.response(tmf)
      xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]
    } else {
      if(all(miss.xy) && !is.null(bws$call)){
        xdat <- data.frame(eval(bws$call[["xdat"]], environment(bws$call)))
        ydat = eval(bws$call[["ydat"]], environment(bws$call))
      }
          
      ## catch and destroy NA's
      xdat = toFrame(xdat)
      goodrows = 1:dim(xdat)[1]
      rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
      goodrows[rows.omit] = 0

      if (all(goodrows==0))
        stop("Data has no rows without NAs")

      xdat = xdat[goodrows,,drop = FALSE]
      ydat = ydat[goodrows]
    }


    if (is.null(plot.errors.bar.num))
      plot.errors.bar.num = min(length(ydat),25)

    if (missing(plot.errors.method) &
        any(!missing(plot.errors.boot.num), !missing(plot.errors.boot.method),
            !missing(plot.errors.boot.blocklen))){
      warning(paste("plot.errors.method must be set to 'bootstrap' to use bootstrapping.",
                    "\nProceeding without bootstrapping."))
    }

    plot.behavior = match.arg(plot.behavior)
    plot.errors.method = match.arg(plot.errors.method)
    .npRmpi_guard_bootstrap_plot_autodispatch(plot.errors.method, where = "plot(...)")
    plot.errors.boot.method = match.arg(plot.errors.boot.method)
    plot.errors.center = match.arg(plot.errors.center)
    plot.errors.type = match.arg(plot.errors.type)

    if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
        is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
      stop("the tail probability plot.errors.alpha must lie in (0,0.5)")
    plot.errors.style = match.arg(plot.errors.style)
    plot.errors.bar = match.arg(plot.errors.bar)

    common.scale = common.scale | (!is.null(ylim))

    if (plot.errors.method == "asymptotic") {
      warning(paste("asymptotic errors are not supported with single index regression.\n",
                    "Proceeding without calculating errors"))
      plot.errors.method = "none"
    }
    
    if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
        is.null(plot.errors.boot.blocklen))
      plot.errors.boot.blocklen = b.star(xdat,round=TRUE)[1,1]
    
    if (plot.errors.method == "none" && plot.errors.type == "all") {
      warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
      plot.errors.method <- "bootstrap"
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
          ymin <- eval(parse(text=paste("min(",
                               ifelse(plot.errors,"na.omit(",""),
                               "c(temp.mean",
                               ifelse(plot.errors,"- temp.err[,1]",""),
                               ", temp.err[,3]",
                               ifelse(plot.errors,"- temp.err[,1]",""),
                               "))",
                               ifelse(plot.errors,")",""))))
          ymax <- eval(parse(text=paste("max(",
                               ifelse(plot.errors,"na.omit(",""),
                               "c(temp.mean",
                               ifelse(plot.errors,"+ temp.err[,2]",""),
                               ", temp.err[,3]",
                               ifelse(plot.errors,"+ temp.err[,2]",""),
                               "))",
                               ifelse(plot.errors,")",""))))
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
