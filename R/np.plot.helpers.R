## the idea is that you have done bandwidth selection
## you just need to supply training data and the bandwidth
## this tool will help you visualize the result

.np_seed_enter <- function(random.seed = 42L) {
  save.seed <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    save.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    exists.seed <- TRUE
  } else {
    exists.seed <- FALSE
  }

  set.seed(random.seed)
  list(exists.seed = exists.seed, save.seed = save.seed)
}

.np_seed_exit <- function(state, remove_if_absent = FALSE) {
  if (isTRUE(state$exists.seed)) {
    assign(".Random.seed", state$save.seed, envir = .GlobalEnv)
  } else if (isTRUE(remove_if_absent) &&
             exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

.np_with_seed <- function(random.seed = 42L, code) {
  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state), add = TRUE)
  force(code)
}

.np_mammen_draws <- function(n, B) {
  a <- (1 - sqrt(5)) / 2
  p.a <- (sqrt(5) + 1) / (2 * sqrt(5))
  u <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  out <- matrix(1 - a, nrow = n, ncol = B)
  out[u <= p.a] <- a
  out
}

.np_rademacher_draws <- function(n, B) {
  u <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  out <- matrix(1.0, nrow = n, ncol = B)
  out[u <= 0.5] <- -1.0
  out
}

.np_wild_draws <- function(n, B, wild = c("mammen", "rademacher")) {
  if (length(wild) > 1L)
    wild <- wild[1L]
  wild <- match.arg(wild, c("mammen", "rademacher"))
  if (identical(wild, "mammen")) {
    return(.np_mammen_draws(n = n, B = B))
  }
  .np_rademacher_draws(n = n, B = B)
}

.np_wildhat_chunk_size <- function(n, B) {
  chunk.opt <- getOption("np.plot.wildhat.chunk.size")
  if (!is.null(chunk.opt)) {
    chunk.opt <- as.integer(chunk.opt)
    if (length(chunk.opt) != 1L || is.na(chunk.opt) || chunk.opt < 1L)
      stop("option 'np.plot.wildhat.chunk.size' must be a positive integer")
    return(min(B, chunk.opt))
  }

  if (n < 1L || B < 1L)
    return(1L)

  # Keep the temporary n x chunk response matrix in a moderate memory range.
  target.bytes <- 64 * 1024 * 1024
  chunk <- as.integer(floor(target.bytes / (8 * n)))
  if (!is.finite(chunk) || is.na(chunk) || chunk < 1L)
    chunk <- 1L
  min(B, chunk)
}

.np_wildhat_boot_t <- function(H, fit.mean, residuals, B, wild = c("mammen", "rademacher")) {
  B <- as.integer(B)
  n <- length(residuals)
  if (length(fit.mean) != n)
    stop("length mismatch between fitted means and residuals for wild-hat bootstrap")
  if (B < 1L)
    stop("argument 'plot.errors.boot.num' must be a positive integer")

  chunk.size <- .np_wildhat_chunk_size(n = n, B = B)
  out <- matrix(NA_real_, nrow = B, ncol = nrow(H))
  fit.mean <- as.double(fit.mean)
  residuals <- as.double(residuals)

  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    draws <- .np_wild_draws(n = n, B = bsz, wild = wild)
    ystar <- matrix(fit.mean, nrow = n, ncol = bsz) +
      matrix(residuals, nrow = n, ncol = bsz) * draws
    out[start:stopi, ] <- t(H %*% ystar)
    start <- stopi + 1L
  }

  out
}

gen.label = function(label, altlabel){
  paste(if (is.null(label)) altlabel else label)
}

gen.tflabel = function(condition, tlabel, flabel){
  paste(if (isTRUE(condition)) tlabel else flabel)
}

draw.error.bands = function(ex, ely, ehy, lty = 2, col = par("col")){
  lines(ex,ely,lty=lty,col=col)
  lines(ex,ehy,lty=lty,col=col)
}

draw.error.bars = function(ex, ely, ehy, hbar = TRUE, hbarscale = 0.3, lty = 2, col = par("col")){
  yy = double(3*length(ex))
  jj = seq_along(ex)*3

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
    
    hdelta = pmin(yg, hbardist)/2
    xx[jj-2] = ex - hdelta
    xx[jj-1] = ex + hdelta
    
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
  dot.args <- list(...)
  dot.names <- names(dot.args)
  has.user.lty <- !is.null(dot.names) && any(dot.names == "lty")

  if (has.user.lty) {
    do.call(plot, c(list(x = f, y = y), dot.args))
  } else {
    plot(x = f, y = y, lty = "blank", ...)
  }

  l.f = rep(f, each=3)
  l.f[3*seq_along(f)] = NA

  l.y = unlist(lapply(y, function (p) { c(0,p,NA) }))

  lines(x = l.f, y = l.y, lty = 2)
  points(x = f, y = y)
}

.np_plot_panel_fun <- function(plot.bootstrap, plot.bxp) {
  if (plot.bootstrap && plot.bxp) bxp else plotFactor
}

.np_plot_resolve_xydat <- function(bws, xdat, ydat, miss.xy) {
  if (any(miss.xy) && !all(miss.xy))
    stop("one of, but not both, xdat and ydat was specified")

  if (all(miss.xy) && !is.null(bws$formula)) {
    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1, m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"), "term.labels"), drop = FALSE]
    return(list(xdat = xdat, ydat = ydat))
  }

  if (all(miss.xy) && !is.null(bws$call)) {
    xdat <- data.frame(.np_eval_bws_call_arg(bws, "xdat"))
    ydat <- .np_eval_bws_call_arg(bws, "ydat")
  }

  xdat <- toFrame(xdat)
  goodrows <- seq_len(nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (!is.null(rows.omit))
    goodrows[rows.omit] <- 0

  if (all(goodrows == 0))
    stop("Data has no rows without NAs")

  goodrows <- goodrows[goodrows != 0]
  list(xdat = xdat[goodrows, , drop = FALSE],
       ydat = ydat[goodrows])
}

.npRmpi_guard_bootstrap_plot_autodispatch <- function(plot.errors.method,
                                                      where = "plot()",
                                                      allow.direct.bootstrap = FALSE) {
  invisible(TRUE)
}

.npRmpi_plot_behavior_for_rank <- function(plot.behavior) {
  if (!.npRmpi_autodispatch_called_from_bcast())
    return(plot.behavior)

  rank <- tryCatch(mpi.comm.rank(), error = function(e) NA_integer_)
  if (is.na(rank) || rank == 0L)
    return(plot.behavior)

  "data"
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

.np_plot_normalize_common_options <- function(plot.behavior,
                                             plot.errors.method,
                                             plot.errors.boot.method,
                                             plot.errors.boot.wild = c("rademacher", "mammen"),
                                             plot.errors.boot.blocklen,
                                             plot.errors.center,
                                             plot.errors.type,
                                             plot.errors.alpha,
                                             plot.errors.style,
                                             plot.errors.bar,
                                             xdat,
                                             common.scale,
                                             ylim,
                                             allow_asymptotic_quantile = TRUE) {
  scalar_choice <- function(value, default) {
    if (is.null(value) || length(value) < 1L || is.na(value[1L])) default else value[1L]
  }

  plot.behavior <- match.arg(
    scalar_choice(plot.behavior, "plot"),
    c("plot", "plot-data", "data")
  )
  plot.errors.method <- match.arg(
    scalar_choice(plot.errors.method, "none"),
    c("none", "bootstrap", "asymptotic")
  )
  plot.errors.boot.method <- match.arg(
    scalar_choice(plot.errors.boot.method, "wild-hat"),
    c("wild-hat", "inid", "fixed", "geom")
  )
  plot.errors.boot.wild <- match.arg(
    scalar_choice(plot.errors.boot.wild, "rademacher"),
    c("rademacher", "mammen")
  )
  plot.errors.center <- match.arg(
    scalar_choice(plot.errors.center, "estimate"),
    c("estimate", "bias-corrected")
  )
  plot.errors.type <- match.arg(
    scalar_choice(plot.errors.type, "simultaneous"),
    c("simultaneous", "pointwise", "bonferroni", "pmzsd", "all")
  )

  if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
      is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
    stop("the tail probability plot.errors.alpha must lie in (0,0.5)")

  plot.errors.style <- match.arg(
    scalar_choice(plot.errors.style, "band"),
    c("band", "bar")
  )
  plot.errors.bar <- match.arg(
    scalar_choice(plot.errors.bar, "|"),
    c("|", "I")
  )

  common.scale <- common.scale | (!is.null(ylim))

  if (plot.errors.method == "none" && plot.errors.type == "all") {
    warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
    plot.errors.method <- "bootstrap"
  }

  if (allow_asymptotic_quantile && plot.errors.method == "asymptotic") {
    if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      warning("bootstrap quantile bands cannot be calculated with asymptotics, calculating pmzsd errors")
      plot.errors.type <- "pmzsd"
    }

    if (plot.errors.center == "bias-corrected") {
      warning("no bias corrections can be calculated with asymptotics, centering on estimate")
      plot.errors.center <- "estimate"
    }
  }

  if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
      is.null(plot.errors.boot.blocklen))
    plot.errors.boot.blocklen <- b.star(xdat, round = TRUE)[1,1]

  list(
    plot.behavior = plot.behavior,
    plot.errors.method = plot.errors.method,
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.wild = plot.errors.boot.wild,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    plot.errors.style = plot.errors.style,
    plot.errors.bar = plot.errors.bar,
    common.scale = common.scale,
    plot.errors = (plot.errors.method != "none")
  )
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
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL
    is.wild.hat <- plot.errors.boot.method == "wild-hat"
    is.inid <- plot.errors.boot.method == "inid"

    if (is.wild.hat && gradients) {
      cont.idx <- which(bws$xdati$icon)
      if (is.na(match(slice.index, cont.idx))) {
        warning("plot.errors.boot.method='wild-hat' supports gradients only for continuous slices; using requested bootstrap method fallback")
        is.wild.hat <- FALSE
      }
    }

    if (is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.train <- suppressWarnings(npreg(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        gradients = FALSE,
        warn.glp.gradient = FALSE
      ))

      s.vec <- NULL
      if (gradients) {
        cont.idx <- which(bws$xdati$icon)
        cpos <- match(slice.index, cont.idx)
        gorder <- if (length(gradient.order) == 1L) {
          rep.int(as.integer(gradient.order), length(cont.idx))
        } else {
          as.integer(gradient.order)
        }
        if (length(gorder) != length(cont.idx))
          gorder <- rep.int(1L, length(cont.idx))
        s.vec <- integer(length(cont.idx))
        s.vec[cpos] <- gorder[cpos]
      }

      H <- suppressWarnings(npreghat(
        bws = bws,
        txdat = xdat,
        exdat = exdat,
        s = s.vec,
        output = "matrix"
      ))

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - fit.train$mean)
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- list(
        t = .np_wildhat_boot_t(
          H = H,
          fit.mean = fit.train$mean,
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
    } else {
      boofun <- if (is.inid) {
        function(data, indices) {
          fit <- suppressWarnings(npreg(
            txdat = xdat[indices, , drop = FALSE],
            tydat = ydat[indices],
            exdat = exdat, bws = bws,
            gradients = gradients,
            gradient.order = gradient.order,
            warn.glp.gradient = FALSE
          ))
          if (gradients) fit$grad[, slice.index] else fit$mean
        }
      } else {
        function(tsb) {
          fit <- suppressWarnings(npreg(
            txdat = tsb[, 1:(ncol(tsb) - 1), drop = FALSE],
            tydat = tsb[, ncol(tsb)],
            exdat = exdat, bws = bws,
            gradients = gradients,
            gradient.order = gradient.order,
            warn.glp.gradient = FALSE
          ))
          if (gradients) fit$grad[, slice.index] else fit$mean
        }
      }

      if (is.inid){
        boot.out <- boot(
          data = data.frame(xdat, ydat),
          statistic = boofun,
          R = plot.errors.boot.num
        )
      } else {
        boot.out <- tsboot(
          tseries = data.frame(xdat, ydat),
          statistic = boofun,
          R = plot.errors.boot.num,
          l = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      }
    }

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
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
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    miss.z <- missing(zdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- plot.errors.boot.method == "wild-hat"
    is.inid <- plot.errors.boot.method == "inid"

    if (is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.args <- list(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        iterate = FALSE
      )
      hat.args <- list(
        bws = bws,
        txdat = xdat,
        exdat = exdat,
        output = "matrix",
        iterate = FALSE
      )
      if (!miss.z) {
        fit.args$tzdat <- zdat
        hat.args$tzdat <- zdat
        hat.args$ezdat <- ezdat
      }

      fit.train <- do.call(npscoef, fit.args)
      H <- do.call(npscoefhat, hat.args)

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - as.vector(fit.train$mean))
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- list(
        t = .np_wildhat_boot_t(
          H = H,
          fit.mean = as.vector(fit.train$mean),
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
    } else {
      xcols <- seq_len(ncol(xdat))
      ycol <- ncol(xdat) + 1L
      zcols <- if (miss.z) integer(0) else (ycol + 1L):(ycol + ncol(zdat))

      boofun <- if (is.inid) {
        function(data, indices) {
          npscoef(
            txdat = xdat[indices, , drop = FALSE],
            tydat = ydat[indices],
            tzdat = if (miss.z) NULL else zdat[indices, , drop = FALSE],
            exdat = exdat,
            ezdat = if (miss.z) NULL else ezdat,
            bws = bws,
            iterate = FALSE
          )$mean
        }
      } else {
        function(tsb) {
          npscoef(
            txdat = tsb[, xcols, drop = FALSE],
            tydat = tsb[, ycol],
            tzdat = if (miss.z) NULL else tsb[, zcols, drop = FALSE],
            exdat = exdat,
            ezdat = if (miss.z) NULL else ezdat,
            bws = bws,
            iterate = FALSE
          )$mean
        }
      }

      boot.data <- if (miss.z) data.frame(xdat, ydat) else data.frame(xdat, ydat, zdat)
      if (is.inid) {
        boot.out <- boot(data = boot.data, statistic = boofun, R = plot.errors.boot.num)
      } else {
        boot.out <- tsboot(
          tseries = boot.data, statistic = boofun, R = plot.errors.boot.num,
          l = plot.errors.boot.blocklen, sim = plot.errors.boot.method
        )
      }
    }

    all.bp <- list()

    if ((slice.index > 0) && (((slice.index <= ncol(xdat)) && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])) ||
                              ((slice.index > ncol(xdat)) && (bws$zdati$iord[slice.index-ncol(xdat)] || bws$zdati$iuno[slice.index-ncol(xdat)])))) {
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

      for (i in seq_along(u.lev)){
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
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- plot.errors.boot.method == "wild-hat"
    is.inid <- plot.errors.boot.method == "inid"

    if (is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.train <- npplreg(
        txdat = xdat,
        tydat = ydat,
        tzdat = zdat,
        bws = bws
      )
      H <- npplreghat(
        bws = bws,
        txdat = xdat,
        tzdat = zdat,
        exdat = exdat,
        ezdat = ezdat,
        output = "matrix"
      )

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - as.vector(fit.train$mean))
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- list(
        t = .np_wildhat_boot_t(
          H = H,
          fit.mean = as.vector(fit.train$mean),
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
    } else {
      boofun <- if (is.inid) {
        function(data, indices) {
          npplreg(
            txdat = xdat[indices, , drop = FALSE],
            tydat = ydat[indices],
            tzdat = zdat[indices, , drop = FALSE],
            exdat = exdat, ezdat = ezdat, bws = bws
          )$mean
        }
      } else {
        function(tsb) {
          npplreg(
            txdat = tsb[, seq_len(ncol(xdat)), drop = FALSE],
            tydat = tsb[, ncol(xdat) + 1L],
            tzdat = tsb[, (ncol(xdat) + 2L):ncol(tsb), drop = FALSE],
            exdat = exdat, ezdat = ezdat, bws = bws
          )$mean
        }
      }

      if (is.inid){
        boot.out = boot(data = data.frame(xdat,ydat,zdat), statistic = boofun,
          R = plot.errors.boot.num)
      } else {
        boot.out = tsboot(tseries = data.frame(xdat,ydat,zdat), statistic = boofun,
          R = plot.errors.boot.num,
          l = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method)
      }
    }

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$zdati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
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
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"

    boofun <- if (is.inid) {
      function(data, indices) {
        fit <- if (cdf) {
          npudist(tdat = xdat[indices, , drop = FALSE], edat = exdat, bws = bws)
        } else {
          npudens(tdat = xdat[indices, , drop = FALSE], edat = exdat, bws = bws)
        }
        if (cdf) fit$dist else fit$dens
      }
    } else {
      function(tsb) {
        fit <- if (cdf) {
          npudist(tdat = tsb, edat = exdat, bws = bws)
        } else {
          npudens(tdat = tsb, edat = exdat, bws = bws)
        }
        if (cdf) fit$dist else fit$dens
      }
    }

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

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
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
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"

    boofun <- if (is.inid) {
      function(data, indices) {
        npudist(tdat = xdat[indices, , drop = FALSE], edat = exdat, bws = bws)$dist
      }
    } else {
      function(tsb) {
        npudist(tdat = tsb, edat = exdat, bws = bws)$dist
      }
    }

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

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
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
    exdat = toFrame(exdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    tboo =
      if(quantreg) "quant"
      else if (cdf) "dist"
      else "dens"

    is.inid = plot.errors.boot.method=="inid"

    fit.cond <- function(tx, ty) {
      switch(
        tboo,
        quant = npqreg(txdat = tx, tydat = ty, exdat = exdat, tau = tau, bws = bws, gradients = gradients),
        dist = npcdist(txdat = tx, tydat = ty, exdat = exdat, eydat = eydat, bws = bws, gradients = gradients),
        dens = npcdens(txdat = tx, tydat = ty, exdat = exdat, eydat = eydat, bws = bws, gradients = gradients)
      )
    }
    out.cond <- function(fit) {
      switch(
        tboo,
        quant = if (gradients) fit$yqgrad[, gradient.index] else fit$quantile,
        dist = if (gradients) fit$congrad[, gradient.index] else fit$condist,
        dens = if (gradients) fit$congrad[, gradient.index] else fit$condens
      )
    }
    boofun <- if (is.inid) {
      function(data, indices) out.cond(fit.cond(
        tx = xdat[indices, , drop = FALSE],
        ty = ydat[indices, , drop = FALSE]
      ))
    } else {
      function(tsb) out.cond(fit.cond(
        tx = tsb[, seq_len(ncol(xdat)), drop = FALSE],
        ty = tsb[, (ncol(xdat) + 1L):ncol(tsb), drop = FALSE]
      ))
    }
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
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
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
    exdat = toFrame(exdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    tboo =
      if(quantreg) "quant"
      else if (cdf) "dist"
      else "dens"

    is.inid = plot.errors.boot.method=="inid"

    fit.cond <- function(tx, ty) {
      switch(
        tboo,
        quant = npqreg(txdat = tx, tydat = ty, exdat = exdat, tau = tau, bws = bws, gradients = gradients),
        dist = npcdist(txdat = tx, tydat = ty, exdat = exdat, eydat = eydat, bws = bws, gradients = gradients),
        dens = npcdens(txdat = tx, tydat = ty, exdat = exdat, eydat = eydat, bws = bws, gradients = gradients)
      )
    }
    out.cond <- function(fit) {
      switch(
        tboo,
        quant = if (gradients) fit$yqgrad[, gradient.index] else fit$quantile,
        dist = if (gradients) fit$congrad[, gradient.index] else fit$condist,
        dens = if (gradients) fit$congrad[, gradient.index] else fit$condens
      )
    }
    boofun <- if (is.inid) {
      function(data, indices) out.cond(fit.cond(
        tx = xdat[indices, , drop = FALSE],
        ty = ydat[indices, , drop = FALSE]
      ))
    } else {
      function(tsb) out.cond(fit.cond(
        tx = tsb[, seq_len(ncol(xdat)), drop = FALSE],
        ty = tsb[, (ncol(xdat) + 1L):ncol(tsb), drop = FALSE]
      ))
    }
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
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
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
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){

    boot.err = matrix(data = NA, nrow = nrow(xdat), ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- plot.errors.boot.method == "wild-hat"
    is.inid <- plot.errors.boot.method=="inid"

    if (is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.train <- npindex(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        gradients = FALSE
      )
      H <- npindexhat(
        bws = bws,
        txdat = xdat,
        exdat = xdat,
        output = "matrix",
        s = if (gradients) 1L else 0L
      )

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - as.vector(fit.train$mean))
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- list(
        t = .np_wildhat_boot_t(
          H = H,
          fit.mean = as.vector(fit.train$mean),
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
    } else {
      ## beta[1] is always 1.0, so use first column of gradients matrix ... 
      boofun <- if (is.inid) {
        function(data, indices) {
          fit <- npindex(
            txdat = xdat[indices, , drop = FALSE],
            tydat = ydat[indices],
            exdat = xdat, bws = bws,
            gradients = gradients
          )
          if (gradients) fit$grad[,1] else fit$mean
        }
      } else {
        function(tsb) {
          fit <- npindex(
            txdat = tsb[, 1:(ncol(tsb)-1), drop = FALSE],
            tydat = tsb[, ncol(tsb)],
            exdat = xdat, bws = bws,
            gradients = gradients
          )
          if (gradients) fit$grad[,1] else fit$mean
        }
      }

      if (is.inid){
        boot.out = boot(data = data.frame(xdat,ydat), statistic = boofun,
          R = plot.errors.boot.num)
      } else {
        boot.out = tsboot(tseries = data.frame(xdat,ydat), statistic = boofun,
          R = plot.errors.boot.num,
          l = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method)
      }
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
  if(anyNA(prob)) stop("'prob' contains missing values")
  if(any(prob < 0 | prob > 1, na.rm = TRUE)) stop("'prob' outside [0,1]")
  if(anyNA(x)) stop("missing values and NaN's not allowed")
  if (is.ordered(x)){
    x <- droplevels(x)
    tq = unclass(table(x))
    tq = tq / sum(tq)
    tq[length(tq)] <- 1.0
    bscape <- levels(x)
    tq <- cumsum(tq)
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
