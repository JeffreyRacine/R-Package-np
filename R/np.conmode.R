npconmode <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npconmode",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npconmode",bws$call)
        else if (!is.call(bws))
          UseMethod("npconmode",bws)
        else
          UseMethod("npconmode",NULL)
      } else {
        UseMethod("npconmode", NULL)
      }
    } else {
      UseMethod("npconmode", NULL)
    }
  }

.npConmodeValidateNewdataTerms <- function(newdata, xnames) {
  nd <- toFrame(newdata)
  missing.names <- setdiff(xnames, names(nd))
  if (length(missing.names))
    stop(sprintf(
      "newdata must contain columns: %s",
      paste(shQuote(xnames), collapse = ", ")
    ), call. = FALSE)
  invisible(TRUE)
}

.npConmodeValidateCategoricalResponse <- function(tydat) {
  tydat <- toFrame(tydat)
  if (NCOL(tydat) != 1L)
    stop("'tydat' must consist of one (1) factor or ordered factor response",
         call. = FALSE)
  y <- tydat[[1L]]
  if (!is.factor(y))
    stop("npconmode requires a categorical response: supply 'tydat' as a factor or ordered factor",
         call. = FALSE)
  invisible(TRUE)
}

npconmode.formula <-
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    tydat <- tmf[, bws$variableNames[["response"]], drop = FALSE]
    txdat <- tmf[, bws$variableNames[["terms"]], drop = FALSE]
    .npConmodeValidateCategoricalResponse(tydat)

    has.eval <- !is.null(newdata)
    if (has.eval) {
      .npConmodeValidateNewdataTerms(newdata, bws$variableNames[["terms"]])
      has.ey <- bws$variableNames[["response"]] %in% names(newdata)

      if (has.ey){
        umf.args <- list(formula = tt, data = newdata)
        umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
        emf <- umf
        eydat <- emf[, bws$variableNames[["response"]], drop = FALSE]
      } else {
        umf.args <- list(formula = formula(bws)[-2], data = newdata)
        umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
        emf <- umf
      }

      exdat <- emf[, bws$variableNames[["terms"]], drop = FALSE]
    }
    
    cm.args <- list(txdat = txdat, tydat = tydat)
    if (has.eval) {
      cm.args$exdat <- exdat
      if (has.ey)
        cm.args$eydat <- eydat
    }
    cm.args$bws <- bws
    ev <- do.call(npconmode, c(cm.args, list(...)))

    ev$omit <- attr(umf,"na.action")
    ev$rows.omit <- as.vector(ev$omit)
    ev$nobs.omit <- length(ev$rows.omit)

    ev$conmode <- napredict(ev$omit, ev$conmode)
    ev$condens <- napredict(ev$omit, ev$condens)
    ev$conderr <- napredict(ev$omit, ev$conderr)
    if (!is.null(ev$probabilities))
      ev$probabilities <- napredict(ev$omit, ev$probabilities)
    if (!is.null(ev$probability.errors))
      ev$probability.errors <- napredict(ev$omit, ev$probability.errors)
    if (!is.null(ev$probability.repaired.rows))
      ev$probability.repaired.rows <- napredict(ev$omit, ev$probability.repaired.rows)
    if (!is.null(ev$probability.gradients))
      ev$probability.gradients <- .npConmodeNapredictArray(ev$omit, ev$probability.gradients)

    return(ev)
  }

npconmode.call <-
  function(bws, ...) {
    npconmode(txdat = .np_eval_bws_call_arg(bws, "xdat"),
              tydat = .np_eval_bws_call_arg(bws, "ydat"),
              bws = bws, ...)
  }

.npConmodeEffectiveProper <- function(bws, proper) {
  if (!is.null(proper)) {
    if (!is.logical(proper) || length(proper) != 1L || is.na(proper))
      stop("'proper' must be TRUE, FALSE, or NULL")
    return(proper)
  }

  regtype <- bws$regtype
  if (!is.null(bws$regtype.engine))
    regtype <- bws$regtype.engine
  !identical(as.character(regtype)[1L], "lc")
}

.npConmodeProjectSimplexVector <- function(x) {
  u <- sort(x, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  ind <- seq_along(u)
  keep <- which(u - cssv / ind > 0)
  if (!length(keep))
    return(rep(1 / length(x), length(x)))
  theta <- cssv[max(keep)] / max(keep)
  pmax(x - theta, 0)
}

.npConmodeNormalizeProperControl <- function(proper.control) {
  if (is.null(proper.control))
    return(list())
  if (!is.list(proper.control) || is.data.frame(proper.control))
    stop("'proper.control' must be a list")

  nms <- names(proper.control)
  if (is.null(nms))
    nms <- rep("", length(proper.control))
  bad <- nms == "" | !(nms %in% "tol")
  if (any(bad))
    stop(sprintf("unused argument%s in proper.control: %s",
                 if (sum(bad) == 1L) "" else "s",
                 paste(sQuote(nms[bad]), collapse = ", ")))

  proper.control
}

.npConmodeProperProbabilities <- function(pmat,
                                          levels,
                                          proper = TRUE,
                                          proper.control = list()) {
  if (!is.matrix(pmat) || !is.numeric(pmat))
    stop("internal error: conditional mode probabilities must be a numeric matrix")
  if (ncol(pmat) != length(levels))
    stop("internal error: probability matrix does not match the response support")

  proper.control <- .npConmodeNormalizeProperControl(proper.control)
  tol <- proper.control$tol
  if (is.null(tol))
    tol <- sqrt(.Machine$double.eps)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0)
    stop("'proper.control$tol' must be a non-negative scalar")

  raw <- pmat
  finite.row <- rowSums(is.finite(raw)) == ncol(raw)
  raw.rowsum <- rep(NA_real_, nrow(raw))
  raw.rowsum[finite.row] <- rowSums(raw[finite.row,, drop = FALSE])
  raw.min <- rep(NA_real_, nrow(raw))
  raw.min[finite.row] <- apply(raw[finite.row,, drop = FALSE], 1L, min)
  raw.neg.violation <- ifelse(is.na(raw.min), NA_real_, pmax(-raw.min, 0))

  repaired <- raw
  projection.distance <- rep(NA_real_, nrow(raw))
  invalid.rows <- !finite.row

  if (isTRUE(proper)) {
    for (i in which(finite.row)) {
      x <- raw[i,]
      if (any(x < -tol) || abs(sum(x) - 1) > tol) {
        repaired[i,] <- .npConmodeProjectSimplexVector(x)
      } else {
        repaired[i,] <- pmax(x, 0)
        repaired[i,] <- repaired[i,] / sum(repaired[i,])
      }
      projection.distance[i] <- max(abs(repaired[i,] - raw[i,]))
    }
  } else {
    projection.distance[finite.row] <- 0
  }

  repaired.rows <- isTRUE(proper) & is.finite(projection.distance) & projection.distance > tol
  repaired.rows[is.na(repaired.rows)] <- FALSE

  binary.complement.max <- NA_real_
  if (ncol(repaired) == 2L) {
    ok <- rowSums(is.finite(repaired)) == 2L
    if (any(ok))
      binary.complement.max <- max(abs(repaired[ok, 2L] - (1 - repaired[ok, 1L])))
  }

  colnames(repaired) <- as.character(levels)
  info <- list(
    reason = if (!isTRUE(proper)) "not_requested" else if (any(repaired.rows)) "projected" else "already_proper",
    tol = tol,
    nlevels = ncol(raw),
    invalid.rows = sum(invalid.rows),
    repaired.rows = sum(repaired.rows),
    max.negative.violation.raw = suppressWarnings(max(raw.neg.violation, na.rm = TRUE)),
    max.row.sum.deviation.raw = suppressWarnings(max(abs(raw.rowsum - 1), na.rm = TRUE)),
    max.projection.distance = suppressWarnings(max(projection.distance, na.rm = TRUE)),
    binary.complement.max.deviation = binary.complement.max
  )
  for (nm in c("max.negative.violation.raw", "max.row.sum.deviation.raw", "max.projection.distance")) {
    if (!is.finite(info[[nm]]))
      info[[nm]] <- NA_real_
  }

  list(
    probabilities = repaired,
    raw.probabilities = raw,
    repaired.rows = repaired.rows,
    proper.requested = isTRUE(proper),
    proper.applied = any(repaired.rows),
    proper.info = info
  )
}

.npConmodeSelect <- function(pmat, perr) {
  enrow <- nrow(pmat)
  nlev <- ncol(pmat)
  mdens <- double(enrow)
  mderr <- rep(NA_real_, enrow)
  indices <- integer(enrow)
  for (i in seq_len(nlev)) {
    tf <- pmat[, i] >= mdens
    tf[is.na(tf)] <- FALSE
    indices[tf] <- i
    mdens[tf] <- pmat[tf, i]
    mderr[tf] <- perr[tf, i]
  }
  list(indices = indices, condens = mdens, conderr = mderr)
}

.npConmodeNapredictArray <- function(omit, x) {
  if (is.null(x) || !length(omit))
    return(x)
  dx <- dim(x)
  if (is.null(dx) || length(dx) < 3L)
    return(napredict(omit, x))

  keep <- seq_len(dx[1L] + length(omit))[-as.integer(omit)]
  dn <- dimnames(x)
  if (!is.null(dn))
    dn[[1L]] <- NULL
  out <- array(NA_real_,
               dim = c(dx[1L] + length(omit), dx[-1L]),
               dimnames = dn)
  out[keep,,] <- x
  out
}


npconmode.conbandwidth <-
  function (bws,
            txdat = stop("invoked without training data 'txdat'"),
            tydat = stop("invoked without training data 'tydat'"),
            exdat, eydat,
            proper = NULL,
            proper.control = list(),
            probabilities = FALSE,
            gradients = FALSE,
            level = NULL,
            ...){

    probabilities <- npValidateScalarLogical(probabilities, "probabilities")
    gradients <- npValidateScalarLogical(gradients, "gradients")
    txdat = toFrame(txdat)
    tydat = toFrame(tydat)

    no.ex = missing(exdat)
    no.ey = missing(eydat)

    if (!no.ex)
      exdat = toFrame(exdat)

    if (!no.ey)
      eydat = toFrame(eydat)

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(data.frame(txdat, tydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]
    tydat <- tydat[keep.rows,,drop = FALSE]

    if (!no.ex){
      keep.eval <- rep_len(TRUE, nrow(exdat))
      eval.df <- data.frame(exdat)
      if (!no.ey)
        eval.df <- data.frame(eval.df, eydat)
      rows.omit <- attr(na.omit(eval.df), "na.action")
      if (length(rows.omit) > 0L)
        keep.eval[as.integer(rows.omit)] <- FALSE

      exdat <- exdat[keep.eval,,drop = FALSE]

      if (!no.ey)
        eydat <- eydat[keep.eval,,drop = FALSE]

      if (!any(keep.eval))
        stop("Evaluation data has no rows without NAs")
    }


    tnrow = dim(txdat)[1]
    enrow = if (no.ex) tnrow else dim(exdat)[1]

    if (!no.ey && no.ex)
      stop("npconmode: invalid invocation: 'eydat' provided but not 'exdat'")

    if (bws$yndim != 1 || bws$yncon > 0)
      stop("'tydat' must consist of one (1) discrete variable")

    if(no.ey)
      efac <- factor(bws$ydati$all.lev[[1]],levels = bws$ydati$all.lev[[1]], ordered = is.ordered(tydat[,1]))
    else
      efac <- factor(union(bws$ydati$all.lev[[1]], levels(eydat[,1])),
                     levels = union(bws$ydati$all.lev[[1]], levels(eydat[,1])), ordered = is.ordered(tydat[,1]))

    nlev <- nlevels(efac)
    level.values <- levels(efac)
    if (is.null(level)) {
      gradient.level <- level.values[1L]
    } else {
      if (length(level) != 1L || is.na(level) ||
          !(as.character(level) %in% level.values))
        stop("'level' must identify one response level in the fitted conmode object")
      gradient.level <- as.character(level)
    }
    gradient.level.index <- match(gradient.level, level.values)
    if (isTRUE(gradients)) {
      if (bws$xndim < 1L)
        stop("npconmode class-probability gradients/effects require at least one conditioning variable")
    }
    pmat <- matrix(NA_real_, enrow, nlev, dimnames = list(NULL, level.values))
    perr <- matrix(NA_real_, enrow, nlev, dimnames = list(NULL, level.values))
    pgrad <- if (isTRUE(gradients)) {
      matrix(NA_real_,
             nrow = enrow,
             ncol = bws$xndim,
             dimnames = list(NULL, bws$xnames))
    } else {
      NULL
    }

    for (i in seq_len(nlevels(efac))) {
        dens.obj <- npcdens(
          txdat = txdat,
          tydat = tydat,
          exdat = if (no.ex) txdat else exdat,
          eydat = rep(efac[i], enrow),
          bws = bws,
          gradients = isTRUE(gradients) && identical(i, gradient.level.index)
        )
        pmat[, i] <- dens.obj$condens
        perr[, i] <- dens.obj$conderr
        if (isTRUE(gradients) && identical(i, gradient.level.index)) {
          if (is.null(dens.obj$congrad))
            stop("internal error: conditional-density gradient was not returned")
          pgrad[,] <- dens.obj$congrad
        }
    }

    proper.effective <- .npConmodeEffectiveProper(bws, proper)
    proper.out <- .npConmodeProperProbabilities(
      pmat,
      levels = levels(efac),
      proper = proper.effective,
      proper.control = proper.control
    )
    select <- .npConmodeSelect(proper.out$probabilities, perr)
    indices <- select$indices
    mdens <- select$condens
    mderr <- select$conderr
    mderr[proper.out$repaired.rows] <- NA_real_
    cm.args <- list(
      bws = bws,
      xeval = if (no.ex) txdat else exdat,
      conmode = efac[indices],
      condens = mdens,
      conderr = mderr,
      ntrain = nrow(txdat),
      trainiseval = no.ex,
      proper.requested = proper.out$proper.requested,
      proper.applied = proper.out$proper.applied,
      proper.info = proper.out$proper.info,
      gradients = gradients
    )
    if (isTRUE(probabilities) || isTRUE(gradients)) {
      cm.args$probabilities <- proper.out$probabilities
      cm.args$probability.levels <- efac
      cm.args$probability.errors <- perr
      cm.args$probability.errors[proper.out$repaired.rows, ] <- NA_real_
      cm.args$probability.repaired.rows <- proper.out$repaired.rows
      if (!no.ex) {
        cm.args$xtrain <- txdat
        cm.args$ytrain <- tydat
      }
    }
    if (isTRUE(gradients)) {
      cm.args$probability.gradients <- pgrad
      cm.args$probability.gradient.level <- gradient.level
      cm.args$probability.gradient.names <- bws$xnames
      cm.args$probability.gradient.info <- list(
        target = "class probability",
        response = if (nlev == 2L) "binary" else "multinomial",
        level = gradient.level,
        default.level = level.values[1L],
        semantics = "smooth class-probability gradients/effects for one selected response level"
      )
    }
    if (!(no.ey && !no.ex))
      cm.args$yeval <- if (no.ey) tydat else eydat
    con.mode <- do.call(conmode, cm.args)
    
    if (!(no.ey && !no.ex)){
      confusion.matrix <- 
        table(factor(if (no.ex) tydat[,1] else eydat[,1], exclude = NULL),
              factor(con.mode$conmode,exclude = NULL), dnn=c("Actual", "Predicted"))

      cj <- match(levels(factor(if (no.ex) tydat[,1] else eydat[,1], exclude = NULL)),
                  levels(factor(con.mode$conmode,exclude = NULL)), nomatch = 0)
      rj <- cj > 0

      t.diag <- cj
      t.diag[rj] <-  diag(confusion.matrix[rj,cj,drop=FALSE])
      
      CCR.overall <- sum(t.diag)/enrow
      
      CCR.byoutcome <- t.diag/rowSums(confusion.matrix)
      names(CCR.byoutcome) <- levels(factor(if (no.ex) tydat[,1] else eydat[,1], exclude = NULL))

      con.mode$confusion.matrix <- confusion.matrix
      con.mode$CCR.overall <- CCR.overall
      con.mode$CCR.byoutcome <- CCR.byoutcome

      confusion.matrix <- confusion.matrix/enrow
      t.diag <- t.diag/enrow

      fit.mcfadden <- sum(t.diag) - (sum(confusion.matrix^2)-sum(t.diag^2))
      con.mode$fit.mcfadden <- fit.mcfadden
    }
    con.mode
  }

npconmode.default <- function(bws, txdat, tydat,
                              nomad = FALSE,
                              proper = NULL,
                              proper.control = list(),
                              probabilities = FALSE,
                              gradients = FALSE,
                              level = NULL,
                              ...){
  nomad <- npValidateScalarLogical(nomad, "nomad")
  probabilities <- npValidateScalarLogical(probabilities, "probabilities")
  gradients <- npValidateScalarLogical(gradients, "gradients")
  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  txdat.named <- any(sc.names == "txdat")
  tydat.named <- any(sc.names == "tydat")

  no.bws <- missing(bws)
  no.txdat <- missing(txdat)
  no.tydat <- missing(tydat)
  has.explicit.bws <- (!no.bws) && isa(bws, "conbandwidth")

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  if(!no.tydat) {
    if(!tydat.named)
      tydat <- toFrame(tydat)
    .npConmodeValidateCategoricalResponse(tydat)
  }

  sc.bw <- sc
  
  sc.bw[[1]] <- quote(npcdensbw)
  sc.bw$proper <- NULL
  sc.bw$proper.control <- NULL
  sc.bw$probabilities <- NULL
  sc.bw$gradients <- NULL
  sc.bw$level <- NULL
  sc.bw$newdata <- NULL
  sc.bw$exdat <- NULL
  sc.bw$eydat <- NULL

  bws.formula <- (!no.bws) && inherits(bws, "formula")
  if (bws.formula) {
    ib <- match("bws", names(sc.bw), nomatch = 0L)
    if (ib > 0L) names(sc.bw)[ib] <- "formula"
  }

  if(bws.named && !bws.formula){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('txdat','tydat')
  nstxy <- c('xdat','ydat')
  
  m.txy <- match(ostxy, names(sc.bw), nomatch = 0)

  if(any(m.txy > 0)) {
    names(sc.bw)[m.txy] <- nstxy[m.txy > 0]
  }

  if (bws.formula && no.tydat) {
    mf.call <- sc.bw
    mf.call[[1]] <- quote(stats::model.frame)
    keep <- match(c("formula", "data", "subset", "na.action"),
                  names(mf.call), nomatch = 0L)
    mf.call <- mf.call[c(1L, keep)]
    if (!("formula" %in% names(mf.call)))
      mf.call$formula <- bws
    mf <- eval(mf.call, parent.frame())
    y <- stats::model.response(mf)
    if (is.null(y))
      stop("npconmode requires a categorical response in the formula",
           call. = FALSE)
    .npConmodeValidateCategoricalResponse(data.frame(y))
  }

  use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
    sc.bw,
    caller_env = parent.frame()
  )

  tbw <- if (!has.explicit.bws) {
    if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting conditional density bandwidth",
        .np_eval_bw_call(sc.bw, caller_env = parent.frame())
      )
    } else {
      .np_eval_bw_call(sc.bw, caller_env = parent.frame())
    }
  } else {
    .np_eval_bw_call(sc.bw, caller_env = parent.frame())
  }

  call.args <- list(bws = tbw)
  if (no.bws) {
    call.args$txdat <- txdat
    call.args$tydat <- tydat
  } else {
    if (txdat.named) call.args$txdat <- txdat
    if (tydat.named) call.args$tydat <- tydat
    if ((!bws.named) && (!txdat.named) && (!no.tydat) && (!tydat.named)) {
      call.args <- c(call.args, list(tydat))
    }
  }
  call.args$proper <- proper
  call.args$proper.control <- proper.control
  call.args$probabilities <- probabilities
  call.args$gradients <- gradients
  call.args$level <- level
  do.call(npconmode, c(call.args, list(...)))
}
