.npConmodeMetadataFromBws <- function(bws) {
  fields <- c(
    "regtype",
    "pregtype",
    "basis",
    "degree",
    "bernstein.basis",
    "regtype.engine",
    "basis.engine",
    "degree.engine",
    "bernstein.basis.engine",
    "degree.search",
    "nomad.shortcut",
    "nomad.time",
    "powell.time"
  )
  out <- setNames(vector("list", length(fields)), fields)
  for (field in fields) {
    if (!is.null(bws[[field]]))
      out[[field]] <- bws[[field]]
  }
  out$nomad <- if (!is.null(bws$nomad)) {
    isTRUE(bws$nomad)
  } else if (!is.null(bws$nomad.shortcut)) {
    isTRUE(bws$nomad.shortcut$enabled)
  } else {
    NULL
  }
  out$search.engine <- if (!is.null(bws$search.engine)) {
    as.character(bws$search.engine)[1L]
  } else if (!is.null(bws$degree.search$mode)) {
    as.character(bws$degree.search$mode)[1L]
  } else {
    NULL
  }
  out
}

conmode =
  function(bws, xeval, yeval = NA,
           conmode,
           condens, conderr = NA,
           confusion.matrix = NA,
           CCR.overall = NA,
           CCR.byoutcome = NA,
           fit.mcfadden = NA,
           ntrain, trainiseval = FALSE,
           proper.requested = FALSE,
           proper.applied = FALSE,
           proper.info = NULL,
           probabilities = NULL,
           probability.levels = NULL,
           probability.gradients = NULL,
           probability.gradient.level = NULL,
           probability.gradient.names = NULL,
           probability.gradient.info = NULL,
           gradients = FALSE){

    if (missing(bws) || missing(xeval) || missing(conmode) || missing(condens) || missing(ntrain))
      stop("improper invocation of conmode constructor")

    d = list(
      xbw = bws$xbw,
      ybw = bws$ybw,
      bws = bws,
      xnames = bws$xnames,
      ynames = bws$ynames,
      nobs = nrow(xeval),
      xndim = bws$xndim,
      yndim = bws$yndim,
      xnord = bws$xnord,
      xnuno = bws$xnuno,
      xncon = bws$xncon,
      ynord = bws$ynord,
      ynuno = bws$ynuno,
      yncon = bws$yncon,
      pscaling = bws$pscaling,
      ptype = bws$ptype,
      pcxkertype = bws$pcxkertype,
      puxkertype = bws$puxkertype,
      poxkertype = bws$poxkertype,
      pcykertype = bws$pcykertype,
      puykertype = bws$puykertype,
      poykertype = bws$poykertype,
      xeval = xeval,
      yeval = yeval,
      conmode = conmode,
      condens = condens,
      conderr = conderr,
      proper.requested = proper.requested,
      proper.applied = proper.applied,
      proper.info = proper.info,
      confusion.matrix = confusion.matrix,
      CCR.overall = CCR.overall,
      CCR.byoutcome = CCR.byoutcome,
      fit.mcfadden = fit.mcfadden,
      ntrain = ntrain,
      trainiseval = trainiseval,
      gradients = gradients)

    if (!is.null(probabilities)) {
      d$probabilities <- probabilities
      d$probability.levels <- probability.levels
    }
    if (!is.null(probability.gradients)) {
      d$probability.gradients <- probability.gradients
      d$probability.gradient.level <- probability.gradient.level
      d$probability.gradient.names <- probability.gradient.names
      d$probability.gradient.info <- probability.gradient.info
    }

    metadata <- .npConmodeMetadataFromBws(bws)
    d[names(metadata)] <- metadata

    class(d) = "conmode"

    d
  }


print.conmode <- function(x, ...){
  cat("\nConditional Mode data: ", x$ntrain, " training points,",
      if (x$trainiseval) "" else paste(" and ", x$nobs, " evaluation points,\n", sep=""),
      " in ", x$xndim + x$yndim, " variable(s)",
      "\n(", x$yndim, " dependent variable(s), and ", x$xndim, " explanatory variable(s))\n\n",
      sep="")
  print(matrix(x$ybw,ncol=x$yndim,dimnames=list(paste("Dep. Var. ",x$pscaling,":",sep=""),x$ynames)))

  print(matrix(x$xbw,ncol=x$xndim,dimnames=list(paste("Exp. Var. ",x$pscaling,":",sep=""),x$xnames)))

  cat(genBwSelStr(x$bws))

  cat(genBwKerStrs(x$bws))
  cat(genTimingStr(x$bws))
  if (!is.null(x$proper.info) && isTRUE(x$proper.requested)) {
    cat("Proper conditional probabilities: ",
        if (isTRUE(x$proper.applied)) "projected" else "already proper",
        "\n", sep="")
  }
  if (isTRUE(x$gradients)) {
    cat("Class-probability gradients/effects: available",
        if (!is.null(x$probability.gradient.level))
          paste0(" for level ", as.character(x$probability.gradient.level))
        else "",
        "\n", sep="")
  }
  
  cat("\n\n")
  if(!missing(...))
    print(...)
  invisible(x)
}

fitted.conmode <- function(object, ...) {
 object$condens 
}
mode.conmode <- function(x) { x$conmode }
predict.conmode <- function(object,
                            newdata = NULL,
                            type = c("class", "prob"),
                            ...) {
  type <- match.arg(type)
  if (!is.null(newdata))
    stop("newdata prediction for conmode objects is not yet available; refit or evaluate npconmode() with the desired evaluation data")
  if (identical(type, "class"))
    return(object$conmode)
  probs <- object$probabilities
  if (is.null(probs))
    stop("class probabilities are not stored: fit with probabilities=TRUE")
  probs
}
gradients.conmode <- function(x, level = NULL, errors = FALSE, ...) {
  if (isTRUE(errors))
    stop("gradient standard errors are not available for conmode objects")
  gout <- x$probability.gradients
  if (is.null(gout))
    stop("class-probability gradients/effects are not available: fit with gradients=TRUE")
  stored.level <- x$probability.gradient.level
  if (!is.null(level) && !identical(as.character(level), as.character(stored.level)))
    stop(sprintf("stored class-probability gradients are for level %s; refit with level=%s to obtain that level",
                 sQuote(as.character(stored.level)),
                 sQuote(as.character(level))))
  gout
}

summary.conmode <- function(object, ...){
  cat("\nConditional Mode data: ", object$ntrain, " training points,",
      if (object$trainiseval) "" else paste(" and ", object$nobs, " evaluation points,\n", sep=""),
      " in ", object$xndim + object$yndim, " variable(s)",
      "\n(", object$yndim, " dependent variable(s), and ", object$xndim, " explanatory variable(s))\n\n",
      sep="")

  cat(genOmitStr(object))
  print(matrix(object$ybw,ncol=object$yndim,dimnames=list(paste("Dep. Var. ",object$pscaling,":",sep=""),object$ynames)))

  print(matrix(object$xbw,ncol=object$xndim,dimnames=list(paste("Exp. Var. ",object$pscaling,":",sep=""),object$xnames)))

  cat(genBwSelStr(object$bws))
  cat('\n')
  pCatGofStr(object)

  cat(genBwKerStrs(object$bws))
  cat(genTimingStr(object$bws))
  if (!is.null(object$proper.info)) {
    cat("Proper conditional probabilities: ",
        if (isTRUE(object$proper.requested)) {
          if (isTRUE(object$proper.applied)) "projected" else "already proper"
        } else {
          "not requested"
        },
        "\n", sep="")
    if (!is.null(object$proper.info$max.negative.violation.raw))
      cat("  Max negative violation (raw): ",
          signif(object$proper.info$max.negative.violation.raw, 6), "\n", sep="")
    if (!is.null(object$proper.info$max.row.sum.deviation.raw))
      cat("  Max row-sum deviation (raw): ",
          signif(object$proper.info$max.row.sum.deviation.raw, 6), "\n", sep="")
    if (!is.null(object$proper.info$repaired.rows))
      cat("  Repaired rows: ", object$proper.info$repaired.rows, "\n", sep="")
  }
  if (isTRUE(object$gradients)) {
    cat("Class-probability gradients/effects: available",
        if (!is.null(object$probability.gradient.level))
          paste0(" for level ", as.character(object$probability.gradient.level))
        else "",
        "\n", sep="")
  }
  cat('\n\n')
  
}
