npcdisthat <- function(bws,
                       txdat = stop("training data 'txdat' missing"),
                       tydat = stop("training data 'tydat' missing"),
                       exdat,
                       eydat,
                       y = NULL,
                       output = c("matrix", "apply")) {
  if (!inherits(bws, "condbandwidth")) {
    stop("argument 'bws' must inherit from class 'condbandwidth' in npcdisthat()")
  }

  exdat.arg <- NULL
  eydat.arg <- NULL
  if (!missing(exdat))
    exdat.arg <- exdat
  if (!missing(eydat))
    eydat.arg <- eydat

  .npcdhat_core(
    bws = bws,
    txdat = txdat,
    tydat = tydat,
    exdat = exdat.arg,
    eydat = eydat.arg,
    y = y,
    output = output,
    operator = "integral",
    class_name = "npcdisthat",
    where = "npcdisthat"
  )
}
