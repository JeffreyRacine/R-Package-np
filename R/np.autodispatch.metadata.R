# Unleased unconditional results are finalized once, on the coordinator. Data
# roles establish which labels were generated; a transport-looking user column
# name is not evidence that a field should be rewritten.
.npRmpi_autodispatch_finalize_unleased_names <- function(result, mc, definition,
                                                        owner, caller) {
  fit <- inherits(result, c("npdensity", "npdistribution"))
  bandwidth <- if (fit) result[["bws", exact = TRUE]] else result
  if (!inherits(bandwidth, c("bandwidth", "dbandwidth")) ||
      !is.function(definition)) return(result)
  # External bandwidths own their names/provenance, even if a name resembles a
  # private token. Do not manufacture a new selection signature for such input.
  if (inherits(mc[["bws", exact = TRUE]], c("bandwidth", "dbandwidth")))
    return(result)
  roles <- intersect(c("dat", "tdat", "edat"),
                     intersect(names(formals(definition)), names(mc)))
  labels <- list()
  for (role in roles) {
    value <- get(role, envir = owner, inherits = FALSE)
    if (!is.atomic(value) || !is.null(dim(value))) next
    expr <- eval(substitute(substitute(ARG), list(ARG = as.name(role))), owner)
    expr <- .np_formula_dot_expression(expr, caller)
    labels[[role]] <- paste(deparse(expr), collapse = "")
  }
  training <- intersect(c("dat", "tdat"), names(labels))
  if (!length(training)) return(result)
  old <- bandwidth[["xnames", exact = TRUE]]
  role <- training[[1L]]
  if (!is.character(old) || length(old) != 1L || is.na(old) ||
      !grepl(paste0("^\\.__npRmpi_autod_", role, "_[0-9]+$"), old))
    return(result)
  descriptor <- bandwidth[["selection.procedure", exact = TRUE]]
  if (!is.null(descriptor) &&
      (!is.list(descriptor) ||
       !identical(descriptor[["signature", exact = TRUE]],
                  .np_bw_procedure_signature(bandwidth))))
    return(result)
  label <- labels[[role]]
  replace <- function(x) {
    if (is.character(x)) x[!is.na(x) & x == old] <- label
    x
  }
  index <- function(x) { names(x) <- replace(names(x)); x }
  frame <- function(x) {
    if (is.data.frame(x)) names(x) <- replace(names(x))
    else if (is.matrix(x)) colnames(x) <- replace(colnames(x))
    x
  }
  info <- function(x) {
    if (!is.list(x)) return(x)
    for (field in intersect(c("iord", "iuno", "icon", "inumord", "all.lev",
      "all.ulev", "all.dlev", "all.nlev", "all.min", "all.max",
      "all.min.next", "all.max.prev"), names(x))) x[[field]] <- index(x[[field]])
    x
  }
  bandwidth <- updateBwNameMetadata(list(xnames = label), bandwidth)
  for (field in c("icon", "iuno", "iord"))
    if (!is.null(bandwidth[[field, exact = TRUE]]))
      bandwidth[[field]] <- index(bandwidth[[field]])
  bandwidth$xdati <- info(bandwidth[["xdati", exact = TRUE]])
  if (is.list(bandwidth[["dati", exact = TRUE]]))
    bandwidth$dati$x <- info(bandwidth$dati[["x", exact = TRUE]])
  if (is.list(bandwidth[[".np.native.training", exact = TRUE]]))
    bandwidth$.np.native.training$dat <- frame(bandwidth$.np.native.training[["dat", exact = TRUE]])
  if (!is.null(descriptor)) {
    descriptor$signature$xnames <- replace(descriptor$signature$xnames)
    for (field in c("icon", "iuno", "iord"))
      descriptor$signature[[field]] <- index(descriptor$signature[[field]])
    bandwidth$selection.procedure <- descriptor
  }
  if (!fit) return(bandwidth)
  result$bws <- bandwidth
  result$xnames <- replace(result[["xnames", exact = TRUE]])
  # Only the default training evaluation inherits the training role. Explicit
  # evaluation data retain their own column names, including token collisions.
  if (isTRUE(result[["trainiseval", exact = TRUE]]))
    result$eval <- frame(result[["eval", exact = TRUE]])
  result
}
