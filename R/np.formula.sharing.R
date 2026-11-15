# Reuse exact ordinary columns after the existing model-frame transaction.
# Retain the frame header, NA-policy result and copy-on-modify ownership.
.np_formula_share_prepared_columns <- function(frame, aligned) {
  if (isS4(frame) || !identical(class(frame), "data.frame") ||
      typeof(aligned) != "list" || is.object(aligned) || isS4(aligned))
    return(frame)
  columns <- frame
  attributes(columns) <- NULL
  plain <- function(x) {
    if (isS4(x) || !typeof(x) %in% c("double", "integer", "logical", "raw"))
      return(FALSE)
    a <- attributes(x)
    cl <- attr(x, "class", exact = TRUE)
    if (is.null(cl))
      return(all(names(a) %in% "names") &&
               all(vapply(a, function(v)
                 is.character(v) && is.null(attributes(v)), logical(1L))))
    if (typeof(x) != "integer" ||
        (!identical(cl, "factor") && !identical(cl, c("ordered", "factor"))))
      return(FALSE)
    all(names(a) %in% c("names", "levels", "class")) &&
      all(vapply(a, function(v)
        is.character(v) && is.null(attributes(v)), logical(1L)))
  }
  encoding <- function(x)
    lapply(attributes(x), function(v) if (is.character(v)) Encoding(v) else NULL)
  selected <- integer()
  for (i in seq_len(min(length(columns), length(aligned)))) {
    a <- .subset2(frame, i)
    b <- .subset2(aligned, i)
    if (plain(a) && plain(b) &&
        identical(a, b, num.eq = FALSE, single.NA = FALSE,
                  attrib.as.set = FALSE) &&
        identical(encoding(a), encoding(b)))
      selected <- c(selected, i)
  }
  if (!length(selected)) return(frame)
  a <- attributes(frame)
  a[["row.names"]] <- .row_names_info(frame, 0L)
  value <- frame
  attributes(value) <- NULL
  for (i in selected) value[[i]] <- .subset2(aligned, i)
  attributes(value) <- a
  value
}
