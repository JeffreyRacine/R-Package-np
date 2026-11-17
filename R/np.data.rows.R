# Retain column storage for an ordinary all-row extraction. R's data-frame
# extraction changes the automatic row-name record even when no row is dropped,
# so only the columns may be shared; the result still needs its own header.
.np_data_plain_column <- function(x) {
  typeof(x) %in% c("double", "integer", "logical", "character", "raw", "complex") &&
    !is.object(x) && !isS4(x) && all(names(attributes(x)) %in% "names")
}

.np_data_plain_frame <- function(x) {
  !isS4(x) && identical(class(x), "data.frame") &&
    all(names(attributes(x)) %in% c("names", "row.names", "class")) &&
    all(vapply(x, function(column) {
      .np_data_plain_column(column) ||
        (!isS4(column) && typeof(column) == "integer" &&
         identical(names(attributes(column)), c("levels", "class")) &&
         is.character(attr(column, "levels")) &&
         (identical(class(column), "factor") ||
          identical(class(column), c("ordered", "factor"))))
    }, logical(1L)))
}

# Test for a workspace method without forcing an active binding. The original
# extraction owns method lookup and any user side effects on the fallback path.
.np_data_custom_method <- function(generic, class) {
  exists(paste0(generic, ".", class), envir = globalenv(), inherits = FALSE) ||
    !is.null(getS3method(generic, class, optional = TRUE, envir = emptyenv()))
}

.np_data_rows <- function(x, rows, frame = TRUE) {
  plain <- if (frame) .np_data_plain_frame(x) else .np_data_plain_column(x)
  if (plain && typeof(rows) == "logical" && is.null(attributes(rows))) {
    n <- if (frame) .row_names_info(x, 2L) else length(x)
    if (length(rows) == n && isTRUE(all(rows))) {
      if (!frame) return(x)
      factor.columns <- any(vapply(x, is.factor, logical(1L)))
      ordered.columns <- any(vapply(x, is.ordered, logical(1L)))
      if (factor.columns &&
          .np_data_custom_method("dim", "factor")) return(x[rows, , drop = FALSE])
      if (ordered.columns &&
          (.np_data_custom_method("[", "ordered") ||
           .np_data_custom_method("dim", "ordered"))) return(x[rows, , drop = FALSE])
      rn <- attr(x, "row.names")
      if (identical(getS3method("[", "data.frame", envir = emptyenv()),
                    get("[.data.frame", baseenv())) &&
          !.np_data_custom_method("names", "data.frame") &&
          !.np_data_custom_method("length", "data.frame") &&
          all(vapply(x, function(column) length(unclass(column)) == n,
                     logical(1L)))) {
        record <- .row_names_info(x, 0L)
        compact <- is.integer(record) && length(record) == 2L &&
          is.na(record[1L]) && !is.na(record[2L])
        valid <- is.null(attributes(record)) &&
          (compact || (!anyNA(record) && !anyDuplicated(record)))
        if (!valid) return(x[rows, , drop = FALSE])
        value <- x
        oldClass(value) <- NULL
        attr(value, "row.names") <- NULL
        attr(value, "row.names") <- rn
        oldClass(value) <- "data.frame"
        return(value)
      }
    }
  }
  if (frame) x[rows, , drop = FALSE] else x[rows]
}
