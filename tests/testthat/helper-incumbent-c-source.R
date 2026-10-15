# Static contract view of the incumbent definitions at their actual include sites.
# This is deliberately not a C preprocessor or proof of conditional-clone behavior.
np_test_read_incumbent_c_source <- function(con, ...) {
  lines <- base::readLines(con, ...)
  if (!is.character(con) || length(con) != 1L ||
      !basename(con) %in% c("np.c", "jksum.c"))
    return(lines)

  boundary <- if (basename(con) == "np.c")
    "#define NP_CF316_LITERAL_ENTRIES 1" else "#define NP_CF316_LITERAL 1"
  cut <- which(lines == boundary)
  if (length(cut) != 1L)
    stop("Expected one explicit conditional clone boundary", call. = FALSE)
  pattern <- '^#include "(conditional_literal/[A-Za-z0-9_]+[.]inc)"$'
  include <- grep(pattern, lines)
  before <- include[include < cut]
  after <- include[include > cut]
  names.before <- sub(pattern, "\\1", lines[before])
  names.after <- sub(pattern, "\\1", lines[after])
  extra <- if (basename(con) == "jksum.c") paste0(
    "conditional_literal/jksum_np_plot_gate_",
    c("prepare", "body", "dispatch_native"), ".inc") else character()
  if (!length(before) || anyDuplicated(names.before) || anyDuplicated(names.after) ||
      !setequal(names.after, c(names.before, extra)) ||
      any(grepl('^[[:space:]]*#[[:space:]]*include.*conditional_literal/', lines[-include])))
    stop("Conditional fragment include contract changed", call. = FALSE)

  # Keep the clone wiring visible, without counting macro-renamed bodies as
  # duplicate incumbent definitions. Existing raw-source wiring tests stay raw.
  out <- as.list(lines)
  source.dir <- normalizePath(dirname(con), mustWork = TRUE)
  for (i in seq_along(before)) {
    fragment <- normalizePath(file.path(source.dir, names.before[[i]]),
                              mustWork = TRUE)
    if (!startsWith(fragment, paste0(source.dir, "/conditional_literal/")))
      stop("Conditional fragment escaped source directory", call. = FALSE)
    body <- base::readLines(fragment, ...)
    if (any(grepl('^[[:space:]]*#[[:space:]]*include', body)))
      stop("Nested fragment include requires explicit review", call. = FALSE)
    out[[before[[i]]]] <- body
  }
  unlist(out, use.names = FALSE)
}
