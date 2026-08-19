npRmpi_namespace_hygiene_root <- function() {
  source_root <- Sys.getenv("NP_SOURCE_ROOT", unset = "")
  if (nzchar(source_root)) {
    source_root <- tryCatch(
      normalizePath(source_root, mustWork = TRUE),
      error = function(e) ""
    )
  }
  starts <- unique(Filter(
    nzchar,
    c(
      source_root,
      tryCatch(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), error = function(e) ""),
      tryCatch(normalizePath(getwd(), mustWork = TRUE), error = function(e) ""),
      tryCatch(normalizePath(file.path(getwd(), ".."), mustWork = TRUE), error = function(e) ""),
      tryCatch(normalizePath(file.path(getwd(), "..", ".."), mustWork = TRUE), error = function(e) "")
    )
  ))

  seen <- character()
  for (start in starts) {
    candidate <- start
    repeat {
      if (!nzchar(candidate) || candidate %in% seen)
        break
      search_paths <- unique(c(candidate, file.path(candidate, "00_pkg_src", "npRmpi")))

      for (path in search_paths) {
        if (!nzchar(path) || path %in% seen)
          next
        seen <- c(seen, path)

        desc <- file.path(path, "DESCRIPTION")
        if (!file.exists(desc))
          next
        dcf <- tryCatch(read.dcf(desc), error = function(e) NULL)
        if (!is.null(dcf) &&
            nrow(dcf) &&
            identical(unname(dcf[1L, "Package"]), "npRmpi")) {
          return(path)
        }
      }

      parent <- dirname(candidate)
      if (identical(parent, candidate))
        break
      candidate <- parent
    }
  }

  stop("Could not locate npRmpi package root for namespace hygiene checks", call. = FALSE)
}

npRmpi_test_source_path <- function(...) {
  root <- tryCatch(
    npRmpi_namespace_hygiene_root(),
    error = function(e) NULL
  )
  testthat::skip_if(
    is.null(root),
    "npRmpi package source is unavailable in this test context"
  )

  path <- file.path(root, ...)
  testthat::skip_if_not(
    file.exists(path),
    paste("npRmpi package source file is unavailable:", path)
  )
  path
}

npRmpi_test_extract_c_function <- function(lines, name) {
  stopifnot(is.character(lines), length(name) == 1L, nzchar(name))
  candidates <- grep(
    paste0("\\b", name, "[[:space:]]*\\("),
    lines,
    perl = TRUE
  )

  for (start in candidates) {
    opened <- FALSE
    depth <- 0L
    rejected <- FALSE

    for (line_index in seq.int(start, length(lines))) {
      chars <- strsplit(lines[[line_index]], "", fixed = TRUE)[[1L]]
      open_at <- match("{", chars, nomatch = 0L)
      semicolon_at <- match(";", chars, nomatch = 0L)

      if (!opened) {
        if (semicolon_at > 0L && (open_at == 0L || semicolon_at < open_at)) {
          rejected <- TRUE
          break
        }
        if (open_at == 0L)
          next
        opened <- TRUE
        chars <- chars[seq.int(open_at, length(chars))]
      }

      depth <- depth + sum(chars == "{") - sum(chars == "}")
      if (depth == 0L)
        return(paste(lines[start:line_index], collapse = "\n"))
    }

    if (opened && !rejected)
      break
  }

  stop(sprintf("Could not extract C function definition: %s", name),
       call. = FALSE)
}

npRmpi_test_compact_source <- function(source) {
  gsub("[[:space:]]+", " ", source)
}

npRmpi_namespace_hygiene_scan <- function(root = npRmpi_namespace_hygiene_root(), pkg = "npRmpi") {
  root <- normalizePath(root, mustWork = TRUE)

  surface_for <- function(rel) {
    if (grepl("^tests/", rel))
      return("tests")
    if (grepl("^(R|inst|demo)/", rel))
      return("runtime")
    "other"
  }

  is_r_like_file <- function(rel) {
    grepl("\\.R$", rel) | basename(rel) %in% c("Rprofile", ".Rprofile")
  }

  all <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  all <- all[file.info(all)$isdir %in% FALSE]
  rel <- sub(paste0("^", root, "/?"), "", normalizePath(all, mustWork = TRUE))
  keep <- is_r_like_file(rel) & grepl("^(R|inst|demo|tests)/", rel)

  files <- data.frame(
    path = normalizePath(all[keep], mustWork = TRUE),
    rel = rel[keep],
    surface = vapply(rel[keep], surface_for, character(1L)),
    stringsAsFactors = FALSE
  )

  parse_ops <- function(path) {
    expr <- parse(path, keep.source = TRUE)
    pd <- getParseData(expr, includeText = TRUE)
    if (is.null(pd) || !nrow(pd))
      return(pd[FALSE, , drop = FALSE])
    pd
  }

  out <- lapply(seq_len(nrow(files)), function(i) {
    info <- files[i, , drop = FALSE]
    pd <- parse_ops(info$path)
    if (!nrow(pd))
      return(list(ns = NULL, str = NULL))

    ns_rows <- pd[pd$token %in% c("NS_GET", "NS_GET_INT"), , drop = FALSE]
    str_rows <- pd[pd$token == "STR_CONST" & grepl("::", pd$text, fixed = TRUE), , drop = FALSE]

    ns_out <- lapply(seq_len(nrow(ns_rows)), function(j) {
      row <- ns_rows[j, , drop = FALSE]
      siblings <- pd[pd$parent == row$parent, , drop = FALSE]
      pkg_row <- siblings[siblings$token == "SYMBOL_PACKAGE", , drop = FALSE]
      sym_row <- siblings[siblings$token %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL", "STR_CONST"), , drop = FALSE]
      if (!nrow(pkg_row) || !nrow(sym_row))
        return(NULL)
      data.frame(
        rel = info$rel,
        surface = info$surface,
        token = row$token,
        pkg = pkg_row$text[[1L]],
        sym = sym_row$text[[1L]],
        line = row$line1[[1L]],
        col = row$col1[[1L]],
        stringsAsFactors = FALSE
      )
    })
    ns_out <- do.call(rbind, ns_out[!vapply(ns_out, is.null, logical(1L))])

    str_out <- if (nrow(str_rows)) {
      data.frame(
        rel = info$rel,
        surface = info$surface,
        token = str_rows$token,
        pkg = "",
        sym = str_rows$text,
        line = str_rows$line1,
        col = str_rows$col1,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }

    list(ns = ns_out, str = str_out)
  })

  ns_df <- do.call(rbind, lapply(out, `[[`, "ns"))
  str_df <- do.call(rbind, lapply(out, `[[`, "str"))

  if (is.null(ns_df)) {
    ns_df <- data.frame(
      rel = character(),
      surface = character(),
      token = character(),
      pkg = character(),
      sym = character(),
      line = integer(),
      col = integer(),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(str_df)) {
    str_df <- data.frame(
      rel = character(),
      surface = character(),
      token = character(),
      pkg = character(),
      sym = character(),
      line = integer(),
      col = integer(),
      stringsAsFactors = FALSE
    )
  }

  runtime_ns <- ns_df[ns_df$surface == "runtime", , drop = FALSE]
  tests_ns <- ns_df[ns_df$surface == "tests", , drop = FALSE]

  list(
    files = files,
    runtime_same_package = runtime_ns[runtime_ns$pkg == pkg, , drop = FALSE],
    runtime_external_triple = runtime_ns[runtime_ns$pkg != pkg & runtime_ns$token == "NS_GET_INT", , drop = FALSE],
    runtime_external_double = runtime_ns[runtime_ns$pkg != pkg & runtime_ns$token == "NS_GET", , drop = FALSE],
    runtime_parser_literals = str_df[str_df$surface == "runtime", , drop = FALSE],
    tests_same_package = tests_ns[tests_ns$pkg == pkg, , drop = FALSE],
    tests_external_triple = tests_ns[tests_ns$pkg != pkg & tests_ns$token == "NS_GET_INT", , drop = FALSE]
  )
}

npRmpi_namespace_hygiene_format <- function(df) {
  if (!nrow(df))
    return("<none>")
  paste(
    sprintf(
      "%s:%d:%d %s%s%s",
      df$rel, df$line, df$col,
      ifelse(nzchar(df$pkg), paste0(df$pkg, ifelse(df$token == "NS_GET_INT", ":::", "::")), ""),
      ifelse(df$token == "STR_CONST", "", df$sym),
      ifelse(df$token == "STR_CONST", df$sym, "")
    ),
    collapse = "\n"
  )
}
