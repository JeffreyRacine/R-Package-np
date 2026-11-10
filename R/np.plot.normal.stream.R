# Explicit seeded plot normals. Never run caller code, load this package, or
# initialize MPI in the child. Box-Muller's hidden cache stays in the caller.
.np_plot_normal_request <- function(seed, kind, count) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max)
    stop("invalid plot normal seed", call. = FALSE)
  uniform <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister", "Knuth-TAOCP", "Knuth-TAOCP-2002", "L'Ecuyer-CMRG")
  if (!.np_plot_rng_kind_shape(kind, 2L) ||
      !kind[1L] %in% uniform || kind[2L] != "Box-Muller")
    stop("seeded plot Box-Muller draws require a built-in uniform generator and supported RNG kinds",
         call. = FALSE)
  if (!is.numeric(count) || length(count) != 1L || !is.finite(count) ||
      count < 0 || count != floor(count) || count > floor((2^53 - 1)/8))
    stop("invalid plot normal draw count", call. = FALSE)
  list(seed = as.integer(seed), kind = kind, count = as.double(count),
       version = R.version$version.string)
}

.np_plot_normal_child <- function() {
  args <- commandArgs(TRUE)
  request <- readRDS(args[1L])
  stopifnot(identical(request$version, R.version$version.string))
  do.call(RNGkind, as.list(request$kind))
  set.seed(request$seed)
  entry <- .Random.seed
  output <- file(args[2L], "wb")
  on.exit(close(output), add = TRUE)
  remaining <- request$count
  while (remaining > 0) {
    count <- min(remaining, 65536L)
    writeBin(stats::rnorm(count), output, size = 8L, endian = "little")
    remaining <- remaining - count
  }
  flush(output)
  saveRDS(list(kind = RNGkind(), count = request$count,
    version = R.version$version.string, entry = entry, end = .Random.seed,
    pid = Sys.getpid(), namespaces = loadedNamespaces()), args[3L])
}

.np_plot_normal_run <- function(script, request, output, receipt, log, timeout = 600L) {
  command <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
    "Rscript.exe" else "Rscript")
  if (!file.exists(command)) stop("Rscript is unavailable for seeded plot normals", call. = FALSE)
  system2(command, c("--vanilla", "--default-packages=stats",
    shQuote(c(script, request, output, receipt))), stdout = log, stderr = log,
    timeout = timeout)
}

.np_plot_normal_stream <- function(seed, kind, count) {
  request <- .np_plot_normal_request(seed, kind, count)
  stream <- new.env(parent = emptyenv())
  stream$remaining <- request$count; stream$closed <- FALSE
  stream$path <- NULL; stream$connection <- NULL
  stream$info <- list(kind = kind, count = request$count, pid = NULL)
  if (request$count == 0) return(stream)
  path <- tempfile("np-plot-normals-")
  if (!dir.create(path, mode = "0700")) stop("cannot create plot normal workspace", call. = FALSE)
  stream$path <- path
  complete <- FALSE
  on.exit(if (!complete) .np_plot_normal_close(stream), add = TRUE)
  files <- file.path(path, c("sampler.R", "request.rds", "normals.bin", "receipt.rds", "sampler.log"))
  writeLines(c("(", deparse(.np_plot_normal_child, width.cutoff = 500L), ")()"), files[1L])
  saveRDS(request, files[2L])
  status <- .np_plot_normal_run(files[1L], files[2L], files[3L], files[4L], files[5L])
  if (!identical(as.integer(status), 0L)) {
    detail <- if (file.exists(files[5L])) tail(readLines(files[5L], warn = FALSE), 8L) else character()
    stop(paste(c(sprintf("seeded plot normal sampler failed (status %s)", status), detail),
               collapse = "\n"), call. = FALSE)
  }
  if (!file.exists(files[4L]) || !file.exists(files[3L]) ||
      !identical(as.double(file.info(files[3L])$size), 8 * request$count))
    stop("incomplete seeded plot normal output", call. = FALSE)
  info <- readRDS(files[4L])
  if (!identical(info$kind, request$kind) || !identical(info$count, request$count) ||
      !identical(info$version, request$version) ||
      !is.integer(info$pid) || length(info$pid) != 1L || info$pid == Sys.getpid())
    stop("invalid seeded plot normal receipt", call. = FALSE)
  stream$connection <- file(files[3L], "rb")
  stream$info <- info
  complete <- TRUE
  stream
}

.np_plot_normal_read <- function(stream, count) {
  if (stream$closed || !is.numeric(count) || length(count) != 1L ||
      !is.finite(count) || count < 0 || count != floor(count) ||
      count > stream$remaining)
    stop("invalid plot normal stream read", call. = FALSE)
  if (count == 0) return(numeric())
  values <- readBin(stream$connection, "double", n = count, size = 8L, endian = "little")
  if (length(values) != count || any(!is.finite(values)))
    stop("incomplete or nonfinite plot normal stream", call. = FALSE)
  stream$remaining <- stream$remaining - count
  values
}

.np_plot_normal_close <- function(stream) {
  if (!stream$closed) {
    stream$closed <- TRUE
    tryCatch(if (!is.null(stream$connection)) close(stream$connection),
      finally = if (!is.null(stream$path)) unlink(stream$path, recursive = TRUE))
  }
  invisible(NULL)
}
