npRmpi_full_test_record_status <- function(witness, shard, shard_count, status) {
  scalar_integer <- function(value, minimum) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= minimum && value <= .Machine$integer.max &&
      value == floor(value)
  }
  if (!is.character(witness) || length(witness) != 1L ||
      is.na(witness) || !nzchar(witness) ||
      !scalar_integer(shard, 1L) || !scalar_integer(shard_count, 1L) ||
      shard > shard_count || !scalar_integer(status, 0L)) {
    stop("invalid npRmpi full-suite raw-status identity")
  }
  receipt <- paste0(witness, ".raw-status")
  temporary <- paste0(receipt, ".tmp")
  line <- sprintf("NP_RMPI_FULL_SHARD_RAW_STATUS %d/%d status=%d",
                  shard, shard_count, status)
  writeLines(line, temporary, useBytes = TRUE)
  if (!file.rename(temporary, receipt))
    stop("could not publish npRmpi full-suite raw status")
  cat(line, "\n", sep = "")
  flush.console()
  invisible(NULL)
}
