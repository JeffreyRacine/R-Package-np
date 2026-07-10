.npRmpi_protocol_rank_limit <- 2047L
.npRmpi_protocol_tag_ub_min <- 32767L

.npRmpi_protocol_tag_table <- function() {
  c(
    manual_bcast_base = 4096L,
    plreg_task = 6144L,
    plreg_result = 6145L,
    scoef_req_base = 8192L,
    scoef_res_base = 10240L,
    scoef_ack_base = 12288L,
    attach_ack_base = 14336L,
    attach_release_base = 16384L
  )
}

.npRmpi_protocol_tag <- function(name) {
  table <- .npRmpi_protocol_tag_table()
  if (!is.character(name) || length(name) != 1L || !nzchar(name) ||
      !(name %in% names(table)))
    stop("unknown npRmpi MPI protocol tag name", call. = FALSE)
  as.integer(table[[name]])
}

.npRmpi_protocol_rank_tag <- function(name,
                                      rank,
                                      min_rank = 0L,
                                      where = "npRmpi MPI protocol") {
  base <- .npRmpi_protocol_tag(name)
  rank <- suppressWarnings(as.integer(rank)[1L])
  min_rank <- suppressWarnings(as.integer(min_rank)[1L])
  if (!is.finite(min_rank))
    min_rank <- 0L
  if (!is.finite(rank) || rank < min_rank) {
    stop(sprintf("%s received invalid MPI rank for protocol tag", where),
         call. = FALSE)
  }
  if (rank > .npRmpi_protocol_rank_limit) {
    stop(sprintf(
      "%s requires MPI rank <= %d for portable protocol tags",
      where,
      .npRmpi_protocol_rank_limit
    ), call. = FALSE)
  }
  tag <- base + rank
  if (!is.finite(tag) || tag > .npRmpi_protocol_tag_ub_min) {
    stop(sprintf("%s generated non-portable MPI tag %s", where, tag),
         call. = FALSE)
  }
  as.integer(tag)
}
