test_that("npRmpi protocol tags are portable and collision-free", {
  table.fun <- getFromNamespace(".npRmpi_protocol_tag_table", "npRmpi")
  tag.fun <- getFromNamespace(".npRmpi_protocol_tag", "npRmpi")
  rank.tag.fun <- getFromNamespace(".npRmpi_protocol_rank_tag", "npRmpi")
  rank.limit <- getFromNamespace(".npRmpi_protocol_rank_limit", "npRmpi")
  tag.ub.min <- getFromNamespace(".npRmpi_protocol_tag_ub_min", "npRmpi")

  tags <- table.fun()
  expect_true(is.integer(tags) || is.numeric(tags))
  expect_true(all(!is.na(tags)))
  expect_true(all(tags < tag.ub.min))
  expect_identical(tag.fun("plreg_task"), as.integer(tags[["plreg_task"]]))

  ranks0 <- 0L:rank.limit
  ranks1 <- 1L:rank.limit
  active.tags <- c(
    vapply(ranks0, function(rank)
      rank.tag.fun("manual_bcast_base", rank, min_rank = 0L),
      integer(1L)),
    tag.fun("plreg_task"),
    tag.fun("plreg_result"),
    vapply(ranks1, function(rank)
      rank.tag.fun("scoef_req_base", rank, min_rank = 1L),
      integer(1L)),
    vapply(ranks1, function(rank)
      rank.tag.fun("scoef_res_base", rank, min_rank = 1L),
      integer(1L)),
    vapply(ranks1, function(rank)
      rank.tag.fun("scoef_ack_base", rank, min_rank = 1L),
      integer(1L)),
    vapply(ranks1, function(rank)
      rank.tag.fun("attach_ack_base", rank, min_rank = 1L),
      integer(1L)),
    vapply(ranks1, function(rank)
      rank.tag.fun("attach_release_base", rank, min_rank = 1L),
      integer(1L))
  )

  expect_true(all(active.tags < tag.ub.min))
  expect_equal(length(unique(active.tags)), length(active.tags))
  expect_false(any(active.tags >= 20000L & active.tags <= 30000L))
  expect_error(
    rank.tag.fun("manual_bcast_base", rank.limit + 1L, min_rank = 0L),
    "rank <="
  )
})

test_that("legacy high explicit MPI protocol tags are absent from active routes", {
  root <- testthat::test_path("..", "..")
  paths <- file.path(
    root,
    "R",
    c("Rcoll.R", "session.R", "np.plregression.bw.R", "np.smoothcoef.bw.R")
  )
  skip_if_not(all(file.exists(paths)))

  src <- unlist(lapply(paths, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl("50000\\s*\\+", src)))
  expect_false(any(grepl("\\b57101L\\b|\\b57102L\\b", src)))
  expect_false(any(grepl("\\b61100L\\b|\\b61200L\\b", src)))
  expect_false(any(grepl("\\b62000L\\b|\\b62100L\\b", src)))
})
