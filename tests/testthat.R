library(testthat)
library(npRmpi)

np_check_filter <- Sys.getenv("NP_CHECK_FILTER", unset = "check-minimal")
np_check_full <- identical(Sys.getenv("NP_CHECK_FULL", unset = ""), "1")
np_full_nslaves <- suppressWarnings(as.integer(
  Sys.getenv("NP_RMPI_TEST_NSLAVES", unset = "1")
))
np_full_shard_size <- suppressWarnings(as.integer(
  Sys.getenv("NP_RMPI_TEST_SHARD_SIZE", unset = "40")
))

local({
  test_error <- NULL

  if (isTRUE(np_check_full)) {
    if (length(np_full_nslaves) != 1L || is.na(np_full_nslaves) ||
        np_full_nslaves < 1L) {
      test_error <- simpleError(
        "NP_RMPI_TEST_NSLAVES must be a positive integer for the full suite"
      )
    } else if (length(np_full_shard_size) != 1L ||
               is.na(np_full_shard_size) || np_full_shard_size < 1L) {
      test_error <- simpleError(
        "NP_RMPI_TEST_SHARD_SIZE must be a positive integer for the full suite"
      )
    } else {
      test_error <- tryCatch({
        test_files <- sort(list.files(
          "testthat", pattern = "^test-.*\\.[rR]$", full.names = FALSE
        ))
        if (!length(test_files) || anyDuplicated(test_files)) {
          stop("npRmpi full-suite runner found no unique test files")
        }

        shard_count <- as.integer(ceiling(length(test_files) / np_full_shard_size))
        witness_dir <- normalizePath(
          file.path("testthat", "_npRmpi_full_shard_witnesses"),
          mustWork = FALSE
        )
        dir.create(witness_dir, recursive = TRUE, showWarnings = FALSE)
        shard_runner <- normalizePath(
          file.path("validation", "run_full_mpi_test_shard.R"),
          mustWork = TRUE
        )
        test_dir <- normalizePath("testthat", mustWork = TRUE)
        rscript <- file.path(R.home("bin"), "Rscript")
        setsid <- unname(Sys.which("setsid"))
        statuses <- integer(shard_count)

        for (shard in seq_len(shard_count)) {
          witness <- file.path(
            witness_dir,
            sprintf("shard-%03d-of-%03d.ok", shard, shard_count)
          )
          if (file.exists(witness)) unlink(witness)

          env <- c(
            paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)),
            "_R_CHECK_PACKAGE_NAME_=npRmpi",
            "NP_RMPI_TEST_SUITE_POOL=1",
            paste0("NOT_CRAN=", Sys.getenv("NOT_CRAN", unset = "true"))
          )
          shard_command <- rscript
          shard_args <- c(
            "--vanilla", shard_runner,
            as.character(shard), as.character(shard_count),
            as.character(np_full_shard_size), as.character(np_full_nslaves),
            shQuote(witness), shQuote(test_dir)
          )
          if (nzchar(setsid)) {
            # Some MPI launchers terminate the shard's process group during
            # their known teardown-137 path.  Give each shard its own session
            # so that teardown cannot terminate the suite orchestrator.
            shard_command <- setsid
            shard_args <- c(shQuote(rscript), shard_args)
          }
          statuses[[shard]] <- system2(
            shard_command,
            args = shard_args,
            stdout = "", stderr = "", env = env
          )

          expected <- sprintf(
            "NP_RMPI_FULL_SHARD_OK %d/%d", shard, shard_count
          )
          witnessed <- file.exists(witness) && identical(
            readLines(witness, warn = FALSE), expected
          )
          if (!witnessed || !statuses[[shard]] %in% c(0L, 137L)) {
            statuses[[shard]] <- 1L
          }
        }

        if (any(statuses != 0L & statuses != 137L)) {
          failed <- which(statuses != 0L & statuses != 137L)
          stop(
            "npRmpi full-suite shard failure(s): ",
            paste(failed, collapse = ", ")
          )
        }
        message(
          "NP_RMPI_FULL_SUITE_OK ", length(test_files), " files in ",
          shard_count, " fresh-pool shard(s)"
        )
        NULL
      }, error = identity)
    }
  }

  if (is.null(test_error) && !isTRUE(np_check_full)) {
    tryCatch({
      if (identical(np_check_filter, "")) {
        test_check("npRmpi")
      } else {
        test_check("npRmpi", filter = np_check_filter)
      }
    }, error = function(e) {
      test_error <<- e
    })
  }

  if (!is.null(test_error)) {
    message(conditionMessage(test_error))
    quit(save = "no", status = 1, runLast = FALSE)
  }
})
