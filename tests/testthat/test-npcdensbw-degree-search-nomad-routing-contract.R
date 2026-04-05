test_that("npcdensbw NOMAD degree search fails fast when crs is unavailable", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = runif(16), y = rnorm(16))

  expect_error(
    with_nprmpi_npcdens_degree_bindings(
      list(.np_nomad_require_crs = function() stop("crs missing", call. = FALSE)),
      npcdensbw(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "coordinate",
        search.engine = "nomad",
        degree.min = 0L,
        degree.max = 1L,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L
      )
    ),
    "crs missing"
  )
})

test_that("npcdensbw NOMAD shadow search keeps collective routing in manual-broadcast context", {
  body_text <- paste(
    deparse(body(get(".npcdensbw_nomad_search", envir = asNamespace("npRmpi"), inherits = FALSE))),
    collapse = "\n"
  )

  expect_true(grepl(".npRmpi_has_active_slave_pool(comm = 1L)", body_text, fixed = TRUE))
  expect_match(body_text, "getOption\\(\"npRmpi\\.local\\.regression\\.mode\",\\s*FALSE\\)")
  expect_true(grepl("if (isTRUE(.npRmpi_autodispatch_called_from_bcast()))", body_text, fixed = TRUE))
  expect_true(grepl("search.result <- eval(mc, envir = environment())", body_text, fixed = TRUE))
  expect_true(grepl("search.result <- .npRmpi_bcast_cmd_expr(mc, comm = 1L,", body_text, fixed = TRUE))
  expect_true(grepl("caller.execute = TRUE)", body_text, fixed = TRUE))
  expect_false(grepl("!isTRUE(.npRmpi_autodispatch_called_from_bcast())", body_text, fixed = TRUE))
})
