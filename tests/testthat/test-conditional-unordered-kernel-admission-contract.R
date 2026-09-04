test_that("conditional unordered admission uses independent Y and X kernels", {
  pkg <- if ("package:npRmpi" %in% search()) "npRmpi" else "np"
  if (pkg == "npRmpi")
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")

  options(np.messages = FALSE, np.tree = FALSE)
  ns <- asNamespace(pkg)
  n <- 36L
  x <- data.frame(xu = factor(rep(letters[1:3], length.out = n)))
  y <- data.frame(yu = factor(rep(LETTERS[1:3], each = 3L,
                              length.out = n)))

  native <- function(prep, point, package) {
    on.exit(.Call("C_np_density_conditional_prepared_destroy",
                  PACKAGE = package))
    ok <- do.call(.Call,
                  c(list("C_np_density_conditional_prepared_prepare"),
                    unname(prep), list(PACKAGE = package)))
    stopifnot(isTRUE(as.logical(ok)))
    .Call("C_np_density_conditional_prepared_eval_raw",
          as.double(point), prep$degree, PACKAGE = package)[1L]
  }
  if (pkg == "npRmpi")
    mpi.bcast.Robj2slave(native)

  evaluate <- function(uykertype, uxkertype, point, bwmethod, bwtype) {
    bw <- do.call(
      get("npcdensbw", ns),
      list(xdat = x, ydat = y, bws = c(.2, .2),
           bandwidth.compute = FALSE, bwmethod = bwmethod,
           bwtype = bwtype, regtype = "lc",
           uykertype = uykertype, uxkertype = uxkertype)
    )
    prep <- get(".npcdensbw_prepared_prepare_args", ns)(
      x, y, bw, invalid.penalty = "dbmax"
    )
    prep$rbw <- as.double(point)
    prep$myopti[27L] <- 0L

    if (pkg == "np")
      return(native(prep, point, pkg))
    mpi.bcast.Robj2slave(prep)
    mpi.bcast.Robj2slave(point)
    mpi.bcast.cmd(native(prep, point, "npRmpi"), caller.execute = TRUE)
  }

  for (bwmethod in c("cv.ml", "cv.ls"))
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      valid <- evaluate("liracine", "aitchisonaitken", c(.75, .30),
                        bwmethod, bwtype)
      invalid <- evaluate("aitchisonaitken", "liracine", c(.75, .30),
                          bwmethod, bwtype)
      expect_true(is.finite(valid) && abs(valid) < 1e150)
      expect_true(!is.finite(invalid) || abs(invalid) >= 1e150)
    }
})
