conditional_layout_fixture <- function() {
  n <- 24L
  list(x = data.frame(x = seq(.05, .95, length.out = n),
                      u = factor(rep(letters[1:2], length.out = n)),
                      o = ordered(rep(1:4, each = 2L, length.out = n))),
       y = data.frame(y = cos(seq_len(n) / 5),
                      u = factor(rep(LETTERS[1:3], each = 3L, length.out = n)),
                      o = ordered(rep(1:5, length.out = n))))
}

conditional_layout_native <- function(prep, point) {
  on.exit(.Call("C_np_density_conditional_prepared_destroy", PACKAGE = "npRmpi"))
  ok <- do.call(.Call, c(list("C_np_density_conditional_prepared_prepare"),
                        unname(prep), list(PACKAGE = "npRmpi")))
  stopifnot(isTRUE(as.logical(ok)))
  .Call("C_np_density_conditional_prepared_eval_raw",
        as.double(point), prep$degree, PACKAGE = "npRmpi")[[1L]]
}

test_that("conditional transforms preserve the independent Yu Yo Xu Xo map", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- conditional_layout_fixture()
  for (full in c(FALSE, TRUE)) {
    x <- if (full) d$x else d$x[c("x", "u")]
    y <- if (full) d$y else d$y["o"]
    for (ky in c("aitchisonaitken", "liracine"))
      for (kx in c("aitchisonaitken", "liracine")) {
        bw <- npcdensbw(xdat = x, ydat = y, bws = rep(.25, ncol(x) + ncol(y)),
          bandwidth.compute = FALSE, regtype = "lc",
          uykertype = ky, uxkertype = kx)
        prep <- npRmpi:::.npcdensbw_prepared_prepare_args(
          x, y, bw, invalid.penalty = "dbmax")
        ## Independent caps in raw decoder order, not the transform's loop order.
        caps <- c(rep(1, bw$ncon),
          if (any(bw$iyuno)) vapply(y[bw$iyuno], function(z)
            if (ky == "liracine") 1 else (nlevels(z) - 1) / nlevels(z), 0),
          rep(1, sum(bw$iyord)),
          if (any(bw$ixuno)) vapply(x[bw$ixuno], function(z)
            if (kx == "liracine") 1 else (nlevels(z) - 1) / nlevels(z), 0),
          rep(1, sum(bw$ixord)))
        for (fraction in c(.2, .6)) {
          raw <- fraction * caps
          if (bw$ncon) raw[seq_len(bw$ncon)] <- .6
          transformed <- qlogis(raw / caps)
          if (bw$ncon) transformed[seq_len(bw$ncon)] <- log(raw[seq_len(bw$ncon)])
          prep$rbw <- as.double(raw)
          ## Private native-ingress contract (CBW_TBNDI), not a public NOMAD option.
          prep$myopti[27L] <- 0L
          direct <- {
            mpi.bcast.Robj2slave(conditional_layout_native)
            mpi.bcast.Robj2slave(prep); mpi.bcast.Robj2slave(raw)
            mpi.bcast.cmd(conditional_layout_native(prep, raw), caller.execute = TRUE)
          }
          prep$myopti[27L] <- 1L
          actual <- {
            mpi.bcast.Robj2slave(conditional_layout_native)
            mpi.bcast.Robj2slave(prep); mpi.bcast.Robj2slave(transformed)
            mpi.bcast.cmd(conditional_layout_native(prep, transformed), caller.execute = TRUE)
          }
          expect_true(is.finite(direct) && abs(direct) < 1e150)
          expect_equal(actual, direct, tolerance = 1e-11)
        }
      }
  }
})

test_that("conditional ordinary transformed searches replay on the raw scale", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- conditional_layout_fixture()
  x <- d$x[c("x", "u")]; y <- d$y["o"]
  for (family in c("npcdensbw", "npcdistbw")) {
    npseed(9404)
    args <- list(xdat = x, ydat = y, regtype = "lc", bws = c(.9, .3, .2),
                 nmulti = 1L, itmax = 1L, powell.remin = FALSE,
                 transform.bounds = TRUE)
    if (family == "npcdistbw") args$ngrid <- 5L
    bw <- do.call(getFromNamespace(family, "npRmpi"), args)
    expect_true(bw$ybw[[1L]] >= 0 && bw$ybw[[1L]] <= 1)
    expect_true(bw$xbw[[2L]] >= 0 && bw$xbw[[2L]] <= .5)
    replay <- list(xdat = x, ydat = y, bws = bw, invalid.penalty = "dbmax")
    if (family == "npcdistbw") replay$ngrid <- 5L
    raw <- do.call(getFromNamespace(paste0(".", family, "_eval_only"), "npRmpi"),
                   replay)
    expect_true(is.finite(bw$fval) && abs(bw$fval) < 1e150)
    expect_equal(raw$objective, bw$fval, tolerance = 5e-7)
  }
})
