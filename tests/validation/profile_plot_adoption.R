# Opt-in external lifecycle sentinel; not a default testthat/CRAN test.
# Run each mode in a fresh process. Create a fresh output directory first.
# Set R_LIBS to the exact private installed library and TMPDIR to scratch.
# Profile (one/two workers = two/three ranks):
# R_PROFILE='' R_PROFILE_USER="$R_LIBS/npRmpi/Rprofile" \
#   mpiexec -n 2 Rscript --no-save profile_plot_adoption.R profile 1 OUTPUT
# Session: R_PROFILE='' R_PROFILE_USER='' Rscript --vanilla \
#   profile_plot_adoption.R session 1 OUTPUT
# Attach: R_PROFILE='' R_PROFILE_USER='' mpiexec -n 2 Rscript --vanilla \
#   profile_plot_adoption.R attach 1 OUTPUT
# One-rank refusal: the profile launch above with -n 1 and arguments
#   single 0 OUTPUT. It must refuse pool-dependent plotting, not adopt/spawn.
# Acceptance requires RESULT.rds, PROFILE_PLOT_SENTINEL_PASS, raw zero,
# no timeout/forced cleanup, and independently verified owned-process absence.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[1L] %in% c("profile", "session", "attach", "single"))
  stop("usage: profile_plot_adoption.R <profile|session|attach|single> <workers> <fresh-output-dir>")
mode <- args[1L]
workers <- suppressWarnings(as.integer(args[2L]))
if (is.na(workers) || !workers %in% if (mode == "single") 0L else c(1L, 2L))
  stop("workers must be 1 or 2 (0 only for the one-rank refusal control)")
output <- normalizePath(args[3L], mustWork = TRUE)
if (file.exists(file.path(output, "RESULT.rds")))
  stop("refusing to overwrite an existing result")
suppressPackageStartupMessages(library(npRmpi))
ns <- asNamespace("npRmpi")
options(np.messages = FALSE, np.plot.progress = FALSE)
master <- TRUE
if (mode == "attach") master <- isTRUE(npRmpi.init(mode = "attach", quiet = TRUE))
if (mode == "session") stopifnot(isTRUE(npRmpi.init(nslaves = workers, quiet = TRUE)))

if (master) {
  refusals <- list()
  expect_refusal <- function(expr, key, message) {
    answer <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(answer, "error"), grepl(message, conditionMessage(answer), fixed = TRUE))
    refusals[[key]] <<- conditionMessage(answer)
  }
  require.pool <- get(".npRmpi_require_active_slave_pool", ns)
  if (mode %in% c("profile", "single")) {
    stopifnot(isTRUE(getOption("npRmpi.profile.active")),
              isTRUE(getOption("npRmpi.mpi.initialized")),
              !isTRUE(getOption("npRmpi.pool.active")))
    expect_refusal(require.pool(where = "profile plot sentinel"), "unadopted",
                   "requires an active MPI slave pool")
  }

  if (mode == "single") {
    stopifnot(mpi.comm.size(0L) == 1L)
    saveRDS(list(mode = mode, workers = workers, refusals = refusals,
                 outputs = list(), library = find.package("npRmpi")),
            file.path(output, "RESULT.rds"))
    cat("PROFILE_PLOT_SENTINEL_PASS single expected-refusal\n")
    flush.console()
    mpi.quit()
  } else {
    stopifnot(mpi.comm.rank(1L) == 0L, mpi.comm.size(1L) == workers + 1L)
    # Existing transport helper publishes this test artifact directory only.
    mpi.bcast.Robj2slave(output)
    mpi.bcast.cmd({
      np.mpi.initialize()
      options(np.messages = FALSE, np.plot.progress = FALSE)
      RNGkind("Mersenne-Twister", "Inversion"); set.seed(42)
      saveRDS(list(pid = Sys.getpid(), rank = mpi.comm.rank(), size = mpi.comm.size(),
                   library = find.package("npRmpi")),
              file.path(output, paste0("before-", mpi.comm.rank(), ".rds")))
      mpi.barrier()
    }, caller.execute = TRUE)

    if (mode == "profile") {
      # A deliberately unleased negative fixture reproduces the original call
      # before adoption. It is discarded before the positive user recipe.
      mpi.bcast.cmd({
        profile.x <- data.frame(x = seq(.08, .92, length.out = 20L))
        profile.fixture <- npudens(tdat = profile.x, bws = .15,
          bandwidth.compute = FALSE, bwscaling = FALSE, ckertype = "beta",
          ckerbound = "fixed", ckerlb = 0, ckerub = 1, se = FALSE)
      }, caller.execute = TRUE)
      expect_refusal(plot(profile.fixture, errors = "none", output = "data", neval = 3L),
        "profile_master", "requires an active MPI slave pool")
      expect_refusal(plot(profile.fixture, errors = "asymptotic", center = "bias-corrected",
        band = "pmzsd", output = "data", neval = 3L),
        "profile_master_pmzsd", "requires an active MPI slave pool")
      expect_refusal(mpi.bcast.cmd(plot(profile.fixture), caller.execute = TRUE),
        "direct_broadcast", "plot(...) inside mpi.bcast.cmd(...) is unsupported")
      mpi.bcast.cmd(rm(profile.fixture, profile.x), caller.execute = TRUE)
      # Explicit public adoption; do not set internal options or use auto/attach.
      stopifnot(isTRUE(npRmpi.init(mode = "spawn", nslaves = workers, quiet = TRUE)))
    }
    stopifnot(isTRUE(getOption("npRmpi.pool.active")), mpi.comm.size(1L) == workers + 1L)
    mpi.bcast.cmd({
      RNGkind("Mersenne-Twister", "Inversion"); set.seed(42)
      saveRDS(list(pid = Sys.getpid(), rank = mpi.comm.rank(), size = mpi.comm.size(),
                   library = find.package("npRmpi")),
              file.path(output, paste0("after-", mpi.comm.rank(), ".rds")))
      mpi.barrier()
    }, caller.execute = TRUE)
    identities <- lapply(0:workers, function(rank) {
      before <- readRDS(file.path(output, paste0("before-", rank, ".rds")))
      after <- readRDS(file.path(output, paste0("after-", rank, ".rds")))
      stopifnot(identical(before, after), identical(after$library, find.package("npRmpi")))
      after
    })

    values <- list()
    x <- data.frame(x = seq(.08, .92, length.out = 20L))
    for (kernel in c("beta", "gaussian")) {
      model <- npudens(tdat = x, bws = .15, bandwidth.compute = FALSE,
        bwscaling = FALSE, ckertype = kernel, ckerbound = "fixed",
        ckerlb = 0, ckerub = 1, se = FALSE)
      original <- serialize(model, NULL)
      for (kind in c("none", "pmzsd")) {
        seed <- .Random.seed
        value <- if (kind == "none") {
          plot(model, errors = "none", output = "data", neval = 3L)
        } else {
          plot(model, errors = "asymptotic", center = "bias-corrected",
               band = "pmzsd", output = "data", neval = 3L)
        }
        frame <- as.data.frame(value)
        stopifnot(nrow(frame) == 3L, identical(original, serialize(model, NULL)),
                  identical(seed, .Random.seed))
        if (kind == "pmzsd") stopifnot(all(is.finite(frame$lower)),
          all(is.finite(frame$upper)), all(frame$lower <= frame$upper))
        values[[paste(kernel, kind, sep = "_")]] <- frame
      }
    }

    # A fresh real collective fit proves reuse. All ranks enter the barrier
    # within the same broadcast; a master-only barrier would deadlock here.
    mpi.bcast.cmd({
      reuse.x <- data.frame(x = seq(.08, .92, length.out = 20L))
      reuse.fit <- npudens(tdat = reuse.x, bws = .15, bandwidth.compute = FALSE,
        bwscaling = FALSE, ckertype = "gaussian", ckerbound = "fixed",
        ckerlb = 0, ckerub = 1, se = FALSE)
      stopifnot(all(is.finite(fitted(reuse.fit))))
      saveRDS(list(pid = Sys.getpid(), rank = mpi.comm.rank(), value = fitted(reuse.fit)),
              file.path(output, paste0("reuse-", mpi.comm.rank(), ".rds")))
      mpi.barrier()
    }, caller.execute = TRUE)
    reuse <- lapply(0:workers, function(rank) readRDS(file.path(output, paste0("reuse-", rank, ".rds"))))
    for (i in seq_along(reuse)) stopifnot(identical(reuse[[i]]$pid, identities[[i]]$pid),
      isTRUE(all.equal(reuse[[i]]$value, reuse[[1L]]$value, tolerance = 1e-12)))
    if (mode == "session") {
      npRmpi.quit(force = TRUE)
      expect_refusal(plot(model, errors = "none", output = "data", neval = 3L),
        "closed_pool", "requires an active MPI slave pool")
    }
    saveRDS(list(mode = mode, workers = workers, identities = identities,
      refusals = refusals, outputs = values, reuse = reuse,
      library = find.package("npRmpi")), file.path(output, "RESULT.rds"))
    cat("PROFILE_PLOT_SENTINEL_PASS", mode, "workers", workers, "plots", length(values), "\n")
    flush.console()
    if (mode == "profile") mpi.bcast.cmd(mpi.quit(), caller.execute = TRUE)
    if (mode == "attach") npRmpi.quit(mode = "attach")
  }
}
