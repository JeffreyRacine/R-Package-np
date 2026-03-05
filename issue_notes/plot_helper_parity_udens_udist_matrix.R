np_parity_load_package <- function() {
  use.installed <- tolower(Sys.getenv("NP_PLOT_PARITY_USE_INSTALLED", "0")) %in% c("1", "true", "yes")
  if (use.installed) {
    suppressPackageStartupMessages(library(npRmpi))
    return(invisible(NULL))
  }

  file.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script.dir <- if (length(file.arg)) {
    dirname(normalizePath(sub("^--file=", "", file.arg[1L]), mustWork = FALSE))
  } else {
    getwd()
  }
  repo.root <- normalizePath(file.path(script.dir, ".."), mustWork = FALSE)
  suppressPackageStartupMessages(library(devtools))
  load_all(repo.root, quiet = TRUE)
}

np_parity_load_package()

options(np.messages = FALSE)
options(npRmpi.autodispatch = FALSE)
nslaves <- as.integer(Sys.getenv("NP_PLOT_PARITY_NSLAVES", "1"))
npRmpi.init(nslaves = nslaves, autodispatch = FALSE)
options(npRmpi.autodispatch = FALSE)
on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)

set.seed(20260305)
n <- 280L
x <- runif(n)
tx <- data.frame(x = x)

tol <- as.numeric(Sys.getenv('NP_PLOT_PARITY_TOL', '1e-6'))
out.dir <- Sys.getenv('NP_PLOT_PARITY_OUTDIR', '')
if (!nzchar(out.dir)) out.dir <- file.path('/tmp', sprintf('plot_helper_parity_udens_udist_serial_%s', format(Sys.time(), '%Y%m%d_%H%M%S')))
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

results <- list()
add_result <- function(row) results[[length(results) + 1L]] <<- row

run_udens <- function(case, bw.args) {
  out <- tryCatch({
    bw <- do.call(npudensbw, c(list(dat = tx), bw.args))
    fit <- npudens(tdat = tx, bws = bw)
    pdata <- suppressWarnings(plot(
      fit,
      plot.behavior = 'data',
      neval = 80,
      plot.errors.method = 'bootstrap',
      plot.errors.boot.method = 'inid',
      plot.errors.boot.num = 49,
      plot.errors.type = 'pointwise'
    ))

    eval.df <- as.data.frame(pdata$d1$eval)
    pred <- as.numeric(predict(fit, newdata = eval.df))
    plotv <- as.numeric(pdata$d1$dens)
    d <- max(abs(plotv - pred))
    list(ok = TRUE, d = d, err = "")
  }, error = function(e) {
    list(ok = FALSE, d = NA_real_, err = conditionMessage(e))
  })

  add_result(data.frame(
    family = 'npudens',
    case = case,
    max_abs_plot_minus_predict = out$d,
    pass = isTRUE(out$ok) && is.finite(out$d) && out$d <= tol,
    error = out$err,
    stringsAsFactors = FALSE
  ))

  if (isTRUE(out$ok)) {
    cat(sprintf('CASE %-34s d(plot,predict)=%.3e\n', paste('npudens', case), out$d))
  } else {
    cat(sprintf('CASE %-34s ERROR=%s\n', paste('npudens', case), out$err))
  }
}

run_udist <- function(case, bw.args) {
  out <- tryCatch({
    bw <- do.call(npudistbw, c(list(dat = tx), bw.args))
    fit <- npudist(tdat = tx, bws = bw)
    pdata <- suppressWarnings(plot(
      fit,
      plot.behavior = 'data',
      neval = 80,
      plot.errors.method = 'bootstrap',
      plot.errors.boot.method = 'inid',
      plot.errors.boot.num = 49,
      plot.errors.type = 'pointwise'
    ))

    eval.df <- as.data.frame(pdata$d1$eval)
    pred <- as.numeric(predict(fit, edat = eval.df))
    plotv <- as.numeric(pdata$d1$dist)
    d <- max(abs(plotv - pred))
    list(ok = TRUE, d = d, err = "")
  }, error = function(e) {
    list(ok = FALSE, d = NA_real_, err = conditionMessage(e))
  })

  add_result(data.frame(
    family = 'npudist',
    case = case,
    max_abs_plot_minus_predict = out$d,
    pass = isTRUE(out$ok) && is.finite(out$d) && out$d <= tol,
    error = out$err,
    stringsAsFactors = FALSE
  ))

  if (isTRUE(out$ok)) {
    cat(sprintf('CASE %-34s d(plot,predict)=%.3e\n', paste('npudist', case), out$d))
  } else {
    cat(sprintf('CASE %-34s ERROR=%s\n', paste('npudist', case), out$err))
  }
}

# npudens cases: fixed+bounds and adaptive_nn without fixed bounds (valid combinations).
run_udens('fixed-epan4-bounded', list(
  bws = 0.15,
  bandwidth.compute = FALSE,
  bwtype = 'fixed',
  ckertype = 'epanechnikov',
  ckerorder = 4,
  ckerbound = 'fixed',
  ckerlb = 0,
  ckerub = 1
))

run_udens('fixed-gaussian', list(
  bws = 0.16,
  bandwidth.compute = FALSE,
  bwtype = 'fixed',
  ckertype = 'gaussian',
  ckerorder = 2
))

run_udens('adaptive-nn-gaussian', list(
  nmulti = 1,
  bandwidth.compute = TRUE,
  bwtype = 'adaptive_nn',
  ckertype = 'gaussian',
  ckerorder = 2
))

# npudist cases
run_udist('fixed-epan4-bounded', list(
  bws = 0.15,
  bandwidth.compute = FALSE,
  bwtype = 'fixed',
  ckertype = 'epanechnikov',
  ckerorder = 4,
  ckerbound = 'fixed',
  ckerlb = 0,
  ckerub = 1
))

run_udist('fixed-gaussian', list(
  bws = 0.16,
  bandwidth.compute = FALSE,
  bwtype = 'fixed',
  ckertype = 'gaussian',
  ckerorder = 2
))

run_udist('adaptive-nn-gaussian', list(
  nmulti = 1,
  bandwidth.compute = TRUE,
  bwtype = 'adaptive_nn',
  ckertype = 'gaussian',
  ckerorder = 2
))

res <- do.call(rbind, results)
out.csv <- file.path(out.dir, 'udens_udist_plot_helper_parity_summary.csv')
write.csv(res, out.csv, row.names = FALSE)
cat(sprintf('SUMMARY_CSV=%s\n', out.csv))
if (any(!res$pass)) {
  cat('PARITY_STATUS=FAIL\n')
  q(status = 1)
}
cat('PARITY_STATUS=PASS\n')
