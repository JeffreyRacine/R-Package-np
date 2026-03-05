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
nslaves <- as.integer(Sys.getenv('NP_PLOT_PARITY_NSLAVES', '1'))
npRmpi.init(nslaves = nslaves, autodispatch = FALSE)
options(npRmpi.autodispatch = FALSE)
on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)

set.seed(20260305)
n <- 320L
x <- sort(runif(n))
y <- sin(2 * pi * x) + rnorm(n, sd = 0.08)
tx <- data.frame(x = x)

tol <- as.numeric(Sys.getenv('NP_PLOT_PARITY_TOL', '1e-6'))
out.dir <- Sys.getenv('NP_PLOT_PARITY_OUTDIR', '')
if (!nzchar(out.dir)) out.dir <- file.path('/tmp', sprintf('plot_helper_parity_npreg_mpi_%s', format(Sys.time(), '%Y%m%d_%H%M%S')))
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

mk_case <- function(name, regtype, ckertype, basis = NA_character_, degree = NA_integer_, bernstein = NA) {
  data.frame(
    name = name,
    regtype = regtype,
    ckertype = ckertype,
    basis = basis,
    degree = degree,
    bernstein = bernstein,
    stringsAsFactors = FALSE
  )
}

cases <- do.call(rbind, list(
  mk_case('lc-gaussian', 'lc', 'gaussian'),
  mk_case('lc-epanechnikov', 'lc', 'epanechnikov'),
  mk_case('ll-gaussian', 'll', 'gaussian'),
  mk_case('ll-epanechnikov', 'll', 'epanechnikov'),
  mk_case('lp-glp-raw-gaussian', 'lp', 'gaussian', basis = 'glp', degree = 3L, bernstein = FALSE),
  mk_case('lp-glp-bern-gaussian', 'lp', 'gaussian', basis = 'glp', degree = 3L, bernstein = TRUE),
  mk_case('lp-add-raw-gaussian', 'lp', 'gaussian', basis = 'additive', degree = 3L, bernstein = FALSE),
  mk_case('lp-add-bern-gaussian', 'lp', 'gaussian', basis = 'additive', degree = 3L, bernstein = TRUE),
  mk_case('lp-tensor-raw-epan', 'lp', 'epanechnikov', basis = 'tensor', degree = 3L, bernstein = FALSE),
  mk_case('lp-tensor-bern-epan', 'lp', 'epanechnikov', basis = 'tensor', degree = 3L, bernstein = TRUE)
))

results <- vector('list', nrow(cases))
for (i in seq_len(nrow(cases))) {
  cs <- cases[i, , drop = FALSE]

  bw.args <- list(
    xdat = tx,
    ydat = y,
    bws = 0.18,
    regtype = cs$regtype,
    ckertype = cs$ckertype,
    bandwidth.compute = FALSE
  )
  if (identical(cs$regtype, 'lp')) {
    bw.args$basis <- cs$basis
    bw.args$degree <- as.integer(cs$degree)
    bw.args$bernstein.basis <- isTRUE(cs$bernstein)
  }

  bw <- do.call(npregbw, bw.args)
  fit <- npreg(txdat = tx, tydat = y, bws = bw, gradients = FALSE, warn.glp.gradient = FALSE)

  pdata <- suppressWarnings(plot(
    fit,
    neval = 80,
    perspective = FALSE,
    plot.behavior = 'data',
    plot.errors.method = 'bootstrap',
    plot.errors.boot.method = 'inid',
    plot.errors.boot.num = 49,
    plot.errors.type = 'pointwise'
  ))

  eval.df <- as.data.frame(pdata$r1$eval)
  names(eval.df) <- names(tx)

  pred <- as.numeric(predict(fit, newdata = eval.df))
  hat <- as.numeric(npreghat(bws = bw, txdat = tx, y = y, exdat = eval.df, s = 0L, output = 'apply'))
  plot.mean <- as.numeric(pdata$r1$mean)

  d.plot.pred <- max(abs(plot.mean - pred))
  d.hat.pred <- max(abs(hat - pred))
  d.plot.hat <- max(abs(plot.mean - hat))

  results[[i]] <- data.frame(
    case = cs$name,
    regtype = cs$regtype,
    ckertype = cs$ckertype,
    basis = if (is.na(cs$basis)) '' else cs$basis,
    degree = if (is.na(cs$degree)) '' else as.character(cs$degree),
    bernstein = if (is.na(cs$bernstein)) '' else as.character(cs$bernstein),
    max_abs_plot_minus_predict = d.plot.pred,
    max_abs_hat_minus_predict = d.hat.pred,
    max_abs_plot_minus_hat = d.plot.hat,
    pass = (d.plot.pred <= tol) && (d.hat.pred <= tol) && (d.plot.hat <= tol),
    stringsAsFactors = FALSE
  )

  cat(sprintf('CASE %-28s d(plot,pred)=%.3e d(hat,pred)=%.3e d(plot,hat)=%.3e\n',
              cs$name, d.plot.pred, d.hat.pred, d.plot.hat))
}

res <- do.call(rbind, results)
out.csv <- file.path(out.dir, 'npreg_plot_helper_parity_summary.csv')
write.csv(res, out.csv, row.names = FALSE)

cat(sprintf('SUMMARY_CSV=%s\n', out.csv))
if (any(!res$pass)) {
  cat('PARITY_STATUS=FAIL\n')
  q(status = 1)
}
cat('PARITY_STATUS=PASS\n')
