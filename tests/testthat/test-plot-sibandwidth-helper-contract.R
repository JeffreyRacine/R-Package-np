align_singleindex_plot_mean <- function(plot_obj, fit) {
  idx <- as.vector(fit$index)
  plot_idx <- as.vector(plot_obj$index)
  ord <- match(plot_idx, idx)
  if (anyNA(ord))
    stop("failed to align single-index plot payload to fit index", call. = FALSE)
  fit$mean[ord]
}

align_singleindex_plot_grad <- function(plot_obj, fit) {
  idx <- as.vector(fit$index)
  plot_idx <- as.vector(plot_obj$index)
  ord <- match(plot_idx, idx)
  if (anyNA(ord))
    stop("failed to align single-index plot payload to fit index", call. = FALSE)
  as.matrix(fit$grad)[ord, , drop = FALSE]
}

with_plot_device <- function(expr) {
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(force(expr))
}

test_that("sibandwidth plot/helper route is free of public core re-entry", {
  skip_if_not_installed("np")

  touched <- c(
    ".np_indexhat_rbw",
    ".np_indexhat_kbw",
    ".np_indexhat_core",
    ".np_indexhat_exact",
    "npindexhat",
    ".np_plot_sibandwidth_engine",
    ".np_inid_boot_from_index",
    ".np_inid_boot_from_index_exact",
    "compute.bootstrap.errors.sibandwidth"
  )
  forbidden <- c(
    "\\bnpindex\\(",
    "\\bnpindexbw\\(",
    "\\bnpreg\\(",
    "\\bnpregbw\\(",
    "\\bnpplreg\\(",
    "\\bnpplregbw\\(",
    "\\bnpscoef\\(",
    "\\bnpscoefbw\\(",
    "\\bnpcdens\\(",
    "\\bnpcdensbw\\(",
    "\\bnpksum\\(",
    "\\bnpregbw\\(",
    "\\brbandwidth\\(",
    "\\bkbandwidth\\("
  )

  for (nm in touched) {
    fn <- getFromNamespace(nm, "np")
    txt <- paste(deparse(body(fn), width.cutoff = 500L), collapse = "\n")
    for (pat in forbidden) {
      expect_false(
        grepl(pat, txt, perl = TRUE),
        info = paste(nm, "contains forbidden public call pattern", pat)
      )
    }
  }
})

test_that("npindex helpers match fitted and gradient outputs across regtype and bwtype", {
  skip_if_not_installed("np")

  indexhat_rbw <- getFromNamespace(".np_indexhat_rbw", "np")
  regression_direct <- getFromNamespace(".np_regression_direct", "np")

  set.seed(20260307)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.06)
  tx <- data.frame(x1 = x1, x2 = x2)

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL, label = "lc"),
    list(regtype = "ll", basis = NULL, degree = NULL, label = "ll"),
    list(regtype = "lp", basis = "tensor", degree = 2L, label = "lp")
  )

  for (bt in c("fixed", "generalized_nn", "adaptive_nn")) {
    h <- if (identical(bt, "fixed")) 0.25 else 0.85
    for (cfg in cfgs) {
      bw.args <- list(
        xdat = tx,
        ydat = y,
        bws = c(1, 1, h),
        bandwidth.compute = FALSE,
        regtype = cfg$regtype,
        bwtype = bt
      )
      if (!is.null(cfg$basis)) {
        bw.args$basis <- cfg$basis
        bw.args$degree <- cfg$degree
      }

      bw <- do.call(npindexbw, bw.args)
      fit.mean <- npindex(
        bws = bw,
        txdat = tx,
        tydat = y,
        gradients = FALSE
      )
      fit.grad <- npindex(
        bws = bw,
        txdat = tx,
        tydat = y,
        gradients = TRUE
      )
      helper.ref <- list(index = as.vector(as.matrix(tx) %*% bw$beta))
      fit.mean.aligned <- align_singleindex_plot_mean(helper.ref, fit.mean)
      fit.grad.mean.aligned <- align_singleindex_plot_mean(helper.ref, fit.grad)
      fit.grad.aligned <- align_singleindex_plot_grad(helper.ref, fit.grad)
      idx.train <- data.frame(index = helper.ref$index)
      rbw <- indexhat_rbw(bws = bw, idx.train = idx.train)
      helper.grad <- regression_direct(
        bws = rbw,
        txdat = idx.train,
        tydat = y,
        exdat = idx.train,
        gradients = TRUE,
        gradient.order = 1L
      )

      hat.mean <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = tx,
        y = y,
        output = "apply",
        s = 0L
      )
      hat.mean.matrix <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = tx,
        output = "matrix",
        s = 0L
      )

      expect_equal(
        hat.mean,
        fit.mean.aligned,
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "hat mean")
      )
      expect_equal(
        as.vector(hat.mean.matrix %*% y),
        hat.mean,
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "hat mean matrix")
      )

      expect_equal(
        helper.grad$mean,
        fit.grad.mean.aligned,
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "helper gradient mean")
      )
      expect_equal(
        helper.grad$grad[, 1L],
        fit.grad.aligned[, 1L] / bw$beta[1L],
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "helper gradient apply")
      )
      expect_equal(
        fit.grad.aligned,
        helper.grad$grad[, 1L] %o% as.vector(bw$beta),
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "full gradient reconstruction")
      )
    }
  }
})

test_that("sibandwidth plot payload matches single-index fits across regtype and bwtype", {
  skip_if_not_installed("np")

  set.seed(20260307)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.06)
  tx <- data.frame(x1 = x1, x2 = x2)

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL, label = "lc"),
    list(regtype = "ll", basis = NULL, degree = NULL, label = "ll"),
    list(regtype = "lp", basis = "tensor", degree = 2L, label = "lp")
  )

  for (bt in c("fixed", "generalized_nn", "adaptive_nn")) {
    h <- if (identical(bt, "fixed")) 0.25 else 0.85
    for (cfg in cfgs) {
      bw.args <- list(
        xdat = tx,
        ydat = y,
        bws = c(1, 1, h),
        bandwidth.compute = FALSE,
        regtype = cfg$regtype,
        bwtype = bt
      )
      if (!is.null(cfg$basis)) {
        bw.args$basis <- cfg$basis
        bw.args$degree <- cfg$degree
      }

      bw <- do.call(npindexbw, bw.args)

      out.data.mean <- suppressWarnings(plot(
        bw,
        xdat = tx,
        ydat = y,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = FALSE
      ))[[1]]
      out.plot.mean <- with_plot_device(plot(
        bw,
        xdat = tx,
        ydat = y,
        plot.behavior = "plot-data",
        perspective = FALSE,
        gradients = FALSE
      ))[[1]]
      fit.mean <- npindex(
        bws = bw,
        txdat = tx,
        tydat = y,
        gradients = FALSE
      )

      expect_equal(
        out.data.mean$mean,
        align_singleindex_plot_mean(out.data.mean, fit.mean),
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "mean")
      )
      expect_equal(
        out.data.mean$mean,
        out.plot.mean$mean,
        tolerance = 1e-12,
        info = paste(bt, cfg$label, "mean plot-data")
      )

      out.data.grad <- suppressWarnings(plot(
        bw,
        xdat = tx,
        ydat = y,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = TRUE
      ))[[1]]
      out.plot.grad <- with_plot_device(plot(
        bw,
        xdat = tx,
        ydat = y,
        plot.behavior = "plot-data",
        perspective = FALSE,
        gradients = TRUE
      ))[[1]]
      fit.grad <- npindex(
        bws = bw,
        txdat = tx,
        tydat = y,
        gradients = TRUE
      )

      expect_equal(
        out.data.grad$mean,
        align_singleindex_plot_mean(out.data.grad, fit.grad),
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "grad mean")
      )
      expect_equal(
        out.data.grad$grad,
        align_singleindex_plot_grad(out.data.grad, fit.grad),
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "grad")
      )
      expect_equal(
        out.data.grad$grad,
        out.plot.grad$grad,
        tolerance = 1e-12,
        info = paste(bt, cfg$label, "grad plot-data")
      )
    }
  }
})
