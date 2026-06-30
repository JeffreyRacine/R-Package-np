expected_singleindex_plot <- function(bw, tx, y, neval = 50L, gradients = FALSE) {
  eval_grid <- getFromNamespace(".np_plot_singleindex_eval_grid", "np")
  local_eval <- getFromNamespace(".np_plot_singleindex_local_eval", "np")

  info <- eval_grid(bws = bw, xdat = tx, neval = neval, where = "test")
  out <- local_eval(
    bws = bw,
    idx.train = info$idx.train,
    idx.eval = info$idx.eval,
    ydat = y,
    gradients = gradients
  )
  out$trainiseval <- info$trainiseval
  out
}

with_plot_device <- function(expr) {
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(force(expr))
}

run_singleindex_bootstrap_plot <- function(bw,
                                           xdat,
                                           ydat,
                                           boot_method,
                                           gradients = FALSE,
                                           output = "data",
                                           boot_num = 9L,
                                           neval = 50L) {
  suppressWarnings(plot(
    bw,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    output = behavior,
    perspective = FALSE,
    gradients = gradients,
    errors = "bootstrap",
    bootstrap = boot_method,
    B = boot_num
  ))[[1]]
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
    h <- if (identical(bt, "fixed")) 0.25 else 5L
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
      fit.mean.aligned <- fit.mean$mean
      fit.grad.mean.aligned <- fit.grad$mean
      fit.grad.aligned <- fit.grad$grad
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
    h <- if (identical(bt, "fixed")) 0.25 else 5L
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
      neval <- 17L
      expected.mean <- expected_singleindex_plot(
        bw = bw,
        tx = tx,
        y = y,
        neval = neval,
        gradients = FALSE
      )

      out.data.mean <- suppressWarnings(plot(
        bw,
        xdat = tx,
        ydat = y,
        neval = neval,
        output = "data",
        perspective = FALSE,
        gradients = FALSE
      ))[[1]]
      out.plot.mean <- with_plot_device(plot(
        bw,
        xdat = tx,
        ydat = y,
        neval = neval,
        output = "plot-data",
        perspective = FALSE,
        gradients = FALSE
      ))[[1]]

      expect_equal(
        out.data.mean$index,
        expected.mean$index,
        tolerance = 1e-12,
        info = paste(bt, cfg$label, "index")
      )
      expect_equal(
        out.data.mean$mean,
        expected.mean$mean,
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "mean")
      )
      expect_equal(
        out.data.mean$mean,
        out.plot.mean$mean,
        tolerance = 1e-12,
        info = paste(bt, cfg$label, "mean plot-data")
      )
      expect_false(out.data.mean$trainiseval)
      expect_equal(length(out.data.mean$mean), neval)

      expected.grad <- expected_singleindex_plot(
        bw = bw,
        tx = tx,
        y = y,
        neval = neval,
        gradients = TRUE
      )

      out.data.grad <- suppressWarnings(plot(
        bw,
        xdat = tx,
        ydat = y,
        neval = neval,
        output = "data",
        perspective = FALSE,
        gradients = TRUE
      ))[[1]]
      out.plot.grad <- with_plot_device(plot(
        bw,
        xdat = tx,
        ydat = y,
        neval = neval,
        output = "plot-data",
        perspective = FALSE,
        gradients = TRUE
      ))[[1]]

      expect_equal(
        out.data.grad$index,
        expected.grad$index,
        tolerance = 1e-12,
        info = paste(bt, cfg$label, "grad index")
      )
      expect_equal(
        out.data.grad$mean,
        expected.grad$mean,
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "grad mean")
      )
      expect_equal(
        out.data.grad$grad,
        expected.grad$grad,
        tolerance = 1e-10,
        info = paste(bt, cfg$label, "grad")
      )
      expect_equal(
        out.data.grad$grad,
        out.plot.grad$grad,
        tolerance = 1e-12,
        info = paste(bt, cfg$label, "grad plot-data")
      )
      expect_false(out.data.grad$trainiseval)
      expect_equal(nrow(out.data.grad$grad), neval)
    }
  }
})

test_that("sibandwidth bootstrap helper routes work across bwtype, regtype, and method", {
  skip_if_not_installed("np")

  set.seed(20260307)
  n <- 35
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL, label = "lc"),
    list(regtype = "ll", basis = NULL, degree = NULL, label = "ll"),
    list(regtype = "lp", basis = "tensor", degree = 2L, label = "lp")
  )

  for (bt in c("fixed", "generalized_nn", "adaptive_nn")) {
    h <- if (identical(bt, "fixed")) 0.25 else 5L
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
      neval <- 19L
      for (boot_method in c("wild", "inid", "fixed", "geom")) {
        out.data <- run_singleindex_bootstrap_plot(
          bw = bw,
          xdat = tx,
          ydat = y,
          boot_method = boot_method,
          gradients = FALSE,
          output = "data",
          neval = neval
        )
        out.plot <- with_plot_device(run_singleindex_bootstrap_plot(
          bw = bw,
          xdat = tx,
          ydat = y,
          boot_method = boot_method,
          gradients = FALSE,
          output = "plot-data",
          neval = neval
        ))

        expect_equal(
          length(out.data$mean),
          neval,
          info = paste(bt, cfg$label, boot_method, "mean length")
        )
        expect_equal(
          dim(out.data$merr),
          c(neval, 2L),
          info = paste(bt, cfg$label, boot_method, "merr shape")
        )
        expect_true(
          all(is.finite(out.data$mean)),
          info = paste(bt, cfg$label, boot_method, "mean finite")
        )
        expect_true(
          all(is.finite(out.data$merr)),
          info = paste(bt, cfg$label, boot_method, "merr finite")
        )
        expect_equal(
          out.data$mean,
          out.plot$mean,
          tolerance = 1e-12,
          info = paste(bt, cfg$label, boot_method, "plot-data center")
        )
      }
    }
  }
})

test_that("sibandwidth nonfixed bootstrap helpers fail fast for gradients", {
  skip_if_not_installed("np")

  set.seed(20260307)
  n <- 30
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)

  for (bt in c("generalized_nn", "adaptive_nn")) {
    h <- 5L
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, h),
      bandwidth.compute = FALSE,
      regtype = "ll",
      bwtype = bt
    )

    for (boot_method in c("inid", "fixed", "geom")) {
      expect_error(
        run_singleindex_bootstrap_plot(
          bw = bw,
          xdat = tx,
          ydat = y,
          boot_method = boot_method,
          gradients = TRUE
        ),
        "requires helper mode with gradients=FALSE",
        info = paste(bt, boot_method, "gradient bootstrap guard")
      )
    }
  }
})
