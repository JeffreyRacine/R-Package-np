test_that("retired plot spellings are rejected before their values are forced", {
  old_names <- c('data.overlay','data.rug','factor.boxplot','boxplot.outliers',
               'coef.index','common.scale','proper.method','proper.control',
               'boot.control','grid.control','render.control','bias.control',
               'gaussian.draws','random.seed','gradient.order','plot.vars','fit.mean.train')
  current <- old_names
  retired <- chartr('.', '_', current)
  for(i in seq_along(current)) {
    bad <- setNames(list(quote(stop('must not be evaluated'))),retired[i])
    expect_error(.np_plot_validate_public_dots(bad),paste0('retired plot argument: ', retired[i], '; use ', current[i]),fixed=TRUE)
    both <- c(bad,setNames(list(NULL),current[i]))
    expect_error(.np_plot_validate_public_dots(both),retired[i],fixed=TRUE)
  }
  expect_error(plot.npcopula(NULL,boot_control=stop('forced')), 'boot.control',fixed=TRUE)
  expect_error(np.pairs(y.vars='x',y.dat=data.frame(x=1:3),y_vars=stop('forced')),
               'y.vars',fixed=TRUE)
  expect_error(np_render_control(bar_num=NULL),'unused argument',fixed=TRUE)
  expect_error(np.pairs.plot(pair_list=NULL),'unused argument',fixed=TRUE)
})

test_that("ordinary plot typo suggestions use dotted canonical names", {
  for (name in c('data.ovvverlay', 'data.ruug', 'boot.controll')) {
    expected <- c('data.ovvverlay'='data.overlay', 'data.ruug'='data.rug',
                  'boot.controll'='boot.control')[[name]]
    bad <- setNames(list(quote(stop('must not be evaluated'))), name)
    expect_error(.np_plot_validate_public_dots(bad),
                 paste0('did you mean ', expected, '?'), fixed=TRUE)
  }
  expect_error(plot(structure(list(), class='npregression'),
                    data.ovvverlay=stop('must not be evaluated')),
               'did you mean data.overlay?', fixed=TRUE)
})

test_that("control constructors retain their contract while plot inputs use dotted names", {
  ctrl<-np_render_control(bar.num=3L)
  expect_identical(ctrl,structure(list(style='band',bar='|',bar_num=3L),class='np_render_control'))
  expect_identical(do.call(np_render_control,list(bar.num=3L)),ctrl)
  expect_identical(np_render_control('band','|',3L),ctrl)
  expect_error(np_render_control(bar.num=1L,bar.num=2L),'matched by multiple')
  expect_error(np_render_control(bar.nu=2.5),'bar.num',fixed=TRUE)
  opts<-list(errors='bootstrap',data.overlay=FALSE,data.rug=TRUE,
    factor.boxplot=TRUE,boxplot.outliers=FALSE,common.scale=FALSE,coef.index=2L,
    proper.method='isotonic',proper.control=list(mode='slice'),
    boot.control=np_boot_control(blocklen=3L),grid.control=np_grid_control(xtrim=c(.1,.9)),
    render.control=ctrl)
  out<-.np_plot_normalize_public_dots(opts)
  expect_false(out$plot.data.overlay);expect_true(out$plot.rug)
  expect_true(out$plot.bxp);expect_false(out$plot.bxp.out)
  expect_identical(out$coef.index,2L);expect_false(out$common.scale)
  expect_identical(out$proper.method,'isotonic')
  expect_identical(out$proper.control,list(mode='slice'))
  expect_identical(out$plot.errors.bar.num,3L)
  expect_identical(out$plot.errors.boot.blocklen,3L)
  expect_identical(out$xtrim,structure(list(probs=c(.1,.9)),class="np_plot_trim_range"))
  # The canonical dotted key already admits NULL at its existing family owner.
  expect_identical(.np_plot_normalize_public_dots(list(proper.control=NULL)),list(proper.control=NULL))
  expect_error(.np_plot_normalize_public_dots(list(data.overlay=NA)), 'data.overlay',fixed=TRUE)
  expect_error(.np_plot_normalize_public_dots(list(grid.control=list())), 'np_grid_control',fixed=TRUE)
})


test_that("retired spellings are rejected by every registered plot boundary", {
  pkg <- getNamespaceName(environment(npudens))
  classes <- c("npdensity", "npdistribution", "npregression", "condensity",
    "condistribution", "qregression", "lsqregression", "plregression",
    "smoothcoefficient", "singleindex", "conmode", "bandwidth", "dbandwidth",
    "rbandwidth", "conbandwidth", "condbandwidth", "scbandwidth", "sibandwidth",
    "plbandwidth", "npcopula", "npregiv", "npregivderiv", "np_plot_data")
  for (cl in classes) {
    method <- getS3method("plot", cl, envir = asNamespace(pkg))
    dummy <- structure(list(bws = structure(list(), class = "dbandwidth")), class = cl)
    expect_error(method(dummy, data_overlay = stop("must not be forced")),
      if (cl == "np_plot_data") "retired plot argument: data_overlay; plot-data replay accepts graphical overrides, not estimation controls" else
        "retired plot argument: data_overlay; use data.overlay", fixed = TRUE, info = cl)
  }
  expect_error(.np_plot_validate_public_dots(alist(data.rug = TRUE, data.rug = FALSE)),
    "duplicate plot argument", fixed = TRUE)
  expect_error(.np_plot_normalize_public_dots(list(data.rug = TRUE, data.rug = FALSE)),
    "duplicate plot argument", fixed = TRUE)
})


test_that("conditional-mode rugs use shared argument validation", {
  dummy <- structure(list(), class = "conmode")
  expect_error(plot(dummy, data.rug = NA),
               "data.rug must be TRUE or FALSE.", fixed = TRUE)
  expect_error(plot(dummy, data.rug = NULL),
               "data.rug must be TRUE or FALSE.", fixed = TRUE)
  expect_error(plot(dummy, data.rug = TRUE, data.rug = FALSE),
               "duplicate plot argument", fixed = TRUE)
  expect_error(plot(dummy, data.rug = ), "argument is missing")
  expect_error(plot(dummy, data_rug = stop("must not be forced")),
               "retired plot argument: data_rug; use data.rug", fixed = TRUE)
})


test_that("replay admission preserves its graphical-override contract", {
  expect_null(.np_plot_check_public_names(c("col", "col"), replay = TRUE))
  expect_null(.np_plot_check_public_names(c("persp", "", ""), replay = TRUE))
  expect_error(.np_plot_check_public_names("data_overlay", replay = TRUE),
    "plot-data replay accepts graphical overrides, not estimation controls", fixed = TRUE)
  expect_error(.np_plot_check_public_names(c("col", "col")),
    "duplicate plot argument", fixed = TRUE)
})

test_that("conditional-mode rug defaults and explicit selection reach drawing", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages = FALSE, plot.par.mfrow = FALSE)
  on.exit(options(old), add = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- data.frame(x = seq(.05, .95, length.out = 24L),
                  y = factor(rep(c("a", "b"), 12L)))
  fit <- npconmode(y ~ x, data = d, bws = c(.2, .3),
                  bandwidth.compute = FALSE, probabilities = TRUE)
  rug.calls <- 0L
  testthat::local_mocked_bindings(rug = function(x, ...) {
    rug.calls <<- rug.calls + 1L
  }, .package = "graphics")
  omitted <- plot(fit, output = "plot-data", neval = 3L)
  expect_identical(rug.calls, 0L)
  disabled <- plot(fit, output = "plot-data", neval = 3L, data.rug = FALSE)
  expect_identical(rug.calls, 0L)
  evaluated <- 0L
  enabled <- plot(fit, output = "plot-data", neval = 3L,
                  data.rug = { evaluated <- evaluated + 1L; TRUE })
  expect_gt(rug.calls, 0L)
  expect_identical(evaluated, 1L)
  expect_identical(omitted, disabled)
  expect_identical(enabled, disabled)
})
