test_that("retired plot spellings are rejected before their values are forced", {
  current <- c('data.overlay','data.rug','factor.boxplot','boxplot.outliers',
               'coef.index','common.scale','proper.method','proper.control',
               'boot.control','grid.control','render.control')
  retired <- chartr('.', '_',current)
  for(i in seq_along(current)) {
    bad <- setNames(list(quote(stop('must not be evaluated'))),retired[i])
    expect_error(.np_plot_validate_public_dots(bad),paste0('did you mean ',current[i]),fixed=TRUE)
    both <- c(bad,setNames(list(NULL),current[i]))
    expect_error(.np_plot_validate_public_dots(both),retired[i],fixed=TRUE)
  }
  expect_error(plot.npcopula(NULL,boot_control=stop('forced')), 'boot.control',fixed=TRUE)
  expect_error(np.pairs(y.vars='x',y.dat=data.frame(x=1:3),y_vars=stop('forced')),
               'y.vars',fixed=TRUE)
  expect_error(np_render_control(bar_num=NULL),'unused argument',fixed=TRUE)
  expect_error(np.pairs.plot(pair_list=NULL),'unused argument',fixed=TRUE)
})

test_that("ordinary typo suggestions use the dotted canonical names", {
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

test_that("dotted controls preserve saved schemas and normalization", {
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
  expect_identical(out$xtrim,c(.1,.9))
  # The canonical dotted key already admits NULL at its existing family owner.
  expect_identical(.np_plot_normalize_public_dots(list(proper.control=NULL)),list(proper.control=NULL))
  expect_error(.np_plot_normalize_public_dots(list(data.overlay=NA)), 'data.overlay',fixed=TRUE)
  expect_error(.np_plot_normalize_public_dots(list(grid.control=list())), 'np_grid_control',fixed=TRUE)
})
