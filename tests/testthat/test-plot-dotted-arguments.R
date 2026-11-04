test_that("retired plot spellings are rejected before their values are forced", {
  old_names <- c('data.overlay','data.rug','factor.boxplot','boxplot.outliers',
               'coef.index','common.scale','proper.method','proper.control',
               'boot.control','grid.control','render.control')
  current <- chartr('.', '_', old_names)
  retired <- old_names
  for(i in seq_along(current)) {
    bad <- setNames(list(quote(stop('must not be evaluated'))),retired[i])
    expect_error(.np_plot_validate_public_dots(bad),paste0('retired plot argument: ', retired[i], '; use ', current[i]),fixed=TRUE)
    both <- c(bad,setNames(list(NULL),current[i]))
    expect_error(.np_plot_validate_public_dots(both),retired[i],fixed=TRUE)
  }
  expect_error(plot.npcopula(NULL,boot.control=stop('forced')), 'boot_control',fixed=TRUE)
  expect_error(np.pairs(y.vars='x',y.dat=data.frame(x=1:3),y_vars=stop('forced')),
               'y.vars',fixed=TRUE)
  expect_error(np_render_control(bar_num=NULL),'unused argument',fixed=TRUE)
  expect_error(np.pairs.plot(pair_list=NULL),'unused argument',fixed=TRUE)
})

test_that("ordinary plot typo suggestions use snake_case canonical names", {
  for (name in c('data_ovvverlay', 'data_ruug', 'boot_controll')) {
    expected <- c('data_ovvverlay'='data_overlay', 'data_ruug'='data_rug',
                  'boot_controll'='boot_control')[[name]]
    bad <- setNames(list(quote(stop('must not be evaluated'))), name)
    expect_error(.np_plot_validate_public_dots(bad),
                 paste0('did you mean ', expected, '?'), fixed=TRUE)
  }
  expect_error(plot(structure(list(), class='npregression'),
                    data_ovvverlay=stop('must not be evaluated')),
               'did you mean data_overlay?', fixed=TRUE)
})

test_that("control constructors retain their contract while plot inputs use snake_case", {
  ctrl<-np_render_control(bar.num=3L)
  expect_identical(ctrl,structure(list(style='band',bar='|',bar_num=3L),class='np_render_control'))
  expect_identical(do.call(np_render_control,list(bar.num=3L)),ctrl)
  expect_identical(np_render_control('band','|',3L),ctrl)
  expect_error(np_render_control(bar.num=1L,bar.num=2L),'matched by multiple')
  expect_error(np_render_control(bar.nu=2.5),'bar.num',fixed=TRUE)
  opts<-list(errors='bootstrap',data_overlay=FALSE,data_rug=TRUE,
    factor_boxplot=TRUE,boxplot_outliers=FALSE,common_scale=FALSE,coef_index=2L,
    proper_method='isotonic',proper_control=list(mode='slice'),
    boot_control=np_boot_control(blocklen=3L),grid_control=np_grid_control(xtrim=c(.1,.9)),
    render_control=ctrl)
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
  expect_identical(.np_plot_normalize_public_dots(list(proper_control=NULL)),list(proper.control=NULL))
  expect_error(.np_plot_normalize_public_dots(list(data_overlay=NA)), 'data_overlay',fixed=TRUE)
  expect_error(.np_plot_normalize_public_dots(list(grid_control=list())), 'np_grid_control',fixed=TRUE)
})
