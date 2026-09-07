test_that("smooth-coefficient demo rows do not claim unsupported gradients", {
  env <- new.env(parent = baseenv())
  path <- system.file('demo_family_npscoef.R', package = 'npRmpi')
  sys.source(path, envir = env)
  for (tier in c('smoke', 'sentinel')) {
    rows <- read.csv(system.file('demo_matrices', paste0('npscoef-',tier,'.csv'),
                                package='npRmpi'), stringsAsFactors=FALSE,
                     na.strings=character(), colClasses='character')
    expect_identical(rows$gradients, rep('FALSE',nrow(rows)))
    for (i in seq_len(nrow(rows)))
      expect_true(env$npscoef_demo_validate_row(rows[i,,drop=FALSE]))
    bad <- rows[1,,drop=FALSE]; bad$gradients <- 'TRUE'
    expect_error(env$npscoef_demo_validate_row(bad),
                 'demo rows require gradients=FALSE', fixed=TRUE)
  }
})
