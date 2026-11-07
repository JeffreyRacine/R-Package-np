test_that("bounded companions allocate residual risk over distinct targets only", {
  r <- c(-.2,.4); region <- list(lower=c(.3,.2),upper=c(.8,1.2))
  scales <- list(c(.1,.3),c(.2,.1))
  for(alpha in c(.1,.05,.01)) {
    got <- .np_udens_bounded_analytic(r,region,scales,alpha,2L)
    for(band in names(got)) {
      tail <- alpha / if(band=="pointwise")4 else 8
      z <- uniroot(function(x)pnorm(x,lower.tail=FALSE)-tail,c(0,10),tol=1e-13)$root
      expect_equal(got[[band]]$lower,pmin(r+region$lower-z*scales[[1]],r+region$lower-z*scales[[2]]),tolerance=1e-12)
      expect_equal(got[[band]]$upper,pmax(r+region$upper+z*scales[[1]],r+region$upper+z*scales[[2]]),tolerance=1e-12)
    }
    expect_true(all(got$bonferroni$lower<=got$pointwise$lower))
    expect_true(all(got$bonferroni$upper>=got$pointwise$upper))
    one <- .np_udens_bounded_analytic(r[1],lapply(region,`[`,1),list(.1),alpha,1L)
    expect_identical(one$pointwise,one$bonferroni)
  }
  expect_error(.np_udens_bounded_analytic(r,region,scales,.05,3L),'request')
  expect_error(.np_udens_bounded_analytic(r,region,list(c(-1,1)),.05,2L),'request')
})

test_that("bounded all assembly uses original units and signed offsets", {
  band <- list(lower=c(-1,2),upper=c(3,4),critical=2)
  result <- list(center=c(.5,1),family_size=2L,bands=list(pointwise=band,simultaneous=band,bonferroni=band))
  out <- list(inference=list(settings=list(),resolution=list(),center=result$center,
    pilot.residual.scale=c(1,2),empirical.residual.scale=c(2,3)))
  got <- .np_udens_bounded_public_companions(out,result,'all',2,.05)
  expect_identical(names(got$boot.all.err),c('pointwise','simultaneous','bonferroni'))
  expect_equal(got$inference$bands$pointwise[,1,1],c(-.5,1))
  expect_equal(got$boot.err[,1:2],cbind(c(1,0),c(1,1)))
  expect_identical(got$inference$settings$family.unique.targets,2L)
  sub <- .np_udens_beta_plot_subset(got,c(2L,1L,2L))
  expect_equal(sub$boot.all.err$bonferroni,got$boot.all.err$bonferroni[c(2,1,2),])
  expect_identical(sub$inference$settings$family.unique.targets,2L)
  expect_equal(sub$inference$resolution$bonferroni.signed.tail.probability,.05/8)
})
