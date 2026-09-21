test_that("bootstrap labels distinguish an NA category from missing observations", {
  f <- factor(c("a",NA,"b"),exclude=NULL)
  dat <- untangle(data.frame(f=f))
  draws <- matrix(seq_len(27),9,3)
  out <- .np_plot_boot_factor_boxplots(draws,dat,1L,9L,f)
  expect_identical(out$names,as.character(f))
  for(j in 1:3) expect_equal(out$stats[,j],boxplot.stats(draws[,j])$stats)
  order <- c(3L,1L,2L)
  perm <- .np_plot_boot_factor_boxplots(draws[,order],dat,1L,9L,f[order])
  expect_equal(perm$stats,out$stats[,order])
  missing <- f; is.na(missing)[2L] <- TRUE
  expect_error(.np_plot_boot_factor_boxplots(draws,dat,1L,9L,missing),"do not match")
  expect_error(.np_plot_boot_factor_boxplots(draws,dat,1L,9L,c("a","bad","b")),"do not match")
  expect_error(.np_plot_boot_factor_boxplots(draws,dat,1L,9L,f[1:2]),"do not match")
})

test_that("NA-labelled factor bootstrap plots match a relabelled oracle", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  d <- data.frame(x=factor(rep(c("a",NA,"b"),8),exclude=NULL),
                  y=sin(seq_len(24)))
  a <- npreg(y~x,data=d,bws=.2)
  levels(d$x)[is.na(levels(d$x))] <- "missing-label"
  z <- npreg(y~x,data=d,bws=.2)
  set.seed(815)
  p <- plot(a,output="data",errors="bootstrap",bootstrap="inid",B=9)
  seed <- .Random.seed
  set.seed(815)
  q <- plot(z,output="data",errors="bootstrap",bootstrap="inid",B=9)
  expect_identical(.Random.seed,seed)
  expect_equal(fitted(p[[1L]]),fitted(q[[1L]]),tolerance=0)
  expect_equal(p[[1L]]$merr,q[[1L]]$merr,tolerance=0)
})
