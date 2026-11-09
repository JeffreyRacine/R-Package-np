test_that("bounded normal public bands are annotated as asymptotic SE intervals", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  ns <- asNamespace("npRmpi")
  original <- get(".np_plot_variability_descriptor", ns)
  seen <- list()
  capture <- function(...) { z <- original(...); seen[[length(seen)+1L]] <<- z; z }
  utils::assignInNamespace(".np_plot_variability_descriptor", capture, ns="npRmpi")
  on.exit(utils::assignInNamespace(".np_plot_variability_descriptor", original, ns="npRmpi"), add=TRUE)
  path <- tempfile(fileext=".pdf"); grDevices::pdf(path)
  on.exit({grDevices::dev.off();unlink(path)}, add=TRUE)
  x <- qbeta((seq_len(40L)-.5)/40L, 2, 3)
  for (kernel in c("beta", "gaussian")) {
    bw <- npudensbw(dat=data.frame(x=x),bws=.12,bandwidth.compute=FALSE,
      ckertype=kernel,ckerbound="fixed",ckerlb=0,ckerub=1)
    seen <- list()
    z <- plot(bw,xdat=data.frame(x=x),errors="asymptotic",center="bias-corrected",
      band="pmzsd",neval=3L,output="plot-data")
    expect_s3_class(z,"np_plot_data")
    expect_gt(length(seen),0L)
    labels <- unlist(lapply(seen,`[[`,"candidates"))
    expect_true(all(grepl("asymptotic.*SE",labels)))
    expect_false(any(grepl("guide|residual envelope",labels)))
  }
})
