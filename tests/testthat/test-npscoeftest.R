test_that("smooth-coefficient pair panels omit self-pairs and preserve signed scores", {
  ns <- asNamespace("npRmpi")
  expect_identical(getExportedValue("npRmpi", "npscoeftest"), ns$npscoeftest)
  for (generic in c("print","summary"))
    expect_identical(getS3method(generic, "scoeftest", envir=ns),
                     get(paste0(generic, ".scoeftest"), envir=ns))
  w <- data.frame(x=c(.1,.8,.5,.4,.2,.9),z=c(.5,.7,.3,.1,.4,.2))
  h <- c(.4,.9)
  e <- cbind(c(-1,1,0,0,0,0), c(2,.8,1,-1,3,-.6))
  k <- exp(-as.matrix(dist(sweep(as.matrix(w),2,h,"/")))^2/2)
  diag(k) <- 0
  q <- colSums(e*(k %*% e))
  expect_lt(q[1],0)
  for (student in c(FALSE,TRUE)) {
    s <- list(w=w,h=h,student=student,width=4L)
    bank <- ns$.sct_score_bank(e,s$width,student)
    # The final panel is partial; the diagonal square occurs only once.
    got <- ns$.sct_panel_score(s,bank,1L)+ns$.sct_panel_score(s,bank,2L)
    expected <- if(student) c(q,colSums(e^2*((k^2) %*% e^2))) else q
    expect_equal(got,expected,tolerance=1e-12)
  }
})
