test_that("uocquantile basic functionality works", {
  x_num <- c(1, 2, 3, 4, 5)
  expect_equal(unname(uocquantile(x_num, 0.5)), 3)
  
  x_fact <- factor(c("A", "A", "B", "C"))
  expect_equal(as.character(uocquantile(x_fact, 0.5)), "A")
  
  x_ord <- ordered(c("Low", "Medium", "High"), levels=c("Low", "Medium", "High"))
  expect_equal(as.character(uocquantile(x_ord, 0.5)), "Medium")
})

test_that("nptgauss basic functionality works", {
  # nptgauss is called for its side effect (setting C variables) 
  # or returns NULL/constants.
  expect_silent(nptgauss(3.0))
})

test_that("plot basic functionality works", {
  data("faithful")
  bw <- npudensbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # Use pdf(NULL) to avoid opening a window
  pdf(NULL)
  on.exit(dev.off())
  
  expect_silent(plot(bw))
})

test_that("se and gradients methods work", {
  data("cps71")
  cps71_sub <- cps71[1:50, ]
  bw <- npregbw(logwage~age, data=cps71_sub, bws=1.0, bandwidth.compute=FALSE)
  model <- npreg(bws=bw, gradients=TRUE)
  
  expect_type(se(model), "double")
  expect_length(se(model), 50)
  
  expect_type(gradients(model), "double")
  expect_length(gradients(model), 50)
})
