test_that("density-family ordered mass survives operator conversion", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  x <- data.frame(o = ordered(rep(0:2, 4), levels = 0:2))
  e <- x[c(1, 2), , drop = FALSE]
  y <- seq_len(nrow(x))^2
  delta <- outer(as.numeric(e$o), as.numeric(x$o), "-")
  mass <- .3^abs(delta) * .7 / 1.3 / nrow(x)
  cdf <- ifelse(delta < 0, .3^abs(delta) / 1.3,
                1 - .3^(abs(delta) + 1) / 1.3) / nrow(x)
  for (target in c("density", "distribution")) {
    bw <- if (target == "density")
      npudensbw(dat = x, bws = .3, bandwidth.compute = FALSE) else
      npudistbw(dat = x, bws = .3, bandwidth.compute = FALSE)
    hat <- if (target == "density") npudenshat else npudisthat
    expected <- if (target == "density") mass else cdf
    expect_equal(unname(as.matrix(hat(bw, x, e)))[, ], expected,
                 ignore_attr = TRUE, tolerance = 1e-12)
    expect_equal(hat(bw, x, e, y = y, output = "apply"),
                 drop(expected %*% y), tolerance = 1e-12)
    expect_identical(kbandwidth(bw)$okertype, "nliracine")
  }
  raw <- kbandwidth.numeric(.3, xdati = untangle(x), xnames = "o")
  expect_identical(.np_kbandwidth_okertype(raw), "liracine")
  # Copula marginal adapters deliberately receive plain density-role lists.
  expect_identical(.np_make_kbandwidth_unconditional(
    list(bw = .3, type = "fixed", ckertype = "gaussian", ckerorder = 2,
         ckerbound = "none", ukertype = "aitchisonaitken", okertype = "liracine"),
    x)$okertype, "nliracine")
})

test_that("density profile bootstrap uses the public categorical kernel mass", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add=TRUE)
  pkg <- getNamespaceName(environment(npudens))
  internal <- function(name) getFromNamespace(name,pkg)
  withr::local_options(np.messages=FALSE, np.categorical.compress=TRUE)
  withr::local_preserve_seed()
  set.seed(72019)
  x <- data.frame(u=factor(rep(1:3,20)),o=ordered(rep(0:3,15)))
  e <- x[c(2,5,8,14),,drop=FALSE]
  n <- nrow(x)
  counts <- rmultinom(5L,n,rep(1/n,n))
  for(uk in c("liracine","aitchisonaitken"))
    for(ok in c("liracine","wangvanryzin","racineliyan")) {
      bw <- npudensbw(dat=x,bws=c(.3,.2),bandwidth.compute=FALSE,
                      ukertype=uk,okertype=ok)
      same <- outer(as.integer(e$u),as.integer(x$u),"==")
      wu <- if(uk=="liracine") ifelse(same,1,.3)/1.6 else ifelse(same,.7,.15)
      d <- abs(outer(as.integer(e$o),as.integer(x$o),"-"))
      wo <- switch(ok,liracine=.2^d*.8/1.2,
        wangvanryzin=ifelse(d==0,.8,.4*.2^d),
        racineliyan=t(t(.2^d)/vapply(as.integer(x$o),
          function(a)sum(.2^abs(a-1:4)),numeric(1))))
      W <- wu*wo
      codes <- internal(".np_cat_profile_code_matrix")
      expect_equal(internal(".np_density_cat_profile_kernel_matrix")(
        codes(e),codes(x),x,bw),W,tolerance=1e-13)
      expect_equal(fitted(npudens(bw,tdat=x,edat=e)),rowMeans(W),tolerance=1e-12)
      got <- internal(".np_inid_boot_from_ksum_unconditional")(x,e,bw,
        B=5L,operator="normal",counts=counts)
      expect_equal(got$t,t(W%*%counts/n),tolerance=1e-12)
      expect_equal(got$t0,rowMeans(W),tolerance=1e-12)
      # Disable compression only as a control, not as the production repair.
      reference <- withr::with_options(list(np.categorical.compress=FALSE),
        internal(".np_inid_boot_from_ksum_unconditional")(x,e,bw,
          B=5L,operator="normal",counts=counts))
      expect_equal(got,reference,tolerance=1e-12)
    }
})
