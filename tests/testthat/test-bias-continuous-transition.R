test_that("continuous transition uses a geometry-relative nested projection", {
  x <- matrix(seq(-1,1,length.out=48),ncol=1L,dimnames=list(NULL,"x"))
  Q <- np:::.np_bias_adaptation_basis(x,x,5L)
  geometry <- np:::.np_bias_adaptation_geometry(1L,1+1e-6,1,.4,.35,5L)
  state <- list(x=x,Q=Q,geometry=geometry)
  at <- matrix(.2,1,1,dimnames=list(NULL,"x"))
  z <- np:::.np_bias_adaptation_pilot(state,at,geometry$b1)
  # Independent sum of weighted least-squares projections; the old absolute
  # n*lambda penalty fails this identity even though its output is finite.
  levels <- sort(unique(geometry$lambda)); d <- 1/(1+levels)
  mixture <- d-c(d[-1L],0); expected <- numeric(nrow(x))
  v <- np:::.np_bias_adaptation_basis(x,at,5L)
  for(k in seq_along(levels)) {
    ii <- which(geometry$lambda<=levels[k]); X <- Q[,ii,drop=FALSE]
    A <- solve(crossprod(X,z$g*X),t(X)*rep(z$g,each=length(ii)))
    expected <- expected+mixture[k]*drop(v[,ii,drop=FALSE]%*%A)
  }
  expect_equal(drop(z$H),expected,tolerance=1e-9)
  expect_equal(unname(drop(z$H%*%cbind(1,x))),c(1,.2),tolerance=1e-11)
  local <- state;local$geometry$lambda[] <- 0
  original <- np:::.np_bias_adaptation_pilot(local,at,geometry$b1)
  expect_lt(max(abs(z$H-original$H)),1e-4)
})


test_that("continuous pilot preparation preserves named kernel and derivative rows", {
  x <- cbind(x1=sin(seq_len(72)),x2=cos(seq_len(72)*sqrt(2)))
  rownames(x) <- paste0("row",seq_len(nrow(x)))
  Q <- np:::.np_bias_adaptation_basis(x,x,c(2L,2L))
  geometry <- np:::.np_bias_adaptation_geometry(c(1L,1L),c(1.4,.6),
    c(1,1),c(.4,.4),c(.35,.35),c(2L,2L))
  state <- list(x=x,Q=Q,geometry=geometry)
  at <- x[17,,drop=FALSE]
  prepared <- state
  prepared$geometry$transition <- np:::.np_bias_transition_geometry(Q,geometry$lambda,x)
  for(h in list(geometry$b1,c(Inf,geometry$b1[2]))) {
    expected <- drop(np:::.np_bias_adaptation_kernel(x,at,h))
    actual <- np:::.np_bias_transition_weights(prepared,at,h,prepared$geometry$transition)
    expect_identical(actual,expected)
  }
  original <- np:::.np_bias_adaptation_pilot(state,at,geometry$b1,
    derivatives=1:2,orders=c(1L,1L))
  reused <- np:::.np_bias_adaptation_pilot(prepared,at,geometry$b1,
    derivatives=1:2,orders=c(1L,1L))
  expect_identical(reused,original)
})
