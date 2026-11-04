test_that("D8 fixed pilot moments include the complete corrected polynomial", {
  x <- data.frame(x = c(-1.2, -.3, .1, .8, 1.4))
  g <- data.frame(x = c(-.9, .2, 1.1))
  state <- list(x=x, n=5L, q=1L, r=0L, s=0L, h=.6)
  actual <- .np_udens_d8_pilot_moments(state, g)
  L <- function(y, a) {
    u <- (a-y)/.6
    dnorm(u)/.6 * (35/16-35*u^2/16+7*u^4/16-u^6/48)
  }
  pilot <- function(y) rowMeans(vapply(x$x, function(mu) dnorm(y, mu, .6), numeric(length(y))))
  mu <- vapply(g$x, function(a) integrate(function(y) L(y,a)*pilot(y), -12, 12, rel.tol=1e-9)$value, 0.0)
  second <- vapply(g$x, function(a) integrate(function(y) L(y,a)^2*pilot(y), -12, 12, rel.tol=1e-9)$value, 0.0)
  expect_equal(actual$mean, mu, tolerance=1e-9)
  expect_equal(actual$se^2, (second-mu^2)/5, tolerance=1e-9)
  C <- vapply(seq_len(nrow(g)), function(j)
    .np_udens_d8_pilot_column(state,g,actual$mean,j), numeric(nrow(g)))
  expect_equal(C, t(C), tolerance=1e-12)
  expect_equal(diag(C), actual$se^2, tolerance=1e-12)
  expect_gte(min(eigen(C,symmetric=TRUE,only.values=TRUE)$values), -1e-12)
})

test_that("D8 factor storage is exact, bounded and invocation owned", {
  store <- .np_udens_d8_store(4L)
  on.exit(.np_udens_d8_store_close(store), add=TRUE)
  .np_udens_d8_store_append(store, 1:4)
  .np_udens_d8_store_append(store, 5:8)
  expect_identical(.np_udens_d8_store_panel(store,1L,2L), matrix(as.double(1:8),4))
  expect_identical(.np_udens_d8_store_rows(store,2L,3L), matrix(c(2,3,6,7),2))
  expect_error(.np_udens_d8_store_panel(store,1L,3L), "incomplete")
  path <- store$path
  .np_udens_d8_store_close(store)
  expect_false(file.exists(path))
  expect_silent(.np_udens_d8_store_close(store))
  state <- list(n=100L,q=2L,max.bytes=256*1024^2)
  b <- .np_udens_d8_pilot_budget(state,2500L,9999L,3L,32L,32L,256L)
  expect_lt(b$bytes,state$max.bytes)
  expect_identical(b$maximum.factor.disk.bytes,8*2500^2)
  state$max.bytes <- 1024
  expect_error(.np_udens_d8_pilot_budget(state,2500L,9999L,3L,32L,32L,256L), "byte budget")
})
