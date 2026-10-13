test_that("copula display quotient preserves original rows and rejects conflicts", {
  grid.eval <- getFromNamespace(".npcopula_grid_eval", "np")
  grid.matrix <- getFromNamespace(".npcopula_grid_matrix", "np")
  u <- expand.grid(u1 = c(.5, .2, .5), u2 = c(.25, .5, .75))
  value <- u$u1 + 2 * u$u2
  dat <- data.frame(copula = value, u, x = u$u1 - 1, y = u$u2 + 1)
  fit <- structure(list(eval = dat, copula = value, xnames = c("x", "y"),
                         evaluation = "grid", grid.dim = c(3L, 3L), u.grid = u),
                    class = "npcopula")
  original <- serialize(fit, NULL)
  grid <- grid.eval(fit)
  expect_identical(grid$xgrid, dat[c("x", "y")])
  expect_length(grid$row.order, 9L)
  expect_identical(dim(grid$z), c(2L, 3L))
  expect_identical(grid$z, outer(c(.2, .5), c(.25, .5, .75), function(x,y) x + 2*y))
  expect_identical(serialize(fit, NULL), original)
  lower <- value - .1
  missing <- u$u1 == .5 & u$u2 == .5
  lower[missing] <- NA_real_
  expect_identical(which(is.na(grid.matrix(lower, grid))), 4L)
  conflict <- lower
  conflict[which(duplicated(grid$display.cell))[1L]] <- 17
  expect_error(grid.matrix(conflict, grid), "conflicting surface values", fixed = TRUE)
  conflict <- lower
  conflict[which(missing)[1L]] <- 0
  expect_error(grid.matrix(conflict, grid), "conflicting surface values", fixed = TRUE)
  incomplete <- fit
  incomplete$u.grid[1L, 1L] <- .9
  expect_error(grid.eval(incomplete), "not rectangular", fixed = TRUE)
})
