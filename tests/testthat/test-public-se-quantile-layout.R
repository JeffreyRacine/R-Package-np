test_that("quantile MPI payload layout follows the uncertainty demand", {
  layout <- getFromNamespace(".npqreg_tau_layout", "npRmpi")
  delta.matrix <- getFromNamespace(".npqreg_quantile_delta_matrix", "npRmpi")
  decode <- getFromNamespace(".npqreg_fit_tau_vector_from_parallel_matrix", "npRmpi")
  for (n in c(0L, 3L)) for (p in c(1L, 2L)) {
    for (tau in list(.3, c(.25, .7))) for (g in c(FALSE, TRUE)) {
      for (s in c(FALSE, TRUE)) {
        mask <- layout(p, g, s)
        expect_identical(mask$width, 1L + as.integer(s) + p * as.integer(g) * (1L + as.integer(s)))
        pieces <- lapply(seq_along(tau), function(j) {
          delta <- list(quanterr = seq_len(n) / 100 + j,
                        quantgrad = matrix(seq_len(n * p) + j, n, p),
                        quantgerr = matrix(seq_len(n * p) / 10 + j, n, p))
          if (!g && !s) matrix(seq_len(n) + j / 10, ncol = 1L)
          else cbind(seq_len(n) + j / 10, delta.matrix(delta, g, s))
        })
        packed <- do.call(cbind, pieces)
        out <- decode(packed, tau, g, se = s, expected.grad.cols = p)
        for (j in seq_along(tau)) {
          actual <- if (length(tau) == 1L) out$yq else out$yq[, j]
          expect_identical(as.double(actual), seq_len(n) + j / 10)
        }
        expect_identical(is.null(out$yqerr), !s)
        expect_identical(is.null(out$yqgerr), !s)
        if (g) {
          expected.dim <- if (length(tau) == 1L) c(n, p) else c(n, p, length(tau))
          expect_identical(dim(out$yqgrad), expected.dim)
          if (s) expect_identical(dim(out$yqgerr), expected.dim)
          for (j in seq_along(tau)) {
            actual <- if (length(tau) == 1L) out$yqgrad else out$yqgrad[, , j]
            expect_identical(as.double(actual), as.double(seq_len(n * p) + j))
          }
        }
      }
    }
  }
  expect_error(decode(matrix(1, 2L, 3L), c(.25, .75), FALSE, se = FALSE),
               "malformed", fixed = TRUE)
  expect_error(decode(matrix(1, 2L, 5L), .5, TRUE, grad.names = "x", se = TRUE),
               "malformed", fixed = TRUE)
  expect_error(decode(matrix(1, 2L, 4L), .5, TRUE, se = FALSE,
                      expected.grad.cols = 2L), "malformed", fixed = TRUE)
})
