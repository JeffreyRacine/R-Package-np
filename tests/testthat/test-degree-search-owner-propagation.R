test_that("all seven cell-search owners propagate an evaluator failure", {

  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 18L; i <- seq_len(n)
  d <- data.frame(x = sin(i * sqrt(2)) + i / 30,
                  z = cos(i * sqrt(3)) + i / 50)
  d$y <- 1 + d$x * (1 + d$z) + sin(i) / 10
  original <- getFromNamespace(".np_degree_search", "np")
  marker <- errorCondition("diagnostic owner evaluator failure",
    class = "np_test_owner_error", call = quote(fixed_degree_owner()), token = 808L)
  calls <- 0L; valid <- 0L
  wrapped <- function(..., eval_fun) {
    calls <<- calls + 1L
    evaluate <- function(degree) {
      if (degree[[1L]] == 1L) stop(marker)
      out <- eval_fun(degree)
      valid <<- valid + 1L
      out
    }
    original(..., eval_fun = evaluate)
  }
  local_mocked_bindings(.np_degree_search = wrapped, .package = "np")
  for (owner in c("npregbw", "npcdensbw", "npcdistbw", "npscoefbw",
                  "npplregbw", "npindexbw", "nplsqregbw")) {
    args <- list(formula = if (owner %in% c("npscoefbw", "npplregbw"))
                   y ~ x | z else y ~ x,
                 data = d, regtype = "lp", nomad = FALSE,
                 degree.select = "exhaustive", search.engine = "cell",
                 degree.min = 0L, degree.max = 2L, degree.start = 0L,
                 nmulti = 1L)
    if (owner %in% c("npindexbw", "npscoefbw")) args$optim.maxit <- 20L
    else {args$itmax <- 2L; args$powell.remin <- FALSE}
    if (owner == "npcdistbw") args$ngrid <- 5L
    if (owner == "npindexbw") args$formula <- y ~ x + z
    if (owner == "nplsqregbw") {
      args$bws <- args$formula
      args$formula <- NULL
      args$nomad <- "auto"
      args$scale <- rep(1, n)
    }
    npseed(934)
    before <- calls
    condition <- tryCatch(do.call(getFromNamespace(owner, "np"), args), error = identity)
    expect_identical(condition, marker, info = owner)
    expect_equal(calls, before + 1L, info = owner)
  }
})
