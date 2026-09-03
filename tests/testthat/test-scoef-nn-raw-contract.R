test_that("smooth-coefficient raw admission excludes finite penalties", {
  valid <- getFromNamespace(".npscoefbw_raw_objective_valid", "npRmpi")
  penalty <- sqrt(.Machine$double.xmax)
  for(value in c(-1, 0, 1, penalty / 2)) expect_true(valid(value))
  for(value in list(penalty, Inf, -Inf, NA_real_, NaN, numeric(), c(1, 2), "1"))
    expect_false(valid(value))
})

test_that("smooth-coefficient radius admission respects typed training coordinates", {
  check <- getFromNamespace(".npscoef_nn_assert_training_radius", "npRmpi")
  z <- data.frame(u = factor(c("a", "b", "a", "b", "a")), x = c(0, 0, 0, 1, 2))
  for(type in c("generalized_nn", "adaptive_nn")) {
    bw <- list(type = type, icon = c(FALSE, TRUE), bw = c(0, 2))
    error <- tryCatch(check(bw, z, "test owner"), np_nn_candidate_invalid = identity)
    expect_s3_class(error, "np_nn_candidate_invalid")
    expect_identical(error$owner, "test owner")
    expect_identical(error$point, c(0, 2))
    bw$bw[2L] <- 3
    expect_true(check(bw, z, "test owner"))
    bw$bw[2L] <- 4
    expect_true(check(bw, z, "test owner"))
    fixed <- bw; fixed$type <- "fixed"; fixed$bw <- c(0, 0)
    expect_true(check(fixed, z, "test owner"))
    categorical <- bw; categorical$icon[] <- FALSE
    expect_true(check(categorical, NULL, "test owner"))
    expect_error(check(bw, z[, 1L, drop = FALSE], "test owner"),
                 "invalid coordinate map", fixed = TRUE)
    # Domain admission, not radius admission, rejects an undecodable k.
    bw$bw[2L] <- 1
    expect_true(check(bw, z, "test owner"))
  }
})

test_that("fixed-owner wrappers penalize only typed NN candidate invalidity", {
  find.definition <- function(e, name) {
    if(!is.call(e)) return(NULL)
    if(length(e) == 3L && identical(e[[1L]], as.name("<-")) &&
       identical(e[[2L]], as.name(name))) return(e[[3L]])
    for(child in as.list(e)[-1L]) {
      found <- find.definition(child, name)
      if(!is.null(found)) return(found)
    }
    NULL
  }
  owner <- body(getFromNamespace("npscoefbw.scbandwidth", "npRmpi"))
  for(name in c("overall.cv.ls", "partial.cv.ls")) {
    definition <- find.definition(owner, name)
    catcher <- definition[[3L]][[2L]]
    expect_identical(catcher[[1L]], as.name("tryCatch"))
    expect_identical(names(as.list(catcher))[-c(1L, 2L)], "np_nn_candidate_invalid")
    scope <- new.env(parent = baseenv())
    scope$maxPenalty <- sqrt(.Machine$double.xmax)
    scope$steps <- 0L
    scope$cv_progress_step <- scope$partial_progress_step <- function(...) {
      scope$steps <- scope$steps + 1L
    }
    raw.name <- paste0(name, ".raw")
    scope[[raw.name]] <- function(...) stop("unrelated evaluator sentinel")
    wrapper <- eval(definition, scope)
    args <- if(name == "overall.cv.ls") list(c(2, 3)) else list(c(2, 3), 1L)
    expect_error(do.call(wrapper, args), "unrelated evaluator sentinel", fixed = TRUE)
    scope[[raw.name]] <- function(...) stop(structure(
      list(message = "typed invalid sentinel", call = NULL),
      class = c("np_nn_candidate_invalid", "error", "condition")))
    expect_identical(do.call(wrapper, args), scope$maxPenalty)
    expect_identical(scope$steps, 1L)
  }
})
