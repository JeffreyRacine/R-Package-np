proper_empty_ns <- function() asNamespace("npRmpi")

proper_empty_object <- function(cdf = FALSE, x = rep(c(1, 2), each = 3L),
                                y = rep(c(0, .5, 1), 2L)) {
  out <- list(gradients = FALSE, trainiseval = FALSE,
    yndim = 1L, yncon = 1L, ynord = 0L, ynuno = 0L,
    xeval = data.frame(x = x), yeval = data.frame(y = y),
    bws = list(iycon = 1L, cykerlb = 0, cykerub = 1))
  out[[if (cdf) "condist" else "condens"]] <- rep(c(-.1, .9, .2), length.out = length(x))
  out
}

test_that("proper projections skip only complete typed-empty X slices", {
  ns <- proper_empty_ns()
  for (cdf in c(FALSE, TRUE)) {
    prefix <- if (cdf) ".np_condist_" else ".np_condens_"
    field <- if (cdf) "condist" else "condens"
    apply <- get(paste0(prefix, "apply_proper_grid"), ns)
    prepare <- get(paste0(prefix, "prepare_proper_plan"), ns)
    project <- get(paste0(prefix, "project_values_with_plan"), ns)
    healthy <- proper_empty_object(cdf)
    before <- apply(healthy)
    plan <- prepare(healthy)
    mixed <- healthy
    mixed[[field]][4:6] <- NA_real_
    attr(mixed, ".np.empty.base.rows") <- c(0L, 0L, 0L, 1L, 1L, 1L)
    expect_identical(prepare(mixed), plan)
    got <- apply(mixed)
    expect_true(got$applied)
    expect_identical(got[[field]][1:3], before[[field]][1:3])
    expect_true(all(is.na(got[[field]][4:6])))
    expect_identical(got$proper.info$empty.slice.count, 1L)
    expect_true(is.na(got$proper.info$projection.distance.l2[2L]))
    if (cdf) {
      expect_true(is.na(got$proper.info$monotone.violations.raw[2L]))
      expect_true(all(is.na(got$proper.info$range.raw[2L, ])))
    } else {
      expect_true(is.na(got$proper.info$negative.count.raw[2L]))
      expect_true(is.na(got$proper.info$integral.raw[2L]))
    }
    expect_null(before$proper.info$empty.slice.count)
    expect_null(attr(plan, ".np.empty.base.rows"))
    expect_null(plan[["empty.slices", exact = TRUE]])
    expect_identical(project(healthy[[field]], plan),
                     project(healthy[[field]], plan, empty.slices = NULL))
    expect_error(project(matrix(mixed[[field]], nrow = 1L), plan,
                         empty.slices = c(0L, 1L)),
                 "matrix/bootstrap projection", fixed = TRUE)
    expect_error(project(matrix(healthy[[field]], nrow = 1L), plan,
                         empty.slices = c(0L, 0L)),
                 "matrix/bootstrap projection", fixed = TRUE)
    expect_error(project(mixed[[field]], plan), "finite", fixed = TRUE)
    expect_error(project(healthy[[field]], plan, empty.slices = c(0L, 1L)),
                 "contradict", fixed = TRUE)

    empty <- mixed
    empty[[field]][] <- NA_real_
    attr(empty, ".np.empty.base.rows")[] <- 1L
    all <- apply(empty)
    expect_false(all$applied)
    expect_true(all$proper.info$supported)
    expect_identical(all$reason, "all_slices_empty")
    expect_identical(all$proper.info$empty.slice.count, 2L)
    expect_true(all(is.na(all$proper.info$projection.distance.l2)))
    env <- new.env(parent = ns)
    dispatcher <- get(paste0(prefix, "apply_proper"), ns)
    environment(dispatcher) <- env
    assign(paste0(prefix, "apply_proper_slice"), function(...) stop("unexpected slice retry"), env)
    expect_identical(dispatcher(empty, proper.control = list(mode = "slice")), all)

    broken <- mixed
    attr(broken, ".np.empty.base.rows") <- as.double(attr(broken, ".np.empty.base.rows"))
    expect_error(apply(broken), "empty-row metadata", fixed = TRUE)
    attr(broken, ".np.empty.base.rows") <- c(0L, 0L, 0L, 1L, 0L, 1L)
    expect_error(apply(broken), "within a conditional proper X slice", fixed = TRUE)
    attr(broken, ".np.empty.base.rows") <- NULL
    attr(broken, ".np.empty.rows") <- c(0L, 0L, 0L, 1L, 1L, 1L)
    expect_error(apply(broken), "finite", fixed = TRUE)
  }
})

test_that("recursive proper slices retain purpose and original query metadata", {
  ns <- proper_empty_ns()
  for (cdf in c(FALSE, TRUE)) {
    prefix <- if (cdf) ".np_condist_" else ".np_condens_"
    field <- if (cdf) "condist" else "condens"
    env <- new.env(parent = ns)
    calls <- list()
    flagged <- TRUE
    child <- function(exdat, eydat, ...) {
      calls[[length(calls) + 1L]] <<- list(...)
      obj <- proper_empty_object(cdf, exdat$x, eydat[[1L]])
      obj[[field]] <- if (cdf) eydat[[1L]] else rep(1, nrow(exdat))
      if (flagged) {
        flags <- as.integer(exdat$x == .7)
        obj[[field]][flags == 1L] <- NA_real_
        attr(obj, ".np.empty.base.rows") <- flags
        attr(obj, ".np.empty.rows") <- flags
      }
      obj
    }
    assign(if (cdf) "npcdist" else "npcdens", child, env)
    assign(paste0(prefix, "slice_dispatch_enabled"), function() TRUE, env)
    owner <- get(paste0(prefix, "apply_proper_slice"), ns)
    environment(owner) <- env
    obj <- proper_empty_object(cdf, c(.7, .3, .7), c(.2, .4, .8))
    obj[[field]][c(1, 3)] <- NA_real_
    attr(obj, ".np.empty.base.rows") <- c(1L, 0L, 1L)
    context <- list(txdat = data.frame(x = c(.3, .5, .7)),
      tydat = data.frame(y = c(.1, .5, .9)), exdat = obj$xeval,
      eydat = obj$yeval, allow.external = TRUE)
    ctrl <- list(mode = "slice", slice.grid.size = 3L)
    got <- owner(obj, slice.context = context, proper.control = ctrl)
    expect_true(got$applied)
    expect_identical(attr(got, ".np.empty.base.rows"), c(1L, 0L, 1L))
    expect_identical(attr(got, ".np.empty.rows"), c(1L, 0L, 1L))
    expect_true(all(is.na(got[[field]][c(1, 3)])))
    expect_true(is.finite(got[[field]][2L]))
    expect_identical(calls[[1L]][[".np.require.complete"]], FALSE)
    expect_identical(calls[[1L]][[".np.defer.empty.rows"]], TRUE)
    expect_identical(got$proper.info$empty.slice.count, 1L)

    bad <- obj
    attr(bad, ".np.empty.base.rows") <- NULL
    expect_error(owner(bad, slice.context = context, proper.control = ctrl),
                 "changed base-X support", fixed = TRUE)

    flagged <- FALSE
    healthy <- proper_empty_object(cdf, c(.7, .3, .7), c(.2, .4, .8))
    context$allow.external <- NULL
    expect_true(owner(healthy, slice.context = context, proper.control = ctrl)$applied)
    expect_identical(tail(calls, 1L)[[1L]][[".np.require.complete"]], TRUE)
  }
})
