with_np_degree_bindings <- function(bindings, code) {
  code <- substitute(code)
  ns <- asNamespace("np")
  old <- lapply(names(bindings), function(name) get(name, envir = ns, inherits = FALSE))
  names(old) <- names(bindings)

  for (name in names(bindings)) {
    was_locked <- bindingIsLocked(name, ns)
    if (was_locked) {
      unlockBinding(name, ns)
    }
    assign(name, bindings[[name]], envir = ns)
    if (was_locked) {
      lockBinding(name, ns)
    }
  }

  on.exit({
    for (name in names(old)) {
      was_locked <- bindingIsLocked(name, ns)
      if (was_locked) {
        unlockBinding(name, ns)
      }
      assign(name, old[[name]], envir = ns)
      if (was_locked) {
        lockBinding(name, ns)
      }
    }
  }, add = TRUE)

  eval(code, envir = parent.frame())
}

capture_degree_messages_only <- function(expr) {
  messages <- character()
  withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  messages
}

degree_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

test_that("npregbw exhaustive degree search matches manual profile minimum", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(28)))
  dat$y <- dat$x^2 + rnorm(nrow(dat), sd = 0.05)

  bw0 <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 0L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  bw1 <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  auto <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expected_fval <- min(bw0$fval, bw1$fval)

  expect_s3_class(auto, "rbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, expected_fval + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))
  expect_identical(nrow(auto$degree.search$trace), auto$degree.search$n.unique)
  expect_identical(auto$degree.search$n.cached, auto$degree.search$n.visits - auto$degree.search$n.unique)

  manual <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  expect_null(manual$degree.search)
})

test_that("npregbw coordinate search can be exhaustively certified on a small grid", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(26),
    x2 = runif(26)
  )
  dat$y <- dat$x1 + dat$x2^2 + rnorm(nrow(dat), sd = 0.05)

  exhaustive <- np::npregbw(
    y ~ x1 + x2,
    data = dat,
    regtype = "lp",
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  coordinate <- np::npregbw(
    y ~ x1 + x2,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = TRUE,
    degree.restarts = 1L,
    degree.max.cycles = 4L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(coordinate$degree.search$mode, "coordinate")
  expect_true(isTRUE(coordinate$degree.search$completed))
  expect_true(isTRUE(coordinate$degree.search$certified))
  expect_equal(as.integer(coordinate$degree), as.integer(exhaustive$degree))
  expect_equal(coordinate$fval, exhaustive$fval, tolerance = 1e-10)
  expect_lte(coordinate$degree.search$best.fval, coordinate$degree.search$baseline.fval + 1e-10)
  expect_identical(nrow(coordinate$degree.search$trace), coordinate$degree.search$n.unique)
  expect_identical(coordinate$degree.search$n.cached, coordinate$degree.search$n.visits - coordinate$degree.search$n.unique)
})

test_that("npregbw automatic degree search enforces pilot guardrails", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(y = rnorm(24), x = runif(24))

  expect_error(
    np::npregbw(
      y ~ x,
      data = dat,
      regtype = "lc",
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "automatic degree search currently requires regtype='lp'"
  )

  expect_error(
    np::npregbw(
      y ~ x,
      data = dat,
      regtype = "lp",
      bernstein.basis = FALSE,
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 4L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "degree.max <= 3"
  )

  bw <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    bernstein.basis = FALSE,
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_false(isTRUE(bw$bernstein.basis))
})

test_that("npreg forwards automatic LP degree search through npregbw", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = runif(24))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.05)

  fit <- local({
    suppressPackageStartupMessages(library(np))
    npreg(
      y ~ x,
      data = dat,
      regtype = "lp",
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    )
  })

  expect_s3_class(fit, "npregression")
  expect_s3_class(fit$bws, "rbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})

test_that("internal degree search returns best-so-far on interrupt", {
  degree_search <- getFromNamespace(".np_degree_search", "np")

  eval_count <- 0L
  result <- degree_search(
    method = "exhaustive",
    candidates = list(0:1),
    baseline_degree = 0L,
    start_degree = 0L,
    eval_fun = function(degree) {
      eval_count <<- eval_count + 1L
      if (eval_count >= 2L) {
        stop(structure(list(message = "interrupt"), class = c("interrupt", "condition")))
      }
      list(
        objective = as.numeric(sum(degree)),
        payload = list(degree = as.integer(degree)),
        num.feval = 1L
      )
    },
    direction = "min",
    trace_level = "full"
  )

  expect_true(isTRUE(result$interrupted))
  expect_false(isTRUE(result$completed))
  expect_identical(as.integer(result$best$degree), 0L)
  expect_identical(as.integer(result$best_payload$degree), 0L)
})

test_that("restart starts are deterministic and RNG-independent", {
  restart_starts <- getFromNamespace(".np_degree_restart_starts", "np")

  set.seed(1)
  junk1 <- runif(10)
  starts1 <- restart_starts(
    candidates = list(0:2, 0:2),
    restarts = 4L,
    exclude = list(c(0L, 0L))
  )

  set.seed(999)
  junk2 <- runif(10)
  starts2 <- restart_starts(
    candidates = list(0:2, 0:2),
    restarts = 4L,
    exclude = list(c(0L, 0L))
  )

  expect_false(identical(junk1, junk2))
  expect_identical(starts1, starts2)
})

test_that("coordinate search skips incumbent cell revisits within a sweep", {
  degree_search <- getFromNamespace(".np_degree_search", "np")

  result <- degree_search(
    method = "coordinate",
    candidates = list(0:2, 0:2),
    baseline_degree = c(0L, 0L),
    start_degree = c(0L, 0L),
    restarts = 0L,
    max_cycles = 1L,
    eval_fun = function(degree) {
      list(
        objective = as.numeric(sum(degree)),
        payload = list(degree = as.integer(degree)),
        num.feval = 1L
      )
    },
    direction = "min",
    trace_level = "full"
  )

  expect_identical(result$n.unique, 5L)
  expect_identical(result$n.visits, 5L)
  expect_identical(result$n.cached, 0L)
  expect_identical(nrow(result$trace), result$n.unique)
  expect_false(any(result$trace$cached))
})

test_that("automatic exhaustive search emits a safety warning on large grids", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.degree.search.warn.grid = 3L)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(20),
    x2 = runif(20)
  )
  dat$y <- dat$x1 + dat$x2 + rnorm(nrow(dat), sd = 0.05)

  expect_warning(
    np::npregbw(
      y ~ x1 + x2,
      data = dat,
      regtype = "lp",
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "exhaustive degree search will evaluate 4 degree combinations"
  )
})

test_that("automatic exhaustive search honors internal hard grid limits", {
  old_opts <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.degree.search.warn.grid = Inf,
    np.degree.search.max.grid = 3L
  )
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(20),
    x2 = runif(20)
  )
  dat$y <- dat$x1 + dat$x2 + rnorm(nrow(dat), sd = 0.05)

  expect_error(
    np::npregbw(
      y ~ x1 + x2,
      data = dat,
      regtype = "lp",
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "exceeding the configured limit of 3"
  )
})

test_that("automatic degree search emits staged progress output", {
  old_opts <- options(np.messages = TRUE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(20)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.05)

  msgs <- with_np_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 20, by = 0.5))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("np"), inherits = FALSE)(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "exhaustive",
        degree.min = 0L,
        degree.max = 1L,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L
      )
    )
  )

  expect_true(any(grepl("Automatic polynomial degree search baseline \\(0\\)", msgs)))
  expect_true(any(grepl("Selecting polynomial degree and bandwidth", msgs)))
  expect_true(any(grepl("exhaustive", msgs)))
  expect_true(any(grepl("best (", msgs, fixed = TRUE)))

  coord_msgs <- with_np_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 40, by = 0.5))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("np"), inherits = FALSE)(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "coordinate",
        degree.min = 0L,
        degree.max = 1L,
        degree.verify = TRUE,
        degree.max.cycles = 2L,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L
      )
    )
  )

  expect_true(any(grepl("Coordinate automatic polynomial degree search over 0:1", coord_msgs)))
  expect_true(any(grepl("max 3 search evaluations", coord_msgs)))
  expect_true(any(grepl("Exhaustively certifying automatic polynomial degree search over 2 degree combinations \\(re-optimizing bandwidths\\)", coord_msgs)))
})
