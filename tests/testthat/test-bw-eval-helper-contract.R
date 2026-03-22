test_that(".np_try_eval_in_frames resolves symbols from eval env first", {
  env <- new.env(parent = emptyenv())
  env$x <- 11L
  out <- .np_try_eval_in_frames(quote(x), eval_env = env)
  expect_true(out$ok)
  expect_identical(out$value, 11L)
})

test_that(".np_try_eval_in_frames preserves NULL symbol values", {
  env <- new.env(parent = emptyenv())
  env$x <- NULL
  out <- .np_try_eval_in_frames(quote(x), eval_env = env)
  expect_true(out$ok)
  expect_null(out$value)
})

test_that(".np_try_eval_in_frames resolves symbols from caller frames", {
  val <- local({
    x <- 17L
    .np_try_eval_in_frames(quote(x), eval_env = new.env(parent = emptyenv()))
  })
  expect_true(val$ok)
  expect_identical(val$value, 17L)
})

test_that(".np_try_eval_in_frames can disable caller frame fallback", {
  out <- local({
    x <- 19L
    .np_try_eval_in_frames(quote(x), eval_env = new.env(parent = emptyenv()), search_frames = FALSE)
  })
  expect_false(out$ok)
  expect_true(inherits(out$error, "error"))
})

test_that(".np_try_eval_in_frames returns an error object when resolution fails", {
  out <- .np_try_eval_in_frames(quote(np_missing_symbol_contract), eval_env = new.env(parent = emptyenv()))
  expect_false(out$ok)
  expect_true(inherits(out$error, "error"))
  expect_match(conditionMessage(out$error), "np_missing_symbol_contract")
})

test_that(".np_try_eval_in_frames evaluates non-symbol expressions in caller frames", {
  out <- local({
    x <- 4L
    .np_try_eval_in_frames(quote(x + 1L), eval_env = new.env(parent = emptyenv()))
  })
  expect_true(out$ok)
  expect_identical(out$value, 5L)
})

test_that(".np_try_eval_in_frames can disable caller fallback for non-symbol expressions", {
  out <- local({
    x <- 4L
    .np_try_eval_in_frames(
      quote(x + 1L),
      eval_env = new.env(parent = emptyenv()),
      search_frames = FALSE
    )
  })
  expect_false(out$ok)
  expect_true(inherits(out$error, "error"))
})

test_that(".np_try_eval_in_frames honors enclos for non-symbol expressions", {
  eval_env <- list()
  enclos <- new.env(parent = baseenv())
  enclos$z <- 8L
  out <- .np_try_eval_in_frames(
    quote(z + 2L),
    eval_env = eval_env,
    enclos = enclos,
    search_frames = FALSE
  )
  expect_true(out$ok)
  expect_identical(out$value, 10L)
})

test_that(".np_eval_bw_call reports underlying evaluation errors", {
  call_obj <- as.call(list(as.name("identity"), as.name("np_missing_bw_arg_contract")))
  environment(call_obj) <- new.env(parent = baseenv())

  expect_error(
    .np_eval_bw_call(call_obj, caller_env = new.env(parent = baseenv())),
    "np_missing_bw_arg_contract"
  )
})

test_that(".np_eval_bw_call evaluates successful calls in caller context", {
  out <- local({
    z <- 23L
    .np_eval_bw_call(quote(identity(z)), caller_env = environment())
  })

  expect_identical(out, 23L)
})

test_that(".np_bw_call_uses_nomad_degree_search tolerates missing dummy formals", {
  out <- local({
    x1 <- rnorm(10)
    x2 <- rnorm(10)
    y <- x1 + rnorm(10)
    .np_bw_call_uses_nomad_degree_search(
      quote(npregbw(y ~ x1 + x2, nmulti = 2)),
      caller_env = environment()
    )
  })

  expect_false(out)
})

test_that(".np_bw_call_uses_nomad_degree_search classifies supported cases narrowly", {
  out <- local({
    x1 <- rnorm(10)
    x2 <- rnorm(10)
    y <- x1 + rnorm(10)

    c(
      lc = .np_bw_call_uses_nomad_degree_search(
        quote(npregbw(y ~ x1 + x2, regtype = "lc", nmulti = 2)),
        caller_env = environment()
      ),
      ll = .np_bw_call_uses_nomad_degree_search(
        quote(npregbw(y ~ x1 + x2, regtype = "ll", nmulti = 2)),
        caller_env = environment()
      ),
      lp_manual = .np_bw_call_uses_nomad_degree_search(
        quote(npregbw(y ~ x1 + x2, regtype = "lp", degree = c(1L, 1L), nmulti = 2)),
        caller_env = environment()
      ),
      lp_cell = .np_bw_call_uses_nomad_degree_search(
        quote(npregbw(
          y ~ x1 + x2,
          regtype = "lp",
          degree.select = "coordinate",
          search.engine = "cell",
          nmulti = 2
        )),
        caller_env = environment()
      ),
      lp_nomad = .np_bw_call_uses_nomad_degree_search(
        quote(npregbw(
          y ~ x1 + x2,
          regtype = "lp",
          degree.select = "coordinate",
          search.engine = "nomad",
          nmulti = 2
        )),
        caller_env = environment()
      ),
      lp_nomad_powell = .np_bw_call_uses_nomad_degree_search(
        quote(npregbw(
          y ~ x1 + x2,
          regtype = "lp",
          degree.select = "exhaustive",
          search.engine = "nomad+powell",
          nmulti = 2
        )),
        caller_env = environment()
      )
    )
  })

  expect_false(out[["lc"]])
  expect_false(out[["ll"]])
  expect_false(out[["lp_manual"]])
  expect_false(out[["lp_cell"]])
  expect_true(out[["lp_nomad"]])
  expect_true(out[["lp_nomad_powell"]])
})

test_that(".np_bw_dispatch_target distinguishes formula, named-data, and bws entry shapes", {
  formula.out <- local({
    y <- rnorm(5)
    x <- rnorm(5)
    .np_bw_dispatch_target(
      dots = list(quote(y ~ x)),
      data_arg_names = c("xdat", "ydat"),
      eval_env = environment()
    )
  })
  expect_true(inherits(formula.out, "formula"))

  data.out <- .np_bw_dispatch_target(
    dots = list(xdat = quote(1:3), ydat = quote(4:6)),
    data_arg_names = c("xdat", "ydat"),
    eval_env = new.env(parent = baseenv())
  )
  expect_null(data.out)

  bw.out <- .np_bw_dispatch_target(
    dots = list(bws = quote(structure(list(), class = "rbandwidth")), xdat = quote(1:3)),
    data_arg_names = c("xdat", "ydat"),
    eval_env = new.env(parent = baseenv())
  )
  expect_s3_class(bw.out, "rbandwidth")
})

test_that(".np_bw_formula_from_call extracts the formula argument call when present", {
  out <- .np_bw_formula_from_call(quote(npregbw(y ~ x1 + x2, nmulti = 2)))
  expect_true(is.call(out))
  expect_identical(deparse(out), "y ~ x1 + x2")

  expect_null(.np_bw_formula_from_call(quote(npregbw(xdat = x1, ydat = y, nmulti = 2))))
})
