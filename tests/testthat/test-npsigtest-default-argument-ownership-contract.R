t2_npsig_default_fixture <- function() {
  n <- 16L
  xdat <- data.frame(x = seq(-1, 1, length.out = n))
  ydat <- sin(1.2 * xdat$x) + seq_len(n) / 100
  list(xdat = xdat, ydat = ydat)
}

test_that("default-route test argument ownership is explicit and signature-locked", {
  pkg <- environmentName(environment(npsigtest))
  owners <- getFromNamespace(".np_npsig_default_test_args", pkg)
  method <- getS3method("npsigtest", "default")

  expect_identical(
    owners,
    c("B", "boot.method", "boot.type", "pivot", "joint", "index", "random.seed")
  )
  expect_identical(names(formals(npsigtest)), c("bws", "...", "B"))
  expect_identical(names(formals(method)), c("bws", "xdat", "ydat", "..."))
})

test_that("default route partitions every single and paired test control", {
  pkg <- environmentName(environment(npsigtest))
  owners <- getFromNamespace(".np_npsig_default_test_args", pkg)
  method <- getS3method("npsigtest", "default")
  fixture <- t2_npsig_default_fixture()
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .np_eval_bw_call = function(call, caller_env) {
      captured$bw.names <- names(call)
      structure(list(marker = "bw"), class = "rbandwidth")
    },
    npsigtest = function(bws, ..., B = 399) {
      captured$test <- c(list(B = B), list(...))
      list(call = NULL, bws = bws)
    },
    .package = pkg
  )

  values <- list(
    B = 9L,
    boot.method = "iid",
    boot.type = "I",
    pivot = FALSE,
    joint = TRUE,
    index = 1L,
    random.seed = 19L
  )
  selections <- c(
    lapply(owners, function(owner) owner),
    combn(owners, 2L, simplify = FALSE),
    list(owners)
  )

  for (selection in selections) {
    controls <- values[selection]
    result <- do.call(
      method,
      c(list(bws = 0.45, xdat = fixture$xdat, ydat = fixture$ydat), controls)
    )
    expect_s3_class(result$bws, "rbandwidth")
    expect_length(intersect(captured$bw.names, owners), 0L)
    expect_true(all(selection %in% names(captured$test)),
                info = paste(selection, collapse = "+"))
    for (owner in selection)
      expect_identical(captured$test[[owner]], values[[owner]],
                       info = paste(selection, collapse = "+"))
  }
})

test_that("default route evaluates bandwidth work before deferred test controls", {
  pkg <- environmentName(environment(npsigtest))
  method <- getS3method("npsigtest", "default")
  fixture <- t2_npsig_default_fixture()
  events <- character()

  local_mocked_bindings(
    .np_eval_bw_call = function(call, caller_env) {
      eval(call$bws, envir = caller_env)
      structure(list(marker = "bw"), class = "rbandwidth")
    },
    npsigtest = function(bws, ..., B = 399) {
      force(B)
      dots <- list(...)
      force(dots$random.seed)
      list(call = NULL, bws = bws, B = B, dots = dots)
    },
    .package = pkg
  )

  result <- method(
    bws = {
      events <- c(events, "bandwidth")
      0.45
    },
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    B = {
      events <- c(events, "B")
      9L
    },
    random.seed = {
      events <- c(events, "random.seed")
      23L
    }
  )

  expect_identical(events, c("bandwidth", "B", "random.seed"))
  expect_identical(result$B, 9L)
  expect_identical(result$dots$random.seed, 23L)
})

test_that("default route executes numeric-start and computed-bandwidth calls", {
  fixture <- t2_npsig_default_fixture()

  numeric.start <- npsigtest(
    bws = 0.45,
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    B = 9L,
    boot.method = "iid",
    boot.type = "I",
    pivot = FALSE,
    joint = TRUE,
    index = 1L,
    random.seed = 31L
  )
  expect_s3_class(numeric.start, "sigtest")
  expect_identical(numeric.start$boot.num, 9L)
  expect_identical(numeric.start$ixvar, 1L)
  expect_identical(numeric.start$joint, TRUE)

  computed <- npsigtest(
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    B = 9L,
    boot.method = "iid",
    boot.type = "I",
    pivot = FALSE,
    joint = TRUE,
    index = 1L,
    random.seed = 31L,
    nmulti = 1L,
    itmax = 1L
  )
  expect_s3_class(computed, "sigtest")
  expect_s3_class(computed$bws, "rbandwidth")
  expect_identical(computed$boot.num, 9L)
})

test_that("default route preserves regression bandwidth arguments", {
  pkg <- environmentName(environment(npsigtest))
  method <- getS3method("npsigtest", "default")
  fixture <- t2_npsig_default_fixture()
  captured <- NULL

  local_mocked_bindings(
    .np_eval_bw_call = function(call, caller_env) {
      captured <<- names(call)
      structure(list(marker = "bw"), class = "rbandwidth")
    },
    npsigtest = function(bws, ..., B = 399) list(call = NULL, bws = bws),
    .package = pkg
  )
  method(
    bws = 0.45,
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    B = 9L,
    regtype = "ll",
    ckertype = "epanechnikov"
  )
  expect_true(all(c("regtype", "ckertype") %in% captured))
})

test_that("default route preserves strict diagnostics", {
  method <- getS3method("npsigtest", "default")
  fixture <- t2_npsig_default_fixture()
  expect_error(
    method(
      bws = 0.45, xdat = fixture$xdat, ydat = fixture$ydat,
      B = 9L, definitely.not.an.np.argument = TRUE
    ),
    "npregbw(): unused argument 'definitely.not.an.np.argument'",
    fixed = TRUE
  )
  expect_error(
    npsigtest(
      bws = 0.45, xdat = fixture$xdat, ydat = fixture$ydat,
      boot.num = 9L
    ),
    "use 'B' instead",
    fixed = TRUE
  )
})
