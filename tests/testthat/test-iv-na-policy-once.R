test_that("IV formula NA factories run once while policies handle each frame", {
  pkg <- environmentName(environment(npregiv))
  methods <- lapply(c("npregiv", "npregivderiv"), getS3method, class = "formula")
  names(methods) <- c("npregiv", "npregivderiv")
  d <- data.frame(y = seq_len(10L), z = seq_len(10L)/10, w = cos(seq_len(10L)))
  for (name in c("npregiv", "npregivderiv")) for (has.newdata in c(FALSE, TRUE)) {
    method <- methods[[name]]
    hits <- c(data = 0L, factory = 0L, policy = 0L, newdata = 0L)
    events <- character()
    data <- function() {hits["data"] <<- hits["data"] + 1L; events <<- c(events, "data"); d}
    newdata <- function() {hits["newdata"] <<- hits["newdata"] + 1L; events <<- c(events, "newdata"); d[1:4, ]}
    factory <- function() {
      hits["factory"] <<- hits["factory"] + 1L; events <<- c(events, "factory")
      function(object) {hits["policy"] <<- hits["policy"] + 1L; events <<- c(events, "policy"); stats::na.omit(object)}
    }
    native <- function(y, z, w, ...) list(y = y, z = z, w = w, args = list(...))
    local_mocked_bindings(npregiv = native, npregivderiv = native, .package = pkg)
    out <- if (has.newdata) method(y ~ z | w, data = data(), newdata = newdata(), na.action = factory()) else
      method(y ~ z | w, data = data(), na.action = factory())
    expect_identical(unname(hits), c(1L, 1L, if (has.newdata) 2L else 1L, as.integer(has.newdata)))
    expect_identical(events, c("data", "factory", "policy", if (has.newdata) c("newdata", "policy")))
    expect_identical(unname(out$y), d$y)
    expect_identical(out$call$na.action, quote(factory()))
    expect_false(grepl("evalq|<environment", paste(deparse(out$call), collapse = " ")))
  }
})

test_that("IV NA policy reuse retains missing, character and function policies", {
  pkg <- environmentName(environment(npregiv))
  methods <- lapply(c("npregiv", "npregivderiv"), getS3method, class = "formula")
  names(methods) <- c("npregiv", "npregivderiv")
  d <- data.frame(y = seq_len(10L), z = seq_len(10L)/10, w = cos(seq_len(10L)))
  d$z[3L] <- NA_real_
  for (name in c("npregiv", "npregivderiv")) {
    method <- methods[[name]]
    native <- function(y, z, w, ...) list(y = y, z = z, w = w, args = list(...))
    local_mocked_bindings(npregiv = native, npregivderiv = native, .package = pkg)
    a <- method(y ~ z | w, data = d, newdata = d[1:4, ], na.action = stats::na.omit)
    b <- method(y ~ z | w, data = d, newdata = d[1:4, ], na.action = "na.omit")
    old <- options(na.action = "na.omit"); on.exit(options(old), add = TRUE)
    c <- method(y ~ z | w, data = d, newdata = d[1:4, ])
    expect_identical(a$y, b$y); expect_identical(a$y, c$y)
    expect_identical(a$args, b$args); expect_identical(a$args, c$args)
    expect_identical(as.integer(a$omit), 3L)
    expect_identical(as.integer(a$eval.omit), 3L)
  }
})


test_that("IV missing NA promises retain model-frame default policy through wrappers", {
  pkg <- environmentName(environment(npregiv))
  methods <- lapply(c("npregiv", "npregivderiv"), getS3method, class = "formula")
  d <- data.frame(y = seq_len(10L), z = seq_len(10L)/10, w = cos(seq_len(10L)))
  d$z[3L] <- NA_real_
  old <- options(na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  native <- function(y, z, w, ...) list(y = y, z = z, w = w, args = list(...))
  local_mocked_bindings(npregiv = native, npregivderiv = native, .package = pkg)
  for (method in methods) for (with.newdata in c(FALSE, TRUE)) {
    expected <- if (with.newdata) method(y ~ z | w, data = d, newdata = d[1:4, ]) else
      method(y ~ z | w, data = d)
    wrap <- function(na.action) {
      if (with.newdata) method(y ~ z | w, data = d, newdata = d[1:4, ], na.action = na.action) else
        method(y ~ z | w, data = d, na.action = na.action)
    }
    outer <- function(na.action) wrap(na.action = na.action)
    for (actual in list(wrap(), outer())) {
      expect_identical(actual$y, expected$y)
      expect_identical(actual$args, expected$args)
      expect_identical(actual$omit, expected$omit)
      expect_identical(actual$eval.omit, expected$eval.omit)
    }
  }
})
