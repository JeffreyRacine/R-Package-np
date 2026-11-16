test_that("new conditional publication requires all literal slots to be reconstructible", {
  complete <- getFromNamespace(".npRmpi_bw_literals_complete", "npRmpi")
  empty <- list(formula = list(call = quote(npcdensbw(xdat = NULL, ydat = NULL))),
                selector.native = list(xdat = NULL, ydat = NULL),
                literal = c(f.xdat = "C0", s.ydat = "C0"))
  expect_true(complete(NULL))
  expect_true(complete(empty))
  literal.call <- as.call(list(quote(npcdensbw), xdat = 1:2))
  expect_true(complete(list(formula.cache = list(call = literal.call))))
  expect_true(complete(list(formula = list(call.cache = literal.call))))
  for (target in c("f", "s")) for (role in c("xdat", "ydat")) {
    for (value in list(seq_len(240L), sort(sin(seq_len(240L))),
                       data.frame(x = seq_len(240L)),
                       structure(c(.1, .2), class = "caller_owned"))) {
      metadata <- empty
      if (target == "f") metadata$formula$call[role] <- list(value)
      else metadata$selector.native[role] <- list(value)
      before <- serialize(metadata, NULL)
      expect_false(complete(metadata))
      expect_identical(serialize(metadata, NULL), before)
    }
    metadata <- empty
    if (target == "f") metadata$formula$call[role] <- list(quote(stop("not evaluated")))
    else metadata$selector.native[role] <- list(quote(stop("not evaluated")))
    expect_true(complete(metadata))
  }
})

test_that("failed literal probes preserve compact and sorted input representations", {
  ns <- asNamespace("npRmpi")
  code <- get(".npRmpi_bw_literal_code", ns)
  for (kind in c("integer", "double", "sorted", "sorted.frame", "deferred")) {
    x <- switch(kind, integer = seq_len(240L), double = as.double(seq_len(240L)),
      sorted = sort(sin(seq_len(240L) * .7)),
      sorted.frame = data.frame(x = sort(sin(seq_len(240L) * .7))),
      deferred = as.character(seq_len(240L)))
    role <- get("toFrame", ns)(x)
    # Incumbent untagging can wrap a compact column without expanding it.
    role <- get(".npRmpi_autodispatch_untag", ns)(role)
    literal.before <- serialize(x, NULL)
    role.before <- serialize(role, NULL)
    expect_null(code(x, role), info = kind)
    expect_identical(serialize(x, NULL), literal.before, info = kind)
    expect_identical(serialize(role, NULL), role.before, info = kind)
  }
})

test_that("matrix reconstruction preserves a certified compact wrapped role", {
  jit <- compiler::enableJIT(0L)
  withr::defer(compiler::enableJIT(jit))
  wrap <- function(x) { w <- x; names(w) <- NULL; w }
  keep <- list(as.double(seq_len(240L)), as.double(seq.int(241L, 480L)))
  role <- list(V1 = wrap(keep[[1L]]), V2 = wrap(keep[[2L]]))
  attr(role, "row.names") <- .set_row_names(240L)
  attr(role, "class") <- "data.frame"
  before <- serialize(role, NULL)
  caller.before <- serialize(keep, NULL)
  expected <- matrix(0, 240L, 2L)
  for (i in seq_len(480L)) expected[i] <- i
  expect_lt(length(before), length(serialize(expected, NULL)) / 2)
  actual <- getFromNamespace(".npRmpi_bw_literal_rebuild", "npRmpi")(role, "M0")
  expect_identical(actual, expected)
  expect_identical(serialize(role, NULL), before)
  expect_identical(serialize(keep, NULL), caller.before)
})

test_that("an empty matrix role declines literal reconstruction", {
  ns <- asNamespace("npRmpi")
  x <- matrix(numeric(), 3L, 0L)
  expect_null(get(".npRmpi_bw_literal_code", ns)(x, get("toFrame", ns)(x)))
})

test_that("conditional name capture does not call custom name accessors again", {
  table <- get(".__S3MethodsTable__.", asNamespace("base"))
  key <- "names.np_literal_test_frame"
  old <- get0(key, table, inherits = FALSE)
  withr::defer({
    if (is.null(old)) rm(list = key, envir = table) else assign(key, old, table)
  })
  registerS3method("names", "np_literal_test_frame", function(x)
    stop("name accessor should not run during metadata capture"),
    envir = asNamespace("base"))
  owner <- new.env(parent = emptyenv())
  owner$xdat <- structure(data.frame(x = 1:3),
                         class = c("np_literal_test_frame", "data.frame"))
  owner$bw.args <- list(xnames = "x", ynames = "y")
  out <- getFromNamespace(".npRmpi_bw_selector_name_plan", "npRmpi")(
    list(xdat = quote(user.x)), owner)
  expect_null(out)
  owner$bw.args$xnames <- ".np_native_xdat"
  owner$bw.args$ynames <- ".np_native_ydat"
  expr <- quote(c(123456789, 234567890, 345678901, 456789012, 567890123))
  out <- getFromNamespace(".npRmpi_bw_selector_name_plan", "npRmpi")(
    list(xdat = quote(user.x), ydat = expr), owner)
  expect_identical(out, list(xdat = "user.x", ydat = deparse(expr)))
  owner$bw.args$xnames <- c(".np_native_xdat", "second")
  owner$bw.args$ynames <- "actual.y"
  expect_null(getFromNamespace(".npRmpi_bw_selector_name_plan", "npRmpi")(
    list(xdat = quote(user.x), ydat = expr), owner))
})

test_that("literal metadata reconstructs original shapes from invocation data", {
  ns <- asNamespace("npRmpi")
  encode <- get(".npRmpi_bw_literals_encode", ns)
  decode <- get(".npRmpi_bw_literals_decode", ns)
  wire <- get(".npRmpi_bw_metadata_wire", ns)
  fingerprint <- get(".npRmpi_autodispatch_fingerprint", ns)
  to.frame <- get("toFrame", ns)
  untag <- get(".npRmpi_autodispatch_untag", ns)
  x <- sin(seq_len(24L))
  inputs <- list(
    vector = x,
    named = setNames(x, paste0("r", seq_along(x))),
    factor = factor(rep(c("a", "b"), 12L)),
    ordered = ordered(rep(c("a", "b"), 12L)),
    matrix = matrix(x, ncol = 2L),
    named.matrix = matrix(x, ncol = 2L,
      dimnames = list(paste0("r", 1:12), c("a", "b"))),
    frame = data.frame(x = x),
    attributed.frame = structure(data.frame(x = x), annotation = rep("kept", 24L)))
  materialized <- quote(npcdensbw(xdat = .x, ydat = .y))
  for (value in inputs) {
    original <- list(formula = list(producer = "npcdensbw",
      call = as.call(list(quote(npcdensbw), xdat = value, ydat = value))))
    values <- list(.x = untag(to.frame(value)), .y = untag(to.frame(value)))
    encoded <- encode(wire(original), original, materialized, values)
    expect_length(encoded$literal, 2L)
    expect_null(encoded$formula$call$xdat)
    restored <- decode(encoded, materialized, values)
    expect_null(restored$literal)
    expect_identical(fingerprint(original), fingerprint(restored))
    expect_true(identical(value, restored$formula$call$xdat, attrib.as.set = FALSE))
    remote <- unserialize(serialize(values, NULL))
    invisible(lapply(remote, function(d) lapply(d, function(z)
      if (is.numeric(z)) sum(z) else length(as.character(z)))))
    expect_identical(fingerprint(original), fingerprint(decode(
      unserialize(serialize(encoded, NULL)), materialized, remote)))
  }
})

test_that("literal references do not evaluate syntax or certify wrong data", {
  ns <- asNamespace("npRmpi")
  encode <- get(".npRmpi_bw_literals_encode", ns)
  code <- get(".npRmpi_bw_literal_code", ns)
  metadata <- list(selector.native = list(xdat = quote(stop("never evaluate")),
                                        ydat = quote(y)))
  call <- quote(npcdensbw(xdat = .x, ydat = .y))
  values <- list(.x = data.frame(x = c(.1, .2)), .y = data.frame(y = c(.3, .4)))
  expect_identical(encode(metadata, metadata, call, values), metadata)
  expect_null(code(c(.1, .3), values$.x))
  expect_null(code(structure(c(.1, .2), class = "caller_owned"), values$.x))
  expect_null(code(setNames(c(.1, .2), c("same", "same")), values$.x))
})

test_that("malformed private literal references fail before reconstruction", {
  ns <- asNamespace("npRmpi")
  decode <- get(".npRmpi_bw_literals_decode", ns)
  metadata <- list(formula = list(call = quote(npcdensbw(xdat = NULL))),
                   literal = c(f.xdat = "C0"))
  call <- quote(npcdensbw(xdat = .x))
  values <- list(.x = data.frame(x = c(.1, .2)))
  for (refs in list(character(), c(f.xdat = "M9"), c(f.xdat = "F1:9"),
      c(f.xdat = "F1:1,1"), c(f.xdat = "system('never')"),
      c(f.xdat = NA_character_), c(unknown.xdat = "C0"),
      c(f.xdat = "C0", f.xdat = "C0"))) {
    bad <- metadata
    bad$literal <- refs
    expect_error(decode(bad, call, values), "invalid internal literal")
  }
  expect_error(decode(metadata, call, list()), "reference is unavailable")
  bad <- metadata
  bad$formula$call <- NULL
  expect_error(decode(bad, call, values), "invalid internal literal target")
})
