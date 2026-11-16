test_that("formula role descriptors preserve values and ordinary headers", {
  ns <- asNamespace("npRmpi")
  owner <- environment(get(".npRmpi_autodispatch_derive_native_roles", ns))
  describe <- get(".npRmpi_formula_role_descriptor", owner)
  extract <- get(".npRmpi_formula_role_extract", owner)
  d <- data.frame(x = c(.1, .3, .7), y = c(3, 1, 2),
                  f = factor(c("a", "b", "a")), check.names = FALSE)
  form <- y ~ x + f
  environment(form) <- baseenv()
  frame <- model.frame(form, d)
  before <- serialize(frame, NULL)
  targets <- list(d["x"], d[c("f", "x")], d$y,
                  d[c(1L, 2L, 3L), "x", drop = FALSE])
  for (target in targets) {
    descriptor <- describe(frame, target, "y")
    expect_false(is.null(descriptor))
    expect_identical(extract(frame, descriptor), target)
    if (is.data.frame(target))
      expect_identical(.row_names_info(extract(frame, descriptor), 0L),
                       .row_names_info(target, 0L))
    expect_identical(extract(unserialize(serialize(frame, NULL)), descriptor), target)
  }
  expect_identical(serialize(frame, NULL), before)
  expect_null(describe(frame, d$y))
  expect_null(describe(frame, d[3:1, "x", drop = FALSE]))
  expect_null(describe(frame, transform(d["x"], x = x + 1)))
  expect_null(describe(frame, structure(d["x"], class = c("custom", "data.frame"))))
  expect_null(describe(frame, structure(d["x"], terms = terms(y ~ x))))
  expect_null(describe(frame, data.frame(unknown = d$x)))
  dup <- frame
  names(dup) <- rep("x", length(dup))
  expect_null(describe(dup, d["x"]))
  expect_null(describe(frame, as.matrix(d["x"])))
  matrix.frame <- data.frame(y = d$y)
  matrix.frame$x <- I(matrix(1:6, ncol = 2))
  expect_null(describe(matrix.frame, d["x"]))
  expect_error(extract(frame, list(kind = "frame", columns = 0L)), "columns are unavailable")
  expect_error(extract(frame, list(kind = "unknown", columns = 1L)), "invalid.*descriptor")
})

test_that("formula training sharing requires the already admitted canonical lease", {
  ns <- asNamespace("npRmpi")
  materialize <- get(".npRmpi_autodispatch_materialize_call", ns)
  compact <- get(".npRmpi_autodispatch_derive_native_roles", ns)
  receive <- get(".npRmpi_spmd_eval_payload", ns)
  untag <- get(".npRmpi_autodispatch_untag", ns)
  state <- get(".npRmpi_lease_state", ns)
  id <- "formula-role-test-only"
  on.exit({
    for (env in list(state$store, state$registry))
      if (exists(id, env, inherits = FALSE)) rm(list = id, envir = env)
  }, add = TRUE)
  d <- data.frame(y = c(1, 3, 2), x = c(.1, .3, .7))
  form <- y ~ x
  environment(form) <- baseenv()
  frame <- model.frame(form, d)
  b <- list(ynames = "y", .np.formula.training = list(frame = frame))
  x <- d["x"]
  y <- d$y
  mc <- quote(npreg(bws = b, txdat = x, tydat = y))
  attr(mc, ".npRmpi.bound.call.owner") <- TRUE
  p <- materialize(mc, environment())
  p$tmpvals <- lapply(p$tmpvals, untag)
  p$prepublish <- lapply(p$prepublish, untag)
  expect_identical(compact(p, mc), p) # An inline formula frame is not authority.
  source <- as.character(p$call[["bws"]])
  p$tmpvals[source] <- NULL
  p$prepublish[source] <- NULL
  p$lease.bindings <- setNames(list(id), source)
  p$lease.replacements <- setNames(list(b), source)
  assign(id, b, state$store)
  assign(id, list(state = "live"), state$registry)
  before <- serialize(get(id, state$store), NULL)
  q <- compact(p, mc)
  expect_length(q$derived.bindings, 2L)
  expect_identical(q$lease.bindings, p$lease.bindings)
  expect_identical(q$lease.replacements, p$lease.replacements)
  expect_identical(serialize(get(id, state$store), NULL), before)
  consumer <- function(bws, txdat, tydat) list(x = txdat, y = tydat, fitted = txdat$x + tydat)
  payload <- function(z) {
    call <- z$call
    call[[1L]] <- consumer
    list(call = call, tmpvals = c(z$tmpvals, z$prepublish), tmpnames = z$tmpnames,
         lease.bindings = z$lease.bindings, derived.bindings = z$derived.bindings)
  }
  expect_identical(receive(payload(q), list())$result, receive(payload(p), list())$result)
  for (field in c("lease.bindings", "lease.replacements")) {
    missing <- p
    missing[[field]] <- NULL
    expect_identical(compact(missing, mc), missing)
  }
  changed <- p
  response <- as.character(p$call[["tydat"]])
  changed$tmpvals[[response]] <- y + 1
  expect_length(compact(changed, mc)$derived.bindings, 1L)
  expect_identical(compact(changed, mc)$tmpvals[[response]], y + 1)
  changed <- p
  names(changed$call)[names(changed$call) == "txdat"] <- "exdat"
  expect_length(compact(changed, mc)$derived.bindings, 1L)
  for (generic in c("npsigtest", "npcopula", "nplsqreg")) {
    excluded <- mc
    excluded[[1L]] <- as.name(generic)
    expect_identical(compact(p, excluded), p)
  }
  assign(id, list(state = "prepared"), state$registry)
  expect_error(receive(payload(q), list()), "is not live")
  assign(id, list(state = "live"), state$registry)
  changed <- b
  changed[[".np.formula.training"]][["frame"]] <- NULL
  assign(id, changed, state$store)
  expect_identical(compact(p, mc), p)
  expect_error(receive(payload(q), list()), "columns are unavailable")
})

test_that("public formula fits reuse certified roles and restored objects stay ordinary", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  old <- options(l1.test.formula.derived = 0L)
  on.exit(options(old), add = TRUE)
  trace(".npRmpi_autodispatch_derive_native_roles", exit = quote({
    options(l1.test.formula.derived = getOption("l1.test.formula.derived", 0L) +
      length(Filter(is.list, returnValue()$derived.bindings)))
  }), where = ns, print = FALSE)
  on.exit(untrace(".npRmpi_autodispatch_derive_native_roles", where = ns), add = TRUE)
  i <- seq_len(24L); d <- data.frame(y = cos(i*.3), x = sin(i*.7))
  b <- npregbw(y ~ x, data = d, bws = .5, bandwidth.compute = FALSE)
  resolve <- get(".npRmpi_lease_resolve", ns)
  lease <- resolve(b, "npreg", "bws")
  expect_false(is.null(lease))
  fit <- npreg(b, se = FALSE)
  expect_identical(getOption("l1.test.formula.derived"), 2L)
  expect_identical(resolve(b, "npreg", "bws"), lease)
  options(l1.test.formula.derived = 0L)
  restored <- unserialize(serialize(b, NULL))
  expect_null(resolve(restored, "npreg", "bws"))
  other <- npreg(restored, se = FALSE)
  expect_identical(getOption("l1.test.formula.derived"), 0L)
  expect_identical(fitted(other), fitted(fit))
  options(l1.test.formula.derived = 0L)
  changed <- npreg(b, txdat = transform(d["x"], x = x + .1), tydat = d$y, se = FALSE)
  expect_identical(getOption("l1.test.formula.derived"), 1L)
  expect_identical(resolve(b, "npreg", "bws"), lease)
  expect_false(any(grepl(".__npRmpi", all.names(changed$call), fixed = TRUE)))
})
