test_that("native roles use only an existing validated lease binding", {
  ns <- asNamespace("npRmpi")
  materialize <- get(".npRmpi_autodispatch_materialize_call", ns)
  compact <- get(".npRmpi_autodispatch_derive_native_roles", ns)
  x <- data.frame(x = c(.1, .3, .7))
  y <- c(1, 3, 2)
  b <- list(.np.native.training = list(xdat = x, ydat = y))
  mc <- quote(npreg(bws = b, txdat = x, tydat = y))
  attr(mc, ".npRmpi.bound.call.owner") <- TRUE
  p <- materialize(mc, environment())
  untag <- get(".npRmpi_autodispatch_untag", ns)
  p$tmpvals <- lapply(p$tmpvals, untag)
  p$prepublish <- lapply(p$prepublish, untag)
  state <- get(".npRmpi_lease_state", ns)
  id <- "validated-model-lease"
  assign(id, b, state$store)
  on.exit(rm(list = id, envir = state$store), add = TRUE)
  source <- as.character(p$call[["bws"]])
  p$tmpvals[source] <- NULL
  p$prepublish[source] <- NULL
  p$lease.bindings <- setNames(list("validated-model-lease"), source)
  p$lease.replacements <- setNames(list(b), source)
  result <- compact(p, mc)
  expect_length(result$derived.bindings, 2L)
  expect_identical(result$lease.bindings, p$lease.bindings)
  expect_identical(result$lease.replacements, p$lease.replacements)
  for (name in names(result$derived.bindings)) {
    descriptor <- result$derived.bindings[[name]]
    expect_identical(descriptor[["source"]], source)
    expect_identical(result$derived.replacements[[name]],
      b[[".np.native.training"]][[descriptor[["component"]]]])
    expect_false(name %in% names(result$tmpvals))
  }
  no.binding <- p
  no.binding$lease.bindings <- NULL
  expect_identical(compact(no.binding, mc), no.binding)
  no.value <- p
  no.value$lease.replacements <- NULL
  expect_identical(compact(no.value, mc), no.value)
  unequal <- p
  response <- as.character(p$call[["tydat"]])
  unequal$tmpvals[[response]] <- y + 1
  expect_length(compact(unequal, mc)$derived.bindings, 1L)
  expect_identical(compact(unequal, mc)$tmpvals[[response]], y + 1)
  classed <- p
  changed <- b
  changed[[".np.native.training"]][["xdat"]] <-
    structure(x, class = c("custom", "data.frame"))
  assign(id, changed, state$store)
  expect_length(compact(classed, mc)$derived.bindings, 1L)
})


test_that("leased roles use canonical metadata and receiving bindings", {
  ns <- asNamespace("npRmpi")
  materialize <- get(".npRmpi_autodispatch_materialize_call", ns)
  compact <- get(".npRmpi_autodispatch_derive_native_roles", ns)
  receive <- get(".npRmpi_spmd_eval_payload", ns)
  untag <- get(".npRmpi_autodispatch_untag", ns)
  sanitize <- get(".npRmpi_autodispatch_sanitize_object", ns)
  state <- get(".npRmpi_lease_state", ns)
  id <- "retained-role-receive-model"
  on.exit({
    if (exists(id, state$store, inherits = FALSE)) rm(list = id, envir = state$store)
    if (exists(id, state$registry, inherits = FALSE)) rm(list = id, envir = state$registry)
  }, add = TRUE)
  x <- data.frame(x = c(.1, .3, .7))
  y <- c(1, 3, 2)
  b <- list(.np.native.training = list(xdat = x, ydat = y))
  canonical <- b
  assign(id, canonical, state$store)
  assign(id, list(state = "live"), state$registry)
  prepare <- function(value = x) {
    mc <- quote(npreg(bws = b, txdat = value, tydat = y))
    attr(mc, ".npRmpi.bound.call.owner") <- TRUE
    p <- materialize(mc, environment())
    p$tmpvals <- lapply(p$tmpvals, untag)
    p$prepublish <- lapply(p$prepublish, untag)
    source <- as.character(p$call[["bws"]])
    p$tmpvals[source] <- NULL
    p$prepublish[source] <- NULL
    p$lease.bindings <- setNames(list(id), source)
    p$lease.replacements <- setNames(list(untag(b)), source)
    list(original = p, reduced = compact(p, mc))
  }
  consumer <- function(bws, txdat, tydat) {
    list(call = match.call(), fitted = txdat$x + tydat,
         rows = .row_names_info(txdat, 0L))
  }
  payload <- function(p) {
    call <- p$call
    call[[1L]] <- consumer
    list(call = call, tmpvals = c(p$tmpvals, p$prepublish),
         tmpnames = p$tmpnames, lease.bindings = p$lease.bindings,
         derived.bindings = p$derived.bindings)
  }
  q <- prepare()
  expect_length(q$reduced$derived.bindings, 2L)
  expected <- receive(payload(q$original), list())$result
  actual <- receive(payload(q$reduced), list())$result
  expect_identical(actual, expected)
  public <- sanitize(actual, q$reduced$lease.replacements)
  expect_false(any(grepl(".__npRmpi", all.names(public$call), fixed = TRUE)))
  object <- list(call = q$original$call)
  expect_identical(sanitize(object, c(q$original$tmpvals, q$original$lease.replacements)),
    sanitize(object, c(q$reduced$tmpvals, q$reduced$lease.replacements,
                       q$reduced$derived.replacements)))
  # A public object's normalized fingerprint cannot certify canonical row flags.
  explicit <- x[seq_len(nrow(x)), , drop = FALSE]
  b[[".np.native.training"]][["xdat"]] <- explicit
  q <- prepare(explicit)
  expect_length(q$reduced$derived.bindings, 1L)
  expected <- receive(payload(q$original), list())$result
  actual <- receive(payload(q$reduced), list())$result
  expect_identical(actual$rows, .row_names_info(explicit, 0L))
  expect_identical(actual, expected)
  b <- canonical
  q <- prepare()
  assign(id, list(state = "prepared"), state$registry)
  expect_error(receive(payload(q$reduced), list()), "is not live")
  assign(id, list(state = "live"), state$registry)
  assign(id, list(.np.native.training = list(ydat = y)), state$store)
  expect_error(receive(payload(q$reduced), list()), "retained training component is unavailable")
  rm(list = id, envir = state$store)
  expect_identical(compact(q$original, structure(quote(npreg(bws = b)),
    .npRmpi.bound.call.owner = TRUE)), q$original)
  expect_error(receive(payload(q$reduced), list()), "is unavailable")
})

test_that("live leases preserve explicit row metadata and later reuse", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  resolve <- get(".npRmpi_lease_resolve", ns)
  x <- data.frame(x = sort(sin(seq_len(20L))))
  y <- sort(cos(seq_len(20L)))
  b <- npregbw(xdat = x, ydat = y, bws = .4, bandwidth.compute = FALSE)
  lease <- resolve(b, "npreg", "bws")
  expect_false(is.null(lease))
  explicit <- x[seq_len(nrow(x)), , drop = FALSE]
  actual <- npreg(b, txdat = explicit, tydat = y, se = FALSE)
  expect_identical(.row_names_info(actual$bws[[".np.native.training"]][["xdat"]], 0L),
                   .row_names_info(explicit, 0L))
  expect_identical(resolve(b, "npreg", "bws"), lease)
  copy <- unserialize(serialize(x, NULL))
  copied <- npreg(b, txdat = copy, tydat = y, se = FALSE)
  expect_identical(fitted(copied), fitted(actual))
  expect_identical(resolve(b, "npreg", "bws"), lease)
  copy$x[10L] <- copy$x[10L] + .01
  changed <- npreg(b, txdat = copy, tydat = y, se = FALSE)
  expect_identical(resolve(b, "npreg", "bws"), lease)
  expect_false(any(grepl(".__npRmpi", all.names(changed$call), fixed = TRUE)))
  expect_identical(fitted(eval(actual$call)), fitted(actual))
})
