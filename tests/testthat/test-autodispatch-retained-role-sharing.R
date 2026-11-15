test_that("same-call native roles preserve exact values and explicit overrides", {
  scope <- environment(getFromNamespace(".npRmpi_distributed_call_impl", "npRmpi"))
  compact <- get(".npRmpi_autodispatch_derive_native_roles", scope)
  plain <- get(".npRmpi_transport_plain_value", scope)
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")
  untag <- getFromNamespace(".npRmpi_autodispatch_untag", "npRmpi")
  x <- data.frame(x = c(-0, 1, NA_real_, NaN), g = factor(c("a","b","a","b")))
  y <- c(1, 2, 3, 4)
  b <- list(.np.native.training = list(xdat=x, ydat=y))
  prep <- function(value=x, response=y, owner=TRUE, family="npreg") {
    mc <- as.call(list(as.name(family), bws=quote(b), txdat=quote(value), tydat=quote(response)))
    attr(mc,".npRmpi.bound.call.owner") <- owner
    p <- materialize(mc,environment())
    p$tmpvals <- lapply(p$tmpvals,untag);p$prepublish<-lapply(p$prepublish,untag)
    list(p=p,mc=mc)
  }
  q <- prep(); ans <- compact(q$p,q$mc)
  expect_length(ans$derived.bindings,2L)
  for(nm in names(ans$derived.bindings)) {
    desc <- ans$derived.bindings[[nm]]
    values <- c(ans$tmpvals,ans$prepublish)
    recovered <- values[[desc[["source"]]]][[".np.native.training"]][[desc[["component"]]]]
    expect_identical(recovered,ans$derived.replacements[[nm]])
    expect_false(nm %in% names(values))
    expect_true(nm %in% ans$tmpnames)
  }
  override <- y+1
  q <- prep(response=override); expect_length(compact(q$p,q$mc)$derived.bindings,1L)
  x2 <- x; x2$x[1] <- 0
  q <- prep(value=x2); expect_length(compact(q$p,q$mc)$derived.bindings,1L)
  q <- prep(owner=FALSE);expect_identical(compact(q$p,q$mc),q$p)
  q <- prep(family="npksum");expect_identical(compact(q$p,q$mc),q$p)
  # Same row-name values do not imply the same automatic-row-name state.
  rows <- x[seq_len(nrow(x)), , drop = FALSE]
  q <- prep(value = rows)
  expect_length(compact(q$p, q$mc)$derived.bindings, 1L)
  expect_identical(.row_names_info(q$p$tmpvals[[as.character(q$p$call$txdat)]], 0L),
                   .row_names_info(rows, 0L))
  expect_true(plain(x));expect_true(plain(matrix(1:6,2)))
  expect_true(plain(ordered(c("a","b"))));expect_false(plain(as.Date("2020-01-01")))
  expect_false(plain(structure(1:4, hidden=new.env())))
  expect_false(plain(structure(x,class=c("custom","data.frame"))))
  expect_false(plain(list(x=1:4)))

  q <- prep(value = unserialize(serialize(x, NULL)))
  expect_length(compact(q$p, q$mc)$derived.bindings, 2L)
  q <- prep(response = NULL)
  expect_length(compact(q$p, q$mc)$derived.bindings, 1L)
  q <- prep()
  source <- as.character(q$p$call[["bws"]])
  q$p$tmpvals[source] <- NULL
  q$p$prepublish[source] <- NULL
  q$p$lease.bindings[source] <- list("existing-lease")
  expect_identical(compact(q$p, q$mc), q$p)
  b <- list()
  q <- prep()
  expect_identical(compact(q$p, q$mc), q$p)
  x$day <- as.Date("2020-01-01") + seq_len(nrow(x))
  b <- list(.np.native.training = list(xdat = x, ydat = y))
  q <- prep()
  expect_length(compact(q$p, q$mc)$derived.bindings, 1L)
})

test_that("receiving ranks reconstruct roles and sanitize their returned calls", {
  ns <- asNamespace("npRmpi")
  materialize <- get(".npRmpi_autodispatch_materialize_call", ns)
  compact <- get(".npRmpi_autodispatch_derive_native_roles", ns)
  receive <- get(".npRmpi_spmd_eval_payload", ns)
  sanitize <- get(".npRmpi_autodispatch_sanitize_object", ns)
  x <- data.frame(x = c(1, 3, 5))
  y <- c(2, 4, 6)
  b <- list(.np.native.training = list(xdat = x, ydat = y))
  mc <- quote(npreg(bws = b, txdat = x, tydat = y))
  attr(mc, ".npRmpi.bound.call.owner") <- TRUE
  original <- materialize(mc, environment())
  reduced <- compact(original, mc)
  expect_length(reduced$derived.bindings, 2L)
  consumer <- function(bws, txdat, tydat) {
    list(call = match.call(), fitted = txdat$x + tydat)
  }
  payload <- function(prepared) {
    call <- prepared$call
    call[[1L]] <- consumer
    list(call = call, tmpvals = c(prepared$tmpvals, prepared$prepublish),
         tmpnames = prepared$tmpnames,
         derived.bindings = prepared$derived.bindings)
  }
  expected <- receive(payload(original), list())$result
  actual <- receive(payload(reduced), list())$result
  expect_identical(actual, expected)
  expect_identical(actual$fitted, c(3, 7, 11))
  # A missing receiving step fails before evaluation; missing sanitation
  # leaves temporary argument names in the returned call.
  expect_identical(actual$call[["txdat"]], expected$call[["txdat"]])
  expect_false(any(grepl(".__npRmpi", all.names(actual$call), fixed = TRUE)))
  object <- list(call = original$call)
  expect_identical(
    sanitize(object, c(original$tmpvals, original$prepublish)),
    sanitize(object, c(reduced$tmpvals, reduced$prepublish,
                       reduced$derived.replacements)))
})

test_that("automatic native dispatch sends and restores retained roles", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 20L))
  y <- x$x^2
  b <- npregbw(xdat = x, ydat = y, bws = .4,
               bandwidth.compute = FALSE)
  b <- unserialize(serialize(b, NULL)) # Restored objects have no valid lease.
  original <- npreg(b, txdat = x, tydat = y, se = FALSE)
  derived <- npreg(b, se = FALSE)
  expect_identical(fitted(derived), fitted(original))
  expect_false(any(grepl(".__npRmpi", all.names(derived$call), fixed = TRUE)))
  replay <- eval(derived$call)
  expect_identical(fitted(replay), fitted(derived))
})


test_that("transport certificates preserve text encoding in data and attributes", {
  scope <- environment(getFromNamespace(".npRmpi_autodispatch_derive_native_roles", "npRmpi"))
  same <- get(".npRmpi_transport_same_value", scope)
  utf8 <- "caf\u00e9"
  latin <- iconv(utf8, from = "UTF-8", to = "latin1")
  expect_false(is.na(latin))
  Encoding(utf8) <- "UTF-8"
  Encoding(latin) <- "latin1"
  expect_identical(utf8, latin)
  expect_false(same(utf8, latin))
  expect_true(same(utf8, utf8))
  a <- structure(c(1L, 1L), levels = utf8, class = "factor")
  b <- structure(c(1L, 1L), levels = latin, class = "factor")
  expect_identical(a, b)
  expect_false(same(a, b))
  expect_false(same(data.frame(g = a), data.frame(g = b)))
  expect_true(same(data.frame(g = a), data.frame(g = a)))
  expect_false(same(setNames(1:2, c(utf8, "x")),
                    setNames(1:2, c(latin, "x"))))
  expect_false(same(matrix(1:2, dimnames = list(c(utf8, "x"), NULL)),
                    matrix(1:2, dimnames = list(c(latin, "x"), NULL))))
})
