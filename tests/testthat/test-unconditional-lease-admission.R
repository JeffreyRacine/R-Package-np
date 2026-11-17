test_that("density and distribution publication retains unique existing capabilities", {
  tab <- getFromNamespace(".npRmpi_lease_capability_table", "npRmpi")()
  expect_false(anyDuplicated(tab$producer) > 0L)
  expect_false(anyDuplicated(tab$capability) > 0L)
  expect_identical(tab$producer[1:8], c("npregbw", "nplsqregbw", "npplregbw",
    "npindexbw", "npscoefbw", ".npscoefbw_nomad_context_prepare", "npcdensbw", "npcdistbw"))
  expect_identical(tab$kind[1:8], c(rep("leased_public", 5L), "scoped_internal", rep("leased_public",2L)))
  expect_identical(tab$argument[1:8], c(rep("bws", 5L), "ctx", rep("bws",2L)))
  expect_identical(tab$fingerprint_policy[1:8], c(rep("bandwidth.v1", 5L), "context.v1", rep("bandwidth.v1",2L)))
  expect_identical(tab$capability[1:8], c("bw.npreg", "bw.nplsqreg", "bw.npplreg",
    "bw.npindex", "bw.npscoef", "internal.npscoef.nomad_context", "bw.npcdens", "bw.npcdist"))
  expect_identical(tab$consumer[1:8], c("npreg", "nplsqreg", "npplreg", "npindex",
    "npscoef", ".npscoefbw_nomad_pool_start", "npcdens", "npcdist"))
  expect_identical(tab$expected_class[1:8], c("rbandwidth", "lsqregressionbandwidth",
    "plbandwidth", "sibandwidth", "scbandwidth", "list", "conbandwidth", "condbandwidth"))
  families <- c(npudens = "bandwidth", npudist = "dbandwidth",
                npcdens = "conbandwidth", npcdist = "condbandwidth")
  validate <- getFromNamespace(".npRmpi_lease_validate_plan", "npRmpi")
  for (family in names(families)) {
    row <- tab[tab$producer == paste0(family, "bw"), , drop = FALSE]
    expect_equal(nrow(row), 1L)
    expect_identical(row$consumer, family)
    expect_identical(row$expected_class, unname(families[[family]]))
    expect_identical(row$argument, "bws")
    expect_identical(row$fingerprint_policy, "bandwidth.v1")
    plan <- as.list(row)
    plan$session <- "test-session"; plan$id <- "test-id"; plan$remote_name <- "test-id"
    expect_true(validate(plan))
    plan$expected_class <- "rbandwidth"
    expect_error(validate(plan), "not capability-authorized", fixed = TRUE)
  }
  expect_false(any(tab$consumer %in% c("npqreg", "npconmode", "npcopula")))
})

test_that("prepublication vector names use the full existing schema contract", {
  apply <- getFromNamespace(".npRmpi_native_names_apply", "npRmpi")
  signature <- getFromNamespace(".np_bw_procedure_signature", "npRmpi")
  old <- ".__npRmpi_autod_dat_1"
  for (cls in c("bandwidth", "dbandwidth")) {
    b <- structure(list(bw = .4, xnames = old, varnames = list(x = old),
      icon = setNames(TRUE, old), iuno = setNames(FALSE, old),
      iord = setNames(FALSE, old), type = "fixed",
      xdati = list(all.lev = setNames(list(old), old)),
      .np.native.training = list(dat = setNames(data.frame(x = 1:3), old))), class = cls)
    b$selection.procedure <- list(signature = signature(b),
      arguments = list(unrelated.literal = old))
    before <- serialize(b, NULL)
    out <- apply(b, list(xnames = "d$x"))
    expect_identical(b, unserialize(before))
    expect_identical(out$xnames, "d$x")
    expect_identical(names(out$.np.native.training$dat), "d$x")
    expect_identical(names(out$xdati$all.lev), "d$x")
    expect_identical(out$xdati$all.lev[[1L]], old)
    expect_identical(out$selection.procedure$arguments, b$selection.procedure$arguments)
    expect_identical(out$selection.procedure$signature, signature(out))
    expect_identical(apply(out, list(xnames = "other")), out)
    stale <- b; stale$selection.procedure$signature$bw <- .5
    expect_identical(apply(stale, list(xnames = "d$x")), stale)
    expect_identical(apply(b, list()), b)
  }
})


test_that("supplied unconditional bandwidths retain their own vector labels", {
  plan <- getFromNamespace(".npRmpi_native_name_plan", "npRmpi")
  capture <- function(dat, bws = NULL) {
    plan(sys.function(), environment(), parent.frame(), list(dat = dat, bws = bws))
  }
  d <- data.frame(x = 1:3)
  expect_identical(capture(d$x), list(xnames = "d$x"))
  for (cls in c("bandwidth", "dbandwidth")) {
    b <- structure(list(xnames = ".__npRmpi_autod_dat_1"), class = cls)
    expect_identical(capture(d$x, b), list())
  }
})

test_that("new family formula finishes preserve order and caller ownership", {
  apply <- getFromNamespace(".npRmpi_bw_formula_apply", "npRmpi")
  project <- getFromNamespace(".npRmpi_bw_metadata_wire", "npRmpi")
  restore <- getFromNamespace(".npRmpi_bw_restore_public_owner", "npRmpi")
  retain <- getFromNamespace(".np_bws_retain_formula_training", "npRmpi")
  classes <- c(npudensbw = "bandwidth", npudistbw = "dbandwidth",
    npcdensbw = "conbandwidth", npcdistbw = "condbandwidth")
  d <- data.frame(x = 1:4, y = c(1, NA, 3, 4))
  owner <- new.env(parent = baseenv())
  for (producer in names(classes)) {
    uncond <- producer %in% c("npudensbw", "npudistbw")
    f <- if (uncond) ~y else y ~ x
    environment(f) <- owner
    frame <- stats::model.frame(f, d, na.action = stats::na.omit)
    cl <- quote(bandwidth_call()); environment(cl) <- owner
    plan <- list(producer = producer, call = cl, formula = f, frame = frame,
      terms = attr(frame, "terms"), na.action = stats::na.omit,
      variableNames = list(x = "x", y = "y"))
    b <- structure(list(bw = .4), class = unname(classes[[producer]]))
    expected <- b; expected$call <- cl; expected$formula <- f
    if (uncond) {
      expected$terms <- plan$terms
      expected <- retain(expected, frame, plan$na.action)
      expected$rows.omit <- as.vector(attr(frame, "na.action"))
      expected$nobs.omit <- length(expected$rows.omit)
    } else {
      expected$rows.omit <- as.vector(attr(frame, "na.action"))
      expected$nobs.omit <- length(expected$rows.omit)
      expected$terms <- plan$terms
      expected <- retain(expected, frame, plan$na.action)
      expected$variableNames <- plan$variableNames
    }
    out <- apply(b, plan)
    expect_identical(out, expected)
    expect_identical(names(out), names(expected))
    expect_identical(apply(out, plan), out)
    stored <- apply(b, project(list(formula = plan))$formula)
    restored <- restore(list(bws = stored), list(bws = out))
    expect_identical(restored$bws, out)
    expect_identical(environment(restored$bws$formula), owner)
  }
})


test_that("prepared-frame references preserve exact retained columns and header", {
  skip_if(getRversion() < "4.2.0", "formula reuse requires R 4.2.0 condition classes")
  owner <- environment(getFromNamespace(".npRmpi_autodispatch_call", "npRmpi"))
  admit <- get(".npRmpi_bw_frame_admit", owner, inherits = TRUE)
  decode <- get(".npRmpi_bw_frame_decode", owner, inherits = TRUE)
  project <- getFromNamespace(".npRmpi_bw_metadata_wire", "npRmpi")
  d <- data.frame(x = sin(seq_len(12L)), g = factor(rep(c("a", "b"), 6L)))
  for (producer in c("npudensbw", "npudistbw")) {
    formula <- stats::as.formula("~x + g", env = globalenv())
    frame <- stats::model.frame(formula, data = d)
    dat <- frame[, c("x", "g"), drop = FALSE]
    plan <- list(producer = producer, call = quote(npudensbw(~x + g, data = d)),
      formula = formula, frame = frame, na.action = stats::na.omit,
      terms = attr(frame, "terms"), xterms = NULL, chromoly = NULL,
      ynames = character(), omit = NULL, scalar.children = FALSE)
    metadata <- list(names = list(), formula = plan)
    call <- quote(npudensbw(dat = prepared_dat))
    values <- list(prepared_dat = dat)
    before <- serialize(metadata, NULL)
    # The dispatcher binds an empty native-identity vector for dat-only calls.
    # Its presence must not exclude an otherwise representable formula plan.
    metadata$identity <- character()
    before <- serialize(metadata, NULL)
    encoded <- admit(metadata, call, values)
    expect_identical(encoded$status, "admit")
    wire <- unserialize(serialize(list(m = encoded$metadata, v = values), NULL))
    restored <- decode(wire$m, call, wire$v)
    expected <- project(metadata)
    attr(expected$formula$frame, "row.names") <- .row_names_info(frame, 0L)
    expect_identical(restored, expected)
    expect_identical(.row_names_info(restored$formula$frame, 0L), .row_names_info(frame, 0L))
    expect_identical(serialize(metadata, NULL), before)
    expect_null(restored[["frame.reference", exact = TRUE]])
    expect_identical(admit(list(names = list(), formula = NULL), call, values)$status, "none")

    invalid <- encoded$metadata
    invalid$frame.reference$role <- "xdat"
    expect_error(decode(invalid, call, values), "invalid prepared-frame reference", fixed = TRUE)
    invalid <- encoded$metadata
    invalid$frame.reference$extra <- seq_len(12L)
    expect_error(decode(invalid, call, values), "invalid prepared-frame reference", fixed = TRUE)
    invalid <- encoded$metadata
    invalid$frame.reference$order <- setdiff(invalid$frame.reference$order, "terms")
    expect_error(decode(invalid, call, values), "invalid prepared-frame reference", fixed = TRUE)

    # A same-shaped replay is not a shape error. The existing publication
    # certificate must distinguish its changed retained training content.
    changed <- wire$v
    changed$prepared_dat$x[1L] <- changed$prepared_dat$x[1L] + 0.125
    replay <- decode(wire$m, call, changed)
    fingerprint <- getFromNamespace(".npRmpi_autodispatch_fingerprint", "npRmpi")
    expect_false(identical(fingerprint(restored$formula$frame), fingerprint(replay$formula$frame)))
  }
})

test_that("unbounded formula metadata declines only the new optimization", {
  owner <- environment(getFromNamespace(".npRmpi_autodispatch_call", "npRmpi"))
  admit <- get(".npRmpi_bw_frame_admit", owner, inherits = TRUE)
  d <- data.frame(x = seq_len(12L) / 13)
  formula <- stats::as.formula("~x", env = globalenv())
  frame <- stats::model.frame(formula, data = d)
  dat <- frame[, "x", drop = FALSE]
  plan <- list(producer = "npudensbw", call = quote(npudensbw(~x, data = d)),
    formula = formula, frame = frame, na.action = stats::na.omit,
    terms = attr(frame, "terms"), xterms = NULL, chromoly = NULL,
    ynames = character(), omit = NULL, scalar.children = FALSE)
  metadata <- list(names = list(), formula = plan)
  call <- quote(npudensbw(dat = prepared_dat))
  values <- list(prepared_dat = dat)
  mutations <- list(
    nonempty_identity = function(m) { m$identity <- c(prepared_dat = "x"); m },
    attributed_identity = function(m) {
      m$identity <- structure(character(), payload = 1L); m
    },
    literal = function(m) { m$formula$call[["data"]] <- d; m },
    local_policy = function(m) {
      policy <- function(x) stats::na.omit(x)
      environment(policy) <- baseenv()
      m$formula$na.action <- policy; m
    },
    call_attribute = function(m) { attr(m$formula$call, "payload") <- d; m },
    frame_attribute = function(m) { attr(m$formula$frame, "payload") <- d; m },
    omission = function(m) {
      omit <- structure(13L, class = "omit")
      m$formula$omit <- omit; attr(m$formula$frame, "na.action") <- omit; m
    })
  for (mutate in mutations) {
    original <- mutate(metadata); saved <- serialize(original, NULL)
    result <- admit(original, call, values)
    expect_identical(result$status, "decline")
    expect_null(result$metadata)
    expect_identical(serialize(original, NULL), saved)
  }
})


test_that("prepared-frame syntax admission is bounded without evaluating owners", {
  skip_if(getRversion() < "4.2.0", "formula reuse requires R 4.2.0 condition classes")
  owner <- environment(getFromNamespace(".npRmpi_autodispatch_call", "npRmpi"))
  admit <- get(".npRmpi_bw_frame_admit", owner, inherits = TRUE)
  preflight <- get(".npRmpi_bw_frame_preflight", owner, inherits = TRUE)
  d <- data.frame(x = sin(seq_len(12L)))
  formula <- stats::as.formula("~x", env = globalenv())
  frame <- stats::model.frame(formula, d)
  plan <- list(producer = "npudensbw", call = quote(npudensbw(~x, data = d)),
    formula = formula, frame = frame, na.action = stats::na.omit,
    terms = attr(frame, "terms"), xterms = NULL, chromoly = NULL,
    ynames = character(), omit = NULL, scalar.children = FALSE)
  metadata <- list(names = list(), formula = plan)
  call <- quote(npudensbw(dat = prepared_dat))
  values <- list(prepared_dat = frame[, "x", drop = FALSE])
  for (producer in c("npudensbw", "npudistbw")) {
    metadata$formula$producer <- producer
    missing.string <- metadata
    missing.string$formula$call["call"] <- list(NA_character_)
    expect_identical(admit(missing.string, call, values)$status, "admit")

    local.owner <- new.env(parent = baseenv())
    local.owner$unused <- rep(3.14159, 5000L)
    nested <- metadata
    f <- formula; environment(f) <- local.owner
    nested$formula$call[[2L]] <- f
    expect_identical(admit(nested, call, values)$status, "admit")
    long.tag <- metadata
    long.tag$formula$call[[strrep("t", 9000L)]] <- 1
    expect_identical(admit(long.tag, call, values)$status, "decline")
    wide <- metadata
    wide$formula$call["call"] <- list(structure(
      as.call(c(list(quote(`~`)), rep(list(quote(x)), 5000L))), class = "formula"))
    expect_identical(admit(wide, call, values)$status, "decline")
    deep <- metadata; expr <- quote(x)
    for (i in seq_len(5000L)) expr <- call("-", expr)
    deep$formula$call["call"] <- list(expr)
    expect_identical(admit(deep, call, values)$status, "decline")

    # Deferred character scalars can serialize larger than their printed value.
    # The work estimate is not the exact runtime-format byte certificate.
    deferred <- metadata
    deferred$formula$call["call"] <- list(as.call(c(list(quote(c)),
      lapply(seq_len(100L), function(i) as.character(i + .12345)))))
    budget <- new.env(parent = emptyenv()); budget$left <- 6144
    expect_true(preflight(deferred$formula, frame, budget))
    result <- admit(deferred, call, values)
    projected <- get(".npRmpi_bw_metadata_wire", owner, inherits = TRUE)(deferred,
      plain.language = TRUE)
    syntax <- list(call = projected$formula$call, formula = projected$formula$formula,
                   subset = NULL)
    # Some runtimes materialize these scalars. Both decisions follow the actual
    # bytes; this host's frozen private qualification also proves a real decline.
    expect_identical(result$status, if (length(serialize(syntax, NULL)) > 6144)
      "decline" else "admit")
    if (identical(result$status, "decline"))
      expect_identical(result$reason, "projected call syntax exceeds the lease's metadata bound")
  }
})


test_that("older coordinators decline new formula reuse before syntax inspection", {
  owner <- environment(getFromNamespace(".npRmpi_autodispatch_call", "npRmpi"))
  admit <- get(".npRmpi_bw_frame_admit", owner, inherits = TRUE)
  native <- list(names = list(), formula = NULL)
  for (producer in c("npudensbw", "npudistbw")) {
    metadata <- list(formula = list(producer = producer, frame = data.frame(x = 1:3)))
    for (version in c("3.5.0", "4.1.3", "4.2.0")) {
      scope <- new.env(parent = environment(admit))
      scope$getRversion <- local({ v <- version; function() numeric_version(v) })
      scope$.npRmpi_bw_frame_preflight <- function(...) stop("inspection sentinel")
      tested <- admit; environment(tested) <- scope
      if (version != "4.2.0") {
        out <- tested(metadata, quote(npudensbw(dat = d)), list(d = data.frame(x = 1:3)))
        expect_identical(out$status, "decline")
        expect_identical(out$reason, "prepared-frame reuse requires R 4.2.0 condition classes")
      } else {
        expect_error(tested(metadata, quote(npudensbw(dat = d)), list(d = data.frame(x = 1:3))),
          "inspection sentinel", fixed = TRUE)
      }
      expect_identical(tested(native, quote(npudensbw(dat = d)), list(d = data.frame(x = 1:3))),
        admit(native, quote(npudensbw(dat = d)), list(d = data.frame(x = 1:3))))
    }
  }
})

test_that("subset formula objects decline admission and invalid references fail", {
  skip_if(getRversion() < "4.2.0", "formula reuse requires R 4.2.0 condition classes")
  owner <- environment(getFromNamespace(".npRmpi_autodispatch_call", "npRmpi"))
  admit <- get(".npRmpi_bw_frame_admit", owner, inherits = TRUE)
  decode <- get(".npRmpi_bw_frame_decode", owner, inherits = TRUE)
  d <- data.frame(x = sin(seq_len(12L)))
  formula <- stats::as.formula("~x", env = globalenv())
  frame <- stats::model.frame(formula, d)
  plan <- list(producer = "npudensbw", call = quote(npudensbw(~x, data = d)),
    formula = formula, frame = frame, na.action = stats::na.omit,
    terms = attr(frame, "terms"), xterms = NULL, chromoly = NULL,
    ynames = character(), omit = NULL, scalar.children = FALSE)
  attr(plan$frame, ".np.subset") <- list(expression = quote(x < .8), envir = globalenv())
  metadata <- list(names = list(), formula = plan)
  call <- quote(npudensbw(dat = prepared_dat))
  values <- list(prepared_dat = frame[, "x", drop = FALSE])
  for (producer in c("npudensbw", "npudistbw")) {
    metadata$formula$producer <- producer
    valid <- admit(metadata, call, values)
    expect_identical(valid$status, "admit")
    for (expression in list(formula, as.call(list(quote(inherits), formula, "formula")))) {
      m <- metadata; attr(m$formula$frame, ".np.subset")$expression <- expression
      expect_identical(admit(m, call, values)$status, "decline")
      ref <- valid$metadata; ref$frame.reference$subset$expression <- expression
      expect_error(decode(ref, call, values), "invalid prepared-frame reference", fixed = TRUE)
    }
  }
})
