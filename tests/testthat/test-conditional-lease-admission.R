test_that("density and distribution publication retains unique existing capabilities", {
  tab <- getFromNamespace(".npRmpi_lease_capability_table", "npRmpi")()
  expect_false(anyDuplicated(tab$producer) > 0L)
  expect_false(anyDuplicated(tab$capability) > 0L)
  expect_identical(tab$producer[1:6], c("npregbw", "nplsqregbw", "npplregbw",
    "npindexbw", "npscoefbw", ".npscoefbw_nomad_context_prepare"))
  expect_identical(tab$kind[1:6], c(rep("leased_public", 5L), "scoped_internal"))
  expect_identical(tab$argument[1:6], c(rep("bws", 5L), "ctx"))
  expect_identical(tab$fingerprint_policy[1:6], c(rep("bandwidth.v1", 5L), "context.v1"))
  expect_identical(tab$capability[1:6], c("bw.npreg", "bw.nplsqreg", "bw.npplreg",
    "bw.npindex", "bw.npscoef", "internal.npscoef.nomad_context"))
  expect_identical(tab$consumer[1:6], c("npreg", "nplsqreg", "npplreg", "npindex",
    "npscoef", ".npscoefbw_nomad_pool_start"))
  expect_identical(tab$expected_class[1:6], c("rbandwidth", "lsqregressionbandwidth",
    "plbandwidth", "sibandwidth", "scbandwidth", "list"))
  families <- c(npcdens = "conbandwidth", npcdist = "condbandwidth")
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

test_that("new family formula finishes preserve order and caller ownership", {
  apply <- getFromNamespace(".npRmpi_bw_formula_apply", "npRmpi")
  project <- getFromNamespace(".npRmpi_bw_metadata_wire", "npRmpi")
  restore <- getFromNamespace(".npRmpi_bw_restore_public_owner", "npRmpi")
  retain <- getFromNamespace(".np_bws_retain_formula_training", "npRmpi")
  classes <- c(npcdensbw = "conbandwidth", npcdistbw = "condbandwidth")
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


test_that("only completed conditional default-method calls gain publication", {
  plan <- getFromNamespace(".npRmpi_lease_publication_plan", "npRmpi")
  producer_for <- getFromNamespace(".npRmpi_lease_producer", "npRmpi")
  for (producer in c("npcdensbw", "npcdistbw")) {
    complete <- as.call(list(as.name(paste0(producer, ".default"))))
    expect_identical(producer_for(complete), "")
    attr(complete, ".npRmpi.bound.call.owner") <- TRUE
    expect_identical(producer_for(complete), producer)
    for (method in c("", ".formula", ".NULL", ".conbandwidth", ".condbandwidth")) {
      call <- as.call(list(as.name(paste0(producer, method))))
      attr(call, ".npRmpi.bound.call.owner") <- TRUE
      expect_identical(producer_for(call), "")
      expect_identical(plan(call)$kind, "none")
    }
  }
  for (producer in c("npregbw", "nplsqregbw", "npplregbw", "npindexbw", "npscoefbw"))
    for (method in c("", ".default", ".rbandwidth"))
      expect_identical(producer_for(as.call(list(as.name(paste0(producer, method))))), producer)
  for (producer in c("npudensbw", "npudistbw", "npconmode", "npqreg"))
    expect_identical(plan(as.call(list(as.name(paste0(producer, ".default")))))$kind, "none")
})


test_that("conditional call-only metadata preserves native state and owner", {
  apply <- getFromNamespace(".npRmpi_bw_formula_apply", "npRmpi")
  for (producer in c("npcdensbw", "npcdistbw")) {
    owner <- new.env(parent = baseenv())
    call <- quote(bandwidth_call(xdat = x, ydat = y))
    environment(call) <- owner
    plan <- list(producer = producer, call = call)
    before <- list(bw = .4, xnames = "x", ynames = "y", call = quote(old()))
    expected <- before; expected$call <- call
    after <- apply(before, plan)
    expect_identical(after, expected)
    expect_identical(apply(after, plan), after)
    expect_identical(environment(after$call), owner)
    expect_null(after[["formula", exact = TRUE]])
  }
})
