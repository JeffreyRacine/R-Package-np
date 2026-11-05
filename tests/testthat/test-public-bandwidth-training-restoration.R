test_that("owner restoration retains effective frames and immutable policies", {
  ns <- asNamespace("npRmpi")
  restore <- get(".npRmpi_bw_restore_public_owner", ns)
  wire <- get(".npRmpi_bw_metadata_wire", ns)
  owner <- new.env(parent = baseenv())
  owner$offset <- 3
  f <- y ~ x; environment(f) <- owner
  frame <- model.frame(f, data.frame(y = 1:6, x = 7:12))
  call <- quote(npregbw(y ~ x, bws = .3, bandwidth.compute = FALSE))
  environment(call) <- owner
  policy <- function(x) stats::na.omit(x)
  environment(policy) <- owner
  reference <- structure(list(call = call, formula = f, terms = attr(frame, "terms"),
    .np.formula.training = list(frame = frame, na.action = policy,
      subset = list(expression = quote(x > 0), environment = owner)),
    bw = .3), class = "rbandwidth")
  value <- wire(reference)
  new.frame <- frame[c(2L, 4L, 6L), , drop = FALSE]
  new.frame$y <- new.frame$y + 10
  tt <- attr(frame, "terms"); environment(tt) <- emptyenv()
  attr(tt, "dataClasses") <- NULL
  attr(new.frame, "terms") <- tt
  value[[".np.formula.training"]]$frame <- new.frame
  result <- restore(list(bws = value, mean = 1:3), list(reference))
  expected <- new.frame
  tt <- attr(expected, "terms"); environment(tt) <- owner
  attr(expected, "terms") <- tt
  expect_identical(result$bws[[".np.formula.training"]]$frame, expected)
  expect_identical(result$bws[[".np.formula.training"]]$na.action, policy)
  expect_identical(result$bws[[".np.formula.training"]]$subset,
                   reference[[".np.formula.training"]]$subset)
  expect_identical(result$bws$formula, reference$formula)
  expect_identical(result$bws$call, reference$call)
  expect_identical(result$mean, 1:3)
  expect_identical(result$bws$bw, .3)
  expect_identical(restore(result, list(reference)), result)
  expect_identical(reference[[".np.formula.training"]]$frame, frame)
  for (field in c("formula", "terms", "call")) {
    bad <- value
    bad[[field]] <- if (field == "terms") terms(y ~ z) else if (field == "formula") y ~ z else quote(npregbw(y ~ z))
    expect_error(restore(list(bws = bad), list(reference)),
                 "retained formula metadata changed", info = field)
  }
  bad <- value; bad[[".np.formula.training"]]["na.action"] <- list(NULL)
  expect_error(restore(list(bws = bad), list(reference)), "retained formula metadata changed")
  bad <- value; bad[[".np.formula.training"]]$subset$expression <- quote(x < 0)
  expect_error(restore(list(bws = bad), list(reference)), "retained formula metadata changed")
  bad <- value; bad[[".np.formula.training"]]$frame <- NULL
  expect_error(restore(list(bws = bad), list(reference)), "retained formula metadata changed")
  for (shape in c("absent", "null", "frame-null", "policy-null")) {
    ref <- reference
    if (shape == "absent") ref[[".np.formula.training"]] <- NULL
    if (shape == "null") ref[".np.formula.training"] <- list(NULL)
    if (shape == "frame-null") ref[[".np.formula.training"]]["frame"] <- list(NULL)
    if (shape == "policy-null") ref[[".np.formula.training"]]["na.action"] <- list(NULL)
    got <- restore(list(bws = wire(ref)), list(ref))
    expect_identical(got$bws, ref, info = shape)
  }
  parent <- reference; class(parent) <- "lsqregressionbandwidth"
  parent$tau.bws <- list(reference, reference)
  child.value <- value; class(child.value) <- "lsqregressionbandwidth"
  child.value$tau.bws <- list(value, value)
  nested <- restore(list(bws = child.value), list(parent))
  expect_identical(lapply(nested$bws$tau.bws, function(x) x[[".np.formula.training"]]$frame),
                   list(expected, expected))
})
