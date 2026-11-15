test_that("LSQ wire metadata omits only its unused projected frame", {
  wire <- getFromNamespace(".npRmpi_bw_metadata_wire", "npRmpi")
  f <- y ~ x
  environment(f) <- baseenv()
  frame <- model.frame(f, data.frame(x = 1:6, y = 6:1))
  for (child in c(FALSE, TRUE)) {
    plan <- list(producer = "nplsqregbw", formula = f, terms = terms(f),
      frame = frame, omit = structure(2L, class = "exclude"), ynames = "y",
      scalar.children = TRUE, scalar.child = child,
      call = quote(nplsqregbw(y ~ x, tau = c(.25, .75))))
    metadata <- list(identity = list(type = "native"), formula = plan)
    before <- serialize(metadata, NULL)
    got <- wire(metadata)
    control <- metadata
    control$formula$producer <- "npregbw"
    expected <- wire(control)
    expected$formula$producer <- "nplsqregbw"
    expected$formula$frame <- NULL
    expect_identical(got, expected)
    expect_false("frame" %in% names(got$formula))
    expect_identical(got$formula$omit, plan$omit)
    expect_identical(got$formula$scalar.child, child)
    expect_identical(serialize(metadata, NULL), before)
    expect_lt(length(serialize(got, NULL)), length(serialize(wire(control), NULL)))
  }
  expect_null(wire(NULL))
  expect_identical(wire(list(identity = list(a = 1L), formula = NULL)),
                   list(identity = list(a = 1L), formula = NULL))
  for (producer in c("npregbw", "npindexbw", "npplregbw", "npscoefbw")) {
    got <- wire(list(formula = list(producer = producer, frame = frame)))
    expect_identical(got$formula$frame, wire(list(frame = frame))$frame)
  }
})

test_that("LSQ wire omission retains prior projection callbacks", {
  wire <- getFromNamespace(".npRmpi_bw_metadata_wire", "npRmpi")
  calls <- 0L
  method <- "[[.np_wire_frame_probe"
  prior <- get0(method, envir = globalenv(), inherits = FALSE)
  assign(method, function(x, ...) {
    calls <<- calls + 1L
    NextMethod("[[")
  }, envir = globalenv())
  on.exit(if (is.null(prior)) rm(list = method, envir = globalenv()) else
    assign(method, prior, envir = globalenv()), add = TRUE)
  frame <- structure(list(y = 1:5, x = 6:10), class = "np_wire_frame_probe")
  result <- wire(list(formula = list(producer = "nplsqregbw", frame = frame)))
  expect_identical(calls, 2L)
  expect_false("frame" %in% names(result$formula))
  # Unusual outer containers keep their existing shape and payload.
  outer <- structure(list(formula = list(producer = "nplsqregbw", frame = 1:5)),
                     class = "np_wire_outer_probe")
  expect_identical(wire(outer), outer)
  plan <- structure(list(producer = "nplsqregbw", frame = 1:5),
                    class = "np_wire_plan_probe")
  expect_identical(wire(list(formula = plan)), list(formula = plan))
})
