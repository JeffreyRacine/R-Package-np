test_that("wire-only subset owners leave public metadata and callbacks intact", {
  ns <- asNamespace("npRmpi")
  wire <- get(".npRmpi_bw_metadata_wire", ns)
  hash <- get(".npRmpi_autodispatch_fingerprint", ns)
  owner <- new.env(parent = baseenv())
  owner$unused <- 1:8
  reads <- 0L
  makeActiveBinding("tripwire", function() { reads <<- reads + 1L; stop("evaluated owner") }, owner)
  frame <- data.frame(x = 1:4)
  attr(frame, ".np.subset") <- list(expression = quote(x < cutoff), envir = owner)
  for (producer in c("npregbw", "npplregbw", "npindexbw", "npscoefbw", "npcdensbw", "npcdistbw")) {
    input <- list(formula = list(producer = producer, frame = frame))
    output <- wire(input)
    actual <- attr(output$formula$frame, ".np.subset")
    expect_identical(actual$expression, quote(x < cutoff), info = producer)
    expect_identical(environmentName(actual$envir), "", info = producer)
    expect_identical(parent.env(actual$envir), emptyenv(), info = producer)
    expect_length(ls(actual$envir, all.names = TRUE), 0L)
    expect_identical(attr(input$formula$frame, ".np.subset")$envir, owner)
    expect_identical(output$formula$frame$x, frame$x)
    expect_identical(hash(output), hash(input))
  }
  for (producer in c("npudensbw", "npudistbw", "other")) {
    output <- wire(list(formula = list(producer = producer, frame = frame)))
    expect_identical(attr(output$formula$frame, ".np.subset")$envir, owner)
  }
  for (reference in list(globalenv(), asNamespace("stats"), structure(new.env(), name = "named-owner"),
                         structure(new.env(), class = "owner-class"))) {
    other <- frame
    attr(other, ".np.subset")$envir <- reference
    output <- wire(list(formula = list(producer = "npregbw", frame = other)))
    expect_identical(attr(output$formula$frame, ".np.subset")$envir, reference)
  }
  callback <- function(x) x
  environment(callback) <- owner
  output <- wire(list(formula = list(producer = "npregbw", frame = frame, na.action = callback)))
  expect_identical(environment(output$formula$na.action), owner)
  expect_identical(reads, 0L)
})


test_that("projected subset owners are restored only on public output", {
  ns <- asNamespace("npRmpi")
  wire <- get(".npRmpi_bw_metadata_wire", ns)
  restore <- get(".npRmpi_bw_restore_public_owner", ns)
  owner <- new.env(parent = baseenv())
  owner$cutoff <- 3L
  formula <- ~x
  environment(formula) <- owner
  frame <- model.frame(formula, data.frame(x = 1:4))
  subset <- list(expression = quote(x < cutoff), envir = owner)
  attr(frame, ".np.subset") <- subset
  original <- structure(list(bw = .4, call = quote(npregbw(~x)),
    formula = formula, terms = terms(frame),
    .np.formula.training = list(frame = frame, na.action = stats::na.omit,
                               subset = subset)), class = "rbandwidth")
  projected <- wire(list(formula = list(producer = "npregbw", frame = frame)))
  worker <- wire(original)
  worker$.np.formula.training$subset <- attr(projected$formula$frame, ".np.subset")
  expect_identical(parent.env(worker$.np.formula.training$subset$envir), emptyenv())
  result <- restore(list(bws = worker), list(original))
  expect_identical(result$bws$.np.formula.training$subset$envir, owner)
  expect_identical(environment(result$bws$formula), owner)
  expect_identical(original$.np.formula.training$subset$envir, owner)
})
