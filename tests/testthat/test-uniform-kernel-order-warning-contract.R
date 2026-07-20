test_that("uniform kernel order warnings require an explicit user order", {
  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)

  collect_warnings <- function(expr) {
    messages <- character()
    value <- withCallingHandlers(
      force(expr),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(value = value, messages = messages)
  }

  x <- data.frame(x = seq(0.1, 0.9, length.out = 12L))
  y <- sin(2 * pi * x$x)
  y.frame <- data.frame(y = y)
  expected <- "[npRmpi] ignoring kernel order specified with uniform kernel type"

  calls <- list(
    regression = list(
      fun = npregbw,
      args = list(
        xdat = x, ydat = y, bws = 0.2, bandwidth.compute = FALSE,
        ckertype = "uniform"
      ),
      order = list(ckerorder = 4L)
    ),
    density = list(
      fun = npudensbw,
      args = list(
        dat = x, bws = 0.2, bandwidth.compute = FALSE,
        ckertype = "uniform"
      ),
      order = list(ckerorder = 4L)
    ),
    distribution = list(
      fun = npudistbw,
      args = list(
        dat = x, bws = 0.2, bandwidth.compute = FALSE,
        ckertype = "uniform"
      ),
      order = list(ckerorder = 4L)
    ),
    conditional_density_x = list(
      fun = npcdensbw,
      args = list(
        xdat = x, ydat = y.frame, bws = c(0.3, 0.2),
        bandwidth.compute = FALSE, cxkertype = "uniform"
      ),
      order = list(cxkerorder = 4L)
    ),
    conditional_distribution_x = list(
      fun = npcdistbw,
      args = list(
        xdat = x, ydat = y.frame, bws = c(0.3, 0.2),
        bandwidth.compute = FALSE, cxkertype = "uniform"
      ),
      order = list(cxkerorder = 4L)
    )
  )

  for (case in calls) {
    omitted <- collect_warnings(do.call(case$fun, case$args))
    expect_identical(omitted$messages, character())

    explicit <- collect_warnings(do.call(case$fun, c(case$args, case$order)))
    expect_identical(explicit$messages, expected)
  }

  formula.data <- data.frame(y = y, x = x$x)
  formula.calls <- list(
    regression = list(
      fun = npregbw,
      args = list(
        formula = y ~ x, data = formula.data, bws = 0.2,
        bandwidth.compute = FALSE, ckertype = "uniform"
      ),
      order = list(ckerorder = 4L)
    ),
    density = list(
      fun = npudensbw,
      args = list(
        formula = ~ x, data = formula.data, bws = 0.2,
        bandwidth.compute = FALSE, ckertype = "uniform"
      ),
      order = list(ckerorder = 4L)
    ),
    distribution = list(
      fun = npudistbw,
      args = list(
        formula = ~ x, data = formula.data, bws = 0.2,
        bandwidth.compute = FALSE, ckertype = "uniform"
      ),
      order = list(ckerorder = 4L)
    ),
    conditional_density = list(
      fun = npcdensbw,
      args = list(
        formula = y ~ x, data = formula.data, bws = c(0.3, 0.2),
        bandwidth.compute = FALSE, cxkertype = "uniform"
      ),
      order = list(cxkerorder = 4L)
    ),
    conditional_distribution = list(
      fun = npcdistbw,
      args = list(
        formula = y ~ x, data = formula.data, bws = c(0.3, 0.2),
        bandwidth.compute = FALSE, cxkertype = "uniform"
      ),
      order = list(cxkerorder = 4L)
    )
  )

  for (case in formula.calls) {
    expect_identical(
      collect_warnings(do.call(case$fun, case$args))$messages,
      character()
    )
    expect_identical(
      collect_warnings(do.call(case$fun, c(case$args, case$order)))$messages,
      expected
    )
  }

  estimator.calls <- list(
    regression = list(
      fun = npreg,
      args = list(txdat = x, tydat = y, bws = 0.2, ckertype = "uniform"),
      order = list(ckerorder = 4L)
    ),
    density = list(
      fun = npudens,
      args = list(tdat = x, bws = 0.2, ckertype = "uniform"),
      order = list(ckerorder = 4L)
    ),
    distribution = list(
      fun = npudist,
      args = list(tdat = x, bws = 0.2, ckertype = "uniform"),
      order = list(ckerorder = 4L)
    ),
    conditional_density = list(
      fun = npcdens,
      args = list(
        txdat = x, tydat = y.frame, bws = c(0.3, 0.2),
        cxkertype = "uniform"
      ),
      order = list(cxkerorder = 4L)
    ),
    conditional_distribution = list(
      fun = npcdist,
      args = list(
        txdat = x, tydat = y.frame, bws = c(0.3, 0.2),
        cxkertype = "uniform"
      ),
      order = list(cxkerorder = 4L)
    )
  )

  for (case in estimator.calls) {
    expect_identical(
      collect_warnings(do.call(case$fun, case$args))$messages,
      character()
    )
    expect_identical(
      collect_warnings(do.call(case$fun, c(case$args, case$order)))$messages,
      expected
    )
  }

  conditional_y <- list(
    xdat = x, ydat = y.frame, bws = c(0.3, 0.2),
    bandwidth.compute = FALSE, cykertype = "uniform"
  )
  expect_identical(
    collect_warnings(do.call(npcdensbw, conditional_y))$messages,
    character()
  )
  expect_identical(
    collect_warnings(do.call(
      npcdensbw, c(conditional_y, list(cykerorder = 4L))
    ))$messages,
    expected
  )
})
