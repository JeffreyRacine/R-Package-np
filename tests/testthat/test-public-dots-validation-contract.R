public_dots_routes <- c(
  "npreg", "npregbw", "npudens", "npudensbw", "npcdens", "npcdensbw",
  "npudist", "npudistbw", "npcdist", "npcdistbw", "npscoef", "npscoefbw",
  "npindex", "npindexbw", "npplreg", "npplregbw", "npqreg", "npconmode",
  "npcopula"
)

test_that("principal public dots routes reject unknown names before forcing values", {
  for (route in public_dots_routes) {
    call <- as.call(c(
      list(as.name(route)),
      structure(list(quote(stop("unknown dot was forced"))),
                names = "definitely.not.an.np.argument")
    ))
    expect_error(
      eval(call),
      sprintf("%s\\(\\): unused argument 'definitely.not.an.np.argument'", route),
      info = route
    )
  }
})

test_that("public dots diagnostics make only audited suggestions", {
  expect_error(
    npreg(noqmad = stop("unknown dot was forced")),
    "unused argument 'noqmad'; did you mean 'nomad'\\?"
  )
  expect_error(
    npreg(bernstein = stop("unknown dot was forced")),
    "unused argument 'bernstein'; did you mean 'bernstein.basis'\\?"
  )
  expect_error(
    npconmode(gradient.level = stop("unknown dot was forced")),
    "unused argument 'gradient.level'; did you mean 'level'\\?"
  )

  err <- tryCatch(
    npudensbw(edat = stop("unknown dot was forced")),
    error = conditionMessage
  )
  expect_match(err, "unused argument 'edat'", fixed = TRUE)
  expect_false(grepl("did you mean", err, fixed = TRUE))
})

test_that("unnamed and duplicate valid dots pass through the name validator", {
  validate <- getFromNamespace(".np_validate_public_dots", environmentName(environment(npreg)))
  expect_invisible(validate(pairlist(quote(stop("unnamed dot was forced"))), "npreg"))

  duplicate <- as.pairlist(structure(
    list(quote(stop("first duplicate was forced")),
         quote(stop("second duplicate was forced"))),
    names = c("nomad", "nomad")
  ))
  expect_invisible(validate(duplicate, "npreg"))
})

test_that("npindex admits only its exact internal fixed-lc progress token", {
  validate <- getFromNamespace(".np_validate_public_dots", environmentName(environment(npindex)))
  expect_invisible(validate(
    pairlist(.np_lc_fixed_progress_route = quote(stop("internal token was forced"))),
    "npindex"
  ))
  expect_error(
    validate(pairlist(.np_lc_fixed_progress_routes = TRUE), "npindex"),
    "unused argument '.np_lc_fixed_progress_routes'",
    fixed = TRUE
  )
})

test_that("specialized retired-argument diagnostics retain precedence", {
  expect_error(npreg(errors = TRUE), "use 'se' instead", fixed = TRUE)
  expect_error(npindex(boot.num = 9L), "use 'B' instead", fixed = TRUE)
  expect_error(
    npregbw(cfac.init = 0.5),
    "'cfac.init' has been renamed to 'scale.factor.init'",
    fixed = TRUE
  )
  expect_error(
    npcdensbw(
      xdat = data.frame(x = c(0, 0.5, 1)),
      ydat = data.frame(y = c(0, 0.5, 1)),
      bws = c(0.3, 0.3),
      bandwidth.compute = FALSE,
      cvls.i1.rescue = FALSE
    ),
    "cvls.i1.rescue has been removed",
    fixed = TRUE
  )
})
