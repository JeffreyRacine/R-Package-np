progress_time_counter <- function(start = 0, by = 1.3) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

shadow_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

test_that("distribution and conditional public routes use the generic bandwidth selection label", {
  set.seed(321)
  n <- 24
  x <- runif(n)
  y <- x^2 + rnorm(n, sd = 0.1)

  cases <- list(
    quote(npudistbw(dat = data.frame(x = x), nmulti = 2)),
    quote(npudist(tdat = data.frame(x = x), nmulti = 2)),
    quote(npcdensbw(xdat = data.frame(x = x), ydat = data.frame(y = y), regtype = "lc", nmulti = 2)),
    quote(npcdens(txdat = data.frame(x = x), tydat = data.frame(y = y), regtype = "lc", nmulti = 2)),
    quote(npcdistbw(xdat = data.frame(x = x), ydat = data.frame(y = y), regtype = "lc", nmulti = 2)),
    quote(npcdist(txdat = data.frame(x = x), tydat = data.frame(y = y), regtype = "lc", nmulti = 2)),
    quote(npplregbw(xdat = data.frame(x = x), zdat = data.frame(z = x), ydat = y, nmulti = 2)),
    quote(npplreg(y ~ x | z, data = data.frame(x = x, z = x, y = y), nmulti = 2))
  )

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  for (expr in cases) {
    actual <- capture_progress_shadow_trace(
      eval(expr),
      force_renderer = "single_line",
      now = progress_time_counter()
    )

    lines <- shadow_lines(actual)
    info <- paste(deparse(expr), collapse = "")

    expect_true(any(grepl("^\\[np\\] Bandwidth selection", lines)), info = info)
    expect_true(any(grepl("100\\.0%, eta 0\\.0s\\)$", lines)), info = info)
  }
})
