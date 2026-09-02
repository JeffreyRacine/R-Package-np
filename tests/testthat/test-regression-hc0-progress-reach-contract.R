t3_progress_time_counter <- function(start = 0, by = 0.1) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

t3_with_progress_bindings <- function(pkg, bindings, code) {
  ns <- asNamespace(pkg)
  old <- lapply(names(bindings), get, envir = ns, inherits = FALSE)
  names(old) <- names(bindings)

  for (name in names(bindings)) {
    locked <- bindingIsLocked(name, ns)
    if (locked) unlockBinding(name, ns)
    assign(name, bindings[[name]], envir = ns)
    if (locked) lockBinding(name, ns)
  }

  on.exit({
    for (name in names(old)) {
      locked <- bindingIsLocked(name, ns)
      if (locked) unlockBinding(name, ns)
      assign(name, old[[name]], envir = ns)
      if (locked) lockBinding(name, ns)
    }
  }, add = TRUE)

  code()
}

t3_capture_single_line_output <- function(pkg, ansi, code) {
  path <- tempfile(fileext = ".log")
  con <- file(path, open = "wb")
  closed <- FALSE
  on.exit({
    if (!closed) close(con)
    unlink(path)
  }, add = TRUE)

  output <- t3_with_progress_bindings(
    pkg,
    list(
      .np_progress_single_line_connection = function() con,
      .np_progress_single_line_supports_ansi = function(con) ansi,
      .np_progress_output_width = function() 120L
    ),
    code
  )
  invisible(output)
  close(con)
  closed <- TRUE
  readChar(path, nchars = file.info(path)$size, useBytes = TRUE)
}

t3_terminal_replay <- function(bytes) {
  chars <- strsplit(bytes, "", fixed = TRUE)[[1L]]
  line <- character()
  cursor <- 1L
  visible <- character()
  i <- 1L

  snapshot <- function() {
    if (!length(line)) return("")
    sub("[ ]+$", "", paste0(line, collapse = ""))
  }

  while (i <= length(chars)) {
    if (identical(chars[[i]], "\033") &&
        paste0(chars[i:min(i + 3L, length(chars))], collapse = "") == "\033[2K") {
      line <- character()
      cursor <- 1L
      i <- i + 4L
      next
    }
    if (identical(chars[[i]], "\r")) {
      visible <- c(visible, snapshot())
      cursor <- 1L
      i <- i + 1L
      next
    }
    if (identical(chars[[i]], "\n")) {
      visible <- c(visible, snapshot())
      line <- character()
      cursor <- 1L
      i <- i + 1L
      next
    }
    if (length(line) < cursor) line <- c(line, rep.int(" ", cursor - length(line)))
    line[[cursor]] <- chars[[i]]
    cursor <- cursor + 1L
    i <- i + 1L
  }

  c(visible, snapshot())
}

t3_hc0_fixture <- function(degree) {
  n <- 24L
  x <- data.frame(x = seq(-1, 1, length.out = n))
  y <- sin(1.3 * x$x) + seq_len(n) / 100
  ex <- data.frame(x = c(-0.75, -0.2, 0.35))
  bw <- npregbw(
    xdat = x,
    ydat = y,
    bws = 0.45,
    bandwidth.compute = FALSE,
    bwscaling = FALSE,
    bwtype = "fixed",
    regtype = "lp",
    degree = degree,
    degree.select = "manual"
  )
  list(x = x, y = y, ex = ex, bw = bw)
}

test_that("HC0 fit totals cover scalar and general LP over each bandwidth topology", {
  topologies <- c("fixed", "generalized_nn", "adaptive_nn")
  engines <- c(REGTYPE_LP0, REGTYPE_LP)

  for (topology in topologies) {
    bws <- list(type = topology)
    base <- if (identical(topology, "adaptive_nn")) 24L else 3L
    for (engine in engines) {
      expect_identical(.np_reg_fit_total(bws, 24L, 3L, FALSE, engine), base)
      expect_identical(.np_reg_fit_total(bws, 24L, 3L, TRUE, engine), 24L + base)
    }
    expect_identical(.np_reg_fit_total(bws, 24L, 3L, TRUE, 99L), base)
  }

  expect_identical(
    .np_reg_fit_total(list(type = "fixed"), 1L, 1L, TRUE, REGTYPE_LP),
    2L
  )
})

test_that("HC0 native progress reaches both phases for degrees zero through three", {
  old <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old), add = TRUE)

  for (degree in 0:3) {
    fixture <- t3_hc0_fixture(degree)
    actual <- capture_progress_shadow_trace(
      npreg(
        bws = fixture$bw,
        txdat = fixture$x,
        tydat = fixture$y,
        exdat = fixture$ex,
        gradients = TRUE,
        se = TRUE
      ),
      force_renderer = "single_line",
      now = t3_progress_time_counter()
    )

    trace <- actual$trace
    current <- vapply(trace, function(entry) {
      if (is.null(entry$current)) NA_integer_ else as.integer(entry$current)
    }, integer(1L))
    total <- vapply(trace, function(entry) {
      if (is.null(entry$total)) NA_integer_ else as.integer(entry$total)
    }, integer(1L))
    times <- vapply(trace, `[[`, numeric(1L), "now")
    event <- vapply(trace, `[[`, character(1L), "event")
    info <- paste("degree", degree, paste(current, collapse = ","))

    expect_s3_class(actual$value, "npregression")
    expect_true(length(trace) >= 3L, info = info)
    observed <- current[!is.na(current)]
    expect_true(all(total[!is.na(total)] == 27L), info = info)
    expect_true(all(diff(observed) >= 0L), info = info)
    expect_lte(max(observed), 27L)
    expect_true(any(observed > 0L & observed <= 24L), info = info)
    expect_true(any(observed > 24L), info = info)
    expect_identical(sum(event == "finish"), 1L, info = info)
    expect_identical(tail(observed, 1L), 27L, info = info)
    expect_true(all(diff(times) <= 0.1 + sqrt(.Machine$double.eps)), info = info)
  }
})

test_that("HC0 fit lines replay without ANSI, CR, padding, or newline residue", {
  pkg <- environmentName(environment(npreg))
  render <- getFromNamespace(".np_progress_render_single_line", pkg)
  long <- sprintf("[%s] Fitting regression 24/27 (88.9%%, elapsed 0.2s, eta 0.0s)", pkg)
  done <- sprintf("[%s] Fitting regression 27/27 (100.0%%, elapsed 0.3s, eta 0.0s)", pkg)

  for (ansi in c(FALSE, TRUE)) {
    output <- t3_capture_single_line_output(pkg, ansi, function() {
      render(list(render_line = long, last_width = 0L), event = "render")
      render(list(render_line = done, last_width = nchar(long, type = "width")),
             event = "render")
      render(list(render_line = done, last_width = nchar(done, type = "width")),
             event = "finish")
    })
    visible <- t3_terminal_replay(output)
    expect_true(long %in% visible, info = paste("ansi", ansi))
    expect_true(done %in% visible, info = paste("ansi", ansi))
    expect_identical(tail(visible, 1L), "", info = paste("ansi", ansi))
  }
})
