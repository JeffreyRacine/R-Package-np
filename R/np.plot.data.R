# Plot results describe computed geometry, not fitted estimators. In particular,
# interval endpoints are absolute values and must never be passed through se().
.np_plot_data_empty_intervals <- function() {
  data.frame(row_id = integer(), target_id = character(), band = character(),
             alpha = double(), lower = double(), upper = double(),
             family_id = character(), status = character(),
             stringsAsFactors = FALSE)
}

.np_plot_data_targets <- function() {
  out <- data.frame(target_id = character(), quantity = character(),
    label = character(), variable = character(), level = character(),
    term = character(), inference_id = character(), tau = double(),
    order = integer(), stringsAsFactors = FALSE)
  out$multi_index <- I(list())
  out
}

.np_plot_data_target <- function(id, quantity, label = quantity,
                                 variable = NA_character_, level = NA_character_,
                                 term = NA_character_, inference_id = NA_character_,
                                 tau = NA_real_, order = NA_integer_,
                                 multi_index = integer()) {
  out <- data.frame(target_id = id, quantity = quantity, label = label,
    variable = variable, level = level, term = term, inference_id = inference_id,
    tau = as.double(tau), order = as.integer(order), stringsAsFactors = FALSE)
  out$multi_index <- I(list(as.integer(multi_index)))
  out
}

.np_plot_data_panel <- function(grid, roles, values, intervals, display = list()) {
  stopifnot(is.data.frame(grid), length(roles) == ncol(grid))
  coordinates <- paste0("eval_", seq_len(ncol(grid)))
  variables <- data.frame(coordinate = coordinates, name = names(grid),
    role = roles, label = names(grid), stringsAsFactors = FALSE)
  names(grid) <- coordinates
  list(grid = grid, variables = variables, values = values,
       intervals = intervals, display = display)
}

.np_plot_data_values <- function(estimate, target_id,
                                 bias_corrected = rep(NA_real_, length(estimate)),
                                 std_error = rep(NA_real_, length(estimate))) {
  n <- length(estimate)
  if (length(bias_corrected) != n || length(std_error) != n)
    stop("plot target values have inconsistent row counts", call. = FALSE)
  data.frame(row_id = seq_len(n), target_id = target_id,
    estimate = as.double(estimate), bias_corrected = as.double(bias_corrected),
    std_error = as.double(std_error),
    status = ifelse(is.finite(estimate), "available", "undefined"),
    stringsAsFactors = FALSE)
}

.np_plot_data_intervals <- function(bands, target_id, alpha,
                                    family_id = NA_character_,
                                    status = NULL) {
  if (!length(bands)) return(.np_plot_data_empty_intervals())
  out <- lapply(names(bands), function(band) {
    endpoints <- bands[[band]]
    if (!is.matrix(endpoints) || ncol(endpoints) != 2L)
      stop("plot band endpoints must have exactly two columns", call. = FALSE)
    n <- nrow(endpoints)
    state <- if (is.null(status)) {
      ifelse(is.finite(endpoints[, 1L]) & is.finite(endpoints[, 2L]),
             "available", "undefined")
    } else rep_len(status, n)
    data.frame(row_id = seq_len(n), target_id = target_id, band = band,
      alpha = as.double(alpha), lower = as.double(endpoints[, 1L]),
      upper = as.double(endpoints[, 2L]), family_id = family_id, status = state,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.np_plot_data_validate <- function(x) {
  fail <- function(message) stop(paste0("invalid np_plot_data: ", message), call. = FALSE)
  if (!is.list(x) || !identical(x$schema_version, 1L) ||
      !is.character(x$family) || length(x$family) != 1L || is.na(x$family))
    fail("unknown schema or family")
  if (!is.list(x$panels) || is.null(names(x$panels)) ||
      anyNA(names(x$panels)) || any(!nzchar(names(x$panels))) ||
      anyDuplicated(names(x$panels))) fail("panels require unique identifiers")
  if (!is.data.frame(x$targets) || anyNA(x$targets$target_id) ||
      anyDuplicated(x$targets$target_id)) fail("targets require unique identifiers")
  if (!identical(names(x), c("schema_version", "family", "request", "panels", "targets",
                            "inference", "coverage_families", "diagnostics")))
    fail("unexpected top-level fields")
  if (!identical(names(x$targets), names(.np_plot_data_targets())) ||
      !is.character(x$targets$target_id) || !is.character(x$targets$inference_id) ||
      !is.double(x$targets$tau) || !is.integer(x$targets$order) ||
      !is.list(x$targets$multi_index)) fail("invalid target fields")
  ids <- x$targets$inference_id[!is.na(x$targets$inference_id)]
  if (!is.list(x$inference) || anyDuplicated(names(x$inference)) ||
      any(!ids %in% names(x$inference))) fail("undefined inference record")
  panel.keys <- vector("list", length(x$panels))
  names(panel.keys) <- names(x$panels)
  for (panel_id in names(x$panels)) {
    panel <- x$panels[[panel_id]]
    if (!identical(names(panel), c("grid", "variables", "values", "intervals", "display")))
      fail("unexpected panel fields")
    if (!is.data.frame(panel$grid) || !is.data.frame(panel$variables) ||
        !identical(panel$variables$coordinate, names(panel$grid)))
      fail("coordinate metadata does not match the grid")
    for (table in c("values", "intervals")) {
      d <- panel[[table]]
      if (!is.data.frame(d) || !is.integer(d$row_id) ||
          anyNA(d$row_id) || any(d$row_id < 1L | d$row_id > nrow(panel$grid)) ||
          anyNA(d$target_id) || any(!d$target_id %in% x$targets$target_id))
        fail(paste("invalid row/target keys in", table))
      fields <- if (table == "values")
        c("row_id", "target_id", "estimate", "bias_corrected", "std_error", "status") else
        names(.np_plot_data_empty_intervals())
      numbers <- if (table == "values") c("estimate", "bias_corrected", "std_error") else
        c("alpha", "lower", "upper")
      strings <- setdiff(fields, c("row_id", numbers))
      if (!identical(names(d), fields) ||
          !all(vapply(d[numbers], is.double, logical(1L))) ||
          !all(vapply(d[strings], is.character, logical(1L))) || anyNA(d$status))
        fail(paste("invalid fields in", table))
      keys <- c("row_id", "target_id", if (table == "intervals") c("band", "alpha"))
      if (anyDuplicated(d[keys])) fail(paste("duplicate keys in", table))
    }
    intervals <- panel$intervals
    if (!identical(names(intervals), names(.np_plot_data_empty_intervals())))
      fail("unexpected interval fields")
    finite <- is.finite(intervals$lower) & is.finite(intervals$upper)
    if (any(intervals$lower[finite] > intervals$upper[finite]))
      fail("reversed interval endpoints")
    value.key <- paste(panel$values$row_id, panel$values$target_id, sep = ":")
    panel.keys[[panel_id]] <- value.key
    band.key <- paste(intervals$row_id, intervals$target_id, sep = ":")
    if (any(!band.key %in% value.key)) fail("interval without a target value")
    ids <- intervals$family_id[!is.na(intervals$family_id)]
    if (any(!ids %in% names(x$coverage_families))) fail("undefined coverage family")
  }
  for (family in x$coverage_families) {
    members <- family$members
    if (!is.character(family$event) || length(family$event) != 1L ||
        !is.character(family$scope) || length(family$scope) != 1L ||
        !is.data.frame(members) ||
        !identical(names(members), c("panel_id", "row_id", "target_id")) ||
        (nrow(members) > 1L && anyDuplicated(members)))
      fail("invalid coverage-family definition")
    for (panel_id in unique(members$panel_id)) {
      if (!panel_id %in% names(x$panels)) fail("coverage family has an unknown panel")
      selected <- which(members$panel_id == panel_id)
      if (any(!paste(members$row_id[selected], members$target_id[selected],
                     sep = ":") %in% panel.keys[[panel_id]]))
        fail("coverage family has an unknown target row")
    }
  }
  invisible(x)
}

.np_plot_data <- function(family, request, panels, targets, inference = list(),
                          coverage_families = list(), diagnostics = list()) {
  out <- structure(list(schema_version = 1L, family = family, request = request,
    panels = panels, targets = targets, inference = inference,
    coverage_families = coverage_families, diagnostics = diagnostics),
    class = "np_plot_data")
  out <- .np_plot_data_family_events(out)
  .np_plot_data_validate(out)
  out
}

.np_plot_data_family_events <- function(x) {
  records <- list()
  for (panel_id in names(x$panels)) {
    d <- x$panels[[panel_id]]$intervals
    keep <- which(d$status != "display_only")
    x$panels[[panel_id]]$intervals$family_id[] <- NA_character_
    if (!length(keep)) next
    d <- d[keep, , drop = FALSE]
    scope <- vapply(d$family_id, function(id) {
      prior <- if (!is.na(id)) x$coverage_families[[id]] else NULL
      if (is.null(prior$scope)) "finite grid within panel and scalar target" else prior$scope
    }, character(1L))
    marginal <- d$band %in% c("pointwise", "pmzsd")
    union <- scope == "union of all supplied plot slices" & !marginal
    key <- paste(ifelse(union, "union", paste(panel_id, d$target_id, sep = ":")),
      d$band, d$alpha, ifelse(marginal, d$row_id, "all"), sep = ":")
    records[[panel_id]] <- data.frame(panel_id = panel_id, interval_row = keep,
      row_id = d$row_id, target_id = d$target_id, key = key,
      event = ifelse(marginal, "this scalar target row is covered",
        "every member of this finite target family is covered"),
      scope = ifelse(marginal, "pointwise; no joint coverage claim", scope))
  }
  if (!length(records)) { x$coverage_families <- list(); return(x) }
  d <- do.call(rbind, records)
  group <- match(d$key, unique(d$key))
  rows <- split(seq_len(nrow(d)), factor(group, levels = seq_len(max(group))))
  # These are ordinary, already-typed columns. Avoid dispatching a data-frame
  # slice/unique for every pointwise scalar family while preserving row names.
  members <- d[c("panel_id", "row_id", "target_id")]
  member.names <- row.names(d)
  x$coverage_families <- lapply(rows, function(i) {
    part <- structure(lapply(members, function(column) column[i]),
                      class = "data.frame", row.names = member.names[i])
    if (length(i) > 1L) part <- unique(part)
    list(event = d$event[i[1L]], scope = d$scope[i[1L]], members = part)
  })
  names(x$coverage_families) <- paste0("family_", seq_along(rows))
  ids <- names(x$coverage_families)[group]
  for (panel_id in names(records)) {
    i <- which(d$panel_id == panel_id)
    x$panels[[panel_id]]$intervals$family_id[d$interval_row[i]] <- ids[i]
  }
  x
}

print.np_plot_data <- function(x, ...) {
  cat("Plot data: ", x$family, " (schema ", x$schema_version, ")\n", sep = "")
  cat(length(x$panels), "panel(s),", nrow(x$targets), "target(s),",
      sum(vapply(x$panels, function(p) nrow(p$intervals), integer(1L))),
      "interval row(s)\n")
  print(x$targets[c("target_id", "quantity", "label")], row.names = FALSE)
  invisible(x)
}

.np_plot_data_select <- function(value, choices, what) {
  if (is.null(value)) {
    if (length(choices) != 1L)
      stop(sprintf("select %s explicitly; available: %s", what,
                   paste(choices, collapse = ", ")), call. = FALSE)
    return(choices)
  }
  if (is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == floor(value) && value >= 1L && value <= length(choices))
    return(choices[as.integer(value)])
  if (is.character(value) && length(value) == 1L && !is.na(value) && value %in% choices)
    return(value)
  stop(sprintf("invalid %s selection; available: %s", what,
               paste(choices, collapse = ", ")), call. = FALSE)
}

as.data.frame.np_plot_data <- function(x, row.names = NULL, optional = FALSE,
                                       panel = NULL, target = NULL, band = NULL, ...) {
  if (length(list(...))) stop("unused plot-data extraction arguments", call. = FALSE)
  .np_plot_data_validate(x)
  panel_id <- .np_plot_data_select(panel, names(x$panels), "panel")
  p <- x$panels[[panel_id]]
  target_id <- .np_plot_data_select(target, unique(p$values$target_id), "target")
  values <- p$values[p$values$target_id == target_id, , drop = FALSE]
  intervals <- p$intervals[p$intervals$target_id == target_id, , drop = FALSE]
  out <- cbind(p$grid[values$row_id, , drop = FALSE], values)
  if (nrow(intervals)) {
    band_id <- .np_plot_data_select(band, unique(intervals$band), "band")
    intervals <- intervals[intervals$band == band_id, , drop = FALSE]
    key <- match(values$row_id, intervals$row_id)
    if (anyNA(key) || nrow(intervals) != nrow(values))
      stop("plot band does not cover every retained target row", call. = FALSE)
    intervals <- intervals[key, setdiff(names(intervals), c("row_id", "target_id")), drop = FALSE]
    names(intervals)[names(intervals) == "status"] <- "interval_status"
    out <- cbind(out, intervals)
  } else if (!is.null(band)) {
    stop("no intervals were computed for this target", call. = FALSE)
  }
  rownames(out) <- row.names
  out
}
