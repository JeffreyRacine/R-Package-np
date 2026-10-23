# Bias controls describe auxiliary inference only; they never select a model.
np_bias_control <- function(construction = c("iterated", "legacy"),
                            iterations = 16L, degree = 5L,
                            bernstein.basis = TRUE,
                            p0.scale = 4, p1.scale = 1.5,
                            operator.scale = 4/3, variance.scale = 4/3,
                            p0.bandwidth = NULL, p1.bandwidth = NULL,
                            operator.bandwidth = NULL, variance.bandwidth = NULL) {
  construction <- match.arg(construction)
  supplied <- names(match.call(expand.dots = FALSE))[-1L]
  if (construction == "legacy" && length(setdiff(supplied, "construction")))
    stop("auxiliary bias controls do not apply to construction='legacy'", call. = FALSE)
  for (name in c("iterations", "degree")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max)
      stop(name, " must be a positive integer", call. = FALSE)
  }
  bernstein.basis <- npValidateScalarLogical(bernstein.basis, "bernstein.basis")
  if (degree > 1L && !bernstein.basis)
    stop("auxiliary degree > 1 requires bernstein.basis=TRUE", call. = FALSE)
  requested <- list()
  for (role in c("p0", "p1", "operator", "variance")) {
    scale.name <- paste0(role, ".scale")
    raw.name <- paste0(role, ".bandwidth")
    raw <- get(raw.name)
    if (!is.null(raw) && scale.name %in% supplied)
      stop(role, " scale and bandwidth are mutually exclusive", call. = FALSE)
    name <- if (!is.null(raw)) raw.name else scale.name
    value <- get(name)
    if (!is.numeric(value) || !is.null(dim(value)) || !length(value) || any(!is.finite(value)) || any(value <= 0) ||
        (!is.null(names(value)) && (anyNA(names(value)) || any(!nzchar(names(value))) || anyDuplicated(names(value)))))
      stop(name, " must contain finite positive values with unique names, if named", call. = FALSE)
    requested[[name]] <- value
  }
  structure(list(construction = construction, iterations = as.integer(iterations),
    degree = as.integer(degree), bernstein.basis = bernstein.basis,
    bandwidth.controls = requested, supplied = supplied, version = 1L),
    class = "np_bias_control")
}

.np_bias_vector <- function(value, columns, name) {
  d <- length(columns)
  if (!is.null(names(value))) {
    if (length(value) != d || !setequal(names(value), columns))
      stop(name, " names must match every continuous predictor exactly", call. = FALSE)
    value <- value[columns]
  } else if (length(value) == 1L) value <- rep(value, d)
  else if (length(value) != d)
    stop(name, " must be scalar or have one value per continuous predictor", call. = FALSE)
  setNames(as.double(value), columns)
}

# Check the construction before dispatch, including caller-edited list fields.
.np_bias_validate_control <- function(control) {
  if (!inherits(control, "np_bias_control") ||
      !identical(control[["version"]], 1L) ||
      !identical(length(control[["construction"]]), 1L) ||
      !control[["construction"]] %in% c("iterated", "legacy"))
    stop("invalid np_bias_control object", call. = FALSE)
  if (control[["construction"]] == "legacy") {
    if (!identical(control, np_bias_control(construction="legacy")))
      stop("legacy bias controls cannot contain auxiliary overrides", call. = FALSE)
  } else {
    needed <- c("iterations", "degree", "bernstein.basis", "bandwidth.controls")
    if (!all(needed %in% names(control)) ||
        !is.list(control[["bandwidth.controls"]]) ||
        length(control[["bandwidth.controls"]]) != 4L)
      stop("invalid auxiliary bias controls", call. = FALSE)
    checked <- do.call(np_bias_control, c(list(construction="iterated",
      iterations=control[["iterations"]], degree=control[["degree"]],
      bernstein.basis=control[["bernstein.basis"]]), control[["bandwidth.controls"]]))
    if (!setequal(names(checked[["bandwidth.controls"]]),
                  names(control[["bandwidth.controls"]])))
      stop("incomplete auxiliary bandwidth controls", call. = FALSE)
  }
  control
}

.np_bias_resolve <- function(control, xcon, profile.sizes) {
  control <- .np_bias_validate_control(control)
  if (control[["construction"]] != "iterated")
    stop("iterated bias resolver requires construction='iterated'", call. = FALSE)
  requested <- control
  # Revalidate ordinary list fields if a caller has edited a control object.
  control <- do.call(np_bias_control, c(list(construction="iterated",
    iterations=control[["iterations"]], degree=control[["degree"]],
    bernstein.basis=control[["bernstein.basis"]]), control[["bandwidth.controls"]]))
  d <- ncol(xcon); n <- nrow(xcon); n.min <- min(profile.sizes)
  spread <- setNames(as.double(EssDee(xcon)), names(xcon))
  if (!d || any(!is.finite(spread)) || any(spread <= 0))
    stop("bias correction requires positive continuous predictor spreads", call. = FALSE)
  a <- 1/(4+d)
  bandwidth <- effective <- resolved <- list(); source <- character()
  for (role in c("p0", "p1", "operator", "variance")) {
    raw <- paste0(role, ".bandwidth"); scale <- paste0(role, ".scale")
    is.raw <- raw %in% names(control[["bandwidth.controls"]])
    name <- if (is.raw) raw else scale
    value <- .np_bias_vector(control[["bandwidth.controls"]][[name]], names(xcon), name)
    resolved[[role]] <- value
    source[role] <- if (is.raw) "raw" else "profile_scale"
    bandwidth[[role]] <- if (is.raw) value else value * spread * n.min^(-a)
    effective[[role]] <- bandwidth[[role]]/(spread*n^(-a))
  }
  list(construction = "iterated", version = 1L, policy = "rot-profile-v1",
    requested = requested, resolved_request = resolved, source = source,
    columns = names(xcon), n = n, profile_sizes = profile.sizes, n_min = n.min,
    exponent = a, spread = spread, spread_convention = "EssDee on full retained training sample",
    effective_package_scale = effective, bandwidth = bandwidth,
    auxiliary_degree = control[["degree"]], iterations = control[["iterations"]],
    bernstein.basis = control[["bernstein.basis"]])
}
