.np_nomad_bw_type <- function(template, setup) {
  type <- if (!is.null(setup$type)) setup$type else template$type
  as.character(type)[1L]
}

.np_nomad_bw_cont_index <- function(setup) {
  if (!is.null(setup$cont_idx))
    return(as.integer(setup$cont_idx))
  if (!is.null(setup$cont_flat))
    return(as.integer(setup$cont_flat))
  integer(0L)
}

.np_nomad_bw_cat_index <- function(setup) {
  if (!is.null(setup$cat_idx))
    return(as.integer(setup$cat_idx))
  if (!is.null(setup$cat_flat))
    return(as.integer(setup$cat_flat))
  integer(0L)
}

.np_nomad_bw_ncont <- function(bounds, setup = NULL) {
  if (!is.null(bounds$ncont))
    return(as.integer(bounds$ncont[1L]))
  if (!is.null(bounds$ncon))
    return(as.integer(bounds$ncon[1L]))
  if (!is.null(setup))
    return(length(.np_nomad_bw_cont_index(setup)))
  0L
}

.np_nomad_bw_ncat <- function(bounds, setup = NULL) {
  if (!is.null(bounds$ncat))
    return(as.integer(bounds$ncat[1L]))
  if (!is.null(setup))
    return(length(.np_nomad_bw_cat_index(setup)))
  max(0L, length(bounds$lower) - .np_nomad_bw_ncont(bounds, setup))
}

.np_nomad_bw_cat_scale <- function(setup) {
  if (!is.null(setup$bandwidth.scale.categorical))
    return(as.double(setup$bandwidth.scale.categorical[1L]))
  1e4
}

.np_nomad_bw_ordinary_nn_upper <- function(setup, nn.lower, ncont) {
  if (!is.null(setup$nobs) && length(setup$nobs)) {
    return(rep.int(max(nn.lower, as.integer(setup$nobs[1L]) - 1L), ncont))
  }
  rep.int(nn.lower, ncont)
}

.np_nomad_bw_coordinate_table <- function(lower,
                                          upper,
                                          bbin,
                                          setup,
                                          fixed.lower = NA_real_,
                                          nn.lower = NA_real_,
                                          where = "NOMAD bandwidth search") {
  ncont <- length(.np_nomad_bw_cont_index(setup))
  ncat <- length(.np_nomad_bw_cat_index(setup))
  n <- length(lower)
  if (!n) {
    return(data.frame(
      index = integer(0L),
      class = character(0L),
      solver_lower = numeric(0L),
      solver_upper = numeric(0L),
      bbin = integer(0L),
      decoded_scale = numeric(0L),
      storage_index = integer(0L),
      source = character(0L),
      stringsAsFactors = FALSE
    ))
  }

  class <- c(
    if (ncont > 0L) {
      if (identical(as.character(setup$type)[1L], "fixed")) {
        rep.int("continuous_fixed_scale", ncont)
      } else {
        rep.int("continuous_nn_index", ncont)
      }
    } else {
      character(0L)
    },
    rep.int("categorical_lambda", ncat)
  )
  decoded.scale <- c(
    if (ncont > 0L && !is.null(setup$cont_scale)) {
      as.double(setup$cont_scale)
    } else {
      rep.int(NA_real_, ncont)
    },
    rep.int(.np_nomad_bw_cat_scale(setup), ncat)
  )
  storage.index <- c(.np_nomad_bw_cont_index(setup), .np_nomad_bw_cat_index(setup))
  source <- c(
    if (ncont > 0L) {
      if (identical(as.character(setup$type)[1L], "fixed")) {
        rep.int(sprintf("scale.factor.search.lower=%.15g", as.double(fixed.lower)[1L]), ncont)
      } else {
        rep.int(sprintf("nn.lower=%.15g", as.double(nn.lower)[1L]), ncont)
      }
    } else {
      character(0L)
    },
    rep.int("kernel lambda upper * categorical scale", ncat)
  )

  data.frame(
    index = seq_len(n),
    class = class,
    solver_lower = as.numeric(lower),
    solver_upper = as.numeric(upper),
    bbin = as.integer(bbin),
    decoded_scale = decoded.scale,
    storage_index = storage.index,
    source = source,
    route = as.character(where)[1L],
    stringsAsFactors = FALSE
  )
}

.np_nomad_coordinate_roles <- function(bounds, degree.search = NULL) {
  table <- attr(bounds, "nomad.coordinate.table", exact = TRUE)
  roles <- if (!is.null(table) && "class" %in% names(table)) {
    as.character(table$class)
  } else {
    ncont <- .np_nomad_bw_ncont(bounds)
    ncat <- .np_nomad_bw_ncat(bounds)
    c(
      if (ncont > 0L) {
        ifelse(as.integer(bounds$bbin[seq_len(ncont)]) == 1L,
               "continuous_nn_index",
               "continuous_fixed_scale")
      } else {
        character(0L)
      },
      rep.int("categorical_lambda", ncat)
    )
  }
  if (!is.null(degree.search)) {
    roles <- c(roles, rep.int("degree", length(degree.search$lower)))
  }
  roles
}

.np_nomad_known_coordinate_roles <- function() {
  c("continuous_fixed_scale", "continuous_real", "continuous_nn_index",
    "categorical_lambda", "degree")
}

.np_nomad_validate_coordinate_roles <- function(roles,
                                                n,
                                                where = "NOMAD coordinate geometry") {
  if (is.null(roles))
    return(NULL)

  roles <- as.character(roles)
  n <- as.integer(n)[1L]
  if (length(roles) != n) {
    stop(sprintf("%s requires %d coordinate role(s), got %d",
                 where, n, length(roles)),
         call. = FALSE)
  }
  if (anyNA(roles)) {
    stop(sprintf("%s requires complete coordinate role metadata", where),
         call. = FALSE)
  }
  unknown <- setdiff(unique(roles), .np_nomad_known_coordinate_roles())
  if (length(unknown)) {
    stop(sprintf("%s has unknown coordinate role(s): %s",
                 where, paste(unknown, collapse = ", ")),
         call. = FALSE)
  }
  roles
}

.np_nomad_apply_source_geometry <- function(opts,
                                            user.opts = list(),
                                            roles,
                                            expected.length,
                                            where = "NOMAD source geometry") {
  where <- as.character(where)[1L]
  if (is.na(where) || !nzchar(where))
    where <- "NOMAD source geometry"
  if (missing(expected.length)) {
    stop(sprintf("%s requires an expected coordinate length", where),
         call. = FALSE)
  }
  expected.length <- suppressWarnings(as.numeric(expected.length)[1L])
  if (!is.finite(expected.length) ||
      expected.length < 0 ||
      expected.length != floor(expected.length)) {
    stop(sprintf("%s requires a non-negative integer coordinate length", where),
         call. = FALSE)
  }
  expected.length <- as.integer(expected.length)
  if (is.null(roles)) {
    stop(sprintf("%s requires coordinate role metadata", where),
         call. = FALSE)
  }
  roles <- .np_nomad_validate_coordinate_roles(
    roles,
    expected.length,
    where = where
  )
  n <- expected.length
  if (!n)
    return(opts)

  generated <- list(
    INITIAL_MESH_SIZE = rep.int(1, n),
    MIN_MESH_SIZE = ifelse(roles %in% c("continuous_fixed_scale", "continuous_real"),
                           sqrt(.Machine$double.eps), 1)
  )

  user.names <- names(user.opts)
  if (is.null(user.names))
    user.names <- character()
  keep <- setdiff(names(generated), user.names)
  if (length(keep))
    opts[keep] <- generated[keep]
  opts
}

.np_nomad_bw_bounds <- function(template,
                                setup,
                                fixed.lower,
                                nn.lower = 1L,
                                where = "NOMAD bandwidth search") {
  ncont <- length(.np_nomad_bw_cont_index(setup))
  ncat <- length(.np_nomad_bw_cat_index(setup))
  type <- .np_nomad_bw_type(template, setup)
  fixed.lower <- as.double(fixed.lower)[1L]
  nn.lower <- as.integer(nn.lower[1L])

  if (ncont > 0L) {
    if (identical(type, "fixed")) {
      cont.lower <- rep.int(fixed.lower, ncont)
      cont.upper <- rep.int(1e6, ncont)
      cont.bbin <- rep.int(0L, ncont)
    } else {
      nn.upper <- if (!is.null(setup$cont_extendednn_upper) &&
                      length(setup$cont_extendednn_upper) == ncont) {
        pmax(nn.lower, as.double(setup$cont_extendednn_upper))
      } else {
        .np_nomad_bw_ordinary_nn_upper(setup, nn.lower, ncont)
      }
      cont.lower <- rep.int(nn.lower, ncont)
      cont.upper <- nn.upper
      cont.bbin <- rep.int(1L, ncont)
    }
  } else {
    cont.lower <- cont.upper <- numeric(0L)
    cont.bbin <- integer(0L)
  }

  cat.scale <- .np_nomad_bw_cat_scale(setup)
  cat.lower <- rep.int(0, ncat)
  cat.upper <- as.double(setup$cat_upper) * cat.scale
  cat.bbin <- rep.int(0L, ncat)

  lower <- c(cont.lower, cat.lower)
  upper <- c(cont.upper, cat.upper)
  bbin <- c(cont.bbin, cat.bbin)
  out <- list(
    lower = lower,
    upper = upper,
    bbin = bbin,
    ncon = ncont,
    ncont = ncont,
    ncat = ncat
  )
  attr(out, "nomad.coordinate.table") <- .np_nomad_bw_coordinate_table(
    lower = lower,
    upper = upper,
    bbin = bbin,
    setup = setup,
    fixed.lower = fixed.lower,
    nn.lower = nn.lower,
    where = where
  )
  out
}

.np_nomad_bw_point_to_storage <- function(point,
                                          template,
                                          setup,
                                          storage.length,
                                          clamp.nn = FALSE) {
  point <- as.numeric(point)
  cont.index <- .np_nomad_bw_cont_index(setup)
  cat.index <- .np_nomad_bw_cat_index(setup)
  ncont <- length(cont.index)
  ncat <- length(cat.index)
  bws <- numeric(storage.length)
  type <- .np_nomad_bw_type(template, setup)

  if (ncont > 0L) {
    cont.point <- point[seq_len(ncont)]
    if (identical(type, "fixed")) {
      ext.bw <- cont.point * setup$cont_scale
      bws[cont.index] <- if (isTRUE(template$scaling)) cont.point else ext.bw
    } else {
      nn.bw <- round(cont.point)
      if (isTRUE(clamp.nn)) {
        nn.upper <- if (!is.null(setup$cont_extendednn_upper) &&
                        length(setup$cont_extendednn_upper) == ncont) {
          pmax(1, as.double(setup$cont_extendednn_upper))
        } else if (!is.null(setup$nobs) && length(setup$nobs)) {
          rep.int(max(1L, as.integer(setup$nobs[1L]) - 1L), ncont)
        } else {
          rep.int(Inf, ncont)
        }
        nn.bw[!is.finite(nn.bw)] <- 1
        nn.bw <- pmax(1, pmin(nn.upper, nn.bw))
      }
      bws[cont.index] <- nn.bw
    }
  }

  if (ncat > 0L) {
    lambda.scaled <- point[ncont + seq_len(ncat)]
    ext.bw <- lambda.scaled / .np_nomad_bw_cat_scale(setup)
    bws[cat.index] <- if (isTRUE(template$scaling)) ext.bw / setup$ncatfac else ext.bw
  }

  bws
}

.np_nomad_bw_storage_to_point <- function(bws, template, setup) {
  cont.index <- .np_nomad_bw_cont_index(setup)
  cat.index <- .np_nomad_bw_cat_index(setup)
  point <- numeric(length(cont.index) + length(cat.index))
  type <- .np_nomad_bw_type(template, setup)

  if (length(cont.index) > 0L) {
    raw <- bws[cont.index]
    point[seq_along(cont.index)] <- if (identical(type, "fixed")) {
      if (isTRUE(template$scaling)) raw else raw / setup$cont_scale
    } else {
      round(raw)
    }
  }

  if (length(cat.index) > 0L) {
    raw <- bws[cat.index]
    ext.bw <- if (isTRUE(template$scaling)) raw * setup$ncatfac else raw
    point[length(cont.index) + seq_along(cat.index)] <-
      ext.bw * .np_nomad_bw_cat_scale(setup)
  }

  point
}

.np_nomad_bw_complete_start_point <- function(point,
                                              bounds,
                                              template,
                                              setup = NULL,
                                              initial = NULL,
                                              where = "NOMAD bandwidth search") {
  point <- .np_nomad_explicit_or_initial_start(
    point = point,
    initial = initial,
    n = length(bounds$lower),
    where = where
  )

  type <- .np_nomad_bw_type(template, setup)
  ncont <- .np_nomad_bw_ncont(bounds, setup)
  ncat <- .np_nomad_bw_ncat(bounds, setup)
  if (identical(type, "fixed")) {
    return(.np_nomad_complete_start_point(
      point = point,
      lower = bounds$lower,
      upper = bounds$upper,
      ncont = ncont
    ))
  }

  n <- length(bounds$lower)
  out <- rep(NA_real_, n)
  if (!is.null(point)) {
    point <- as.numeric(point)
    if (length(point) == n)
      out <- point
  }

  if (ncont > 0L) {
    cont.idx <- seq_len(ncont)
    cont.default <- sqrt(bounds$upper[cont.idx])
    cont.default <- pmax(bounds$lower[cont.idx], pmin(bounds$upper[cont.idx], cont.default))
    bad <- !is.finite(out[cont.idx]) |
      out[cont.idx] < bounds$lower[cont.idx] |
      out[cont.idx] > bounds$upper[cont.idx]
    out[cont.idx][bad] <- cont.default[bad]
  }

  if (ncat > 0L) {
    cat.idx <- ncont + seq_len(ncat)
    cat.default <- pmax(bounds$lower[cat.idx], pmin(bounds$upper[cat.idx], 0.5 * bounds$upper[cat.idx]))
    bad <- !is.finite(out[cat.idx]) |
      out[cat.idx] < bounds$lower[cat.idx] |
      out[cat.idx] > bounds$upper[cat.idx]
    out[cat.idx][bad] <- cat.default[bad]
  }

  vapply(
    seq_len(n),
    function(i) .np_nomad_coerce_start_value(
      out[i],
      type = bounds$bbin[i],
      lb = bounds$lower[i],
      ub = bounds$upper[i]
    ),
    numeric(1L)
  )
}
