.np_public_dots_cache <- new.env(parent = emptyenv())

.np_public_dots_owners <- list(
  npregbw = c("npregbw", "rbandwidth"),
  npreg = c("npreg", "npregbw", "rbandwidth"),
  npudensbw = c("npudensbw", "bandwidth"),
  npudens = c("npudens", "npudensbw", "bandwidth"),
  npcdensbw = c("npcdensbw", "conbandwidth"),
  npcdens = c("npcdens", "npcdensbw", "conbandwidth"),
  npudistbw = c("npudistbw", "dbandwidth"),
  npudist = c("npudist", "npudistbw", "dbandwidth"),
  npcdistbw = c("npcdistbw", "condbandwidth"),
  npcdist = c("npcdist", "npcdistbw", "condbandwidth"),
  npscoefbw = c("npscoefbw", "scbandwidth"),
  npscoef = c("npscoef", "npscoefbw", "scbandwidth"),
  npindexbw = c("npindexbw", "sibandwidth", "npregbw", "rbandwidth"),
  npindex = c("npindex", "npindexbw", "sibandwidth", "npregbw", "rbandwidth"),
  npplregbw = c("npplregbw", "plbandwidth", "npregbw", "rbandwidth"),
  npplreg = c("npplreg", "npplregbw", "plbandwidth", "npregbw", "rbandwidth"),
  npqreg = c("npqreg", "npcdistbw", "condbandwidth"),
  npconmode = c("npconmode", "npcdensbw", "conbandwidth"),
  npcopula = c("npcopula", "npudensbw", "bandwidth", "npudistbw", "dbandwidth")
)

.np_public_dots_extras <- list(
  npregbw = c("remin", "random.seed"),
  npreg = c("remin", "random.seed", "warn.glp.gradient",
            "bandwidth.divide", ".np_fit_progress_handoff"),
  npudensbw = c("mads.nmulti", "nomad.nmulti", "nomad.remin",
                ".beta.range.certify"),
  npudens = c("mads.nmulti", "nomad.nmulti", "nomad.remin",
              ".beta.range.certify", ".np_fit_progress_handoff"),
  npcdensbw = c("random.seed", "mads.nmulti"),
  npcdens = c("random.seed", "mads.nmulti", ".np_fit_progress_handoff",
              ".np_categorical_effects"),
  npudistbw = c("mads.nmulti", "nomad.nmulti", "nomad.remin",
                ".beta.range.certify"),
  npudist = c("mads.nmulti", "nomad.nmulti", "nomad.remin",
              ".beta.range.certify", ".np_fit_progress_handoff"),
  npcdistbw = c("random.seed", "mads.nmulti"),
  npcdist = c("random.seed", "mads.nmulti", ".np_fit_progress_handoff"),
  npscoefbw = character(),
  npscoef = c(".np_fit_progress_allow", ".np_fit_progress_handoff"),
  npindexbw = "remin",
  npindex = c("remin", "warn.glp.gradient", ".np_fit_progress_handoff",
              ".np_index_explicit_bws", ".np_lc_fixed_progress_route",
              "gradient.order", "gradient_order"),
  npplregbw = c("random.seed", "remin"),
  npplreg = c("random.seed", "remin", ".np_fit_progress_handoff"),
  npqreg = c("random.seed", "mads.nmulti", ".np_fit_progress_handoff"),
  npconmode = c("random.seed", "mads.nmulti", ".np_fit_progress_handoff",
                ".np_categorical_effects"),
  npcopula = c("mads.nmulti", "nomad.nmulti", "nomad.remin",
               ".beta.range.certify", "u.auto")
)

.np_public_dots_specialized_rejections <- list(
  npregbw = c("glp.degree", "glp.bernstein", "glp.basis"),
  npreg = c("errors", "glp.degree", "glp.bernstein", "glp.basis"),
  npcdensbw = c("cvls.i1.rescue", "cvls.quadrature.adaptive",
                "cvls.quadrature.adaptive.tol",
                "cvls.quadrature.adaptive.max.refine",
                "cvls.quadrature.adaptive.grid.hy.ratio",
                "cvls.quadrature.adaptive.floor.tol"),
  npscoefbw = c("glp.degree", "glp.bernstein", "glp.basis", "remin"),
  npscoef = c("errors", "glp.degree", "glp.bernstein", "glp.basis", "remin"),
  npindex = c("errors", "boot.num"),
  npplregbw = c("glp.degree", "glp.bernstein", "glp.basis", "remin"),
  npplreg = c("glp.degree", "glp.bernstein", "glp.basis", "remin"),
  npqreg = c("gradient.order", "gradient_order")
)

.np_public_dots_policy <- function(route) {
  if (!route %in% names(.np_public_dots_owners))
    stop(sprintf("unknown public dots route '%s'", route), call. = FALSE)

  if (exists(route, envir = .np_public_dots_cache, inherits = FALSE))
    return(get(route, envir = .np_public_dots_cache, inherits = FALSE))

  namespace <- environment(.np_public_dots_policy)
  namespace.names <- ls(namespace, all.names = TRUE)
  owner.names <- unique(unlist(lapply(.np_public_dots_owners[[route]], function(prefix) {
    namespace.names[grepl(paste0("^", prefix, "([.]|$)"), namespace.names)]
  }), use.names = FALSE))
  owner.names <- owner.names[vapply(owner.names, exists, logical(1L),
                                    envir = namespace, mode = "function",
                                    inherits = FALSE)]
  formal.names <- unique(unlist(lapply(owner.names, function(name) {
    setdiff(names(formals(get(name, envir = namespace, inherits = FALSE))), "...")
  }), use.names = FALSE))
  allowed <- sort(unique(c(
    formal.names,
    .np_public_dots_extras[[route]],
    .np_public_dots_specialized_rejections[[route]]
  )))
  policy <- list(
    allowed = allowed,
    suggest = allowed[!startsWith(allowed, ".")]
  )
  assign(route, policy, envir = .np_public_dots_cache)
  policy
}

.np_public_dots_suggestion <- function(argument, candidates) {
  explicit <- if (identical(argument, "bernstein") &&
                  "bernstein.basis" %in% candidates) {
    "bernstein.basis"
  } else if (identical(argument, "gradient.level") &&
             "level" %in% candidates) {
    "level"
  } else {
    NULL
  }
  if (!is.null(explicit))
    return(explicit)
  if (!nzchar(argument) || !length(candidates))
    return(NULL)

  distances <- as.integer(utils::adist(argument, candidates))
  best <- which(distances == min(distances))
  if (length(best) != 1L || distances[best] != 1L)
    return(NULL)

  candidate <- candidates[best]
  if (nchar(argument) > nchar(candidate) && endsWith(argument, candidate))
    return(NULL)
  candidate
}

.np_public_dots_filter_call <- function(call, route) {
  call.names <- names(call)
  if (is.null(call.names) || !length(call.names))
    return(call)
  keep <- is.na(call.names) | !nzchar(call.names) |
    call.names %in% .np_public_dots_policy(route)$allowed
  keep[[1L]] <- TRUE
  call[keep]
}

.np_public_dots_filter_args <- function(args, route) {
  arg.names <- names(args)
  if (is.null(arg.names) || !length(arg.names))
    return(args)
  keep <- is.na(arg.names) | !nzchar(arg.names) |
    arg.names %in% .np_public_dots_policy(route)$allowed
  args[keep]
}

.np_validate_public_dots <- function(dots.call, route) {
  dot.names <- names(dots.call)
  if (is.null(dot.names) || !length(dot.names))
    return(invisible(TRUE))
  dot.names <- unique(dot.names[!is.na(dot.names) & nzchar(dot.names)])
  if (!length(dot.names))
    return(invisible(TRUE))

  policy <- .np_public_dots_policy(route)
  unused <- setdiff(dot.names, policy$allowed)
  if (!length(unused))
    return(invisible(TRUE))

  if (length(unused) > 1L) {
    stop(sprintf(
      "%s(): unused arguments %s",
      route,
      paste(sprintf("'%s'", unused), collapse = ", ")
    ), call. = FALSE)
  }

  suggestion <- .np_public_dots_suggestion(unused, policy$suggest)
  message <- sprintf("%s(): unused argument '%s'", route, unused)
  if (!is.null(suggestion))
    message <- sprintf("%s; did you mean '%s'?", message, suggestion)
  stop(message, call. = FALSE)
}
