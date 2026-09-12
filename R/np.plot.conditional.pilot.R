# Conditional smooth-bootstrap donor mixture. Used only for explicit bias centers.
.np_plot_pilot_continuous <- function(z, kernel, derivative = 0L, cdf = FALSE) {
  if (cdf) {
    if (kernel == "gaussian") {
      return(stats::pnorm(z))
    }
    a <- if (kernel == "epanechnikov") sqrt(5) else 1
    u <- pmax(-1, pmin(1, z / a))
    return(if (kernel == "uniform") {
      (u + 1) / 2
    } else {
      0.5 + 0.75 * u - 0.25 * u^3
    })
  }
  if (kernel == "gaussian") {
    h <- rep.int(1, length(z))
    if (derivative > 0L) {
      previous <- h
      h <- z
      if (derivative > 1L) {
        for (k in 2:derivative) {
          next.h <- z * h - (k - 1) * previous
          previous <- h
          h <- next.h
        }
      }
    }
    return((-1)^derivative * h * stats::dnorm(z))
  }
  a <- if (kernel == "epanechnikov") sqrt(5) else 1
  inside <- abs(z) < a
  value <- numeric(length(z))
  if (kernel == "uniform") {
    if (derivative == 0L) value[inside] <- 0.5
  } else if (derivative == 0L) {
    value[inside] <- 3 / (4 * a) * (1 - z[inside]^2 / 5)
  } else if (derivative == 1L) {
    value[inside] <- -3 * z[inside] / (10 * a)
  } else if (derivative == 2L) {
    value[inside] <- -3 / (10 * a)
  }
  if (derivative > 0L) value[abs(z) == a] <- NA_real_
  value
}

.np_plot_pilot_quantile <- function(p, kernel) {
  if (kernel == "gaussian") {
    return(stats::qnorm(p))
  }
  if (kernel == "uniform") {
    return(2 * p - 1)
  }
  2 * sqrt(5) * sin(asin(2 * p - 1) / 3)
}

.np_plot_pilot_cat_raw <- function(spec, eval, donor) {
  # WVR and normalized LR have a removable common factor at lambda=1.
  # Use the normalized law's exact continuous extension at that endpoint.
  if (is.ordered(spec$data[[1L]]) && spec$lambda == 1 &&
    spec$okernel %in% c("wangvanryzin", "nliracine")) {
    if (spec$okernel == "wangvanryzin") {
      return(ifelse(outer(eval, donor, "=="), 1, 0.5))
    }
    return(matrix(1, length(eval), length(donor)))
  }
  .np_regression_cat_profile_kernel_matrix(
    matrix(eval, ncol = 1L), matrix(donor, ncol = 1L), spec$data,
    list(bw = spec$lambda, ukertype = spec$ukernel, okertype = spec$okernel)
  )
}

.np_plot_pilot_prepare_side <- function(dat, bw, icon, kernel, order,
                                        lb, ub, ukernel, okernel,
                                        smooth.categories) {
  dat <- toFrame(dat)
  codes <- .np_cat_profile_code_matrix(dat)
  specs <- vector("list", ncol(dat))
  for (j in seq_len(ncol(dat))) {
    if (icon[j]) {
      family <- .np_plot_kernel_perturbation_policy(kernel, order)$family
      lower <- if (length(lb)) lb[j] else -Inf
      upper <- if (length(ub)) ub[j] else Inf
      if (is.na(lower)) lower <- -Inf
      if (is.na(upper)) upper <- Inf
      left <- .np_plot_pilot_continuous((lower - codes[, j]) / bw[j], family, cdf = TRUE)
      right <- .np_plot_pilot_continuous((upper - codes[, j]) / bw[j], family, cdf = TRUE)
      mass <- right - left
      if (any(!is.finite(mass) | mass <= 0)) {
        stop("smooth-bootstrap pilot has no representable probability inside the supplied bounds", call. = FALSE)
      }
      specs[[j]] <- list(
        continuous = TRUE, bw = bw[j], kernel = family,
        lower = lower, upper = upper, left = left, mass = mass,
        bounded = is.finite(lower) || is.finite(upper)
      )
    } else {
      support <- .np_cat_profile_code_matrix(data.frame(
        factor(levels(dat[[j]]),
          levels = levels(dat[[j]]),
          ordered = is.ordered(dat[[j]])
        )
      ))[, 1L]
      spec <- list(
        continuous = FALSE, data = dat[j], support = support,
        lambda = bw[j], ukernel = ukernel, okernel = okernel,
        smooth = smooth.categories
      )
      if (smooth.categories) {
        # One category vector at a time, never a C-by-C transition matrix.
        spec$mass <- vapply(support, function(donor) {
          sum(.np_plot_pilot_cat_raw(spec, support, donor))
        }, numeric(1L))
        if (any(!is.finite(spec$mass) | spec$mass <= 0)) {
          stop("categorical smooth-bootstrap pilot has invalid transition probabilities", call. = FALSE)
        }
      }
      specs[[j]] <- spec
    }
  }
  list(data = dat, codes = codes, specs = specs)
}

.np_plot_conditional_pilot_prepare <- function(xdat, ydat, bws, cdf) {
  pilot <- .np_plot_oversmooth_conditional_bws(bws, cdf = cdf)
  prepare <- function(dat, side, smooth) {
    .np_plot_pilot_prepare_side(
      dat, pilot$bandwidth[[side]],
      pilot[[paste0(side, "dati")]]$icon,
      pilot[[paste0("c", side, "kertype")]], pilot[[paste0("c", side, "kerorder")]],
      pilot[[paste0("c", side, "kerlb")]], pilot[[paste0("c", side, "kerub")]],
      pilot[[paste0("u", side, "kertype")]], pilot[[paste0("o", side, "kertype")]],
      smooth
    )
  }
  list(x = prepare(xdat, "x", TRUE), y = prepare(ydat, "y", FALSE), cdf = cdf)
}

.np_plot_pilot_draw_side <- function(side, idx) {
  out <- side$data[idx, , drop = FALSE]
  row.names(out) <- NULL
  for (j in seq_along(side$specs)) {
    spec <- side$specs[[j]]
    if (spec$continuous) {
      noise <- if (spec$bounded) {
        .np_plot_pilot_quantile(
          spec$left[idx] + stats::runif(length(idx)) * spec$mass[idx], spec$kernel
        )
      } else {
        .np_plot_kernel_random(length(idx), spec$kernel, 2L)
      }
      out[[j]] <- side$codes[idx, j] + spec$bw * noise
    } else if (spec$smooth) {
      donor.codes <- side$codes[idx, j]
      result <- integer(length(idx))
      for (donor in unique(donor.codes)) {
        at <- which(donor.codes == donor)
        probs <- as.vector(.np_plot_pilot_cat_raw(spec, spec$support, donor))
        result[at] <- sample.int(length(spec$support), length(at),
          replace = TRUE, prob = probs
        )
      }
      out[[j]] <- factor(levels(side$data[[j]])[result],
        levels = levels(side$data[[j]]), ordered = is.ordered(side$data[[j]])
      )
    }
  }
  out
}

.np_plot_conditional_pilot_boot <- function(xdat, ydat, exdat, eydat, bws, cdf,
                                            plot.errors.boot.method,
                                            plot.errors.boot.blocklen,
                                            plot.errors.boot.num,
                                            progress.label,
                                            gradient.index = NULL,
                                            gradient.order = 1L) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  B <- as.integer(plot.errors.boot.num)
  n <- nrow(xdat)
  neval <- nrow(exdat)
  gradients <- !is.null(gradient.index)
  where <- "conditional smooth-bootstrap pilot"
  if (nrow(ydat) != n || nrow(eydat) != neval) {
    stop("conditional smooth-bootstrap helper requires aligned x/y training and evaluation rows")
  }
  if (n < 1L || neval < 1L || B < 1L) {
    stop("invalid conditional smooth-bootstrap dimensions")
  }
  if (gradients) {
    gradient.index <- .np_plot_resolve_conditional_gradient_index(
      bws, gradient.index, where
    )
  }
  pilot <- .np_plot_conditional_pilot_prepare(xdat, ydat, bws, cdf)
  fit_one <- function(x.train, y.train) {
    fit <- .np_plot_conditional_eval(
      bws = bws, xdat = x.train, ydat = y.train,
      exdat = exdat, eydat = eydat, cdf = cdf, gradients = gradients,
      gradient.order = gradient.order, lp.first.se.demand = FALSE,
      cat.se.demand = FALSE, se = FALSE, gradient.target = gradient.index
    )
    if (gradients) {
      .np_plot_extract_conditional_gradient(
        fit, gradient.index,
        neval, where
      )
    } else {
      as.vector(fit[[if (cdf) "condist" else "condens"]])
    }
  }
  t0 <- fit_one(xdat, ydat)
  reference.order <- 1L
  if (gradients && bws$xdati$icon[gradient.index]) {
    spec <- npConditionalRegEngineSpec(bws, where = where)
    orders <- npConditionalGradientOrder(bws, spec$reg.engine, gradient.order, where)
    reference.order <- orders[match(gradient.index, which(bws$xdati$icon))]
  }
  reference <- .np_plot_conditional_pilot_reference(pilot, exdat, eydat,
    gradient.index = gradient.index, gradient.order = reference.order
  )
  tmat <- matrix(NA_real_, nrow = B, ncol = length(t0))
  is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
  index.drawer <- if (is.block) {
    .np_block_indices_drawer(
      n = n, B = B, blocklen = plot.errors.boot.blocklen, sim = plot.errors.boot.method
    )
  } else {
    NULL
  }
  smooth_one <- function(idx) {
    xstar <- .np_plot_pilot_draw_side(pilot$x, idx)
    ystar <- .np_plot_pilot_draw_side(pilot$y, idx)
    fit_one(xstar, ystar)
  }
  chunk.size <- .np_inid_chunk_size(n = n, B = B, progress_cap = is.block)
  progress <- .np_plot_bootstrap_progress_begin(
    total = B,
    label = if (is.null(progress.label)) {
      if (gradients) "Plot bootstrap gradient smooth" else "Plot bootstrap smooth"
    } else {
      progress.label
    }
  )
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  chunk.controller <- .np_plot_progress_chunk_controller(chunk.size = chunk.size, progress = progress)
  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + chunk.controller$chunk.size - 1L)
    bsz <- stopi - start + 1L
    chunk.started <- .np_progress_now()
    idx.chunk <- if (!is.null(index.drawer)) {
      index.drawer(start, stopi)
    } else {
      matrix(sample.int(n = n, size = n * bsz, replace = TRUE), nrow = n)
    }
    for (jj in seq_len(bsz)) {
      tmat[start + jj - 1L, ] <- smooth_one(idx.chunk[, jj])
    }
    progress <- .np_plot_progress_tick(state = progress, done = stopi)
    chunk.controller <- .np_plot_progress_chunk_observe(
      controller = chunk.controller,
      bsz = bsz, elapsed.sec = .np_progress_now() - chunk.started
    )
    start <- stopi + 1L
  }
  list(t = tmat, t0 = t0, center = reference)
}

.np_plot_pilot_side_row <- function(side, eval, cdf = FALSE,
                                    derivative.index = NULL, derivative = 0L) {
  value <- rep.int(1, nrow(side$data))
  for (j in seq_along(side$specs)) {
    spec <- side$specs[[j]]
    if (spec$continuous) {
      m <- if (identical(j, derivative.index)) derivative else 0L
      z <- (eval[j] - side$codes[, j]) / spec$bw
      part <- .np_plot_pilot_continuous(z, spec$kernel, m, cdf)
      if (cdf) {
        part <- (part - spec$left) / spec$mass
        part <- pmax(0, pmin(1, part))
      } else {
        part <- part / (spec$bw^(m + 1L) * spec$mass)
        if (eval[j] < spec$lower || eval[j] > spec$upper) part[] <- 0
      }
    } else if (spec$smooth) {
      part <- as.vector(.np_plot_pilot_cat_raw(spec, eval[j], side$codes[, j])) /
        spec$mass[match(side$codes[, j], spec$support)]
    } else {
      part <- if (cdf) {
        as.numeric(side$codes[, j] <= eval[j])
      } else {
        as.numeric(side$codes[, j] == eval[j])
      }
    }
    value <- value * part
  }
  value
}

.np_plot_conditional_pilot_reference <- function(pilot, exdat, eydat,
                                                 gradient.index = NULL,
                                                 gradient.order = 1L) {
  if (!is.null(gradient.index) &&
    !pilot$x$specs[[gradient.index]]$continuous) {
    pair <- npCategoricalFirstDifferenceFrames(
      exdat, gradient.index,
      "conditional smooth-bootstrap pilot reference"
    )
    return(.np_plot_conditional_pilot_reference(pilot, pair$upper, eydat) -
      .np_plot_conditional_pilot_reference(pilot, pair$lower, eydat))
  }
  x <- .np_cat_profile_code_matrix(exdat)
  y <- .np_cat_profile_code_matrix(eydat)
  m <- if (is.null(gradient.index)) 0L else as.integer(gradient.order)
  out <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) {
    b <- .np_plot_pilot_side_row(pilot$y, y[i, ], cdf = pilot$cdf)
    den <- num <- truth <- numeric(m + 1L)
    for (k in 0:m) {
      a <- .np_plot_pilot_side_row(pilot$x, x[i, ],
        derivative.index = gradient.index, derivative = k
      )
      den[k + 1L] <- sum(a)
      num[k + 1L] <- sum(a * b)
      truth[k + 1L] <- if (k == 0L) {
        num[1L] / den[1L]
      } else {
        (num[k + 1L] - sum(choose(k, 1:k) * den[2:(k + 1L)] * rev(truth[1:k]))) / den[1L]
      }
    }
    out[i] <- truth[m + 1L]
  }
  if (any(!is.finite(out))) {
    stop("conditional smooth-bootstrap pilot reference is undefined at an evaluation row (no support or a nonsmooth kernel boundary)", call. = FALSE)
  }
  out
}
