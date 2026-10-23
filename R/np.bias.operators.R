# Bounded operators for bias-corrected inference. Regression supplies the
# polynomial/profile adapter; composition and residual application are generic.
# Extracted from the admitted Gate26 arithmetic; no observation-square storage.

.np_bias_solve <-
function (W, k, v)
{
    out <- .npreghat_solve_eval(W, v, k, ridge.base = 0, canonical.lp = TRUE, allow.empty.rows = FALSE)
    if (is.null(out) || !is.finite(out$ridge) || out$ridge != 0)
        stop("UNQUALIFIED: required original/pilot solve uses ridge")
    row <- k * drop(W %*% out$v)
    if (any(!is.finite(row)))
        stop("UNQUALIFIED: nonfinite required row")
    row
}

.np_bias_key <-
function (spec, e)
{
    if (!length(spec$factor.cols))
        return("all")
    paste(vapply(e[spec$factor.cols], as.integer, integer(1)), collapse = ":")
}

.np_bias_projection <-
function (spec, e, degree, nu)
{
    drop(W.lp(spec$xcon, exdat = e[spec$numeric.cols], degree = rep.int(degree, spec$d), gradient.vec = if (sum(nu))
        nu
    else NULL, basis = "glp", bernstein.basis = spec$bernstein))
}

.np_bias_projector <-
function (spec, anchor, y, components = FALSE)
{
    n <- nrow(spec$x)
    stopifnot(length(y) == n, nrow(anchor) == 1L, all(is.finite(y)))
    k <- drop(.np_kernel_weights_direct(spec$pilot.bws, txdat = spec$xcon, exdat = anchor[spec$numeric.cols],
        bandwidth.divide = TRUE))
    if (any(!is.finite(k)) || any(k < 0))
        stop("UNQUALIFIED: HC3 requires finite nonnegative pilot weights")
    residual <- leverage <- numeric(n)
    blocks <- if (components)
        vector("list", length(spec$groups))
    else NULL
    block.index <- 0L
    for (ii in spec$groups) {
        Q <- spec$Q[ii, , drop = FALSE]
        weight <- k[ii]
        G <- crossprod(Q, Q * weight)
        width <- ncol(Q)
        rank.bound <- min(width, sum(weight != 0))
        solve.canonical <- function(rhs) {
            z <- .Call("C_np_lp_adjoint_prepared_ridge", G, as.double(rhs), as.integer(length(ii)), as.integer(rank.bound),
                PACKAGE = "np")
            if (!is.finite(z$ridge) || z$ridge != 0 || any(!is.finite(z$v)))
                stop("UNQUALIFIED: HC3 required canonical solve uses ridge")
            z$v
        }
        beta <- solve.canonical(crossprod(Q, weight * y[ii]))
        inverse <- matrix(0, width, width)
        for (j in seq_len(width)) {
            rhs <- numeric(width)
            rhs[j] <- 1
            inverse[, j] <- solve.canonical(rhs)
        }
        diagonal <- weight * rowSums((Q %*% inverse) * Q)
        denominator <- 1 - diagonal
        if (any(!is.finite(denominator)) || any(denominator <= 0))
            stop("UNQUALIFIED: nonpositive HC3 leverage denominator")
        residual[ii] <- (y[ii] - drop(Q %*% beta))/denominator
        leverage[ii] <- diagonal
        if (components) {
            block.index <- block.index + 1L
            blocks[[block.index]] <- list(indices = ii, Q = Q, coefficient = inverse %*% t(Q * weight),
                denominator = denominator)
        }
    }
    list(residual = residual, leverage = leverage, blocks = blocks)
}

.np_bias_projectors <-
function (spec, grid, map, reference = FALSE)
{
    anchors <- unique(grid[spec$numeric.cols])
    lapply(seq_len(nrow(anchors)), function(a) {
        match <- Reduce(`&`, lapply(seq_len(spec$d), function(j) grid[[spec$numeric.cols[j]]] == anchors[[j]][a]))
        list(rows = which(match[map$row]), e = grid[which(match)[1L], , drop = FALSE], blocks = .np_bias_projector_blocks(spec,
            grid[which(match)[1L], , drop = FALSE], reference))
    })
}

.np_bias_right_residual <-
function (W, projectors)
{
    for (a in projectors) for (b in a$blocks) {
        z <- W[a$rows, b$indices, drop = FALSE]
        W[a$rows, b$indices] <- z - (z %*% b$Q) %*% b$coefficient
    }
    W
}

.np_bias_kernel <-
function (x, z, h)
{
    stopifnot(nrow(z) <= 16L)
    out <- matrix(0, nrow(x), nrow(z))
    for (j in seq_along(h)) out <- out + outer(x[, j], z[, j], "-")^2/h[j]^2
    exp(-out/2)
}

.np_bias_stable_prepare <-
function (spec, grid, map, block = 16L)
{
    stopifnot(block >= 1L, block <= 16L, spec$q >= 1L, length(spec$b) == spec$d, length(spec$h) == spec$d,
        all(is.finite(spec$b)), all(spec$b > 0), all(is.finite(spec$h)), all(spec$h > 0), nrow(map) <=
            500L, ncol(spec$Q) < min(lengths(spec$groups)))
    projectors <- .np_bias_projectors(spec, grid, map)
    stopifnot(length(projectors) <= 25L)
    x <- as.matrix(spec$xcon)
    for (a in seq_along(projectors)) {
        g <- drop(.np_kernel_weights_direct(spec$pilot.bws, txdat = spec$xcon, exdat = projectors[[a]]$e[spec$numeric.cols],
            bandwidth.divide = TRUE))
        if (any(!is.finite(g)) || any(g <= 0))
            stop("UNQUALIFIED: nonpositive stability metric")
        for (b in seq_along(projectors[[a]]$blocks)) {
            p <- projectors[[a]]$blocks[[b]]
            ii <- p$indices
            qrfit <- qr(p$Q * sqrt(g[ii]))
            if (qrfit$rank != ncol(p$Q))
                stop("UNQUALIFIED: metric QR rank")
            p$g <- g[ii]
            p$thin <- qr.Q(qrfit) * sqrt(g[ii])
            p$projected.trace <- 0
            projectors[[a]]$blocks[[b]] <- p
        }
    }
    for (b in seq_along(spec$groups)) {
        ii <- spec$groups[[b]]
        width <- ncol(spec$Q)
        thin <- do.call(cbind, lapply(projectors, function(a) a$blocks[[b]]$thin))
        tr <- numeric(ncol(thin))
        for (first in seq.int(1L, length(ii), by = block)) {
            jj <- seq.int(first, min(length(ii), first + block - 1L))
            K <- .np_bias_kernel(x[ii, , drop = FALSE], x[ii[jj], , drop = FALSE], spec$h)
            tr <- tr + colSums(thin[jj, , drop = FALSE] * crossprod(K, thin))
        }
        for (a in seq_along(projectors)) {
            p <- projectors[[a]]$blocks[[b]]
            p$projected.trace <- sum(tr[(a - 1L) * width + seq_len(width)])
            p$normalizer <- sum(p$g) - p$projected.trace
            if (!is.finite(p$normalizer) || p$normalizer <= 1e-08 * sum(p$g))
                stop("UNQUALIFIED: unresolved positive residual trace")
            p$thin <- NULL
            projectors[[a]]$blocks[[b]] <- p
        }
    }
    list(projectors = projectors, x = x, h = spec$h, groups = spec$groups)
}

.np_bias_stable_product <-
function (W, op, block = 16L)
{
    stopifnot(nrow(W) <= 500L, ncol(W) == nrow(op$x), block >= 1L, block <= 16L)
    out <- matrix(0, nrow(W), ncol(W))
    for (ii in op$groups) for (first in seq.int(1L, length(ii), by = block)) {
        jj <- seq.int(first, min(length(ii), first + block - 1L))
        K <- .np_bias_kernel(op$x[ii, , drop = FALSE], op$x[ii[jj], , drop = FALSE], op$h)
        for (start in seq.int(1L, nrow(W), by = 128L)) {
            rr <- seq.int(start, min(nrow(W), start + 127L))
            out[rr, ii[jj]] <- W[rr, ii, drop = FALSE] %*% K
        }
    }
    for (a in op$projectors) for (b in a$blocks) out[a$rows, b$indices] <- sweep(out[a$rows, b$indices,
        drop = FALSE], 2L, b$g/b$normalizer, "*")
    out
}

.np_bias_stable_target <-
function (spec, grid, map, op)
{
    out <- matrix(0, nrow(map), nrow(spec$x))
    for (a in op$projectors) for (b in a$blocks) for (r in a$rows) {
        f <- .np_bias_frames(spec, grid[map$row[r], , drop = FALSE], map$target[r])
        for (k in seq_along(f$e)) if (.np_bias_key(spec, f$e[[k]]) == spec$keys[b$indices[1L]]) {
            e <- as.matrix(f$e[[k]][spec$numeric.cols])
            ii <- b$indices
            v <- drop(.np_bias_kernel(op$x[ii, , drop = FALSE], e, op$h))
            if (f$nu)
                v <- v * (op$x[ii, f$nu] - e[1L, f$nu])/op$h[f$nu]^2
            out[r, ii] <- out[r, ii] + f$sign[k] * v * b$g/b$normalizer
        }
    }
    out
}

.np_bias_pilot_rows <-
function (spec, at, derivatives = TRUE)
{
    ii <- spec$groups[[.np_bias_key(spec, at)]]
    if (is.null(ii))
        stop("missing pilot profile")
    x <- spec$xcon[ii, , drop = FALSE]
    u <- x
    for (d in seq_len(spec$d)) u[[d]] <- (x[[d]] - at[[spec$numeric.cols[d]]])/spec$b[d]
    zero <- u[1, , drop = FALSE]
    zero[, ] <- 0
    basis <- W.lp
    Q <- basis(u, degree = rep.int(spec$q, spec$d), basis = "glp", bernstein.basis = spec$bernstein)
    v <- drop(basis(u, exdat = zero, degree = rep.int(spec$q, spec$d), basis = "glp", bernstein.basis = spec$bernstein))
    k <- drop(.np_kernel_weights_direct(spec$pilot.bws, txdat = x, exdat = at[spec$numeric.cols], bandwidth.divide = TRUE))
    h <- .np_bias_solve(Q, k, v)
    dh <- list()
    if (derivatives)
        for (d in seq_len(spec$d)) {
            nu <- integer(spec$d)
            nu[d] <- 1L
            dv <- drop(basis(u, exdat = zero, degree = rep.int(spec$q, spec$d), gradient.vec = nu, basis = "glp",
                bernstein.basis = spec$bernstein))/spec$b[d]
            coefficient <- .np_bias_solve(Q, k, dv)
            w <- h * u[[d]]/spec$b[d]
            dh[[d]] <- coefficient + w - .np_bias_solve(Q, k, drop(crossprod(Q, w)))
        }
    list(indices = ii, mean = h, derivatives = dh)
}

.np_bias_right_compose <-
function (W, row_owner, block = 16L)
{
    if (!is.matrix(W) || !is.numeric(W) || nrow(W) < 1L || nrow(W) > 500L || any(!is.finite(W)) || !is.function(row_owner) ||
        length(block) != 1L || !is.finite(block) || block != floor(block) || block < 1L || block > 16L)
        stop("invalid streamed composition")
    n <- ncol(W)
    out <- matrix(0, nrow(W), n)
    for (first in seq.int(1L, n, by = block)) {
        ii <- seq.int(first, min(n, first + block - 1L))
        P <- row_owner(ii)
        if (!is.matrix(P) || !identical(dim(P), c(length(ii), n)) || any(!is.finite(P)))
            stop("invalid streamed row owner")
        out <- out + W[, ii, drop = FALSE] %*% P
    }
    if (any(!is.finite(out)))
        stop("nonfinite streamed composition")
    out
}

.np_bias_pilot_owner <-
function (spec)
{
    force(spec)
    function(ii) {
        if (length(ii) > 16L || anyNA(ii) || any(ii < 1L | ii > nrow(spec$x)) || any(ii != floor(ii)))
            stop("invalid pilot row request")
        ans <- matrix(0, length(ii), nrow(spec$x))
        for (j in seq_along(ii)) {
            row <- .np_bias_pilot_rows(spec, spec$x[ii[j], , drop = FALSE], FALSE)
            ans[j, row$indices] <- row$mean
        }
        ans
    }
}

.np_bias_residual_apply <-
function (owner, responses, block = 16L)
{
    if (!is.function(owner) || !is.matrix(responses) || !is.numeric(responses) || nrow(responses) < 2L ||
        ncol(responses) < 1L || ncol(responses) > 9999L || any(!is.finite(responses)) || length(block) !=
        1L || !is.finite(block) || block != floor(block) || block < 1L || block > 16L)
        stop("invalid residual application")
    n <- nrow(responses)
    raw <- matrix(0, n, ncol(responses))
    leverage <- row_variance <- numeric(n)
    for (first in seq.int(1L, n, by = block)) {
        ii <- seq.int(first, min(n, first + block - 1L))
        H <- owner(ii)
        if (!is.matrix(H) || !identical(dim(H), c(length(ii), n)) || any(!is.finite(H)))
            stop("invalid residual row owner")
        raw[ii, ] <- responses[ii, , drop = FALSE] - H %*% responses
        leverage[ii] <- H[cbind(seq_along(ii), ii)]
        H[cbind(seq_along(ii), ii)] <- H[cbind(seq_along(ii), ii)] - 1
        row_variance[ii] <- rowSums(H^2)
    }
    if (any(!is.finite(row_variance)) || any(row_variance <= 0))
        stop("nonpositive residual row variance")
    list(raw = raw, residual = raw/sqrt(row_variance), leverage = leverage, row_variance = row_variance)
}

.np_bias_p0_residuals <-
function (spec, responses, block = 16L)
{
    if (!is.matrix(responses) || !is.numeric(responses) || nrow(responses) != nrow(spec$x) || ncol(responses) <
        1L || ncol(responses) > 9999L || any(!is.finite(responses)))
        stop("invalid P0 response matrix")
    out <- matrix(0, nrow(responses), ncol(responses))
    raw <- out
    leverage <- row_variance <- numeric(nrow(responses))
    for (ii in spec$groups) {
        x <- spec$xcon[ii, , drop = FALSE]
        owner <- function(jj) {
            H <- npreghat(spec$pilot.bws, txdat = x, exdat = x[jj, , drop = FALSE], output = "matrix",
                .np.require.finite = TRUE)
            ridge <- attr(H, "ridge.used", exact = TRUE)
            if (is.null(ridge) || any(!is.finite(ridge)) || any(ridge != 0))
                stop("unqualified P0 residual ridge")
            H
        }
        v <- .np_bias_residual_apply(owner, responses[ii, , drop = FALSE], block)
        out[ii, ] <- v$residual
        raw[ii, ] <- v$raw
        leverage[ii] <- v$leverage
        row_variance[ii] <- v$row_variance
    }
    list(residual = out, raw = raw, leverage = leverage, row_variance = row_variance)
}

