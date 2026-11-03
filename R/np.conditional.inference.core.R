# Experimental LC conditional inference. Frozen CDF/LCD5/D numerical owners.
# Internal names are isolated; arithmetic and finite-rank rules retain the private references.
# Public validation and bounded top-level orchestration live in np.conditional.inference.R.

.np_conditional_local_gh <- local({
    cache <- new.env(parent = emptyenv())
    function(k) {
        key <- as.character(k)
        if (exists(key, cache, inherits = FALSE))
            return(cache[[key]])
        J <- matrix(0, k, k)
        J[cbind(1:(k - 1), 2:k)] <- sqrt(1:(k - 1))
        J <- J + t(J)
        e <- eigen(J, symmetric = TRUE)
        o <- order(e$values)
        ans <- list(x = e$values[o], w = e$vectors[1, o]^2)
        ans$w <- ans$w/sum(ans$w)
        cache[[key]] <- ans
        ans
    }
})

.np_conditional_context <- function(x, y, grid, h) {
    gx <- sort(unique(grid$x))
    gy <- sort(unique(grid$y))
    nx <- length(gx)
    ny <- length(gy)
    list(x = x, y = y, grid = grid, h = h, gx = gx, gy = gy, nx = nx, ny = ny, ix = match(grid$x, gx),
        iy = match(grid$y, gy), ref = as.numeric(gy > median(y)), n = length(x), m = nrow(grid))
}

.np_conditional_xfeatures <- function(x, z) {
    u <- outer(x, z$gx, function(a, b) (b - a)/z$h[1])
    K <- dnorm(u)/z$h[1]
    cbind(K, K * (u^2 - 1)/2, dnorm(u/sqrt(2))/(sqrt(2) * z$h[1]))
}

.np_conditional_yfeatures <- function(y, z) {
    v <- outer(y, z$gy, function(a, b) (b - a)/z$h[2])
    g <- pnorm(v)
    gg <- pnorm(v/sqrt(2))
    i <- which(z$ref == 1)
    if (length(i)) {
        g[, i] <- -pnorm(v[, i, drop = FALSE], lower.tail = FALSE)
        gg[, i] <- -pnorm(v[, i, drop = FALSE]/sqrt(2), lower.tail = FALSE)
    }
    cbind(g, v * dnorm(v), gg, 1)
}

.np_conditional_params <- function(ex, ey, z, prob = rep(1/nrow(ex), nrow(ex))) {
    a <- z$ix
    b <- z$iy
    fx <- colSums(ex[, a, drop = FALSE] * prob)
    aa <- colSums(ex[, a + z$nx, drop = FALSE] * prob)/fx
    F <- colSums(ex[, a, drop = FALSE] * ey[, b, drop = FALSE] * prob)/fx
    T <- (1 + aa) * F - colSums(ex[, a + z$nx, drop = FALSE] * ey[, b, drop = FALSE] * prob)/fx + 0.5 *
        colSums(ex[, a, drop = FALSE] * ey[, b + z$ny, drop = FALSE] * prob)/fx
    fl <- colSums(ex[, a + 2 * z$nx, drop = FALSE] * prob)
    Fl <- colSums(ex[, a + 2 * z$nx, drop = FALSE] * ey[, b + 2 * z$ny, drop = FALSE] * prob)/fl
    list(fx = fx, aa = aa, F = F, T = T, fl = fl, Fl = Fl)
}

.np_conditional_moments <- function(z, tau = 1, order = 64) {
    q <- .np_conditional_local_gh(order)
    ex <- matrix(0, z$n, 3 * z$nx)
    ey <- matrix(0, z$n, 3 * z$ny + 1)
    cx <- cy <- vector("list", z$n)
    for (i in seq_len(z$n)) {
        X <- .np_conditional_xfeatures(z$x[i] + tau * z$h[1] * q$x, z)
        Y <- .np_conditional_yfeatures(z$y[i] + tau * z$h[2] * q$x, z)
        ex[i, ] <- colSums(X * q$w)
        ey[i, ] <- colSums(Y * q$w)
        cx[[i]] <- crossprod(sweep(X, 2, ex[i, ]) * sqrt(q$w))
        cy[[i]] <- crossprod(sweep(Y, 2, ey[i, ]) * sqrt(q$w))
    }
    list(ex = ex, ey = ey, cx = cx, cy = cy, order = order, tau = tau)
}

.np_conditional_conditional_cov <- function(z, mm, method, full = TRUE) {
    ex <- .np_conditional_xfeatures(z$x, z)
    p <- .np_conditional_params(ex, .np_conditional_yfeatures(z$y, z), z)
    a <- z$ix
    b <- z$iy
    w <- sweep(ex[, a, drop = FALSE], 2, z$n * p$fx, `/`)
    ca <- sweep(ex[, a + z$nx, drop = FALSE]/ex[, a, drop = FALSE], 2, 1 + p$aa, function(x, y) y - x)
    ca <- matrix(1 + p$aa, z$n, z$m, byrow = TRUE) - (outer(z$x, z$grid$x, function(a, b) (b - a)/z$h[1])^2 -
        1)/2
    wl <- sweep(ex[, a + 2 * z$nx, drop = FALSE], 2, z$n * p$fl, `/`)
    out <- if (full)
        matrix(0, z$m, z$m)
    else rep(0, z$m)
    for (i in seq_len(z$n)) {
        pr <- dnorm((z$x[i] - z$x)/(mm$tau * z$h[1]))
        pr <- pr/sum(pr)
        my <- colSums(mm$ey * pr)
        C <- Reduce(`+`, Map(function(c, p) c * p, mm$cy, pr)) + crossprod(sweep(mm$ey, 2, my) * sqrt(pr))
        tt <- switch(method, original = list(list(y = b, c = w[i, ])), analytic = list(list(y = b, c = w[i,
            ] * ca[i, ]), list(y = b + z$ny, c = 0.5 * w[i, ])), jackknife = list(list(y = b, c = 2 *
            w[i, ]), list(y = b + 2 * z$ny, c = -wl[i, ])))
        for (A in tt) for (B in tt) out <- out + if (full)
            C[A$y, B$y, drop = FALSE] * outer(A$c, B$c)
        else C[cbind(A$y, B$y)] * A$c * B$c
    }
    out
}

.np_conditional_contribution_range <- function(c) {
    t <- sqrt(pmax(0, 2 * c + 1))
    h <- function(v) c * pnorm(v) + 0.5 * v * dnorm(v)
    lo <- pmin(0, c)
    hi <- pmax(0, c)
    ok <- c >= -0.5
    lo[ok] <- pmin(lo[ok], h(t)[ok], h(-t)[ok])
    hi[ok] <- pmax(hi[ok], h(t)[ok], h(-t)[ok])
    list(lower = lo, upper = hi)
}

.np_conditional_weighted_range <- function(x, grid, hx) {
    u <- outer(x, grid$x, function(a, b) (b - a)/hx)
    K <- dnorm(u)
    w <- sweep(K, 2, colSums(K), `/`)
    a <- (u * u - 1)/2
    c <- matrix(1 + colSums(w * a), nrow(w), ncol(w), byrow = TRUE) - a
    b <- .np_conditional_contribution_range(c)
    d <- b$upper - b$lower
    A <- colSums(w * d)
    L <- colSums(w * b$lower)
    q <- sweep(w * d, 2, A, `/`)
    stopifnot(all(is.finite(q)), all(A > 0), max(abs(colSums(q) - 1)) < 1e-12)
    list(lower = L, span = A, qmax = apply(q, 2, max), q = q, range = b, weights = w, coefficient = c)
}

.np_conditional_range_interval <- function(center, b, method = "score", alpha = 0.05) {
    p <- (center - b$lower)/b$span
    q <- b$qmax
    if (any(!is.finite(p)) || any(p < 0 | p > 1))
        stop("corrected center outside certified signed range")
    if (method == "score") {
        z <- qnorm(1 - alpha/2)
        den <- 1 + z * z * q
        mid <- (p + z * z * q/2)/den
        rad <- z * sqrt(p * (1 - p) * q + z * z * q * q/4)/den
        lo <- mid - rad
        hi <- mid + rad
    }
    else if (method == "sharp-score") {
        z <- qnorm(1 - alpha/2)
        lo <- hi <- numeric(length(p))
        for (i in seq_along(p)) {
            state <- .np_conditional_sharp_prepare(b$q[, i])
            f <- function(t) (p[i] - t)^2 - z * z * .np_conditional_sharp_variance(t, state)
            lo[i] <- if (p[i] == 0)
                0
            else uniroot(f, c(0, p[i]), tol = 1e-13)$root
            hi[i] <- if (p[i] == 1)
                1
            else uniroot(f, c(p[i], 1), tol = 1e-13)$root
        }
    }
    else if (method == "Bernstein") {
        L <- log(2/alpha)
        lo <- hi <- numeric(length(p))
        for (i in seq_along(p)) {
            f <- function(t) (p[i] - t)^2 - 2 * L * q[i] * t * (1 - t) - (2/3) * L * q[i] * abs(p[i] -
                t)
            lo[i] <- if (f(0) <= 0)
                0
            else uniroot(f, c(0, p[i]), tol = 1e-13)$root
            hi[i] <- if (f(1) <= 0)
                1
            else uniroot(f, c(p[i], 1), tol = 1e-13)$root
        }
    }
    else stop("unknown range envelope")
    cbind(lower = b$lower + b$span * lo, upper = b$lower + b$span * hi)
}

.np_conditional_range_envelope <- function(base, group, b, method = "score", alpha = 0.05) {
    center <- base$center
    bands <- base$bands
    score <- .np_conditional_range_interval(center, b, method, alpha)
    point <- cbind(lower = pmin(bands$pointwise[, 1], score[, 1]), upper = pmax(bands$pointwise[, 2],
        score[, 2]))
    required <- pmax(center - bands$simultaneous[, 1], bands$simultaneous[, 2] - center, center - point[,
        1], point[, 2] - center)
    maximum <- vapply(split(required, group), max, 0)
    width <- unname(maximum[as.character(group)])
    sim <- cbind(lower = center - width, upper = center + width)
    bonf <- cbind(lower = pmin(bands$bonferroni[, 1], sim[, 1]), upper = pmax(bands$bonferroni[, 2],
        sim[, 2]))
    list(center = center, bands = list(pointwise = point, simultaneous = sim, bonferroni = bonf), range_interval = score,
        expansion = width - (bands$simultaneous[, 2] - center), range = list(lower = b$lower, span = b$span,
            qmax = b$qmax), method = method)
}

.np_conditional_sharp_prepare <- function(q) {
    stopifnot(all(is.finite(q)), all(q >= 0), abs(sum(q) - 1) < 1e-12)
    v <- sort(q[q > 0], decreasing = TRUE)
    list(q = v, total = cumsum(v), squares = cumsum(v * v), k = seq_along(v))
}

.np_conditional_sharp_variance <- function(p, s) {
    stopifnot(length(p) == 1, is.finite(p), p >= 0, p <= 1)
    p <- min(p, 1 - p)
    if (p == 0)
        return(0)
    lambda <- (s$total - 2 * p)/s$k
    k <- which(lambda < s$q & lambda >= c(s$q[-1], 0) - 1e-15)[1]
    if (is.na(k))
        stop("sharp variance active-set failure")
    (s$squares[k] - k * lambda[k]^2)/4
}

.np_conditional_bootstrap_slice <- function(center, root, group, mode, alpha = 0.05) {
    stopifnot(is.matrix(root), nrow(root) == length(center), all(is.finite(root)), ncol(root) > 1)
    B <- ncol(root)
    m <- nrow(root)
    se <- sqrt(rowSums((root - rowMeans(root))^2)/(B - 1))
    maximum <- vapply(split(se, group), max, 0)
    sg <- unname(maximum[as.character(group)])
    if (any(!is.finite(sg)) || any(sg <= 0))
        stop("entire conditioning slice has undefined bootstrap scale")
    R <- root/sg
    if (mode == "symmetric-sd") {
        point <- apply(abs(root), 1, quantile, probs = 1 - alpha, type = 1, names = FALSE)
        bonf <- apply(abs(root), 1, quantile, probs = 1 - alpha/m, type = 1, names = FALSE)
        critical <- quantile(apply(abs(R), 2, max), 1 - alpha, type = 1, names = FALSE)
        width <- critical * sg
        bands <- list(pointwise = cbind(center - point, center + point), simultaneous = cbind(center -
            width, center + width), bonferroni = cbind(center - pmax(bonf, width), center + pmax(bonf,
            width)))
    }
    else if (mode == "asymmetric-sd") {
        k <- 1 + floor(B * alpha/2)
        kb <- 1 + floor(B * alpha/(2 * m))
        v <- t(apply(root, 1, sort))
        mins <- sort(apply(R, 2, min))
        maxs <- sort(apply(R, 2, max))
        crit <- c(mins[k], maxs[B + 1 - k])
        critical <- crit
        point <- cbind(center - v[, B + 1 - k], center - v[, k])
        sim <- cbind(center - sg * crit[2], center - sg * crit[1])
        bonf <- cbind(center - v[, B + 1 - kb], center - v[, kb])
        bonf[, 1] <- pmin(bonf[, 1], sim[, 1])
        bonf[, 2] <- pmax(bonf[, 2], sim[, 2])
        bands <- list(pointwise = point, simultaneous = sim, bonferroni = bonf)
    }
    else stop("unknown bootstrap calibration")
    for (b in names(bands)) colnames(bands[[b]]) <- c("lower", "upper")
    list(center = center, se = se, scale = sg, root_mean = rowMeans(root), critical = critical, bands = bands,
        calibration = mode, B = B)
}

.np_conditional_bootstrap_score <- function(base, group, b, mode, alpha = 0.05) {
    score <- .np_conditional_range_interval(base$center, b, "sharp-score", alpha)
    center <- base$center
    bands <- base$bands
    point <- cbind(lower = pmin(bands$pointwise[, 1], score[, 1]), upper = pmax(bands$pointwise[, 2],
        score[, 2]))
    lower <- pmax(center - bands$simultaneous[, 1], center - point[, 1])
    upper <- pmax(bands$simultaneous[, 2] - center, point[, 2] - center)
    aggregate_slice <- function(v) {
        mx <- vapply(split(v, group), max, 0)
        unname(mx[as.character(group)])
    }
    if (mode == "symmetric-sd") {
        lower <- upper <- aggregate_slice(pmax(lower, upper))
    }
    else {
        lower <- aggregate_slice(lower)
        upper <- aggregate_slice(upper)
    }
    sim <- cbind(lower = center - lower, upper = center + upper)
    bonf <- cbind(lower = pmin(bands$bonferroni[, 1], sim[, 1]), upper = pmax(bands$bonferroni[, 2],
        sim[, 2]))
    base$reference_bands <- bands
    base$range_interval <- score
    base$bands <- list(pointwise = point, simultaneous = sim, bonferroni = bonf)
    base
}

.np_conditional_den_context <- function(x, y, grid, h) {
    stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y), length(x) >= 3, is.data.frame(grid),
        identical(names(grid), c("x", "y")), length(h) == 2, all(is.finite(c(x, y, unlist(grid), h))),
        all(h > 0))
    gx <- sort(unique(grid$x))
    gy <- sort(unique(grid$y))
    stopifnot(length(gy) <= 25, nrow(grid) <= 256)
    u <- outer(x, gx, function(a, b) (b - a)/h[1])
    K <- dnorm(u)
    mass <- colSums(K)
    if (any(mass <= 0))
        stop("undefined predictor mass")
    w <- sweep(K, 2, mass, `/`)
    a <- (u * u - 1)/2
    c <- matrix(1 + colSums(w * a), length(x), length(gx), byrow = TRUE) - a
    list(x = x, y = y, grid = grid, h = h, n = length(x), m = nrow(grid), gx = gx, gy = gy, nx = length(gx),
        ny = length(gy), ix = match(grid$x, gx), iy = match(grid$y, gy), w = w, c = c, A = w * (c + 0.5),
        D = -0.5 * w, a = a, effective_n = 1/colSums(w * w))
}

.np_conditional_den_features <- function(y, gy, h) {
    v <- outer(y, gy, function(a, b) (b - a)/h)
    K <- dnorm(v)/h
    cbind(K, v * v * K)
}

.np_conditional_den_apply_features <- function(z, F) {
    stopifnot(nrow(F) == z$n, ncol(F) == 2 * z$ny)
    T <- crossprod(z$A, F[, seq_len(z$ny), drop = FALSE]) + crossprod(z$D, F[, z$ny + seq_len(z$ny),
        drop = FALSE])
    T[cbind(z$ix, z$iy)]
}

.np_conditional_den_center <- function(z) {
    F <- .np_conditional_den_features(z$y, z$gy, z$h[2])
    O <- crossprod(z$w, F[, seq_len(z$ny), drop = FALSE])
    list(corrected = .np_conditional_den_apply_features(z, F), original = O[cbind(z$ix, z$iy)])
}

.np_conditional_den_normal_moments <- function(mu, gy, h, s = sqrt(2) * h) {
    n <- length(mu)
    ny <- length(gy)
    k <- 2L * ny
    mean <- matrix(0, n, k)
    second <- matrix(0, n, k * k)
    for (j in seq_len(ny)) {
        mass <- dnorm(gy[j], mu, sqrt(h * h + s * s))
        M <- (h * h * mu + s * s * gy[j])/(h * h + s * s)
        S2 <- h * h * s * s/(h * h + s * s)
        mean[, j] <- mass
        mean[, ny + j] <- mass * ((gy[j] - M)^2 + S2)/(h * h)
    }
    put <- function(j, k0, v) {
        second[, j + (k0 - 1L) * k] <<- v
    }
    for (j in seq_len(ny)) for (l in j:ny) {
        m <- (gy[j] + gy[l])/2
        r2 <- h * h/2
        J <- dnorm(gy[j], gy[l], sqrt(2) * h) * dnorm(m, mu, sqrt(s * s + r2))
        M <- (r2 * mu + s * s * m)/(s * s + r2)
        S2 <- s * s * r2/(s * s + r2)
        dj <- gy[j] - M
        dl <- gy[l] - M
        jk <- J * (dj * dj + S2)/(h * h)
        lk <- J * (dl * dl + S2)/(h * h)
        both <- J * (dj * dj * dl * dl + (dj * dj + dl * dl + 4 * dj * dl) * S2 + 3 * S2 * S2)/(h^4)
        put(j, l, J)
        put(l, j, J)
        put(ny + j, l, jk)
        put(l, ny + j, jk)
        put(j, ny + l, lk)
        put(ny + l, j, lk)
        put(ny + j, ny + l, both)
        put(ny + l, ny + j, both)
    }
    list(mean = mean, second = second)
}

.np_conditional_den_probabilities <- function(z, ii, tau = sqrt(2)) {
    p <- dnorm(outer(z$x[ii], z$x, `-`)/(tau * z$h[1]))
    s <- rowSums(p)
    if (any(s <= 0) || any(!is.finite(s)))
        stop("undefined pilot probability mass")
    p/s
}

.np_conditional_den_pilot <- function(z, tile = 32L, max_bytes = 256 * 1024^2) {
    k <- 2 * z$ny
    needed <- 8 * (4 * z$n * k * k + 8 * z$n * k + 8 * z$m * z$m + 8 * tile * z$n)
    if (!is.finite(needed) || needed > max_bytes)
        stop("density moment workspace exceeds budget")
    mm <- .np_conditional_den_normal_moments(z$y, z$gy, z$h[2])
    mu <- mm$mean * 0
    second <- mm$second * 0
    for (first in seq.int(1L, z$n, by = tile)) {
        ii <- first:min(z$n, first + tile - 1L)
        p <- .np_conditional_den_probabilities(z, ii)
        mu[ii, ] <- p %*% mm$mean
        second[ii, ] <- p %*% mm$second
    }
    C <- O <- matrix(0, z$m, z$m)
    j <- z$iy
    l <- z$iy + z$ny
    jidx <- outer(j, j, function(a, b) a + (b - 1L) * k)
    jlidx <- outer(j, l, function(a, b) a + (b - 1L) * k)
    ljidx <- outer(l, j, function(a, b) a + (b - 1L) * k)
    lidx <- outer(l, l, function(a, b) a + (b - 1L) * k)
    for (i in seq_len(z$n)) {
        V <- second[i, ] - as.vector(outer(mu[i, ], mu[i, ]))
        Vjj <- matrix(V[jidx], z$m, z$m)
        Vjl <- matrix(V[jlidx], z$m, z$m)
        Vlj <- matrix(V[ljidx], z$m, z$m)
        Vll <- matrix(V[lidx], z$m, z$m)
        A <- z$A[i, z$ix]
        D <- z$D[i, z$ix]
        W <- z$w[i, z$ix]
        C <- C + outer(A, A) * Vjj + outer(A, D) * Vjl + outer(D, A) * Vlj + outer(D, D) * Vll
        O <- O + outer(W, W) * Vjj
    }
    C <- (C + t(C))/2
    O <- (O + t(O))/2
    list(mean = mu, covariance = C, original_covariance = O, expectation = .np_conditional_den_apply_features(z,
        mu), estimated_workspace = needed, tile = tile)
}

.np_conditional_den_range <- function(z) {
    c <- z$c
    t <- sqrt(pmax(0, 2 * c + 3))
    h <- z$h[2]
    value0 <- (c + 0.5) * dnorm(0)/h
    valuet <- (c + 0.5 - 0.5 * t * t) * dnorm(t)/h
    lo <- pmin(value0, 0)
    hi <- pmax(value0, 0)
    ok <- c >= -1.5
    lo[ok] <- pmin(lo[ok], valuet[ok])
    hi[ok] <- pmax(hi[ok], valuet[ok])
    d <- hi - lo
    L <- colSums(z$w * lo)
    A <- colSums(z$w * d)
    q <- sweep(z$w * d, 2, A, `/`)
    list(lower = L[z$ix], span = A[z$ix], q = q[, z$ix, drop = FALSE], qmax = apply(q, 2, max)[z$ix],
        range = list(lower = lo, upper = hi))
}

.np_conditional_jack_context <- function(z) {
    H <- sqrt(2) * z$h
    u <- outer(z$x, z$gx, function(a, b) (b - a)/H[1])
    K <- dnorm(u)
    w <- sweep(K, 2, colSums(K), `/`)
    z$jw <- w
    z$H <- H
    z
}

.np_conditional_jack_moments <- function(mu, gy, h, s) {
    ny <- length(gy)
    target <- rep(gy, 2)
    hh <- rep(c(h, sqrt(2) * h), each = ny)
    k <- length(hh)
    mean <- vapply(seq_len(k), function(j) dnorm(target[j], mu, sqrt(s * s + hh[j]^2)), numeric(length(mu)))
    second <- matrix(0, length(mu), k * k)
    for (j in seq_len(k)) for (l in j:k) {
        a <- hh[j]^2
        b <- hh[l]^2
        m <- (target[j] * b + target[l] * a)/(a + b)
        r <- a * b/(a + b)
        val <- dnorm(target[j], target[l], sqrt(a + b)) * dnorm(m, mu, sqrt(s * s + r))
        second[, j + (l - 1) * k] <- val
        second[, l + (j - 1) * k] <- val
    }
    list(mean = mean, second = second)
}

.np_conditional_jack_center <- function(z) {
    K <- .np_conditional_den_features(z$y, z$gy, z$h[2])[, seq_len(z$ny), drop = FALSE]
    KH <- .np_conditional_den_features(z$y, z$gy, z$H[2])[, seq_len(z$ny), drop = FALSE]
    M <- 2 * crossprod(z$w, K) - crossprod(z$jw, KH)
    M[cbind(z$ix, z$iy)]
}

.np_conditional_jack_pilot <- function(z, tile = 32L) {
    k <- 2 * z$ny
    mm <- .np_conditional_jack_moments(z$y, z$gy, z$h[2], sqrt(2) * z$h[2])
    mu <- mm$mean * 0
    second <- mm$second * 0
    for (first in seq.int(1L, z$n, by = tile)) {
        ii <- first:min(z$n, first + tile - 1L)
        pr <- .np_conditional_den_probabilities(z, ii)
        mu[ii, ] <- pr %*% mm$mean
        second[ii, ] <- pr %*% mm$second
    }
    j <- z$iy
    l <- j + z$ny
    C <- matrix(0, z$m, z$m)
    for (i in seq_len(z$n)) {
        V <- matrix(second[i, ], k, k) - outer(mu[i, ], mu[i, ])
        A <- 2 * z$w[i, z$ix]
        D <- -z$jw[i, z$ix]
        C <- C + outer(A, A) * V[j, j] + outer(A, D) * V[j, l] + outer(D, A) * V[l, j] + outer(D, D) *
            V[l, l]
    }
    M <- 2 * crossprod(z$w, mu[, seq_len(z$ny), drop = FALSE]) - crossprod(z$jw, mu[, z$ny + seq_len(z$ny),
        drop = FALSE])
    list(mean = mu, covariance = (C + t(C))/2, expectation = M[cbind(z$ix, z$iy)])
}

.np_conditional_jack_range <- function(z) {
    A <- 2 * z$w
    B <- z$jw
    h <- z$h[2]
    H <- z$H[2]
    v0 <- (A/h - B/H) * dnorm(0)
    lo <- pmin(v0, 0)
    hi <- pmax(v0, 0)
    both <- A > 0 & B > 0
    s <- matrix(-1, nrow(A), ncol(A))
    s[both] <- 2 * (log(A[both]) - log(B[both]) + 3 * log(H/h))/(1/h^2 - 1/H^2)
    ok <- both & s >= 0
    t <- sqrt(pmax(s, 0))
    v <- A * dnorm(t/h)/h - B * dnorm(t/H)/H
    lo[ok] <- pmin(lo[ok], v[ok])
    hi[ok] <- pmax(hi[ok], v[ok])
    span <- colSums(hi - lo)
    L <- colSums(lo)
    q <- sweep(hi - lo, 2, span, `/`)
    list(lower = L[z$ix], span = span[z$ix], q = q[, z$ix, drop = FALSE], qmax = apply(q, 2, max)[z$ix],
        range = list(lower = lo, upper = hi))
}

.np_conditional_marginal_guard <- function(base, bounds, ref, alpha = 0.05) {
    pointscore <- ref$range_interval(base$center, bounds, "sharp-score", alpha)
    hull <- function(a, b) cbind(lower = pmin(a[, 1], b[, 1]), upper = pmax(a[, 2], b[, 2]))
    point <- hull(base$bands$pointwise, pointscore)
    sim <- hull(base$bands$simultaneous, point)
    bonf <- hull(base$bands$bonferroni, sim)
    list(center = base$center, bands = list(pointwise = point, simultaneous = sim, bonferroni = bonf),
        gaussian_critical = base$critical, marginal_score_critical = qnorm(1 - alpha/2), range_interval = pointscore)
}

.np_conditional_cp_allocation <- function(n, nx, ny, m, B, workspace_bytes) {
    dims <- c(n, nx, ny, m, B)
    if (any(!is.finite(dims)) || any(dims < 1) || any(dims != floor(dims)) || any(dims > .Machine$integer.max) ||
        length(workspace_bytes) != 1L || !is.finite(workspace_bytes) || workspace_bytes <= 0)
        stop("invalid conditional-pilot dimensions or workspace budget")
    dims <- as.double(dims)
    n <- dims[1]
    nx <- dims[2]
    ny <- dims[3]
    m <- dims[4]
    B <- dims[5]
    base_doubles <- m * B + 2 * n * nx + 4 * n * ny + 12 * n * 128 + 4 * m * 128
    base_bytes <- 8 * base_doubles
    row_bytes <- 48 * n
    if (!is.finite(base_bytes) || base_bytes > 2^53 || base_bytes + row_bytes > workspace_bytes)
        stop("conditional-pilot workspace cannot hold the fixed logical block")
    tile <- min(32, n, floor((workspace_bytes - base_bytes)/row_bytes))
    list(tile = as.integer(tile), base_bytes = base_bytes, donor_row_bytes = row_bytes, estimated_bytes = base_bytes +
        row_bytes * tile, logical_block = 128L, workspace_bytes = workspace_bytes)
}

.np_conditional_cp_probabilities <- function(plan, rows) {
    z <- plan$z
    p <- dnorm(outer(z$x[rows], z$x, `-`)/(plan$tau * z$h[1]))
    sums <- rowSums(p)
    if (any(!is.finite(sums)) || any(sums <= 0))
        stop("invalid conditional-pilot probability mass")
    p/sums
}

.np_conditional_cp_plan <- function(z, B, workspace_bytes) {
    if (length(z$x) != z$n || length(z$y) != z$n || length(z$h) != 2L || any(!is.finite(c(z$x, z$y, z$h,
        z$gx, z$gy))) || any(z$h <= 0))
        stop("invalid conditional-pilot input")
    allocation <- .np_conditional_cp_allocation(z$n, z$nx, z$ny, z$m, B, workspace_bytes)
    u <- outer(z$x, z$gx, function(a, b) (b - a)/z$h[1])
    K <- dnorm(u)
    w <- sweep(K, 2, colSums(K), `/`)
    a <- (u * u - 1)/2
    A <- w * (matrix(1 + colSums(w * a), z$n, z$nx, byrow = TRUE) - a)
    plan <- list(z = z, w = w, A = A, tau = sqrt(2), B = as.integer(B), allocation = allocation)
    d <- outer(z$y, z$gy, function(a, b) b - a)
    ss <- sqrt(3) * z$h[2]
    eg <- pnorm(d/ss)
    j <- which(z$ref == 1)
    if (length(j))
        eg[, j] <- -pnorm(d[, j, drop = FALSE]/ss, lower.tail = FALSE)
    eu <- z$h[2]^2 * d * dnorm(d/ss)/ss^3
    mg <- mu <- matrix(0, z$n, z$ny)
    for (first in seq.int(1L, z$n, by = allocation$tile)) {
        ii <- first:min(z$n, first + allocation$tile - 1L)
        p <- .np_conditional_cp_probabilities(plan, ii)
        mg[ii, ] <- p %*% eg
        mu[ii, ] <- p %*% eu
    }
    plan$meanG <- mg
    plan$meanU <- mu
    plan
}

.np_conditional_cp_apply <- function(plan, Y) {
    z <- plan$z
    stopifnot(is.matrix(Y), nrow(Y) == z$n)
    out <- matrix(0, z$m, ncol(Y))
    for (b in seq_len(z$ny)) {
        v <- (z$gy[b] - Y)/z$h[2]
        G <- if (z$ref[b] == 1)
            -pnorm(v, lower.tail = FALSE)
        else pnorm(v)
        U <- v * dnorm(v)
        root <- crossprod(plan$A, G - plan$meanG[, b]) + 0.5 * crossprod(plan$w, U - plan$meanU[, b])
        ii <- which(z$iy == b)
        out[ii, ] <- root[z$ix[ii], , drop = FALSE]
    }
    out
}
