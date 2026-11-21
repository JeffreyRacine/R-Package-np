#include <R.h>
#include <Rinternals.h>
#include <math.h>

/* Use the package's configured Accelerate linkage on Apple Silicon. The
   portable branch has identical clipping, tie and diagonal conventions. */
#if defined(__APPLE__) && defined(__arm64__) && defined(NP_USE_ACCELERATE_GAUSS)
extern void vvacos(double *, const double *, const int *);
#endif

static void angles(double *out, const double *cosines, int count) {
#if defined(__APPLE__) && defined(__arm64__) && defined(NP_USE_ACCELERATE_GAUSS)
    vvacos(out, cosines, &count);
#else
    for (int i = 0; i < count; ++i) out[i] = acos(cosines[i]);
#endif
}

/* The integrated projection weights are symmetric. Accumulate one triangle,
   retaining the original point masses at coincident observations. */
SEXP C_np_pltest_angles(SEXP coordinates) {
    int n = nrows(coordinates), p = ncols(coordinates);
    if (!isReal(coordinates) || !isMatrix(coordinates) || n < 2 || p < 1)
        error("Invalid projection coordinates.");
    const double *x = REAL(coordinates);
    SEXP ans = PROTECT(allocMatrix(REALSXP, n, n));
    double *a = REAL(ans);
    R_xlen_t nn = (R_xlen_t)n * n;
    for (R_xlen_t i = 0; i < nn; ++i) a[i] = 0.0;
    double *u = (double *)R_alloc((size_t)n * p, sizeof(double));
    int *zero = (int *)R_alloc(n, sizeof(int));
    double *cosines = (double *)R_alloc(n, sizeof(double));
    double *angle = (double *)R_alloc(n, sizeof(double));
    for (int k = 0; k < n; ++k) {
        R_CheckUserInterrupt();
        for (int i = 0; i < n; ++i) {
            double *v = u + (size_t)i*p;
            double largest = 0.0, norm = 0.0;
            for (int d = 0; d < p; ++d) {
                const double *col = x + (size_t)d*n;
                v[d] = col[i] - col[k];
                if (!R_FINITE(v[d])) error("Numerically unrepresentable projection differences; rescale the predictors.");
                largest = fmax(largest, fabs(v[d]));
            }
            zero[i] = largest == 0.0;
            if (zero[i]) continue;
            for (int d = 0; d < p; ++d) {
                v[d] /= largest;
                norm += v[d]*v[d];
            }
            norm = sqrt(norm);
            for (int d = 0; d < p; ++d) v[d] /= norm;
        }
        for (int j = 0; j < n; ++j) {
            double *col = a + (size_t)j*n;
            col[j] += zero[j] ? 1.0 : 0.5;
            const double *v = u + (size_t)j*p;
            if (zero[j]) {
                for (int i = j+1; i < n; ++i) col[i] += zero[i] ? 1.0 : 0.5;
                continue;
            }
            int count = n-j-1;
            if (count == 0) continue;
            for (int t = 0; t < count; ++t) {
                int i = j+1+t;
                const double *w = u + (size_t)i*p;
                double dot = 0.0;
                for (int d = 0; d < p; ++d) dot += w[d]*v[d];
                cosines[t] = zero[i] ? -1.0 : -fmax(-1.0, fmin(1.0, dot));
            }
            angles(angle, cosines, count);
            for (int t = 0; t < count; ++t) col[j+1+t] += angle[t] / (2.0*M_PI);
        }
    }
    for (int j = 0; j < n; ++j)
        for (int i = j; i < n; ++i) {
            a[(size_t)j*n+i] /= (double)n*n;
            a[(size_t)i*n+j] = a[(size_t)j*n+i];
        }
    UNPROTECT(1);
    return ans;
}

/* Maistre-Patilea: common scaling of all density-weighted residuals cancels
   on studentization. Pair order matches stats::dist (lower triangle). */
SEXP C_np_indextest_score(SEXP response, SEXP index, SEXP distances,
                          SEXP bandwidth, SEXP widths) {
    int n = LENGTH(response), m = LENGTH(widths);
    if (!isReal(response) || !isReal(index) || LENGTH(index) != n || n < 2 ||
        !isReal(distances) || XLENGTH(distances) != (R_xlen_t)n*(n-1)/2 ||
        !isReal(bandwidth) || XLENGTH(bandwidth) != 1 || !isReal(widths) || m < 1)
        error("Invalid single-index statistic inputs.");
    const double *y = REAL(response), *z = REAL(index), *d2 = REAL(distances);
    double g = asReal(bandwidth);
    if (!R_FINITE(g) || g <= 0.0) error("Invalid residual bandwidth.");
    SEXP ans = PROTECT(allocVector(REALSXP,m));
    double *score = REAL(ans);
    const double *h = REAL(widths);
    double *u = (double *)R_alloc(n,sizeof(double));
    double *squares = (double *)R_alloc(m,sizeof(double));
    for (int i = 0; i < n; ++i) u[i] = 0.0;
    for (int l = 0; l < m; ++l) {
        if (!R_FINITE(h[l]) || h[l] <= 0.0) error("Invalid test bandwidth.");
        score[l] = squares[l] = 0.0;
    }
    for (int i = 0; i < n; ++i) {
        if ((i % 64) == 0) R_CheckUserInterrupt();
        for (int j = i+1; j < n; ++j) {
            double t = (z[i]-z[j])/g;
            double delta = (y[i]-y[j])*exp(-0.5*t*t);
            u[i] += delta;
            u[j] -= delta;
        }
    }
    R_xlen_t pair = 0;
    for (int i = 0; i < n; ++i) {
        if ((i % 64) == 0) R_CheckUserInterrupt();
        for (int j = i+1; j < n; ++j, ++pair) {
            double dz = z[i]-z[j];
            double orthogonal = fmax(d2[pair]-dz*dz,0.0);
            double product = u[i]*u[j];
            for (int l = 0; l < m; ++l) {
                double t = dz/h[l];
                double a = product*exp(-0.5*(orthogonal+t*t));
                score[l] += a;
                squares[l] += a*a;
            }
        }
    }
    for (int l = 0; l < m; ++l) {
        score[l] *= sqrt((double)n/(n-1))/sqrt(squares[l]);
        if (!R_FINITE(score[l])) error("Degenerate single-index pairwise variance.");
    }
    UNPROTECT(1);
    return ans;
}
