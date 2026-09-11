/* Fixed-coordinate copula ratio influence. Kernel weights remain owned by
 * the canonical kernel-sum leaves; only scalar SEs leave this reducer. */
#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <math.h>

SEXP C_np_copula_density_se(SEXP weights)
{
    if (TYPEOF(weights) != VECSXP || XLENGTH(weights) < 1 ||
        XLENGTH(weights) > INT_MAX)
        error("copula SE requires a nonempty marginal weight list");
    const int p = (int) XLENGTH(weights);
    SEXP first = VECTOR_ELT(weights, 0);
    SEXP dims = getAttrib(first, R_DimSymbol);
    if (TYPEOF(first) != REALSXP || TYPEOF(dims) != INTSXP ||
        XLENGTH(dims) != 2)
        error("copula SE requires double weight matrices");
    const int n = INTEGER(dims)[0], m = INTEGER(dims)[1];
    if (n < 1 || m < 1 || XLENGTH(first) != (R_xlen_t)n * m)
        error("copula SE requires nonempty conformable matrices");
    const double **w = (const double **) R_alloc((size_t)p, sizeof(double *));
    for (int j = 0; j < p; ++j) {
        SEXP x = VECTOR_ELT(weights, j);
        SEXP d = getAttrib(x, R_DimSymbol);
        if (TYPEOF(x) != REALSXP || TYPEOF(d) != INTSXP ||
            XLENGTH(d) != 2 || INTEGER(d)[0] != n || INTEGER(d)[1] != m ||
            XLENGTH(x) != (R_xlen_t)n * m)
            error("copula SE marginal matrix dimensions differ");
        w[j] = REAL(x);
    }
    double *product = (double *) R_alloc((size_t)n, sizeof(double));
    double *total = (double *) R_alloc((size_t)n, sizeof(double));
    SEXP out = PROTECT(allocVector(REALSXP, m));
    for (int k = 0; k < m; ++k) {
        R_CheckUserInterrupt();
        REAL(out)[k] = NA_REAL;
        int valid = 1;
        for (int i = 0; i < n; ++i) {
            product[i] = 1.0;
            total[i] = 0.0;
        }
        for (int j = 0; j < p; ++j) {
            const double *column = w[j] + (R_xlen_t)k * n;
            long double sum = 0.0L;
            for (int i = 0; i < n; ++i) sum += column[i];
            const double mean = (double)(sum/n);
            /* No NZD floor: the ratio is undefined at a zero denominator.
             * Signed higher-order kernel weights need not be positive. */
            if (!R_FINITE(mean) || mean == 0.0) {
                valid = 0;
                break;
            }
            for (int i = 0; i < n; ++i) {
                const double u = column[i]/mean;
                product[i] *= u;
                total[i] += u;
            }
        }
        if (!valid) continue;
        long double joint = 0.0L;
        for (int i = 0; i < n; ++i) joint += product[i];
        const double c = (double)(joint/n);
        if (!R_FINITE(c)) continue;
        long double center = 0.0L;
        for (int i = 0; i < n; ++i) {
            product[i] -= c * (total[i] - p + 1.0);
            center += product[i];
        }
        const double mean_phi = (double)(center/n);
        if (!R_FINITE(mean_phi)) continue;
        long double ss = 0.0L;
        for (int i = 0; i < n; ++i) {
            const long double delta = (long double)product[i] - mean_phi;
            ss += delta * delta;
        }
        const double se = (double)(sqrtl(ss)/n);
        if (R_FINITE(se)) REAL(out)[k] = se;
    }
    UNPROTECT(1);
    return out;
}
