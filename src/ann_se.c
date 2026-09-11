/* Unconditional ANN uncertainty only. Kernel values and point estimates stay
 * with the established kernel-sum/density owners. Each query's influence is
 * accumulated on the original joint sample rows, including all radius terms.
 * No n-by-n matrix or cross-query covariance is allocated. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <limits.h>
#include <stdint.h>
#include <math.h>
#include "conditional_ann_direct.h"

static double ann_distance_k(const double *z, int n, int i, int k)
{
    int low = k - (n-i-1);
    if (low < 0) low = 0;
    int high = k < i ? k : i;
    while (low <= high) {
        const int a = low + (high-low)/2, b = k-a;
        const double la = a ? z[i]-z[i-a] : R_NegInf;
        const double rb = b ? z[i+b]-z[i] : R_NegInf;
        const double ln = a < i ? z[i]-z[i-a-1] : R_PosInf;
        const double rn = b < n-i-1 ? z[i+b+1]-z[i] : R_PosInf;
        if (la <= rn && rb <= ln) return la > rb ? la : rb;
        if (la > rn) high = a-1; else low = a+1;
    }
    error("ANN standard-error distance partition failed");
    return NA_REAL;
}

/* Two-sided rank spacing. Correct cube boundaries using integer arithmetic;
 * this pilot does not change the estimator's selected nearest-neighbor rank. */
static int ann_window(int N, int k)
{
    if (k < 1 || k > N) return 0;
    int m = (int)floor(pow((double)N, 2.0/3.0));
    const uint64_t nn = (uint64_t)N * (uint64_t)N;
    while ((uint64_t)(m+1)*(m+1)*(m+1) <= nn) ++m;
    while ((uint64_t)m*m*m > nn) --m;
    const int edge = (k < N-k ? k : N-k)/2;
    return m < edge ? m : edge;
}

SEXP C_np_ann_geometry(SEXP sorted, SEXP index)
{
    if (TYPEOF(sorted) != REALSXP || XLENGTH(sorted) < 2 ||
        XLENGTH(sorted) > INT_MAX || TYPEOF(index) != INTSXP ||
        XLENGTH(index) != 1 || INTEGER(index)[0] == NA_INTEGER)
        error("ANN standard-error geometry has invalid inputs");
    const int n = (int)XLENGTH(sorted), N = n-1, k = INTEGER(index)[0];
    const double *z = REAL(sorted);
    for (int i = 0; i < n; ++i)
        if (!R_FINITE(z[i]) || (i && z[i] < z[i-1]))
            error("ANN standard-error geometry requires finite sorted coordinates");
    const int m = ann_window(N, k);
    SEXP out = PROTECT(allocVector(VECSXP, 5));
    SEXP names = PROTECT(allocVector(STRSXP, 5));
    const char *keys[] = {"radius", "density.sum", "left", "right", "status"};
    for (int j = 0; j < 5; ++j) {
        SET_VECTOR_ELT(out, j, allocVector(j < 2 ? REALSXP : INTSXP, n));
        SET_STRING_ELT(names, j, mkChar(keys[j]));
    }
    setAttrib(out, R_NamesSymbol, names);
    for (int i = 0; i < n; ++i) {
        if ((i & 255) == 0) R_CheckUserInterrupt();
        double r = NA_REAL, b = NA_REAL;
        int status = 0, left = 0, right = 0;
        if (k < 1 || k > N) status = 1;
        else {
            r = ann_distance_k(z, n, i, k);
            if (!R_FINITE(r) || r <= 0) status = 2;
            else if (m < 1) status = 3;
            else {
                const double width = ann_distance_k(z,n,i,k+m) -
                    ann_distance_k(z,n,i,k-m);
                if (!R_FINITE(width) || width <= 0) status = 4;
                else {
                    b = (2.0*m/N)/width;
                    if (!R_FINITE(b) || b <= 0) status = 5;
                }
            }
            if (!status) {
                int lo = 0, hi = i;
                /* Use actual rounded differences, not rounded z +/- r. */
                while (lo < hi) {
                    const int mid = lo + (hi-lo)/2;
                    if (z[i]-z[mid] <= r) hi = mid; else lo = mid+1;
                }
                left = lo; lo = i; hi = n;
                while (lo < hi) {
                    const int mid = lo + (hi-lo)/2;
                    if (z[mid]-z[i] <= r) lo = mid+1; else hi = mid;
                }
                right = lo;
            }
        }
        REAL(VECTOR_ELT(out,0))[i] = r;
        REAL(VECTOR_ELT(out,1))[i] = b;
        INTEGER(VECTOR_ELT(out,2))[i] = left;
        INTEGER(VECTOR_ELT(out,3))[i] = right;
        INTEGER(VECTOR_ELT(out,4))[i] = status;
    }
    UNPROTECT(2);
    return out;
}

static void ann_matrix(SEXP x, int *nr, int *nc)
{
    SEXP dim = getAttrib(x, R_DimSymbol);
    if (TYPEOF(x) != REALSXP || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
        INTEGER(dim)[0] < 0 || INTEGER(dim)[1] < 0 ||
        XLENGTH(x) != (R_xlen_t)INTEGER(dim)[0]*INTEGER(dim)[1])
        error("ANN standard-error input is not a double matrix");
    *nr = INTEGER(dim)[0]; *nc = INTEGER(dim)[1];
}

static SEXP ann_variance(SEXP geometry, SEXP order, SEXP training,
                       SEXP evaluation, SEXP weights, SEXP derivative,
                       SEXP density, SEXP faces)
{
    int n, p, ne, pe, nw, m;
    ann_matrix(training, &n, &p);
    ann_matrix(evaluation, &ne, &pe);
    ann_matrix(weights, &nw, &m);
    SEXP dd = getAttrib(derivative, R_DimSymbol);
    if (n < 2 || p < 1 || pe != p || nw != n || ne != m ||
        TYPEOF(geometry) != VECSXP || XLENGTH(geometry) != p ||
        TYPEOF(order) != VECSXP || XLENGTH(order) != p ||
        TYPEOF(derivative) != REALSXP ||
        (p && XLENGTH(weights) > R_XLEN_T_MAX/p) ||
        XLENGTH(derivative) != XLENGTH(weights)*p ||
        TYPEOF(density) != LGLSXP || XLENGTH(density) != 1 ||
        LOGICAL(density)[0] == NA_LOGICAL)
        error("ANN standard-error dimensions do not conform");
    if (TYPEOF(dd) != INTSXP ||
        !((p == 1 && XLENGTH(dd) == 2) || XLENGTH(dd) == 3) ||
        INTEGER(dd)[0] != n || INTEGER(dd)[1] != m ||
        (XLENGTH(dd) == 3 && INTEGER(dd)[2] != p))
        error("ANN standard-error derivative tensor does not conform");
    const double *face_cut = NULL, *face_coef = NULL;
    if (faces != R_NilValue) {
        int nf, pf, nc, pc;
        if (TYPEOF(faces) != VECSXP || XLENGTH(faces) != 2 ||
            !LOGICAL(density)[0] || p > INT_MAX/2)
            error("ANN standard-error boundary payload is invalid");
        ann_matrix(VECTOR_ELT(faces,0), &nf, &pf);
        ann_matrix(VECTOR_ELT(faces,1), &nc, &pc);
        if (nf != m || nc != m || pf != 2*p || pc != 2*p)
            error("ANN standard-error boundary dimensions do not conform");
        face_cut = REAL(VECTOR_ELT(faces,0));
        face_coef = REAL(VECTOR_ELT(faces,1));
    }
    int *seen = (int *)R_alloc(n, sizeof(int));
    for (int j = 0; j < p; ++j) {
        SEXP g = VECTOR_ELT(geometry,j), ord = VECTOR_ELT(order,j);
        if (TYPEOF(g) != VECSXP || XLENGTH(g) != 5 ||
            TYPEOF(ord) != INTSXP || XLENGTH(ord) != n)
            error("ANN standard-error geometry dimensions do not conform");
        for (int h = 0; h < 5; ++h)
            if (TYPEOF(VECTOR_ELT(g,h)) != (h < 2 ? REALSXP : INTSXP) ||
                XLENGTH(VECTOR_ELT(g,h)) != n)
                error("ANN standard-error geometry buffers do not conform");
        for (int i = 0; i < n; ++i) seen[i] = 0;
        for (int i = 0; i < n; ++i) {
            const int id = INTEGER(ord)[i];
            const double r = REAL(VECTOR_ELT(g,0))[i];
            const double b = REAL(VECTOR_ELT(g,1))[i];
            const int l = INTEGER(VECTOR_ELT(g,2))[i];
            const int u = INTEGER(VECTOR_ELT(g,3))[i];
            if (id < 1 || id > n || seen[id-1]++ ||
                !R_FINITE(r) || r <= 0 || !R_FINITE(b) || b <= 0 ||
                l < 0 || l > i || u <= i || u > n ||
                INTEGER(VECTOR_ELT(g,4))[i] != 0)
                error("ANN standard-error geometry is invalid");
        }
    }
    SEXP out = PROTECT(allocVector(REALSXP,m));
    long double *phi = (long double *)R_alloc(n, sizeof(long double));
    long double *delta = (long double *)R_alloc((size_t)n+1, sizeof(long double));
    const double *x = REAL(training), *e = REAL(evaluation);
    const double *a = REAL(weights), *d = REAL(derivative);
    for (int q = 0; q < m; ++q) {
        R_CheckUserInterrupt();
        int valid = 1;
        for (int i = 0; i < n; ++i) {
            phi[i] = a[i+(R_xlen_t)n*q];
            if (!R_FINITE((double)phi[i])) valid = 0;
        }
        for (int j = 0; j < p && valid; ++j) {
            SEXP g = VECTOR_ELT(geometry,j);
            const int *ord = INTEGER(VECTOR_ELT(order,j));
            const double *r = REAL(VECTOR_ELT(g,0)), *b = REAL(VECTOR_ELT(g,1));
            const int *left = INTEGER(VECTOR_ELT(g,2)), *right = INTEGER(VECTOR_ELT(g,3));
            for (R_xlen_t i = 0; i <= (R_xlen_t)n; ++i) delta[i] = 0;
            for (int s = 0; s < n; ++s) {
                const int i = ord[s]-1;
                const R_xlen_t pos = i+(R_xlen_t)n*q;
                const long double difference = (long double)e[q+(R_xlen_t)m*j]-x[i+(R_xlen_t)n*j];
                long double w = difference*d[pos+(R_xlen_t)n*m*j];
                if (LOGICAL(density)[0]) w += a[pos];
                w /= (long double)r[s]*b[s];
                if (!R_FINITE((double)w)) { valid = 0; break; }
                delta[left[s]] += w;
                delta[right[s]] -= w;
            }
            long double prefix = 0;
            for (int s = 0; s < n && valid; ++s) {
                prefix += delta[s];
                phi[ord[s]-1] += prefix/n;
            }
        }
        if (face_cut && valid) {
            for (int j = 0; j < p && valid; ++j) {
                const double at = e[q+(R_xlen_t)m*j];
                for (int h = 0; h < 2; ++h) {
                    const R_xlen_t slot = q+(R_xlen_t)m*(2*j+h);
                    const double cut = face_cut[slot], coef = face_coef[slot];
                    if (!R_FINITE(cut) || !R_FINITE(coef)) { valid = 0; break; }
                    for (int i = 0; i < n; ++i) {
                        const double z = x[i+(R_xlen_t)n*j];
                        phi[i] += (long double)coef*((z <= at)-(z <= cut));
                    }
                }
            }
        }
        long double mean = 0;
        NPANNConditionalNorm norm = {0.0, 0.0, 0};
        if (valid) {
            for (int i = 0; i < n; ++i) mean += phi[i];
            mean /= n;
            for (int i = 0; i < n; ++i) {
                np_ann_direct_append(&norm, (double)((phi[i]-mean)/n));
            }
        }
        double se = NA_REAL;
        if (valid && !np_ann_direct_finish(&norm, (size_t)n, &se)) se = NA_REAL;
        REAL(out)[q] = se;
    }
    UNPROTECT(1);
    return out;
}

SEXP C_np_ann_variance(SEXP geometry, SEXP order, SEXP training,
                       SEXP evaluation, SEXP weights, SEXP derivative, SEXP density)
{
    return ann_variance(geometry,order,training,evaluation,weights,derivative,density,R_NilValue);
}

SEXP C_np_ann_variance_faces(SEXP geometry, SEXP order, SEXP training,
                       SEXP evaluation, SEXP weights, SEXP derivative, SEXP density, SEXP faces)
{
    return ann_variance(geometry,order,training,evaluation,weights,derivative,density,faces);
}
