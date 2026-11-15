/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 2000--2016 The R Core Team
 *  Copyright (C) 2001--2016 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Reference:
 * R. C. H. Cheng (1978).
 * Generating beta variates with nonintegral shape parameters.
 * Communications of the ACM 21, 317-322.
 * (Algorithms BB and BC)
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <float.h>
/* Adapted for the new plot law only. The returned double and uniform draws
 * retain upstream Cheng arithmetic; paired logs retain the pre-division ratio.
 * Source: R trunk src/nmath/rbeta.c, frozen in the campaign contract.
 * The package does not replace or interpose R's beta generator. */
static double lx_record,lq_record;
static void record_pair(double numerator, double complement) {
 double l=log(numerator),r=log(complement),d=l-r;
 lx_record= d>=0 ? -log1p(exp(-d)) : d-log1p(exp(d));
 lq_record= d>=0 ? -d-log1p(exp(-d)) : -log1p(exp(d));

}

#define expmax	(DBL_MAX_EXP * M_LN2)/* = log(DBL_MAX) */

static double observed_rbeta(double aa, double bb)
{
    if (ISNAN(aa) || ISNAN(bb) || aa < 0. || bb < 0.)
	return R_NaN;
    if (!R_FINITE(aa) && !R_FINITE(bb)) // a = b = Inf : all mass at 1/2
	return 0.5;
    if (aa == 0. && bb == 0.) // point mass 1/2 at each of {0,1} :
	return (unif_rand() < 0.5) ? 0. : 1.;
    // now, at least one of a, b is finite and positive
    if (!R_FINITE(aa) || bb == 0.)
    	return 1.0;
    if (!R_FINITE(bb) || aa == 0.)
    	return 0.0;

    double a, b, alpha;
    double r, s, t, u1, u2, v, w, y, z;
    int qsame;
    /* FIXME:  Keep Globals (properly) for threading */
    /* Uses these GLOBALS to save time when many rv's are generated : */
    static double beta, gamma, delta, k1, k2;
    static double olda = -1.0;
    static double oldb = -1.0;

    /* Test if we need new "initializing" */
    qsame = (olda == aa) && (oldb == bb);
    if (!qsame) { olda = aa; oldb = bb; }

    a = fmin(aa, bb);
    b = fmax(aa, bb); /* a <= b */
    alpha = a + b;

#define v_w_from__u1_bet(AA) 			\
	    v = beta * log(u1 / (1.0 - u1));	\
	    if (v <= expmax) {			\
		w = AA * exp(v);		\
		if(!R_FINITE(w)) w = DBL_MAX;	\
	    } else				\
		w = DBL_MAX


    if (a <= 1.0) {	/* --- Algorithm BC --- */

	/* changed notation, now also a <= b (was reversed) */

	if (!qsame) { /* initialize */
	    beta = 1.0 / a;
	    delta = 1.0 + b - a;
	    k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
	    k2 = 0.25 + (0.5 + 0.25 / delta) * a;
	}
	/* FIXME: "do { } while()", but not trivially because of "continue"s:*/
	for(;;) {
	    u1 = unif_rand();
	    u2 = unif_rand();
	    if (u1 < 0.5) {
		y = u1 * u2;
		z = u1 * y;
		if (0.25 * u2 + z - y >= k1)
		    continue;
	    } else {
		z = u1 * u1 * u2;
		if (z <= 0.25) {
		    v_w_from__u1_bet(b);
		    break;
		}
		if (z >= k2)
		    continue;
	    }

	    v_w_from__u1_bet(b);

	    if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
		break;
	}
	record_pair((aa == a)?a:w, (aa == a)?w:a);
	return (aa == a) ? a / (a + w) : w / (a + w);

    }
    else {		/* Algorithm BB */

	if (!qsame) { /* initialize */
	    beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
	    gamma = a + 1.0 / beta;
	}
	do {
	    u1 = unif_rand();
	    u2 = unif_rand();

	    v_w_from__u1_bet(a);

	    z = u1 * u1 * u2;
	    r = gamma * v - 1.3862944;
	    s = a + r - w;
	    if (s + 2.609438 >= 5.0 * z)
		break;
	    t = log(z);
	    if (s > t)
		break;
	}
	while (r + alpha * log(alpha / (b + w)) < t);

	record_pair((aa != a)?b:w, (aa != a)?w:b);
	return (aa != a) ? b / (b + w) : w / (b + w);
    }
}

SEXP C_np_plot_udens_bb_draw(SEXP aa, SEXP bb) {
 if(TYPEOF(aa)!=REALSXP || TYPEOF(bb)!=REALSXP)
   error("beta-binomial draw requires double shape vectors");
 R_xlen_t n=XLENGTH(aa);
 if(n<1 || n>1000000 || XLENGTH(bb)!=n)
   error("invalid beta-binomial draw dimensions or block budget");
 /* Validate the entire block before touching RNG state. */
 for(R_xlen_t i=0;i<n;i++) {
   double a=REAL(aa)[i],b=REAL(bb)[i];
   if(!R_FINITE(a) || !R_FINITE(b) || a<=0 || b<=0 ||
      !R_FINITE(a+b) || !R_FINITE(1.0/fmin(a,b)) || !R_FINITE(2*a*b))
     error("positive representable finite beta shapes required");
 }
 SEXP ans=PROTECT(allocMatrix(REALSXP,n,3));
 /* No allocating R API or interrupt/error operation occurs while RNG is held.
  * Blocks are bounded by the caller. Check the completed representation after
  * PutRNGstate, so a numerical refusal never loses the consumed RNG state. */
 GetRNGstate();
 for(R_xlen_t i=0;i<n;i++) {
   REAL(ans)[i]=observed_rbeta(REAL(aa)[i],REAL(bb)[i]);
   REAL(ans)[i+n]=lx_record;REAL(ans)[i+2*n]=lq_record;
 }
 PutRNGstate();
 for(R_xlen_t i=0;i<n;i++)
   if(!R_FINITE(REAL(ans)[i]) || !R_FINITE(REAL(ans)[i+n]) ||
      !R_FINITE(REAL(ans)[i+2*n]) || REAL(ans)[i+n]>=0 || REAL(ans)[i+2*n]>=0)
     error("beta-binomial draw has unrepresentable log coordinates");
 UNPROTECT(1);return ans;
}
