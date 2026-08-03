#ifndef NP_JKSUM_LP_BASIS_H
#define NP_JKSUM_LP_BASIS_H

#include <R_ext/Visibility.h>

typedef enum {
  NP_LP_BASIS_OK = 0,
  NP_LP_BASIS_INVALID,
  NP_LP_BASIS_NONFINITE,
  NP_LP_BASIS_RANK_DEFICIENT,
  NP_LP_BASIS_MEMORY,
  NP_LP_BASIS_LAPACK
} NPLPBasisStatus;

/*
 * Construct a numerically conditioned representation of an existing
 * full-column-rank training basis.  Source and destination are column-major
 * n-by-p matrices with a shared leading dimension.  Each source column is
 * first scaled by its maximum magnitude, then a pivoted Householder QR
 * produces an orthonormal basis for exactly the same column space.
 *
 * The operation is blind to the public basis family and never translates
 * coordinates by an evaluation point.  It is intended for basis-invariant
 * influence-row algebra.  Callers own both matrices; this function retains no
 * storage and allocates only O(p) temporary LAPACK work.
 */
attribute_hidden NPLPBasisStatus np_lp_conditioned_basis_prepare(
  double * const *source,
  int n,
  int p,
  int leading_dimension,
  double *destination);

/*
 * Decide once, before row traversal, whether the source basis carries a
 * representation-dependent unweighted Gram condition below min_rcond.
 * Non-finite source arithmetic requests conditioning so a finite canonical
 * representation can still be attempted.  Allocation or layout failures are
 * reported rather than silently selecting either representation.
 */
attribute_hidden NPLPBasisStatus np_lp_basis_requires_conditioning(
  double * const *source,
  int n,
  int p,
  int leading_dimension,
  double min_rcond,
  int *required);

#endif
