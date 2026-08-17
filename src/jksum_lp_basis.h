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
 * Decide once, before row traversal, whether the source basis carries a
 * representation-dependent unweighted Gram condition below min_rcond.
 * Source columns may share contiguous column-major storage or be allocated
 * independently; leading_dimension supplies the contiguous stride when that
 * layout is present.  Both layouts evaluate the same p-by-p Gram diagnostic.
 * Non-finite source arithmetic requests conditioning so a finite canonical
 * representation can still be attempted.  Allocation failures are reported
 * rather than silently selecting either representation.
 */
attribute_hidden NPLPBasisStatus np_lp_basis_requires_conditioning(
  double * const *source,
  int n,
  int p,
  int leading_dimension,
  double min_rcond,
  int *required);

#endif
