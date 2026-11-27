/* Keep the compiler and linker layout of the established fixed-width rows
 * independent of the tiled kernel. Each native row pays one width check. */
#include "jksum_lp_row.h"

#if NP_LP_TILE_WIDE_DONOR
#define NP_LP_COMPILE_TILED_DONOR 1
#include "jksum_lp_wide_donor.inc"

void np_lp_accumulate_dense_resident_row(const NPLPDenseRowContext *ctx)
{
  if (ctx->nterms > 16)
    np_lp_accumulate_dense_row_generic(ctx);
  else
    np_lp_accumulate_standard_dense_row(ctx);
}
#endif
