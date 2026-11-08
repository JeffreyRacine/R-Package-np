/* Public R RNG state access only. No draws or private RNG symbols. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

/* Called only with no .Random.seed binding, under R-side cleanup. Unlike
 * GetRNGstate/RNGkind, the writer does not initialize an absent state. */
SEXP C_np_rng_scope_materialize(void)
{
    PutRNGstate();
    return R_NilValue;
}

/* The R owner has supplied a previously captured, validated built-in state.
 * Load it so a temporary state swap cannot leave its kind active after the
 * owner removes the originally absent binding. No random values are drawn. */
SEXP C_np_rng_scope_restore_absent(void)
{
    GetRNGstate();
    return R_NilValue;
}
