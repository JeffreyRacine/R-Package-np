/* Test-only real interrupt injection through R's cross-platform package API.
 * Delivers R's empty interrupt condition; no private pending flags or signals. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/GraphicsEngine.h>
SEXP probe_interrupt(void)
{
    Rf_onintr();
    return R_NilValue;
}
