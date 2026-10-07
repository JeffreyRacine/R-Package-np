/* Rf_onintr is a package API declared through the graphics headers. Keep
 * those headers out of numerical translation units. Unlike the front-end
 * onintrNoResume entry point, this API is available on Windows too. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/GraphicsEngine.h>
#include "np_interrupt.h"

void np_raise_real_interrupt(void)
{
    Rf_onintr();
    /* A calling handler may select resume, or interrupts may be suspended.
     * The native owner is already released: never continue into its state. */
    Rf_error("the native search was interrupted and cannot be resumed");
}
