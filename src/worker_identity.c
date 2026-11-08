#include <Rinternals.h>
#include "np_worker_identity.h"

/* Literal compiled into the loaded DLL; never read the current disk image. */
SEXP C_np_worker_native_identity(void)
{
    return Rf_mkString(NP_WORKER_NATIVE_SOURCE);
}
