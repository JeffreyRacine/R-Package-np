/* Test-only, Unix event-loop injection; never linked into the package DLL. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/eventloop.h>
static void (*previous)(void);
static SEXP expression = NULL;
static int remaining;
static void poll(void) {
  if (previous) previous();
  if (remaining > 0 && --remaining == 0) Rf_eval(expression, R_GlobalEnv);
}
SEXP probe_arm(SEXP expr) {
  previous = R_PolledEvents;
  expression = expr;
  R_PreserveObject(expression);
  remaining = 10;
  R_PolledEvents = poll;
  return R_NilValue;
}
SEXP probe_off(void) {
  R_PolledEvents = previous;
  remaining = 0;
  if (expression != NULL) R_ReleaseObject(expression);
  expression = NULL;
  return R_NilValue;
}
SEXP probe_interrupt(void) {
  extern int R_interrupts_pending;
  R_interrupts_pending = 1;
  R_CheckUserInterrupt();
  return R_NilValue;
}
