/* Test-only, Unix event-loop injection; never linked into the package DLL. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <execinfo.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/eventloop.h>
static void (*previous)(void);
static SEXP expression = NULL;
static int remaining;
static void *package_base;
static uintptr_t target_lo, target_hi;
static int native_hits;
static int target_poll(void) {
  if (package_base == NULL) return 1;
  void *frames[64];
  int count = backtrace(frames, 64);
  for (int i = 1; i < count; ++i) {
    Dl_info info;
    if (dladdr(frames[i], &info) && info.dli_fbase == package_base) {
      uintptr_t address = (uintptr_t) frames[i] - 1;
      /* Only the first package frame identifies the owner that polled. An
       * outer beta frame while an inner owner polls is not a CF-234 witness. */
      return address >= target_lo && address < target_hi;
    }
  }
  return 0;
}
static void poll(void) {
  if (previous) previous();
  if (remaining > 0 && target_poll()) {
    ++native_hits;
    if (--remaining == 0) Rf_eval(expression, R_GlobalEnv);
  }
}
SEXP probe_arm(SEXP expr) {
  previous = R_PolledEvents;
  expression = expr;
  R_PreserveObject(expression);
  package_base = NULL;
  native_hits = 0;
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

/* A test-only symbol-range gate; no package instrumentation is installed. */
SEXP probe_arm_at(SEXP expr, SEXP anchor, SEXP range) {
  Dl_info info;
  if (!dladdr(R_ExternalPtrAddr(anchor), &info))
    Rf_error("cannot locate the package native library");
  probe_arm(expr);
  package_base = info.dli_fbase;
  target_lo = (uintptr_t) package_base + (uintptr_t) REAL(range)[0];
  target_hi = (uintptr_t) package_base + (uintptr_t) REAL(range)[1];
  remaining = 3;
  return R_NilValue;
}
SEXP probe_native_hits(void) { return Rf_ScalarInteger(native_hits); }
