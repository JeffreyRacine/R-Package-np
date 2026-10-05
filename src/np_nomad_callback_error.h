/* Callback-internal R unwind containment. The caller roots the continuation cell
 * until all native invocation owners have been released. This is deliberately
 * not an MPI error protocol: distributed callbacks retain their existing path. */
#ifndef NP_NOMAD_CALLBACK_ERROR_H
#define NP_NOMAD_CALLBACK_ERROR_H

#include <setjmp.h>

typedef struct {
  crs_nomad_eval_fn eval;
  void *user_data;
  SEXP error_state;
  crs_nomad_observer *observer;
  crs_nomad_observer_poll_fn poll;
  crs_nomad_observe_fn observe;
  void *observe_data;
  char progress_error[256];
} NPNomadCallbackError;

typedef struct {
  NPNomadCallbackError *owner;
  int n, m, status;
  const double *x;
  double *outputs;
  jmp_buf jump;
} NPNomadCallbackCall;

/* Allocate before native resources. Slot zero holds the first failure: either
 * the slot-one objective continuation or the original observer condition.
 * Slot two roots observer catch classes before any native invocation. */
static SEXP np_nomad_unwind_state(void)
{
  SEXP state = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(state, 1, R_MakeUnwindCont());
  SET_VECTOR_ELT(state, 2, Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(VECTOR_ELT(state, 2), 0, Rf_mkChar("interrupt"));
  SET_STRING_ELT(VECTOR_ELT(state, 2), 1, Rf_mkChar("error"));
  UNPROTECT(1);
  return state;
}

static int np_nomad_error_pending(SEXP state)
{
  return state != R_NilValue && VECTOR_ELT(state, 0) != R_NilValue;
}

static SEXP np_nomad_callback_body(void *data)
{
  NPNomadCallbackCall *call = data;
  call->status = call->owner->eval(call->n, call->x, call->m,
                                 call->outputs, call->owner->user_data);
  return R_NilValue;
}

static void np_nomad_callback_unwind(void *data, Rboolean jump)
{
  NPNomadCallbackCall *call = data;
  if (jump) {
    if (!np_nomad_error_pending(call->owner->error_state))
      SET_VECTOR_ELT(call->owner->error_state, 0,
                     VECTOR_ELT(call->owner->error_state, 1));
    /* The target is inside this C callback, below every NOMAD C++ frame.
     * Return normally through NOMAD before the caller resumes R's unwind. */
    longjmp(call->jump, 1);
  }
}

static int np_nomad_callback_contained(int n, const double *x, int m,
                                      double *outputs, void *data)
{
  NPNomadCallbackError *owner = data;
  /* A stop request may still be followed by another callback. Do not enter R
   * again or overwrite the suspended continuation before native cleanup. */
  if (np_nomad_error_pending(owner->error_state))
    return 1;
  NPNomadCallbackCall call = {
    .owner = owner, .n = n, .m = m, .status = 1, .x = x, .outputs = outputs
  };
  if (setjmp(call.jump) == 0)
    R_UnwindProtect(np_nomad_callback_body, &call,
                    np_nomad_callback_unwind, &call,
                    VECTOR_ELT(owner->error_state, 1));
  if (np_nomad_error_pending(owner->error_state)) {
    /* Synchronous polling is essential: otherwise crs treats callback failure
     * as an invalid evaluation and keeps searching. */
    owner->observer->interval_sec = 0.0;
    owner->poll();
    return 1;
  }
  /* call.status can change during evaluation. Read it only on ordinary
   * return, never after longjmp (the pending branch above always returns). */
  return call.status;
}

typedef struct {
  NPNomadCallbackError *owner;
  const crs_nomad_observer_event *event;
  char *message;
  size_t message_size;
  int status;
} NPNomadErrorProgressCall;

static SEXP np_nomad_error_progress_body(void *data)
{
  NPNomadErrorProgressCall *call = data;
  call->status = call->owner->observe(call->event, call->owner->observe_data,
                                     call->message, call->message_size);
  return R_NilValue;
}

static SEXP np_nomad_error_progress_condition(SEXP condition, void *data)
{
  NPNomadErrorProgressCall *call = data;
  /* Root the condition before requesting stop. Do not retain an unwind whose
   * target is crs's observer-only R_ToplevelExec context. */
  if (!np_nomad_error_pending(call->owner->error_state))
    SET_VECTOR_ELT(call->owner->error_state, 0, condition);
  return R_NilValue;
}

static int np_nomad_error_observer(const crs_nomad_observer_event *event,
                                   void *data, char *message, size_t size)
{
  NPNomadCallbackError *owner = data;
  if (np_nomad_error_pending(owner->error_state))
    return CRS_NOMAD_OBSERVER_OUTCOME_INTERRUPT;
  if (owner->observe != NULL) {
    NPNomadErrorProgressCall call = {
      owner, event, message, size, CRS_NOMAD_OBSERVER_OUTCOME_ERROR
    };
    R_tryCatch(np_nomad_error_progress_body, &call,
               VECTOR_ELT(owner->error_state, 2),
               np_nomad_error_progress_condition, &call, NULL, NULL);
    if (np_nomad_error_pending(owner->error_state))
      return CRS_NOMAD_OBSERVER_OUTCOME_INTERRUPT;
    if (call.status == CRS_NOMAD_OBSERVER_OUTCOME_INTERRUPT)
      return call.status;
    if (call.status != CRS_NOMAD_OBSERVER_OUTCOME_OK) {
      /* A broken progress reporter must not disable native error transport. */
      snprintf(owner->progress_error, sizeof(owner->progress_error), "%s",
               message != NULL && message[0] != '\0' ? message :
               "native NOMAD progress dispatcher failed");
      owner->observe = NULL;
      owner->observer->interval_sec = DBL_MAX;
    }
  }
  return CRS_NOMAD_OBSERVER_OUTCOME_OK;
}

static void np_nomad_error_observer_init(NPNomadCallbackError *owner,
                                         crs_nomad_observer *observer,
                                         crs_nomad_observer_poll_fn poll,
                                         crs_nomad_eval_fn eval, void *data,
                                         SEXP error_state)
{
  memset(owner, 0, sizeof(*owner));
  owner->eval = eval;
  owner->user_data = data;
  owner->error_state = error_state;
  owner->observer = observer;
  owner->poll = poll;
  owner->observe = observer->observe;
  owner->observe_data = observer->user_data;
  if (observer->observe == NULL)
    observer->interval_sec = DBL_MAX;
  observer->api_version = CRS_NOMAD_OBSERVER_API_VERSION;
  observer->struct_size = sizeof(*observer);
  observer->observe = np_nomad_error_observer;
  observer->user_data = owner;
}

static void np_nomad_error_raise(SEXP state)
{
  /* All native owners are gone and globals restored. Objective failures keep
   * their original continuation; observer failures re-signal the saved object
   * in the caller context, without claiming the original restart context. */
  if (VECTOR_ELT(state, 0) == VECTOR_ELT(state, 1))
    R_ContinueUnwind(VECTOR_ELT(state, 1));
  if (Rf_inherits(VECTOR_ELT(state, 0), "interrupt")) {
    /* A real R interrupt has no message. signalCondition preserves its class
     * for caller handlers; abort supplies interruption's default control flow
     * when no exiting handler takes it. stop(condition) can turn it into an
     * ordinary error that try() silently catches. */
    SEXP signal = PROTECT(Rf_lang2(Rf_install("signalCondition"),
                                   VECTOR_ELT(state, 0)));
    Rf_eval(signal, R_BaseEnv);
    SEXP restart = PROTECT(Rf_mkString("abort"));
    SEXP abort = PROTECT(Rf_lang2(Rf_install("invokeRestart"), restart));
    Rf_eval(abort, R_BaseEnv);
    UNPROTECT(3);
  }
  SEXP call = PROTECT(Rf_lang2(Rf_install("stop"), VECTOR_ELT(state, 0)));
  Rf_eval(call, R_BaseEnv);
  UNPROTECT(1);
}
#endif
