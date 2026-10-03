/* Callback-internal R error containment. The caller roots the condition cell
 * until all native invocation owners have been released. This is deliberately
 * not an MPI error protocol: distributed callbacks retain their existing path. */
#ifndef NP_NOMAD_CALLBACK_ERROR_H
#define NP_NOMAD_CALLBACK_ERROR_H

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
} NPNomadCallbackCall;

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

static SEXP np_nomad_callback_caught(SEXP condition, void *data)
{
  NPNomadCallbackCall *call = data;
  SET_VECTOR_ELT(call->owner->error_state, 0, condition);
  return R_NilValue;
}

static int np_nomad_callback_contained(int n, const double *x, int m,
                                      double *outputs, void *data)
{
  NPNomadCallbackError *owner = data;
  NPNomadCallbackCall call = {owner, n, m, 1, x, outputs};
  R_tryCatchError(np_nomad_callback_body, &call, np_nomad_callback_caught, &call);
  if (np_nomad_error_pending(owner->error_state)) {
    /* Synchronous polling is essential: otherwise crs treats callback failure
     * as an invalid evaluation and keeps searching. */
    owner->observer->interval_sec = 0.0;
    owner->poll();
    return 1;
  }
  return call.status;
}

typedef struct {
  NPNomadCallbackError *owner;
  const crs_nomad_observer_event *event;
  char *message;
  size_t message_size;
  int status;
} NPNomadErrorProgressCall;

static void np_nomad_error_progress_body(void *data)
{
  NPNomadErrorProgressCall *call = data;
  call->status = call->owner->observe(call->event, call->owner->observe_data,
                                     call->message, call->message_size);
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
    int ok = R_ToplevelExec(np_nomad_error_progress_body, &call);
    if (ok && call.status == CRS_NOMAD_OBSERVER_OUTCOME_INTERRUPT)
      return call.status;
    if (!ok || call.status != CRS_NOMAD_OBSERVER_OUTCOME_OK) {
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
  /* Native R errors carry condition objects; R's stop(condition) is needed to
   * re-signal the original class, call and fields, rather than a new string. */
  SEXP call = PROTECT(Rf_lang2(Rf_install("stop"), VECTOR_ELT(state, 0)));
  Rf_eval(call, R_BaseEnv);
  UNPROTECT(1);
}
#endif
