#ifndef NP_CONDITIONAL_FAILURE_H
#define NP_CONDITIONAL_FAILURE_H
/* Deferred metadata is raised only after the conditional evaluator completes. */
enum { NP_CONDITIONAL_RANK_AMBIGUOUS=1, NP_CONDITIONAL_NUMERICAL_FAILURE=2,
       NP_CONDITIONAL_WORK_EXHAUSTED=3, NP_CONDITIONAL_COEFFICIENT_FAILURE=4 };
typedef struct {
  unsigned long long failure;
  char bandwidth[2048];
} NPConditionalFailure;
/* Invocation-owned copy survives prepared-context and search cleanup. */
void np_conditional_failure_save(NPConditionalFailure *saved);
void np_conditional_failure_restore(const NPConditionalFailure *saved);
void np_conditional_failure_reset(void);
void np_conditional_failure_record(int reason,int row);
int np_conditional_failure_pending(void);
int np_conditional_failure_reduce(int parallel,int failed);
void np_conditional_failure_bandwidth(const double *bandwidth,int count);
void np_conditional_failure_raise(void);
#endif
