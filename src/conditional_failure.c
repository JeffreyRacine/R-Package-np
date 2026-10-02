#include <R.h>
#include <stdint.h>
#include <stdio.h>
#include "conditional_failure.h"
#ifdef MPI2
#include <mpi.h>
extern MPI_Comm *comm;
#endif
static unsigned long long npc_failure;
static char npc_bandwidth[2048];
void np_conditional_failure_reset(void) {
  npc_failure=0; npc_bandwidth[0]='\0';
}
void np_conditional_failure_record(int reason,int row) {
  if(!npc_failure) npc_failure=((unsigned long long)(unsigned)reason<<32) |
    (unsigned)(row>0?row:0);
}
int np_conditional_failure_pending(void) { return npc_failure!=0; }
int np_conditional_failure_reduce(int parallel,int failed) {
  unsigned long long value=npc_failure ? npc_failure : (failed?1:0);
#ifdef MPI2
  if(parallel) {
    unsigned long long all=0;
    MPI_Allreduce(&value,&all,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm[1]);
    value=all;
  }
#else
  (void)parallel;
#endif
  if(value>1) npc_failure=value;
  return value!=0;
}
void np_conditional_failure_bandwidth(const double *bandwidth,int count) {
  if(!npc_failure || npc_bandwidth[0]) return;
  size_t used=0;
  for(int i=1;i<=count && used<sizeof(npc_bandwidth)-1;i++) {
    const int written=snprintf(npc_bandwidth+used,sizeof(npc_bandwidth)-used,
      "%s%.17g",i==1?"":", ",bandwidth[i]);
    if(written<0)break;
    if((size_t)written>=sizeof(npc_bandwidth)-used){
      const size_t last=sizeof(npc_bandwidth)-5;
      snprintf(npc_bandwidth+last,5," ...");break;
    }
    used+=(size_t)written;
  }
}
void np_conditional_failure_raise(void) {
  if(!npc_failure)return;
  const unsigned reason=(unsigned)(npc_failure>>32),row=(unsigned)npc_failure;
  const char *label=reason==NP_CONDITIONAL_RANK_AMBIGUOUS ? "ambiguous numerical rank" :
    reason==NP_CONDITIONAL_WORK_EXHAUSTED ? "integration work budget exhausted" :
    reason==NP_CONDITIONAL_COEFFICIENT_FAILURE ? "polynomial coefficient reconstruction failed accuracy check" :
    "conditional local solve failed";
  error("conditional bandwidth search stopped: %s; row %u; bandwidth/scale factors (native coordinate order) [%s]",
    label,row,npc_bandwidth[0]?npc_bandwidth:"unavailable");
}
