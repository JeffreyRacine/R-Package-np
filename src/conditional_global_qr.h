#ifndef NP_CONDITIONAL_GLOBAL_QR_H
#define NP_CONDITIONAL_GLOBAL_QR_H
/* Unregistered integration component. No selector or package wiring here.
 * Caller must retain the incumbent all-large admission before using this.
 * Canonical original-coordinate policy; global Q only evaluates admitted rows. */
#include "conditional_deleted_qr.h"
#include "conditional_rank_admission.h"
#include "conditional_global_rank_certificate.h"
#include "conditional_local_qr.h"
#include "conditional_failure.h"
typedef struct {
  int n,k,dgemv,nexception;
  double **basis; /* borrowed original-coordinate columns */
  double *q,*h,*lambda,*anchor,*ones,*row;
  int *exception,*indices;
  NPConditionalQRDeleted deleted;
} NPConditionalQRGlobal;
static void np_cqr_clear(NPConditionalQRGlobal *g){
  free(g->q);free(g->h);free(g->lambda);free(g->anchor);free(g->ones);
  free(g->row);free(g->exception);free(g->indices);
  np_cqr_deleted_clear(&g->deleted);memset(g,0,sizeof(*g));
}
static int np_cqr_prepare(NPConditionalQRGlobal *g,int n,int k,double **basis,int dgemv,
  int dimensions,const int *terms,double **x,const int *positions,const int *original_ids){
  NPLPSolveWorkspace policy;np_lp_solve_workspace_init(&policy);
  NPConditionalQRFast local={0};
  NPConditionalQRAcc *gram=NULL;double *tau=NULL,*work=NULL;int *pivot=NULL;
  int one=1,info=0,lwork=-1,status=1;double query=0.;
  if(!g||!basis||n<2||k<1||k>n||k>(INT_MAX-1)/3||
     (size_t)n>SIZE_MAX/(size_t)k/sizeof(double)||
     (size_t)k>SIZE_MAX/(size_t)k/sizeof(NPConditionalQRAcc))return 1;
  np_cqr_clear(g);g->n=n;g->k=k;g->basis=basis;g->dgemv=dgemv;
  g->q=calloc((size_t)n*k,sizeof(double));g->h=calloc(n,sizeof(double));
  g->lambda=calloc(n,sizeof(double));g->anchor=calloc(n,sizeof(double));
  g->ones=calloc(n,sizeof(double));g->row=calloc(n,sizeof(double));
  g->exception=calloc(n,sizeof(int));g->indices=calloc(n,sizeof(int));
  gram=calloc((size_t)k*k,sizeof(NPConditionalQRAcc));tau=calloc(k,sizeof(double));pivot=calloc(k,sizeof(int));
  if(!g->q||!g->h||!g->lambda||!g->anchor||!g->ones||!g->row||!g->exception||!g->indices||!gram||!tau||!pivot)goto cleanup;
  for(int a=0;a<k;a++){
    if(!basis[a])goto cleanup;
    for(int j=0;j<n;j++)if(!isfinite(basis[a][j]))goto cleanup;
    double scale=F77_CALL(dnrm2)(&n,basis[a],&one);
    if(!isfinite(scale)||scale==0)goto cleanup;
    for(int j=0;j<n;j++)g->q[j+(size_t)n*a]=basis[a][j]/scale;
  }
  for(int j=0;j<n;j++){
    g->ones[j]=1.;
    for(int a=0;a<k;a++)for(int b=a;b<k;b++){
      np_cqr_product(gram+a+(size_t)k*b,basis[a][j],basis[b][j]);
      if(!isfinite(gram[a+(size_t)k*b].hi)||!isfinite(gram[a+(size_t)k*b].lo))goto cleanup;
    }
  }
  if(!np_lp_solve_workspace_reserve(&policy,k,1))goto cleanup;
  for(int i=0;i<n;i++){
    for(int a=0;a<k;a++){
      policy.rhs_source[a]=basis[a][i];
      for(int b=a;b<k;b++){
        NPConditionalQRAcc z=gram[a+(size_t)k*b];np_cqr_product(&z,-basis[a][i],basis[b][i]);
        double value=z.hi+z.lo;if(!isfinite(value))goto cleanup;
        policy.gram_source[a+(size_t)k*b]=policy.gram_source[b+(size_t)k*a]=value;
      }
    }
    g->anchor[i]=policy.gram_source[0];
    NPLPSolvePolicyDiagnostics d={0,0.};
    const int position=positions?positions[i]:i;
    const NPConditionalLocalQRStatus rank=np_cqr_local_row(&local,&policy,n,
      dimensions,k,terms,x,g->ones,position,g->row,NULL);
    if(rank==NP_CQR_LOCAL_FULL) {
      g->lambda[i]=0.0;
      continue;
    }
    if(rank!=NP_CQR_LOCAL_DEFICIENT) {
      np_conditional_failure_record(rank==NP_CQR_LOCAL_AMBIGUOUS ?
        NP_CONDITIONAL_RANK_AMBIGUOUS : NP_CONDITIONAL_NUMERICAL_FAILURE,
        1+(original_ids?original_ids[i]:i));
      goto cleanup;
    }
    if(np_conditional_solve_adjoint_ranked(&policy,k,1,1./n,0,&d)
         !=NP_LP_SOLVE_POLICY_OK)goto cleanup;
    g->lambda[i]=d.ridge_total;
  }
  F77_CALL(dgeqp3)(&n,&k,g->q,&n,pivot,tau,&query,&lwork,&info);
  if(info||!isfinite(query)||query<1||query>INT_MAX)goto cleanup;
  lwork=(int)fmax(query,3*k+1);work=calloc(lwork,sizeof(double));if(!work)goto cleanup;
  F77_CALL(dgeqp3)(&n,&k,g->q,&n,pivot,tau,work,&lwork,&info);if(info)goto cleanup;
  for(int a=0;a<k;a++)if(g->q[a+(size_t)n*a]==0)goto cleanup;
  int lq=-1;double qw=0.;
  F77_CALL(dorgqr)(&n,&k,&k,g->q,&n,tau,&qw,&lq,&info);
  if(info||!isfinite(qw)||qw<1||qw>INT_MAX)goto cleanup;
  if(qw>lwork){free(work);lwork=(int)qw;work=calloc(lwork,sizeof(double));if(!work)goto cleanup;}
  F77_CALL(dorgqr)(&n,&k,&k,g->q,&n,tau,work,&lwork,&info);if(info)goto cleanup;
  for(int i=0;i<n;i++){
    g->h[i]=np_cqr_dot4(g->q+i,n,g->q+i,n,k);
    if(!isfinite(g->h[i])||g->h[i]<0)goto cleanup;
    g->exception[i]=(g->h[i]>.5||g->lambda[i]>0);
    if(g->exception[i])g->indices[g->nexception++]=i;
  }
  status=0;
cleanup:
  np_cqr_fast_clear(&local);
  np_lp_solve_workspace_clear(&policy);free(gram);free(tau);free(pivot);free(work);
  if(status)np_cqr_clear(g);return status;
}
static int np_cqr_row(NPConditionalQRGlobal *g,int i,double *row){
  if(!g||i<0||i>=g->n||!row)return 1;
  if(g->exception[i])return np_cqr_deleted_row(&g->deleted,g->n,g->k,g->ones,g->basis,i,g->lambda[i],g->anchor[i],row,NULL,NULL);
  const double den=1.-g->h[i];
  for(int j=0;j<g->n;j++){
    row[j]=i==j?0.:np_cqr_dot4(g->q+i,g->n,g->q+j,g->n,g->k)/den;
    if(!isfinite(row[j]))return 1;
  }
  return 0;
}
static int np_cqr_cross(const NPConditionalQRGlobal *g,const double *rhs,double *cross){
  int one=1;double alpha=1.,beta=0.;
  if(g->dgemv)F77_CALL(dgemv)("T",&g->n,&g->k,&alpha,g->q,&g->n,rhs,&one,&beta,cross,&one FCONE);
  else for(int a=0;a<g->k;a++)cross[a]=F77_CALL(ddot)(&g->n,g->q+(size_t)g->n*a,&one,rhs,&one);
  for(int a=0;a<g->k;a++)if(!isfinite(cross[a]))return 1;return 0;
}
static int np_cqr_linear(NPConditionalQRGlobal *g,int i,const double *rhs,const double *cross,double *fit){
  if(!g||i<0||i>=g->n||!rhs||!cross||!fit)return 1;
  if(g->exception[i]){
    if(np_cqr_row(g,i,g->row))return 1;
    *fit=np_cqr_dot4(g->row,1,rhs,1,g->n);
  }else *fit=(np_cqr_dot4(g->q+i,g->n,cross,1,g->k)-g->h[i]*rhs[i])/(1.-g->h[i]);
  return !isfinite(*fit);
}
typedef int (*NPConditionalQRResponseRow)(void *,int,double *);
/* Tile limit is an administrative memory bound, not an arithmetic selector.
 * Production uses zero; a smaller value is exposed only by the test adapter. */
static int np_cqr_quadratic(NPConditionalQRGlobal *g,NPConditionalQRResponseRow response,void *context,
                        double *values,int test_tile_limit,int *passes_out){
  const int n=g->n,k=g->k;int status=1,passes=0;
  double *rhs=NULL,*cross=NULL,*diag=NULL,*conv=NULL,*tile=NULL;
  NPConditionalQRAcc *moment=NULL,*acc=NULL;
  size_t rowbytes=(size_t)n*sizeof(double),limit=((size_t)1<<20)/rowbytes;
  if(limit<1)limit=1;if(limit>16)limit=16;if(limit>(size_t)k)limit=k;
  if(test_tile_limit>0&&(size_t)test_tile_limit<limit)limit=test_tile_limit;
  int capacity=(int)limit;
  if((size_t)n>SIZE_MAX/limit/sizeof(double))return 1;
  rhs=calloc(n,sizeof(double));cross=calloc(k,sizeof(double));diag=calloc(n,sizeof(double));
  conv=calloc(n,sizeof(double));moment=calloc((size_t)k*k,sizeof(NPConditionalQRAcc));
  tile=calloc((size_t)n*capacity,sizeof(double));acc=calloc(capacity,sizeof(NPConditionalQRAcc));
  if(!rhs||!cross||!diag||!conv||!moment||!tile||!acc)goto cleanup;
  /* First pass supplies ordinary moments and the first exception tile. */
  int first=1;
  for(int start=0;first||start<g->nexception;start+=capacity){
    int count=g->nexception-start;if(count>capacity)count=capacity;
    memset(acc,0,(size_t)capacity*sizeof(NPConditionalQRAcc));
    for(int e=0;e<count;e++)if(np_cqr_row(g,g->indices[start+e],tile+(size_t)n*e))goto cleanup;
    for(int j=0;j<n;j++){
      if(response(context,j,rhs))goto cleanup;
      for(int z=0;z<n;z++)if(!isfinite(rhs[z]))goto cleanup;
      if(first){
        if(np_cqr_cross(g,rhs,cross))goto cleanup;
        conv[j]=np_cqr_dot4(g->q+j,n,cross,1,k);diag[j]=rhs[j];
        for(int a=0;a<k;a++)for(int b=0;b<k;b++)np_cqr_product(moment+a+(size_t)k*b,g->q[j+(size_t)n*a],cross[b]);
      }
      for(int e=0;e<count;e++){
        double *row=tile+(size_t)n*e;
        np_cqr_product(acc+e,row[j],np_cqr_dot4(rhs,1,row,1,n));
      }
    }
    for(int e=0;e<count;e++)values[g->indices[start+e]]=acc[e].hi+acc[e].lo;
    passes++;first=0;
  }
  for(int i=0;i<n;i++)if(!g->exception[i]){
    NPConditionalQRAcc z={0,0};
    for(int a=0;a<k;a++)for(int b=0;b<k;b++)np_cqr_product(&z,g->q[i+(size_t)n*a]*(moment[a+(size_t)k*b].hi+moment[a+(size_t)k*b].lo),g->q[i+(size_t)n*b]);
    const double h=g->h[i],den=1.-h;
    np_cqr_product(&z,-2.*h,conv[i]);np_cqr_product(&z,h*h,diag[i]);
    values[i]=(z.hi+z.lo)/(den*den);
  }
  for(int i=0;i<n;i++)if(!isfinite(values[i]))goto cleanup;
  if(passes_out)*passes_out=passes;status=0;
cleanup:
  free(rhs);free(cross);free(diag);free(conv);free(moment);free(tile);free(acc);return status;
}

/* Each invocation owns exactly one response-grid interval. The callback must
 * not introduce collectives into an independently owned MPI interval. */
typedef int (*NPConditionalQRIndicator)(void *, int, int);
static int np_cqr_grid_loss(NPConditionalQRGlobal *g, NPConditionalQRResponseRow response,
                       NPConditionalQRIndicator indicator, void *context,
                       int start, int count, double *loss,
                       int test_tile_limit, int *passes_out) {
  if (!g || !g->n || !response || !indicator || !loss || start < 0 || count < 0 ||
      start > INT_MAX-count) return 1;
  if (count == 0) { if (passes_out) *passes_out=0; return 0; }
  const int n=g->n, k=g->k;
  size_t limit=(((size_t)1<<20)/sizeof(double))/(size_t)n;
  if (limit<1) limit=1;
  if (limit>16) limit=16;
  if (limit>(size_t)k) limit=k;
  if (test_tile_limit>0 && limit>(size_t)test_tile_limit) limit=test_tile_limit;
  if ((size_t)n>SIZE_MAX/limit/sizeof(double)) return 1;
  double *tile=calloc((size_t)n*limit,sizeof(double));
  double *rhs=calloc(n,sizeof(double)), *cross=calloc(k,sizeof(double));
  NPConditionalQRAcc *sums=calloc(count,sizeof(NPConditionalQRAcc));
  int *slot=calloc(n,sizeof(int));
  int status=1, passes=0, first=1;
  if (!tile || !rhs || !cross || !sums || !slot) goto cleanup;
  for (int begin=0; first || begin<g->nexception; begin+=(int)limit) {
    int entries=g->nexception-begin;
    if (entries>(int)limit) entries=(int)limit;
    for (int i=0;i<n;i++) slot[i]=-1;
    for (int e=0;e<entries;e++) {
      const int i=g->indices[begin+e];slot[i]=e;
      if(np_cqr_row(g,i,tile+(size_t)n*e)) goto cleanup;
    }
    for (int j=0;j<count;j++) {
      if(response(context,start+j,rhs)) goto cleanup;
      for(int i=0;i<n;i++) if(!isfinite(rhs[i])) goto cleanup;
      if(first && np_cqr_cross(g,rhs,cross)) goto cleanup;
      for(int i=0;i<n;i++) {
        if ((g->exception[i] && slot[i]<0) || (!g->exception[i] && !first)) continue;
        const int target=indicator(context,i,start+j);
        if(target<0) continue;
        if(target>1) goto cleanup;
        double fit;
        if(g->exception[i]) fit=np_cqr_dot4(tile+(size_t)n*slot[i],1,rhs,1,n);
        else if(np_cqr_linear(g,i,rhs,cross,&fit)) goto cleanup;
        if(!isfinite(fit)) goto cleanup;
        const double difference=(double)target-fit;
        np_cqr_product(sums+j,difference,difference);
      }
    }
    passes++;first=0;
  }
  for(int j=0;j<count;j++) {
    loss[j]=sums[j].hi+sums[j].lo;
    if(!isfinite(loss[j])) goto cleanup;
  }
  if(passes_out) *passes_out=passes;
  status=0;
cleanup:
  free(tile);free(rhs);free(cross);free(sums);free(slot);return status;
}

#endif
