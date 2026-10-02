/* Conditional adapter to the canonical whole-support GNN integral.
 * Each coordinate overlap uses the deleted response geometry. Products of
 * overlaps and category sums commute with the finite donor-pair contraction;
 * no joint response grid or donor-by-donor matrix is materialized. */
typedef struct {
  NPGNNConditionalIntegralContext adapter;
  NPGNNIntegralOwner integral;
  double *cross;
  double score;
  int status;
} NPGNNConditionalGeneralCall;

static double np_cgnn_general_category(void *raw, int first, int second)
{
  NPGNNConditionalIntegralContext *a = raw;
  return np_gnn_integral_category_overlap(&a->categories, first, second);
}

static int np_cgnn_general_profiles(void *raw, int **ids, int **reps, int *count)
{
  NPGNNConditionalIntegralContext *a = raw;
  return np_gnn_integral_category_profiles(&a->categories, ids, reps, count);
}

static void np_cgnn_general_cleanup(void *raw, Rboolean jump)
{
  NPGNNConditionalGeneralCall *a = raw;
  np_gnn_integral_owner_cleanup(&a->integral, jump);
  free(a->integral.values);
  free(a->integral.bounds);
  free(a->adapter.xrow);
  free(a->adapter.yrow);
  free(a->cross);
}

static SEXP np_cgnn_general_body(void *raw)
{
  NPGNNConditionalGeneralCall *a = raw;
  NPGNNConditionalIntegralContext *b = &a->adapter;
  NPGNNIntegralOwner *c = &a->integral;
  NPConditionalYRowCtx *y = &b->route->legacy_y;
  const int n = num_obs_train_extern;
  const int categorical = num_var_unordered_extern + num_var_ordered_extern;
  int first = 0, count = n, fail;
  b->categories = (NPGNNIntegralUniformContext){
    .n=n, .nuno=num_var_unordered_extern, .nord=num_var_ordered_extern,
    .kernel_u=KERNEL_den_unordered_extern, .kernel_o=KERNEL_den_ordered_extern,
    .unordered=matrix_Y_unordered_train_extern,
    .ordered=matrix_Y_ordered_train_extern,
    .lambda=y->lambday, .num_categories=num_categories_extern_Y,
    .categories=matrix_categorical_vals_extern_Y};
#ifdef MPI2
  const int parallel = np_objective_outer_rows_enabled(1);
  np_objective_outer_owned_rows(0, n, parallel, &first, &count);
#endif
  *c = (NPGNNIntegralOwner){
    .n=n, .dimensions=num_var_continuous_extern, .folded=1,
    .kernel=KERNEL_den_extern, .has_categories=categorical > 0,
    .data=matrix_Y_continuous_train_extern, .counts=y->vsfy,
    .row=np_gnn_conditional_integral_weight,
    .category_overlap=categorical ? np_cgnn_general_category : NULL,
    .category_profiles=categorical ? np_cgnn_general_profiles : NULL,
    .context=b, .budget=NP_CONDITIONAL_LP_TILE_BUDGET_BYTES,
    .pair_ranks=1, .first_fold=first, .end_fold=first+count, .status=1};
  c->values=np_cgnn_calloc(n,sizeof(double));
  c->bounds=np_cgnn_calloc(n,sizeof(double));
  b->xrow=np_cgnn_calloc(n,sizeof(double));
  b->yrow=np_cgnn_calloc(n,sizeof(double));
  a->cross=np_cgnn_calloc(n,sizeof(double));
  fail=!c->values || !c->bounds || !b->xrow || !b->yrow || !a->cross;
  if(!fail) {
    np_gnn_integral_owner_body(c);
    fail=c->status;
#ifdef NP_CF167_TRACE
    REprintf("CF167 general rank=%d first=%d count=%d dims=%d kernel=%d representation=%d status=%d\n",
      np_cgnn_trace_rank(),first,count,c->dimensions,c->kernel,c->representation,fail);
    if(fail && c->geometry) {
      int iw[200];double work[800];int found=0;
      for(int d=0;d<c->dimensions && !found;++d) {
        NPGNNIntegralGeometry *g=c->geometry+d;
        for(int i=0;i<n && !found;++i)for(int j=0;j<=i && !found;++j)
          for(int at=0;at<g->count && !found;++at)for(int slot=0;slot<2;++slot) {
            NPGNNIntegralInterval t=g->intervals[at];double val,err;
            double anchor=slot?t.successor_anchor:t.primary_anchor;
            int rc=np_gnn_integral_pair_interval(g,c->kernel,t.lo,t.hi,anchor,
              c->data[d][i],c->data[d][j],iw,work,&val,&err);
            if(rc){REprintf("CF167 failed piece d=%d i=%d j=%d at=%d slot=%d lo=%.17g hi=%.17g anchor=%.17g rc=%d\n",d,i,j,at,slot,t.lo,t.hi,anchor,rc);found=1;break;}
          }
      }
    }
#endif
  }
#ifdef MPI2
  if(np_conditional_outer_buffer_finish(parallel,n,fail,c->values,NULL,
       "conditional GNN whole-support integral"))return R_NilValue;
  if(parallel)
    np_mpi_allreduce_in_place_double(c->bounds,n,MPI_SUM,
      "conditional GNN whole-support error estimate");
#else
  if(fail)return R_NilValue;
#endif
  /* The integral is in Y-tree order; the canonical cross term is in original
   * query order. Each map is applied exactly once in the final contraction. */
  for(int i=first;i<first+count && !fail;++i) {
    double logscale=0.0, linear;
    np_progress_bandwidth_loop_step();
    if(np_conditional_cvls_provider_x_row(b->route,i,b->xrow) ||
       np_conditional_cvls_provider_y_train_row(b->route,i,b->yrow,&logscale)) {
      fail=1;break;
    }
    linear=np_blas_ddot_int(n,b->xrow,b->yrow);
    if(np_continuous_kernel_scaled_restore(linear,logscale,1,&linear)!=
       NP_CONTINUOUS_ROW_OK) {fail=1;break;}
    a->cross[i]=linear;
  }
#ifdef MPI2
  if(np_conditional_outer_buffer_finish(parallel,n,fail,a->cross,NULL,
       "conditional GNN whole-support I2"))return R_NilValue;
#else
  if(fail)return R_NilValue;
#endif
  double sum=0.0, correction=0.0;
  for(int i=0;i<n;++i) {
    const int yi=int_TREE_Y==NP_TREE_TRUE ? ipt_lookup_extern_Y[i] : i;
    np_cgnn_compensated((c->values[yi]-2*a->cross[i])/n,&sum,&correction);
  }
  a->score=sum+correction;
#ifdef NP_CF167_TRACE
  double i1=0.0,i2=0.0;
  for(int j=0;j<n;++j){i1+=c->values[j]/n;i2+=a->cross[j]/n;}
  REprintf("CF167 general I1=%.17g I2=%.17g score=%.17g\n",i1,i2,a->score);
#endif
  if(R_FINITE(a->score))a->status=0;
  return R_NilValue;
}

static int np_cgnn_general_cvls(NPConditionalCVLSRouteContext *route,double *cv)
{
  NPGNNConditionalGeneralCall call={.adapter={.route=route},.status=1};
  R_UnwindProtect(np_cgnn_general_body,&call,np_cgnn_general_cleanup,&call,NULL);
  if(!call.status)*cv=call.score;
  return call.status;
}
