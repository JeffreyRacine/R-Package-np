static void np_conditional_plot_capture_accepted(
  NPLPSolveWorkspace *workspace, int width, int nrhs)
{
  NPConditionalPlotCapture *out = np_conditional_plot_capture_current;
  if(out == NULL || out->failed) return;
  NPRegressionGeneralLPFitOwner *owner = (NPRegressionGeneralLPFitOwner *)
    ((char *)workspace - offsetof(NPRegressionGeneralLPFitOwner, solve_workspace));
  NPRegressionGeneralLPFitExecution *execution = (NPRegressionGeneralLPFitExecution *)
    ((char *)owner - offsetof(NPRegressionGeneralLPFitExecution, owner));
  const NPRegressionGeneralLPFitCall *call = execution->call;
  if(call == NULL || owner->nterms != width || out->width != width ||
     call->num_obs_train != out->n || call->num_obs_eval != 1 ||
     call->ordinary_response || call->hc0_context != NULL ||
     call->fit_projection == NULL ||
     call->fit_projection->donor_map != out->donor_map ||
     call->fit_projection->eval_map == NULL ||
     call->fit_projection->eval_map[0] != out->row ||
     out->row < 0 || out->row >= out->m || out->donor_map == NULL ||
     out->slots != call->num_reg_continuous + 1 || nrhs > out->slots) {
    out->failed = 1; return;
  }
  /* First-SE-only blocks do not replace the complete level-first block. */
  for(int k=0;k<width;++k)
    if(workspace->rhs_source[k] != owner->eval_basis[k]) return;
  int expected=1;
  if(call->do_grad && call->do_gerr)
    for(int d=0;d<call->num_reg_continuous;++d)
      expected += np_glp_gradient_direction_active(d) != 0;
  if(nrhs != expected) { out->failed=1; return; }
  const size_t offset=(size_t)out->row*out->slots;
  int slot=0;
  out->coordinate[offset]=0; out->order[offset]=0;
  if(call->do_grad && call->do_gerr)
    for(int d=0;d<call->num_reg_continuous;++d)
      if(np_glp_gradient_direction_active(d)) {
        ++slot;
        out->coordinate[offset+slot]=d+1;
        out->order[offset+slot]=vector_glp_gradient_order_extern != NULL ?
          MAX(1,vector_glp_gradient_order_extern[d]) : 1;
      }
  for(int r=0;r<nrhs;++r)
    for(int k=0;k<width;++k) {
      out->directions[(offset+r)*width+k]=workspace->rhs_source[(size_t)r*width+k];
      out->adjoints[(offset+r)*width+k]=workspace->rhs_work[(size_t)r*width+k];
    }
  int status=0;
  if(owner->ann_kernel_row != NULL) {
    memcpy(out->scratch_weights,owner->ann_kernel_row,(size_t)out->n*sizeof(double));
  } else if(call->kernel_route != NULL) {
    /* Match the accepted beta response/basis scaling, using request scratch. */
    const int stride=width+2;
    NPContinuousKernelDerivativeDiagnostics replay_diagnostics={0};
    out->scratch_columns[0]=owner->squared_response;
    out->scratch_columns[1]=call->vector_Y;
    for(int k=0;k<width;++k) out->scratch_columns[k+2]=owner->basis[k];
    status=np_beta_regression_lp_moment_row_canonical(
      call->bandwidth_mode,out->n,call->num_reg_unordered,
      call->num_reg_ordered,call->num_reg_continuous,width,stride,
      call->operator,call->kernel_u,call->kernel_o,
      call->matrix_X_unordered_train,call->matrix_X_ordered_train,
      call->matrix_X_continuous_train,owner->eval_unordered,
      owner->eval_ordered,owner->eval_continuous,out->scratch_columns,
      owner->basis,call->vector_scale_factor,call->matrix_bandwidth,
      call->bandwidth_mode==BW_GEN_NN ? owner->matrix_bandwidth_eval : call->matrix_bandwidth,
      call->lambda,call->num_categories,call->matrix_categorical_vals,
      call->categorical_compress,out->scratch_moments,NULL,NULL,
      out->scratch_weights,call->kernel_route,&replay_diagnostics,NULL);
  } else {
    status=kernel_weighted_sum_np_ctx_ex(
      call->kernel_c,call->kernel_u,call->kernel_o,call->bandwidth_mode,
      out->n,1,call->num_reg_unordered,call->num_reg_ordered,
      call->num_reg_continuous,0,0,1,1,1,0,0,0,0,
      call->operator,OP_NOOP,0,0,NULL,1,0,0,
      call->bandwidth_mode==BW_ADAP_NN ? NP_TREE_FALSE : call->tree_enabled,
      0,call->bandwidth_mode==BW_ADAP_NN ? NULL : kdt_extern_X,
      NULL,NULL,NULL,call->matrix_X_unordered_train,
      call->matrix_X_ordered_train,call->matrix_X_continuous_train,
      owner->eval_unordered,owner->eval_ordered,owner->eval_continuous,
      NULL,NULL,NULL,call->vector_scale_factor,1,call->matrix_bandwidth,
      call->bandwidth_mode==BW_GEN_NN ? owner->matrix_bandwidth_eval : call->matrix_bandwidth,
      call->lambda,call->num_categories,call->matrix_categorical_vals,
      NULL,NULL,NULL,out->scratch_weights,call->gate_context,NULL,
      NULL,NULL,NULL,
#ifdef MPI2
      0,
#endif
      NULL,NULL);
  }
  if(status != 0) { out->failed=1; return; }
  for(int i=0;i<out->n;++i) {
    const int donor=out->donor_map[i];
    if(donor<0 || donor>=out->n || !isfinite(out->scratch_weights[i])) {
      out->failed=1; return;
    }
    out->weights[(size_t)out->row*out->n+donor]=out->scratch_weights[i];
    out->donor_order[i]=donor;
  }
  out->captured[out->row]=nrhs;
}
