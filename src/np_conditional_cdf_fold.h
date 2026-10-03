/* Conditional CDF-only deleted NN response tiles. Included after the
 * conditional route providers. No density or regression owner calls this code.
 * CDF nodes remain fixed; only radii change with held-out observation identity. */
static int np_cdf_fold_radii(NPConditionalCVLSRouteContext *context, int q, double **yc)
{
  const int n = num_obs_train_extern, d = num_var_continuous_extern;
  const int adaptive = BANDWIDTH_den_extern == BW_ADAP_NN;
  const int nr = adaptive ? n : q;
  const NPNNGeometryContext external_geometry = {.mode = NP_NN_QUERY_EXTERNAL};
  np_conditional_cvls_fold_grid_clear(context);
  context->grid_rows = q;
  context->grid_eval = yc;
  if(d == 0) return 0;
  context->grid_primary = alloc_tmatd(nr,d);
  context->grid_successor = alloc_tmatd(nr,d);
  context->grid_scale = alloc_vecd(d);
  if(!context->grid_primary || !context->grid_successor || !context->grid_scale) return 1;
  for(int l = 0; l < d; ++l) {
    if(adaptive) {
      double **primary = context->beta_y ?
        context->route_y.loo_geometry.primary : context->legacy_y.matrix_bandwidth_y;
      double **successor = context->beta_y ?
        context->route_y.loo_geometry.successor :
        context->legacy_y.matrix_bandwidth_y_successor;
      const double *scale = context->beta_y ?
        context->route_y.loo_geometry.fold_scale : context->legacy_y.adaptive_fold_scale_y;
      if(primary == NULL || successor == NULL || scale == NULL) return 1;
      memcpy(context->grid_primary[l], primary[l], (size_t)n*sizeof(double));
      memcpy(context->grid_successor[l], successor[l], (size_t)n*sizeof(double));
      context->grid_scale[l] = scale[l];
    } else {
      const double h = context->beta_y ? context->route_y.scale_factor[l] :
        context->legacy_y.vsfy[l];
      int k;
      double scale;
      if(np_nn_lookup_from_scale(n-1, 1, h, &k, &scale, NULL) != 0 ||
         compute_nn_distance_train_eval_ctx(n, q, 1,
           matrix_Y_continuous_train_extern[l], yc[l], k, &external_geometry,
           context->grid_primary[l]) != NP_NN_GEOMETRY_OK ||
         compute_nn_distance_train_eval_ctx(n, q, 1,
           matrix_Y_continuous_train_extern[l], yc[l], k+1, &external_geometry,
           context->grid_successor[l]) != NP_NN_GEOMETRY_OK) return 1;
      context->grid_scale[l] = scale;
      for(int e = 0; e < q; ++e) {
        context->grid_primary[l][e] *= scale;
        context->grid_successor[l][e] *= scale;
      }
    }
  }
  return 0;
}

static int np_cdf_fold_grid_prepare(
  NPConditionalCVLSRouteContext *context, int q,
  double **yu, double **yo, double **yc, double **rows, double *logs)
{
  const int n = num_obs_train_extern, d = num_var_continuous_extern;
  const int adaptive = BANDWIDTH_den_extern == BW_ADAP_NN;
  const int nr = adaptive ? n : q;
  double **active = NULL;
  NPBetaScaledRowContext beta;
  int status = 1;
  np_beta_scaled_row_context_init(&beta);
  if(!context->fold_geometry || d < 0 || d > 2 || q < 1 || (d > 0 && yc == NULL))
    return 1;
  if(np_cdf_fold_radii(context,q,yc)) return 1;
  context->grid_rows = q;
  context->grid_variants = 1 << d;
  context->grid_eval = yc;
  context->grid_values[0] = rows;
  context->grid_logs[0] = logs;
  if(d == 0) {
    for(int e = 0; e < q; ++e) {
      if(np_conditional_y_eval_from_ctx(&context->legacy_y,e,yu,yo,yc,q,0,rows[e]))
        return 1;
      logs[e] = 0.0;
    }
    return 0;
  }
  active = alloc_tmatd(nr, d);
  if(active == NULL) goto cleanup;
  for(int v = 0; v < context->grid_variants; ++v) {
    if(v > 0) {
      context->grid_values[v] = alloc_tmatd(n, q);
      context->grid_logs[v] = alloc_vecd(q);
      if(context->grid_values[v] == NULL || context->grid_logs[v] == NULL) goto cleanup;
    }
    for(int l = 0; l < d; ++l)
      memcpy(active[l], (v & (1 << l)) ? context->grid_successor[l] :
        context->grid_primary[l], (size_t)nr*sizeof(double));
    if(context->beta_y) {
      if(np_beta_scaled_row_context_prepare(
           &beta, context->execution_context->y_route,
           context->execution_context->y_diagnostics,
           BANDWIDTH_den_extern, n, q, d, num_var_unordered_extern, num_var_ordered_extern,
           matrix_Y_continuous_train_extern, yc,
           matrix_Y_unordered_train_extern, yu, matrix_Y_ordered_train_extern, yo,
           active, active, context->route_y.operator_code,
           context->route_y.kernel_unordered, context->route_y.kernel_ordered,
           context->route_y.lambda, num_categories_extern_Y, matrix_categorical_vals_extern_Y,
           context->execution_context->categorical_compress, context->route_y.row) !=
           NP_CONTINUOUS_ROW_OK) goto cleanup;
      for(int e = 0; e < q; ++e) {
        beta.row_result.row = context->grid_values[v][e];
        if(np_beta_scaled_row_context_fill(&beta, e, NULL, &context->grid_logs[v][e]) !=
           NP_CONTINUOUS_ROW_OK) goto cleanup;
      }
      np_beta_scaled_row_context_clear(&beta);
    } else {
      /* Borrow scratch/kernel metadata; no ownership passes to the view. */
      NPConditionalYRowCtx view = context->legacy_y;
      view.matrix_bandwidth_y = active;
      for(int e = 0; e < q; ++e) {
        if(np_conditional_y_eval_from_ctx(&view,e,yu,yo,yc,q,0,
             context->grid_values[v][e]) != 0) goto cleanup;
        context->grid_logs[v][e] = 0.0;
      }
    }
  }
  status = 0;
cleanup:
  np_beta_scaled_row_context_clear(&beta);
  if(active != NULL) free_tmat(active);
  return status;
}

/* For more than two response coordinates retain two signed-log factors per
 * coordinate, plus the categorical factor. This avoids 2^d response planes. */
typedef struct {
  double **logs;
  signed char *signs;
  double *scratch;
  signed char *scratch_sign;
  unsigned char *selected;
  int capacity;
} NPCDFFoldFactors;

static int np_cdf_fold_factors_prepare(
  NPConditionalCVLSRouteContext *context, int q, double **yu, double **yo,
  double **yc, NPCDFFoldFactors *f)
{
  const int n = num_obs_train_extern, d = num_var_continuous_extern;
  const int adaptive = BANDWIDTH_den_extern == BW_ADAP_NN;
  NPBetaScaledRowContext beta;
  NPBetaScaledRowCategoricalContext cats;
  int status = 1;
  double *categorical_scratch = NULL;
  np_beta_scaled_row_context_init(&beta);
  memset(&cats, 0, sizeof(cats));
  if(np_cdf_fold_radii(context,q,yc)) goto cleanup;
  for(int l = 0; l < d; ++l) {
    for(int v = 0; v < 2; ++v) {
      double *bandwidth = v ? context->grid_successor[l] : context->grid_primary[l];
      const int plane = (2*l+v)*f->capacity;
      if(context->beta_y) {
        NPContinuousKernelRoute route = *context->execution_context->y_route;
        const int op = OP_INTEGRAL;
        route.segment[0].coordinate_count = 1;
        route.segment[0].lower += l;
        route.segment[0].upper += l;
        if(np_beta_scaled_row_context_prepare(&beta,&route,
             context->execution_context->y_diagnostics,BANDWIDTH_den_extern,n,q,1,0,0,
             matrix_Y_continuous_train_extern+l,yc+l,NULL,NULL,NULL,NULL,
             &bandwidth,&bandwidth,&op,NULL,NULL,NULL,NULL,NULL,0,f->scratch) !=
             NP_CONTINUOUS_ROW_OK) goto cleanup;
        for(int e = 0; e < q; ++e) {
          if(np_beta_scaled_row_context_fill(&beta,e,NULL,NULL) != NP_CONTINUOUS_ROW_OK)
            goto cleanup;
          for(int j = 0; j < n; ++j) {
            const int pos = int_TREE_Y == NP_TREE_TRUE ? ipt_lookup_extern_Y[j] : j;
            f->logs[plane+e][j] = beta.row_workspace.primary_log_absolute[pos];
            f->signs[(size_t)(plane+e)*n+j] = beta.row_workspace.primary_sign[pos];
          }
        }
        np_beta_scaled_row_context_clear(&beta);
      } else {
        const double lower = int_cyker_bound_extern ? vector_cykerlb_extern[l] : -DBL_MAX;
        const double upper = int_cyker_bound_extern ? vector_cykerub_extern[l] : DBL_MAX;
        for(int e = 0; e < q; ++e) for(int j = 0; j < n; ++j) {
          const int pos = int_TREE_Y == NP_TREE_TRUE ? ipt_lookup_extern_Y[j] : j;
          int sign;
          if(np_continuous_kernel_scalar_log(NP_CKERNEL_FAMILY_LEGACY,
               KERNEL_den_extern,2,1,yc[l][e],matrix_Y_continuous_train_extern[l][pos],
               bandwidth[adaptive ? pos : e],lower,upper,&f->logs[plane+e][j],&sign) !=
               NP_CONTINUOUS_KERNEL_SCALAR_OK) goto cleanup;
          f->signs[(size_t)(plane+e)*n+j] = (signed char)sign;
        }
      }
    }
  }
  const int plane = 2*d*f->capacity;
  if(num_var_unordered_extern + num_var_ordered_extern > 0) {
    const int *ops = context->beta_y ? context->route_y.operator_code : context->legacy_y.operator_y;
    const int *ku = context->beta_y ? context->route_y.kernel_unordered : context->legacy_y.kernel_uy;
    const int *ko = context->beta_y ? context->route_y.kernel_ordered : context->legacy_y.kernel_oy;
    const double *lambda = context->beta_y ? context->route_y.lambda : context->legacy_y.lambday;
    categorical_scratch = alloc_vecd(n);
    if(!categorical_scratch) goto cleanup;
    if(np_beta_categorical_factor_context_prepare(&cats,n,q,
         num_var_unordered_extern,num_var_ordered_extern,
         matrix_Y_unordered_train_extern,matrix_Y_ordered_train_extern,yu,yo,
         ku,ko,ops+d,lambda,num_categories_extern_Y,matrix_categorical_vals_extern_Y,
         context->execution_context->categorical_compress,categorical_scratch) != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    for(int e = 0; e < q; ++e) {
      if(np_beta_categorical_log_factor(&cats,e,-1,n,f->scratch,f->scratch_sign) !=
           NP_CONTINUOUS_ROW_OK) goto cleanup;
      for(int j = 0; j < n; ++j) {
        const int pos = int_TREE_Y == NP_TREE_TRUE ? ipt_lookup_extern_Y[j] : j;
        f->logs[plane+e][j] = f->scratch[pos];
        f->signs[(size_t)(plane+e)*n+j] = f->scratch_sign[pos];
      }
    }
  } else for(int e = 0; e < q; ++e) for(int j = 0; j < n; ++j) {
    f->logs[plane+e][j] = 0.0;
    f->signs[(size_t)(plane+e)*n+j] = 1;
  }
  status = 0;
cleanup:
  np_beta_categorical_factor_context_release(&cats);
  free(categorical_scratch);
  np_beta_scaled_row_context_clear(&beta);
  return status;
}

static int np_cdf_fold_factors_fit(
  NPConditionalCVLSRouteContext *context, NPCDFFoldFactors *f,
  int first, int ib, double **xrows, int q, double *fits)
{
  const int n = num_obs_train_extern, d = num_var_continuous_extern;
  const int adaptive = BANDWIDTH_den_extern == BW_ADAP_NN;
  const int nr = adaptive ? n : q;
  for(int b = 0; b < ib; ++b) {
    const int fold = first+b;
    const int held = int_TREE_Y == NP_TREE_TRUE ? ipt_lookup_extern_Y[fold] : fold;
    for(int l = 0; l < d; ++l) for(int j = 0; j < nr; ++j) {
      const int pos = adaptive && int_TREE_Y == NP_TREE_TRUE ? ipt_lookup_extern_Y[j] : j;
      const double primary = context->grid_primary[l][pos];
      const double centre = adaptive ? matrix_Y_continuous_train_extern[l][pos] : context->grid_eval[l][j];
      const double distance = fabs(matrix_Y_continuous_train_extern[l][held]-centre)*context->grid_scale[l];
      double radius = primary;
      if(!(adaptive && j == fold) &&
         (np_nn_two_slot_radius_select(primary,context->grid_successor[l][pos],distance,1,&radius) ||
          !(radius > 0.0))) return 1;
      f->selected[(size_t)l*nr+j] = radius != primary;
    }
    for(int e = 0; e < q; ++e) {
      double maximum = -INFINITY, sum = 0.0;
      for(int j = 0; j < n; ++j) {
        const int cats = 2*d*f->capacity+e;
        double value = f->logs[cats][j];
        int sign = f->signs[(size_t)cats*n+j];
        for(int l = 0; l < d; ++l) {
          const int plane = (2*l+f->selected[(size_t)l*nr+(adaptive ? j : e)])*f->capacity+e;
          value += f->logs[plane][j];
          sign *= f->signs[(size_t)plane*n+j];
        }
        f->scratch[j] = value;
        f->scratch_sign[j] = (signed char)sign;
        if(sign && xrows[b][j] != 0.0 && value > maximum) maximum = value;
      }
      if(maximum != -INFINITY) for(int j = 0; j < n; ++j)
        if(f->scratch_sign[j] && xrows[b][j] != 0.0)
          sum += xrows[b][j]*f->scratch_sign[j]*exp(f->scratch[j]-maximum);
      if(np_continuous_kernel_scaled_restore(sum,maximum,1,&fits[b*q+e]) != NP_CONTINUOUS_ROW_OK)
        return 1;
    }
  }
  return 0;
}

static int np_cdf_fold_cvls(NPConditionalCVLSRouteContext *context, double *cv)
{
  const int n = num_obs_train_extern, m = num_obs_eval_extern;
  const int d = num_var_continuous_extern;
  const int b = MIN(n, np_conditional_lp_cvls_block_size(n, (size_t)4*MAX(1,d)+14U, 4U));
  const int dims[3] = {num_var_unordered_extern,num_var_ordered_extern,d};
  double **eval[3] = {matrix_Y_unordered_eval_extern,matrix_Y_ordered_eval_extern,matrix_Y_continuous_eval_extern};
  double **grid[3] = {NULL,NULL,NULL}, **x[4] = {NULL,NULL,NULL,NULL}, **y = NULL;
  double *fits = NULL, *logs = NULL, *contributions = NULL;
  NPCDFFoldFactors factors = {0};
  int status = 1;
  if(n < 2 || m < 1 || b < 1) return 1;
  const int nblocks = n/b+(n%b != 0);
  int first_block = 0, stride = 1;
#ifdef MPI2
  const int parallel = np_objective_outer_rows_enabled(int_conditional_prepared_context_extern);
  if(parallel) { first_block = my_rank; stride = iNum_Processors; }
#endif
  const int groups = MIN(4, nblocks/stride+(nblocks%stride != 0));
  for(int g = 0; g < groups; ++g) if(!(x[g] = alloc_tmatd(n,b))) goto finish;
  y = alloc_tmatd(n,b);
  fits = alloc_vecd(b*b);
  logs = alloc_vecd(b);
  if(!y || !fits || !logs) goto finish;
  for(int a = 0; a < 3; ++a) if(dims[a] && !(grid[a] = alloc_tmatd(b,dims[a]))) goto finish;
  if(d > 2) {
    size_t cells;
    if(d > (INT_MAX/b-1)/2 || !np_size_mul_checked((size_t)(2*d+1)*b,(size_t)n,&cells)) goto finish;
    factors.capacity = b;
    factors.logs = alloc_tmatd(n,(2*d+1)*b);
    factors.signs = np_jksum_malloc_array_or_die(cells,sizeof(signed char),"CDF fold factor signs");
    factors.scratch = alloc_vecd(n);
    factors.scratch_sign = np_jksum_malloc_array_or_die(n,sizeof(signed char),"CDF fold scratch signs");
    factors.selected = np_jksum_malloc_array_or_die((size_t)d*MAX(n,b),sizeof(unsigned char),"CDF fold radius choices");
    if(!factors.logs || !factors.signs || !factors.scratch || !factors.scratch_sign || !factors.selected) goto finish;
  }
  contributions = alloc_vecd(nblocks);
  if(!contributions) goto finish;
  memset(contributions,0,(size_t)nblocks*sizeof(double));
  *cv = 0.0;
  for(int block = first_block; block < nblocks; block += groups*stride) {
    int counts[4] = {0};
    double sums[4] = {0};
    np_progress_bandwidth_loop_step();
    for(int g = 0; g < groups; ++g) {
      const int first = (block+g*stride)*b;
      counts[g] = MIN(b,MAX(0,n-first));
      for(int t = 0; t < counts[g]; ++t)
        if(np_conditional_cvls_provider_x_row(context,first+t,x[g][t])) goto finish;
    }
    for(int j = 0; j < m; j += b) {
      const int q = MIN(b,m-j);
      for(int a = 0; a < 3; ++a) for(int l = 0; l < dims[a]; ++l) for(int e = 0; e < q; ++e) {
        const int pos = cdfontrain_extern && int_TREE_Y == NP_TREE_TRUE ? ipt_lookup_extern_Y[j+e] : j+e;
        grid[a][l][e] = eval[a][l][pos];
      }
      if(d <= 2 ? np_cdf_fold_grid_prepare(context,q,grid[0],grid[1],grid[2],y,logs) :
         np_cdf_fold_factors_prepare(context,q,grid[0],grid[1],grid[2],&factors)) goto finish;
      for(int g = 0; g < groups; ++g) {
        const int ib = counts[g], first = (block+g*stride)*b;
        if(!ib) continue;
        if(d <= 2 ? np_conditional_cvls_fold_fit_block(context,first,ib,x[g],q,fits) :
           np_cdf_fold_factors_fit(context,&factors,first,ib,x[g],q,fits)) goto finish;
        for(int t = 0; t < ib; ++t) for(int e = 0; e < q; ++e) {
          if(cdfontrain_extern && first+t == j+e) continue;
          const double error = np_conditional_indicator_original_order(first+t,j+e)-fits[t*q+e];
          sums[g] += error*error;
        }
      }
    }
    for(int g = 0; g < groups; ++g)
      if(counts[g]) contributions[block+g*stride] = sums[g];
  }
  status = 0;
finish:
#ifdef MPI2
  status = np_conditional_outer_buffer_finish(parallel,nblocks,status,contributions,
    "NP_RMPI_INJECT_CDIST_CVLS_FAIL_RANK","conditional beta CDF fold blocks MPI_Allreduce");
#endif
  if(status) goto cleanup;
  for(int block = 0; block < nblocks; ++block) *cv += contributions[block];
  status = np_distribution_cvls_finalize(*cv,n,m,cdfontrain_extern,cv) != NP_DISTRIBUTION_CVLS_FINALIZE_OK;
cleanup:
  np_conditional_cvls_fold_grid_clear(context);
  for(int g = 0; g < 4; ++g) if(x[g]) free_tmat(x[g]);
  for(int a = 0; a < 3; ++a) if(grid[a]) free_tmat(grid[a]);
  if(y) free_tmat(y);
  if(factors.logs) free_tmat(factors.logs);
  free(factors.signs); free(factors.scratch); free(factors.scratch_sign); free(factors.selected);
  free(fits); free(logs); free(contributions);
  return status;
}
