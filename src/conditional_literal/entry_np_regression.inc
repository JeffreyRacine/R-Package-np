#line 23294 "np.c"
void np_regression(double * tuno, double * tord, double * tcon, double * ty,
                   double * euno, double * eord, double * econ, double * ey,
                   double * rbw,
                   double * mcv, double * padnum,
                   double * nconfac, double * ncatfac, double * mysd,
                   int * myopti,
                   int * glp_degree,
                   int * glp_gradient_order,
                   int * glp_bernstein,
                   int * glp_basis,
                   double * cm, double * cmerr, double * g, double *gerr,
                   double * xtra,
                   double * ckerlb, double * ckerub,
                   const NPContinuousKernelRoute *kernel_route,
                   NPContinuousKernelDerivativeDiagnostics *kernel_route_diagnostics,
                   int categorical_compress,
                   NPRegressionLPEmptyRows *empty_rows,
                   const NPRegressionGradientRequest *gradient_request,
                   int *variance_unavailable, int **gradient_zero_out,
                   NPRegressionFitProjection *fit_projection){
  NPRegressionFittedCall call = {0};
  call.tuno = tuno;
  call.tord = tord;
  call.tcon = tcon;
  call.ty = ty;
  call.euno = euno;
  call.eord = eord;
  call.econ = econ;
  call.ey = ey;
  call.rbw = rbw;
  call.mcv = mcv;
  call.padnum = padnum;
  call.nconfac = nconfac;
  call.ncatfac = ncatfac;
  call.mysd = mysd;
  call.myopti = myopti;
  call.glp_degree = glp_degree;
  call.glp_gradient_order = glp_gradient_order;
  call.glp_bernstein = glp_bernstein;
  call.glp_basis = glp_basis;
  call.cm = cm;
  call.cmerr = cmerr;
  call.g = g;
  call.gerr = gerr;
  call.xtra = xtra;
  call.ckerlb = ckerlb;
  call.ckerub = ckerub;
  call.kernel_route = kernel_route;
  call.kernel_route_diagnostics = kernel_route_diagnostics;
  call.categorical_compress = categorical_compress;
  call.empty_rows = empty_rows;
  call.gradient_request = gradient_request;
  call.variance_unavailable = variance_unavailable;
  call.gradient_zero_out = gradient_zero_out;
  call.fit_projection = fit_projection;
  call.owner.x_columns[0] = myopti[REG_NUNOI];
  call.owner.x_columns[1] = myopti[REG_NORDI];
  call.owner.x_columns[2] = myopti[REG_NCONI];
  call.owner.num_var = call.owner.x_columns[0] + call.owner.x_columns[1] +
    call.owner.x_columns[2];
  R_UnwindProtect(np_regression_fitted_execute, &call,
                  np_regression_fitted_unwind, &call.owner, NULL);
  np_regression_fitted_clear(&call.owner);
  if(empty_rows != NULL) empty_rows->row_map = NULL;
}
