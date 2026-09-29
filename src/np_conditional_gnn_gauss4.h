/* Scalar Gaussian2 Gauss-Legendre remainder envelope, in reciprocal
 * response coordinates. The fixed constants bound Gaussian derivatives of
 * orders 0..8; the panel budget remains 1e-12/(8*n). This only admits four
 * nodes when that rule already meets the existing absolute error contract.
 * No tail truncation, kernel approximation, or tolerance change. */
static double np_cgnn_cert_log_term(double slope, double u, int r) {
  if (slope == 0.0 && r) return R_NegInf;
  const double product = slope*u;
  const double d = fmax(0.0, fabs(1.0+product)-16*DBL_EPSILON*(1+fabs(product)));
  return (r ? r*log(fabs(slope)) : 0.0)-d*d/4;
}
static double np_cgnn_cert_log_rectangle(int r, double lo, double hi,
                                  double smin, double smax) {
  double ans = R_NegInf;
  const double ends[2] = {smin, smax}, us[2] = {lo, hi};
  for (int j=0; j<2; ++j) {
    const double s=ends[j], u=s ? fmin(hi,fmax(lo,-1/s)) : lo;
    ans=fmax(ans,np_cgnn_cert_log_term(s,u,r));
  }
  for (int j=0; j<2; ++j) {
    const double u=us[j];
    if (u==0.0) continue;
    const double d=sqrt(1.0+8*r), roots[2]={(-1-d)/(2*u),4*r/((1+d)*u)};
    for (int k=0; k<2; ++k)
      if (roots[k]>=smin && roots[k]<=smax)
        ans=fmax(ans,np_cgnn_cert_log_term(roots[k],u,r));
  }
  return ans;
}
static double np_cgnn_gauss4_log_bound(double lo,double hi,double smin,double smax,double amp) {
  static const double g[9]={0.39894228044132696,0.342198280348197,
    0.45719550807364207,0.7648386419581483,1.4920035039229589,
    3.2720933887648305,7.888767891935734,20.59363516792374,57.57528256595846};
  static const double choose[9]={1,8,28,56,70,56,28,8,1};
  const double width=hi-lo;
  if (amp==0) return R_NegInf;
  if (!(width>0 && amp>0 && R_FINITE(width) && R_FINITE(smin) && R_FINITE(smax)))
    return R_PosInf;
  const double pad=16*DBL_EPSILON*fmax(1,fmax(fabs(smin),fabs(smax)));
  smin-=pad;smax+=pad;
  double t[9],terms[9],largest=R_NegInf;
  const double logwidth=log(width);
  for(int r=0;r<9;++r)t[r]=r*logwidth+np_cgnn_cert_log_rectangle(r,lo,hi,smin,smax);
  for(int r=0;r<9;++r){terms[r]=log(choose[r]*g[r]*g[8-r])+t[r]+t[8-r];largest=fmax(largest,terms[r]);}
  if(!R_FINITE(largest))return R_PosInf;
  double sum=0;for(int r=0;r<9;++r)sum+=exp(terms[r]-largest);
  const double c4=24.0*24*24*24/(9.0*40320*40320*40320);
  return logwidth+log(c4*amp)+largest+log(sum)+log1p(4096*DBL_EPSILON);
}

#ifdef NP_CF167_TRACE
static unsigned long np_cgnn_trace_tree_nodes, np_cgnn_trace_pruned_nodes;
static unsigned long np_cgnn_trace_projected_nodes, np_cgnn_trace_certificates;
static int np_cgnn_trace_rank(void) {
#ifdef MPI2
  return my_rank;
#else
  return 0;
#endif
}
#endif
