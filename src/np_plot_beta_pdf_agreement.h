/* New simultaneous-request numerical agreement, not an interval enclosure.
 * The canonical PDF remains the point/amplitude owner. This helper compares
 * one canonical component with independent Rmath dbeta on the SAME shape,
 * using the nearer observation tail. No h/tau cutoff and no replacement. */
static double np_plot_beta_pdf_component_discrepancy(const np_beta_shape *s,
  double canonical_log_pdf)
{
  const double a=1.0+s->target_unit*s->concentration;
  const double b=1.0+s->target_complement_unit*s->concentration;
  const double reference_log_pdf=(s->observation_unit<=0.5?
    dbeta(s->observation_unit,a,b,1):
    dbeta(s->observation_complement_unit,b,a,1))-log(s->support_length);
  if(ISNAN(canonical_log_pdf)||ISNAN(reference_log_pdf)||
     canonical_log_pdf==R_PosInf||reference_log_pdf==R_PosInf)
    error("np_plot_beta_numerical_unavailable: component amplitude comparison failed");
  if(canonical_log_pdf==reference_log_pdf)return 0.0;
  double log_error;
  if(canonical_log_pdf==R_NegInf)log_error=reference_log_pdf;
  else if(reference_log_pdf==R_NegInf)log_error=canonical_log_pdf;
  else {
    const double high=fmax(canonical_log_pdf,reference_log_pdf);
    log_error=high+log(-expm1(-fabs(canonical_log_pdf-reference_log_pdf)));
  }
  if(log_error==R_NegInf)return 0.0;
  double value=exp(log_error);
  if(!R_FINITE(value)||value==0.0)
    error("np_plot_beta_numerical_unavailable: component discrepancy is not representable");
  return value;
}
