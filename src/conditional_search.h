#ifndef NP_CONDITIONAL_SEARCH_H
#define NP_CONDITIONAL_SEARCH_H
void np_conditional_powell(int restrict_search,int integer,double *restricted,
  double *point,double **directions,int n,double ftol,double tol,double small,
  int itmax,int *iterations,double *value,double (*objective)(double *));
#endif
