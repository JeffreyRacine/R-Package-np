#ifndef NP_CONDITIONAL_GLOBAL_RANK_CERTIFICATE_H
#define NP_CONDITIONAL_GLOBAL_RANK_CERTIFICATE_H
/* Sufficient full-rank certificate for an already selected global conditional
 * design. Inconclusive bounds retain the canonical per-row SVD. No degree,
 * estimator, ridge, tree or MPI policy is selected here. See the separately
 * bound certificate proof for inverse-residual and rank-one perturbation bounds. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef struct {
 int p,ready;double gamma,r1,ri,m1,mi,err1,erri,inv1,invi;
 double *a,*m,*scale,*u,*v,*w,*d,*t,*c;
} NPCQRCertificate;
static double np_cqr_cert_up(double x){return nextafter(x,INFINITY);}
static double np_cqr_cert_add(double a,double b){return np_cqr_cert_up(a+b);}
static double np_cqr_cert_mul(double a,double b){return np_cqr_cert_up(a*b);}
static double np_cqr_cert_div(double a,double b){return np_cqr_cert_up(a/b);}
static double np_cqr_cert_down(double x){return nextafter(x,-INFINITY);}
static void np_cqr_cert_clear(NPCQRCertificate *s){
 free(s->a);free(s->m);free(s->scale);free(s->u);free(s->v);free(s->w);free(s->d);free(s->t);free(s->c);memset(s,0,sizeof(*s));
}
static void np_cqr_cert_norms(const double *a,int p,double *one,double *inf){
 const double pu=((double)p+2.)*DBL_EPSILON;
 const double den=np_cqr_cert_down(1.-pu/(1.-pu));
 *one=0.;*inf=0.;
 for(int i=0;i<p;i++){double row=0.,col=0.;for(int j=0;j<p;j++){
  const double x=a[i+(size_t)p*j],y=a[j+(size_t)p*i];
  if(!isfinite(x)||!isfinite(y)){*one=*inf=INFINITY;return;}
  row+=fabs(x);col+=fabs(y);
 }row=np_cqr_cert_div(np_cqr_cert_add(row,(double)p*DBL_MIN),den);col=np_cqr_cert_div(np_cqr_cert_add(col,(double)p*DBL_MIN),den);
 if(row>*inf)*inf=row;if(col>*one)*one=col;}
}
/* Inputs are finite full Gram values from the incumbent compensated owner. */
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
static int np_cqr_cert_prepare(NPCQRCertificate *s,int p,const double *g){
 int info=0,lwork=-1,status=0;double query=0.;int *pivot=NULL;double *work=NULL;
 if(!s||!g)return -1;
 np_cqr_cert_clear(s);
 if(p<2||p>INT_MAX/8||(size_t)p>SIZE_MAX/(size_t)p/sizeof(double))return -1;
 s->p=p;s->gamma=np_cqr_cert_up(((double)p+8.)*DBL_EPSILON/(1.-((double)p+8.)*DBL_EPSILON));
 s->a=calloc((size_t)p*p,sizeof(double));s->m=calloc((size_t)p*p,sizeof(double));s->c=calloc((size_t)p*p,sizeof(double));
 s->scale=calloc(p,sizeof(double));s->u=calloc(p,sizeof(double));s->v=calloc(p,sizeof(double));s->w=calloc(p,sizeof(double));s->d=calloc(p,sizeof(double));s->t=calloc(p,sizeof(double));
 pivot=calloc(p,sizeof(int));
 if(!s->a||!s->m||!s->c||!s->scale||!s->u||!s->v||!s->w||!s->d||!s->t||!pivot){status=-1;goto done;}
 for(int i=0;i<p;i++){
  if(!isnormal(g[i+(size_t)p*i])||g[i+(size_t)p*i]<=0.)goto done;
  s->scale[i]=1./sqrt(g[i+(size_t)p*i]);if(!isnormal(s->scale[i]))goto done;
 }
 for(int j=0;j<p;j++)for(int i=0;i<p;i++){
  double a=g[i+(size_t)p*j]*s->scale[i]*s->scale[j];if(!isfinite(a))goto done;s->a[i+(size_t)p*j]=a;
 }
 memcpy(s->m,s->a,(size_t)p*p*sizeof(double));
 F77_CALL(dgetrf)(&p,&p,s->m,&p,pivot,&info);if(info){status=info<0?-1:0;goto done;}
 F77_CALL(dgetri)(&p,s->m,&p,pivot,&query,&lwork,&info);if(info||!isfinite(query)||query<p||query>INT_MAX){status=-1;goto done;}
 lwork=(int)ceil(query);if((size_t)lwork>SIZE_MAX/sizeof(double)){status=-1;goto done;}
 work=calloc(lwork,sizeof(double));if(!work){status=-1;goto done;}
 F77_CALL(dgetri)(&p,s->m,&p,pivot,work,&lwork,&info);if(info){status=info<0?-1:0;goto done;}
 for(size_t j=0;j<(size_t)p*p;j++)if(!isfinite(s->m[j]))goto done;
 np_cqr_cert_norms(s->m,p,&s->m1,&s->mi);
 /* Enclose the right inverse residual including dot/subtraction roundoff. */
 for(int j=0;j<p;j++)for(int i=0;i<p;i++){
  double sum=0.,mag=0.;for(int k=0;k<p;k++){
   const double z=s->a[i+(size_t)p*k]*s->m[k+(size_t)p*j];sum+=z;mag=np_cqr_cert_add(mag,np_cqr_cert_up(fabs(z)));
  }
  double delta=i==j?1.:0.;
  s->c[i+(size_t)p*j]=np_cqr_cert_add(fabs(delta-sum),np_cqr_cert_add(np_cqr_cert_mul(s->gamma,np_cqr_cert_add(mag,delta)),np_cqr_cert_mul((double)p+8.,DBL_MIN)));
 }
 np_cqr_cert_norms(s->c,p,&s->r1,&s->ri);
 if(!(s->r1<1.)||!(s->ri<1.)||!isfinite(s->m1)||!isfinite(s->mi))goto done;
 s->inv1=np_cqr_cert_div(s->m1,np_cqr_cert_down(1.-s->r1));s->invi=np_cqr_cert_div(s->mi,np_cqr_cert_down(1.-s->ri));
 s->err1=np_cqr_cert_mul(s->inv1,s->r1);s->erri=np_cqr_cert_mul(s->invi,s->ri);
 s->ready=isfinite(s->inv1)&&isfinite(s->invi);
done:free(pivot);free(work);return status<0?status:s->ready;
}
/* A false return means no certificate; caller must still use the cold SVD. */
static int np_cqr_cert_deleted(NPCQRCertificate *s,const double *g,const double *x,double *lower_out){
 const int p=s->p;double u1=0.,ui=0.,v1=0.,vi=0.,w1=0.,wi=0.,tmax=0.;
 if(!s->ready)return 0;
 for(int i=0;i<p;i++){
  double diag=fabs(g[i+(size_t)p*i]);if(!isnormal(diag))return 0;
  s->d[i]=1./sqrt(diag);s->t[i]=s->scale[i]/s->d[i];s->u[i]=s->scale[i]*x[i];
  if(!isnormal(s->d[i])||!isnormal(s->t[i])||!isfinite(s->u[i]))return 0;
  tmax=fmax(tmax,fabs(s->t[i]));u1=np_cqr_cert_add(u1,fabs(s->u[i]));ui=fmax(ui,fabs(s->u[i]));
 }
 for(int i=0;i<p;i++){
  s->v[i]=0.;s->w[i]=0.;for(int j=0;j<p;j++){
   s->v[i]+=s->m[i+(size_t)p*j]*s->u[j];s->w[i]+=s->u[j]*s->m[j+(size_t)p*i];
  }
  v1=np_cqr_cert_add(v1,fabs(s->v[i]));vi=fmax(vi,fabs(s->v[i]));w1=np_cqr_cert_add(w1,fabs(s->w[i]));wi=fmax(wi,fabs(s->w[i]));
 }
 const double tiny=np_cqr_cert_mul((double)p+8.,DBL_MIN);
 const double ev=np_cqr_cert_add(np_cqr_cert_mul(np_cqr_cert_add(s->erri,np_cqr_cert_mul(s->gamma,s->mi)),ui),tiny);
 const double ew=np_cqr_cert_add(np_cqr_cert_mul(np_cqr_cert_add(s->err1,np_cqr_cert_mul(s->gamma,s->m1)),ui),tiny);
 double h=0.;for(int i=0;i<p;i++)h+=s->u[i]*s->v[i];
 double eh=np_cqr_cert_add(np_cqr_cert_mul(u1,np_cqr_cert_add(ev,np_cqr_cert_mul(s->gamma,vi))),np_cqr_cert_add(np_cqr_cert_mul(DBL_EPSILON,np_cqr_cert_add(1.,fabs(h))),tiny));
 double den=np_cqr_cert_down(fabs(1.-h)-eh);if(!(den>0.)||!isfinite(den))return 0;
 double bi=np_cqr_cert_add(s->invi,np_cqr_cert_div(np_cqr_cert_mul(np_cqr_cert_add(vi,ev),np_cqr_cert_add(w1,np_cqr_cert_mul(p,ew))),den));
 double b1=np_cqr_cert_add(s->inv1,np_cqr_cert_div(np_cqr_cert_mul(np_cqr_cert_add(v1,np_cqr_cert_mul(p,ev)),np_cqr_cert_add(wi,ew)),den));
 /* C is exactly the incumbent cold matrix, including multiplication order. */
 for(int j=0;j<p;j++)for(int i=0;i<p;i++){
  s->c[i+(size_t)p*j]=g[i+(size_t)p*j]*s->d[i]*s->d[j];if(!isfinite(s->c[i+(size_t)p*j]))return 0;
 }
 double c1,ci;np_cqr_cert_norms(s->c,p,&c1,&ci);
 double ei=0.,e1=0.;
 /* Positive sums are enclosed once per row/column. The two denominators
  * cover entry-bound arithmetic and p-term accumulation, respectively. */
 const double sum_den=np_cqr_cert_down(1.-s->gamma);
 memset(s->w,0,(size_t)p*sizeof(double));
 for(int i=0;i<p;i++){
  double rows=0.;for(int j=0;j<p;j++){
   double z=s->c[i+(size_t)p*j]*s->t[i]*s->t[j],uv=s->u[i]*s->u[j];
   double residual=z-s->a[i+(size_t)p*j]+uv;
   double mag=fabs(z)+fabs(s->a[i+(size_t)p*j])+fabs(uv);
   double err=fabs(residual)+s->gamma*mag+tiny;
   if(!isfinite(err))return 0;
   rows+=err;s->w[j]+=err;
  }ei=fmax(ei,np_cqr_cert_div(np_cqr_cert_div(rows,sum_den),sum_den));
 }
 for(int j=0;j<p;j++)e1=fmax(e1,np_cqr_cert_div(np_cqr_cert_div(s->w[j],sum_den),sum_den));
 double qi=np_cqr_cert_mul(bi,ei),q1=np_cqr_cert_mul(b1,e1);if(!(qi<1.)||!(q1<1.))return 0;
 double tt=np_cqr_cert_mul(tmax,tmax);
 double invi=np_cqr_cert_mul(tt,np_cqr_cert_div(bi,np_cqr_cert_down(1.-qi)));
 double inv1=np_cqr_cert_mul(tt,np_cqr_cert_div(b1,np_cqr_cert_down(1.-q1)));
 double norm=np_cqr_cert_up(sqrt(np_cqr_cert_mul(c1,ci))),inverse=np_cqr_cert_up(sqrt(np_cqr_cert_mul(inv1,invi)));
 double lower=np_cqr_cert_down(1./np_cqr_cert_mul(norm,inverse));
 const double pu=(double)p*DBL_EPSILON,gamma=pu/(1.-pu);
 if(lower_out)*lower_out=lower;
 return isfinite(lower)&&lower>np_cqr_cert_up(sqrt(gamma));
}
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#endif
