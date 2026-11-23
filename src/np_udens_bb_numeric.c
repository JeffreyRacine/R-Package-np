/* Plot-only complete influence arithmetic for the explicit beta-binomial law.
 * The three-moment correction is bounded by the original numerical tolerance;
 * it cannot normalize a materially incorrect law or replace a failed mean. */
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <limits.h>
typedef struct {double h,l;} pair;
static pair val(double x){return (pair){x,0};}
static pair tsum(double a,double b){double s=a+b,v=s-a;return(pair){s,(a-(s-v))+(b-v)};}
static pair add(pair a,pair b){pair s=tsum(a.h,b.h);pair t=tsum(a.l,b.l);pair u=tsum(s.h,s.l+t.h);return tsum(u.h,u.l+t.l);}
static pair neg(pair x){return(pair){-x.h,-x.l};}
static pair sub(pair a,pair b){return add(a,neg(b));}
static pair mul(pair a,pair b){double p=a.h*b.h;double e=fma(a.h,b.h,-p);pair t=tsum(p,e);t=add(t,val(a.h*b.l));t=add(t,val(a.l*b.h));return add(t,val(a.l*b.l));}
static pair divp(pair a,pair b){if(!R_FINITE(b.h)||b.h==0)error("singular pair divisor");pair q=val(a.h/b.h);pair rem=sub(a,mul(b,q));q=add(q,val(rem.h/b.h));rem=sub(a,mul(b,q));return add(q,val(rem.h/b.h));}
static double dval(pair x){return x.h+x.l;}
static pair basis(R_xlen_t i,int j,const double*T,R_xlen_t n,const double*m){return j==0?val(1):sub(val(T[i+(j-1)*n]),val(m[j-1]));}
SEXP C_np_plot_udens_bb_factor(SEXP Q,SEXP T,SEXP A,SEXP HI,SEXP M,SEXP NS){
 if(TYPEOF(Q)!=REALSXP || TYPEOF(T)!=REALSXP || TYPEOF(A)!=REALSXP ||
    TYPEOF(HI)!=REALSXP || TYPEOF(M)!=REALSXP || TYPEOF(NS)!=REALSXP)
   error("beta-binomial factor requires double arrays");
 R_xlen_t n=XLENGTH(Q);
 SEXP da=getAttrib(A,R_DimSymbol), dt=getAttrib(T,R_DimSymbol), dh=getAttrib(HI,R_DimSymbol);
 if(n<3 || n>INT_MAX || TYPEOF(da)!=INTSXP || XLENGTH(da)!=2 ||
    TYPEOF(dt)!=INTSXP || XLENGTH(dt)!=2 || TYPEOF(dh)!=INTSXP || XLENGTH(dh)!=2)
   error("invalid beta-binomial factor dimensions");
 int k=INTEGER(da)[0];
 if(k<1 || INTEGER(da)[1]!=n || INTEGER(dt)[0]!=n || INTEGER(dt)[1]!=2 ||
    INTEGER(dh)[0]!=k || INTEGER(dh)[1]!=2 || XLENGTH(M)!=2 || XLENGTH(NS)!=1 ||
    XLENGTH(A)!=(R_xlen_t)k*n || XLENGTH(T)!=2*n || XLENGTH(HI)!=2*(R_xlen_t)k)
   error("incompatible beta-binomial factor dimensions");
 double nn=REAL(NS)[0];
 if(!R_FINITE(nn) || nn<10 || nn!=floor(nn) || nn>INT_MAX)
   error("invalid beta-binomial sample count");
 if(8.0*(4.0*n+2.0*k*n+6.0*k)>192.0*1024*1024)
   error("beta-binomial factor exceeds its 192 MiB workspace bound");
 for(R_xlen_t i=0;i<XLENGTH(HI);i++) if(!R_FINITE(REAL(HI)[i])) error("nonfinite influence coefficient");
 for(int i=0;i<2;i++) if(!R_FINITE(REAL(M)[i]) || REAL(M)[i]>=0) error("invalid analytic log moment");
 const double*q=REAL(Q),*t=REAL(T),*a=REAL(A),*hi=REAL(HI),*m=REAL(M);pair G[3][3],m0[3],r[3],coef[3];for(int j=0;j<3;j++){m0[j]=val(0);for(int l=0;l<3;l++)G[j][l]=val(0);}
 for(R_xlen_t i=0;i<n;i++){if(!R_FINITE(q[i])||q[i]<=0||!R_FINITE(t[i])||!R_FINITE(t[i+n]))error("invalid positive quadrature input");pair b[3];for(int j=0;j<3;j++){b[j]=basis(i,j,t,n,m);m0[j]=add(m0[j],mul(val(q[i]),b[j]));}for(int j=0;j<3;j++)for(int l=0;l<3;l++)G[j][l]=add(G[j][l],mul(val(q[i]),mul(b[j],b[l])));}
 for(int j=0;j<3;j++){r[j]=sub(val(j==0?1:0),m0[j]);if(fabs(dval(r[j]))>1e-12)error("material original moment error");}
 /* Fixed three-dimensional Gaussian elimination with explicit nonsingularity. */
 for(int j=0;j<3;j++){int pivot=j;for(int l=j+1;l<3;l++)if(fabs(dval(G[l][j]))>fabs(dval(G[pivot][j])))pivot=l;if(pivot!=j){for(int l=0;l<3;l++){pair z=G[j][l];G[j][l]=G[pivot][l];G[pivot][l]=z;}pair z=r[j];r[j]=r[pivot];r[pivot]=z;}if(!R_FINITE(dval(G[j][j]))||dval(G[j][j])==0)error("singular moment Gram");for(int l=j+1;l<3;l++){pair f=divp(G[l][j],G[j][j]);for(int c=j;c<3;c++)G[l][c]=sub(G[l][c],mul(f,G[j][c]));r[l]=sub(r[l],mul(f,r[j]));}}
 for(int j=2;j>=0;j--){pair z=r[j];for(int l=j+1;l<3;l++)z=sub(z,mul(G[j][l],coef[l]));coef[j]=divp(z,G[j][j]);}
 pair *w=(pair*)R_alloc(n,sizeof(pair));double correction=0;pair moments[3]={ {0,0},{0,0},{0,0} };
 for(R_xlen_t i=0;i<n;i++){pair d=val(0);for(int j=0;j<3;j++)d=add(d,mul(basis(i,j,t,n,m),coef[j]));correction=fmax(correction,fabs(dval(d)));if(correction>1e-12)error("material quadrature correction");w[i]=mul(val(q[i]),add(val(1),d));if(!(dval(w[i])>0))error("nonpositive corrected quadrature");for(int j=0;j<3;j++)moments[j]=add(moments[j],mul(w[i],basis(i,j,t,n,m)));}
 SEXP ans=PROTECT(allocMatrix(REALSXP,k,6)); SEXP factor=PROTECT(allocMatrix(REALSXP,k,n));
 for(int row=0;row<k;row++){R_CheckUserInterrupt();pair ea=val(0);for(R_xlen_t i=0;i<n;i++){if(!R_FINITE(a[row+i*k]))error("nonfinite kernel row");ea=add(ea,mul(w[i],val(a[row+i*k])));}pair mean=val(0),variance=val(0);
  for(R_xlen_t i=0;i<n;i++){pair psi=sub(val(a[row+i*k]),ea);for(int j=0;j<2;j++)psi=add(psi,mul(val(hi[row+j*k]),basis(i,j+1,t,n,m)));REAL(factor)[row+i*k]=dval(psi)*sqrt(dval(w[i])/nn);mean=add(mean,mul(w[i],psi));variance=add(variance,mul(w[i],mul(psi,psi)));}
  variance=divp(variance,val(nn));
  if (!R_FINITE(dval(mean)) || fabs(dval(mean)) > 1e-12 ||
      !R_FINITE(dval(variance)) || dval(variance) <= 0) {
    error("beta-binomial complete influence is not centered with positive variance");
  }
  REAL(ans)[row] = dval(mean);
  REAL(ans)[row + k] = dval(variance);
  REAL(ans)[row + 2 * k] = correction;
  REAL(ans)[row + 3 * k] = dval(sub(moments[0], val(1)));
  REAL(ans)[row + 4 * k] = dval(moments[1]);
  REAL(ans)[row + 5 * k] = dval(moments[2]);
 }
 setAttrib(ans,install("factor"),factor);UNPROTECT(2);return ans;
}

/* Plot-only support-polynomial guide expectation; statistical choices stay in R. */
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <limits.h>
#define BB_GUIDE_EVALUATE_LOOP(DEGREE,NODE_SIDE) do { for(R_xlen_t i=0;i<m;i++){R_CheckUserInterrupt();double sum=0,c=0;for(R_xlen_t j=0;j<n;j++){double v=((NODE_SIDE?lv[j]<=-M_LN2:g[i]<=.5)?exp(lv[j])-g[i]:(1-g[i])-exp(lq[j]))/h;if(v<op[i+(DEGREE+1)*m]||v>op[i+(DEGREE+2)*m])continue;double z=(v-op[i+(DEGREE+3)*m])/op[i+(DEGREE+4)*m],prev=1,now=z,poly=op[i]+op[i+m]*z;for(int k=2;k<=DEGREE;k++){double next=((2*k-1)*z*now-(k-1)*prev)/k;poly+=op[i+k*m]*next;prev=now;now=next;}double term=(3/(4*sqrt(5)))*(1-v*v/5)*poly*pr[j];double y=term-c,t=sum+y;c=(t-sum)-y;sum=t;}REAL(ans)[i]=sum/h;} } while(0)
static SEXP bb_guide_evaluate(SEXP G,SEXP H,SEXP LV,SEXP LQ,SEXP PROB,SEXP OP,int degree) {
 if(TYPEOF(G)!=REALSXP || TYPEOF(H)!=REALSXP || TYPEOF(LV)!=REALSXP ||
    TYPEOF(LQ)!=REALSXP || TYPEOF(PROB)!=REALSXP || TYPEOF(OP)!=REALSXP)
   error("beta-binomial guide requires double arrays");
 R_xlen_t m=XLENGTH(G),n=XLENGTH(LV);SEXP d=getAttrib(OP,R_DimSymbol);
 if(m<1 || m>INT_MAX || n<1 || n>192000 || XLENGTH(H)!=1 ||
    XLENGTH(LQ)!=n || XLENGTH(PROB)!=n || TYPEOF(d)!=INTSXP || XLENGTH(d)!=2 ||
    INTEGER(d)[0]!=m || INTEGER(d)[1]!=degree+5 || XLENGTH(OP)!=(degree+5)*m)
   error("invalid beta-binomial guide dimensions");
 double h=REAL(H)[0];if(!R_FINITE(h)||h<=0)error("positive bandwidth required");
 const double*g=REAL(G),*lv=REAL(LV),*lq=REAL(LQ),*pr=REAL(PROB),*op=REAL(OP);
 for(R_xlen_t j=0;j<n;j++)
   if(!R_FINITE(lv[j]) || !R_FINITE(lq[j]) || lv[j]>0 || lq[j]>0 ||
      fabs(exp(lv[j])+exp(lq[j])-1)>1e-14 || !R_FINITE(pr[j]) || pr[j]<0)
     error("invalid guide node or probability");
 for(R_xlen_t i=0;i<m;i++){
   if(!R_FINITE(g[i]) || g[i]<=0 || g[i]>=1)error("guide targets must be strictly interior");
   for(int j=0;j<degree+5;j++)if(!R_FINITE(op[i+j*m]))error("nonfinite guide operator");
   if(op[i+(degree+1)*m]>0 || op[i+(degree+2)*m]<0 || op[i+(degree+1)*m]>=op[i+(degree+2)*m] ||
      op[i+(degree+1)*m]<-sqrt(5) || op[i+(degree+2)*m]>sqrt(5) || op[i+(degree+4)*m]<=0)
     error("invalid guide support operator");
 }
 SEXP ans=PROTECT(allocVector(REALSXP,m));
 if(degree==5) { BB_GUIDE_EVALUATE_LOOP(5,1); }
 else if(degree==7) { BB_GUIDE_EVALUATE_LOOP(7,0); }
 else { BB_GUIDE_EVALUATE_LOOP(9,0); }
 UNPROTECT(1);return ans;
}

/* Retain the original six-argument degree-five entry and operation order. */
SEXP C_np_plot_udens_bb_guide(SEXP G,SEXP H,SEXP LV,SEXP LQ,SEXP PROB,SEXP OP) {
 return bb_guide_evaluate(G,H,LV,LQ,PROB,OP,5);
}
SEXP C_np_plot_udens_bb_guide_order(SEXP G,SEXP H,SEXP LV,SEXP LQ,SEXP PROB,SEXP OP,SEXP DEGREE) {
 if((TYPEOF(DEGREE)!=INTSXP && TYPEOF(DEGREE)!=REALSXP) || XLENGTH(DEGREE)!=1)
  error("beta-binomial guide degree must be 5, 7 or 9");
 double value=TYPEOF(DEGREE)==INTSXP?(double)INTEGER(DEGREE)[0]:REAL(DEGREE)[0];
 if(!R_FINITE(value) || (value!=5 && value!=7 && value!=9))
  error("beta-binomial guide degree must be 5, 7 or 9");
 return bb_guide_evaluate(G,H,LV,LQ,PROB,OP,(int)value);
}

/* Accurate analytic log moments of the same beta-binomial mixture.
 * Binary64 pairs, fma, and bounded series; no long-double precision assumption. */
static pair bb_log(pair x) {
 int e; double f=frexp(x.h,&e); pair y={f,ldexp(x.l,-e)};
 pair z=divp(sub(y,val(1)),add(y,val(1))),z2=mul(z,z),term=z,s=z;
 for(int j=3;j<=65;j+=2){term=mul(term,z2);s=add(s,divp(term,val(j)));}
 pair ln2={0x1.62e42fefa39efp-1,0x1.abc9e3b39803fp-56};
 return add(mul(val(2),s),mul(val(e),ln2));
}
static pair bb_psi(pair x) {
 pair recur=val(0);
 while(x.h<32){recur=add(recur,divp(val(1),x));x=add(x,val(1));}
 pair inv=divp(val(1),x),inv2=mul(inv,inv),power=inv2;
 pair s=sub(bb_log(x),mul(val(.5),inv));
 const double num[]={1,-1,1,-1,1,-691,1,-3617,43867,-174611,77683,-236364091};
 const double den[]={12,120,252,240,132,32760,12,8160,14364,6600,276,65520};
 for(int j=0;j<12;j++){s=sub(s,mul(divp(val(num[j]),val(den[j])),power));power=mul(power,inv2);}
 return sub(s,recur);
}
SEXP C_np_plot_udens_bb_log_moments(SEXP AB,SEXP K,SEXP W) {
 if(TYPEOF(AB)!=REALSXP || XLENGTH(AB)!=2 || TYPEOF(K)!=REALSXP || XLENGTH(K)!=1 || TYPEOF(W)!=REALSXP)
  error("beta-binomial log moments require double shapes, degree and weights");
 double kd=REAL(K)[0],a=REAL(AB)[0],b=REAL(AB)[1];
 if(!R_FINITE(a)||!R_FINITE(b)||a<=0||b<=0||!R_FINITE(kd)||kd<0||kd!=floor(kd)||kd>INT_MAX-1||XLENGTH(W)!=(R_xlen_t)kd+1||!R_FINITE(a+b+kd))
  error("invalid beta-binomial log-moment dimensions or shapes");
 int k=(int)kd;pair mass=val(0),m[2]={ {0,0},{0,0} };
 for(int j=0;j<=k;j++){if(!R_FINITE(REAL(W)[j])||REAL(W)[j]<0)error("invalid beta-binomial log-moment weight");mass=add(mass,val(REAL(W)[j]));}
 if(fabs(dval(sub(mass,val(1))))>1e-12)error("material beta-binomial log-moment mass error");
 pair total=bb_psi(add(add(val(a),val(b)),val(k)));
 SEXP ans=PROTECT(allocMatrix(REALSXP,k+1,2));SEXP mean=PROTECT(allocVector(REALSXP,2));
 for(int j=0;j<=k;j++){
  R_CheckUserInterrupt();pair v[2]={sub(bb_psi(add(val(a),val(j))),total),sub(bb_psi(add(val(b),val(k-j))),total)};
  for(int l=0;l<2;l++){double t=dval(v[l]);if(!R_FINITE(t)||t>=0)error("invalid beta-binomial analytic log moment");REAL(ans)[j+l*(R_xlen_t)(k+1)]=t;m[l]=add(m[l],mul(val(REAL(W)[j]),v[l]));}
 }
 for(int l=0;l<2;l++)REAL(mean)[l]=dval(m[l]);
 setAttrib(ans,install("mean"),mean);UNPROTECT(2);return ans;
}
