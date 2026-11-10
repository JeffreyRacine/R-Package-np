/* Plot-only fixed-bound beta numerical leaves. No estimator, RNG or MPI owner. */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <limits.h>

static void beta_rule(double a, double b, int n, double *u, double *p,
                      double *off, double *vectors, double *work)
{
    if (!R_FINITE(a) || !R_FINITE(b) || a <= 0 || b <= 0 || n < 2 || n > 2048)
        error("invalid beta plot quadrature geometry");
    /* Jacobi exponents b-1, a-1. Stable products also cover their sum=-1. */
    long double alpha = b - 1.0L, beta = a - 1.0L, s = alpha + beta;
    u[0] = (double)((beta - alpha)/(s + 2));
    for (int i=1; i<n; ++i) {
        long double z=2.0L*i+s;
        u[i]=(double)(((beta-alpha)/z)*(s/(z+2)));
        long double ratio=(i==1) ? 1.0L : (i+s)/(z-1);
        long double square=4.0L*(i/z)*((i+alpha)/z)*((i+beta)/(z+1))*ratio;
        if (!(square>0) || !isfinite(square)) error("invalid beta plot Jacobi recurrence");
        off[i-1]=sqrt((double)square);
    }
    int info=0; char job='V';
    F77_CALL(dstev)(&job,&n,u,off,vectors,&n,work,&info FCONE);
    if (info) error("beta plot Jacobi eigensystem failed (%d)",info);
    long double mass=0, mean=0;
    for (int i=0;i<n;++i) {
        u[i]=(u[i]+1)*.5;
        p[i]=vectors[(size_t)n*i]*vectors[(size_t)n*i];
        if (!(u[i]>0 && u[i]<1) || !(p[i]>0) || !R_FINITE(p[i]))
            error("beta plot quadrature is not representable");
        mass+=p[i];
    }
    for (int i=0;i<n;++i) p[i]/=(double)mass;
    /* Equal-shape measures have exact reflection symmetry. Enforce that
       identity before the certificate, not by loosening its tolerance. */
    if (a==b) for (int i=0;i<(n+1)/2;++i) {
        int j=n-1-i;
        double left=.5*(u[i]+(1-u[j])), weight=.5*(p[i]+p[j]);
        u[i]=left;u[j]=1-left;p[i]=weight;p[j]=weight;
    }
    for (int i=0;i<n;++i) mean+=(long double)p[i]*u[i];
    if (fabsl(mean-(long double)a/(a+b))>2e-13L)
        error("beta plot quadrature mean certificate failed");
}

SEXP C_np_plot_udens_beta_rule(SEXP shape, SEXP order)
{
    if (TYPEOF(shape)!=REALSXP || XLENGTH(shape)!=2 || TYPEOF(order)!=INTSXP || XLENGTH(order)!=1)
        error("invalid beta plot quadrature arguments");
    int n=INTEGER(order)[0]; if (n<2 || n>2048) error("invalid beta plot quadrature order");
    double a=REAL(shape)[0], b=REAL(shape)[1]; int reflect=(a>b);
    SEXP ans=PROTECT(allocMatrix(REALSXP,n,2));
    double *v=(double*)R_alloc((size_t)n*n,sizeof(double));
    double *off=(double*)R_alloc(n,sizeof(double)), *work=(double*)R_alloc(2*n,sizeof(double));
    beta_rule(reflect?b:a,reflect?a:b,n,REAL(ans),REAL(ans)+n,off,v,work);
    if (reflect) for (int i=0;i<(n+1)/2;++i) {
        int j=n-1-i; double x=REAL(ans)[i], p=REAL(ans)[n+i];
        REAL(ans)[i]=1-REAL(ans)[j]; REAL(ans)[n+i]=REAL(ans)[n+j];
        REAL(ans)[j]=1-x; REAL(ans)[n+j]=p;
    }
    UNPROTECT(1); return ans;
}

SEXP C_np_plot_udens_beta_axis(SEXP targets, SEXP lv, SEXP l1v, SEXP scales, SEXP order)
{
    if (TYPEOF(targets)!=REALSXP || TYPEOF(lv)!=REALSXP || TYPEOF(l1v)!=REALSXP ||
        TYPEOF(scales)!=REALSXP || XLENGTH(scales)!=2 || XLENGTH(lv)!=XLENGTH(l1v) ||
        TYPEOF(order)!=INTSXP || XLENGTH(order)!=1) error("invalid beta plot composition arguments");
    R_xlen_t nx=XLENGTH(targets), nv=XLENGTH(lv); int n=INTEGER(order)[0];
    double rh=REAL(scales)[0], rb=REAL(scales)[1];
    if (!nx || !nv || nx>INT_MAX || nv>INT_MAX || n<2 || n>2048 ||
        !R_FINITE(rh) || !R_FINITE(rb) || rh<=0 || rb<=0 ||
        (double)nx*(double)nv*24>33554432.0/sizeof(double))
        error("beta plot composition geometry or 32 MiB leaf budget");
    for (R_xlen_t k=0;k<nx;++k) if (!(REAL(targets)[k]>0 && REAL(targets)[k]<1)) error("beta plot requires interior targets");
    /* A finite interior logit quadrature point may round one log coordinate to
       -0 while the other remains finite. Literal public endpoints are rejected
       before this leaf; their infinite logarithms are still inadmissible here. */
    for (R_xlen_t k=0;k<nv;++k) if (!R_FINITE(REAL(lv)[k]) || !R_FINITE(REAL(l1v)[k]) ||
        REAL(lv)[k]>0 || REAL(l1v)[k]>0 || (REAL(lv)[k]==0 && REAL(l1v)[k]==0) ||
        fabs(exp(REAL(lv)[k])+exp(REAL(l1v)[k])-1)>4e-15)
        error("beta plot requires consistent interior log coordinates");
    SEXP ans=PROTECT(allocVector(REALSXP,nx*nv*24));
    double *out=REAL(ans), *u=(double*)R_alloc(n,sizeof(double)), *p=(double*)R_alloc(n,sizeof(double));
    double *off=(double*)R_alloc(n,sizeof(double)), *vec=(double*)R_alloc((size_t)n*n,sizeof(double));
    double *work=(double*)R_alloc(2*n,sizeof(double)), *aa=(double*)R_alloc(4*n,sizeof(double));
    double *bb=(double*)R_alloc(4*n,sizeof(double)), *normal=(double*)R_alloc(4*n,sizeof(double));
    for (R_xlen_t ix=0;ix<nx;++ix) {
        R_CheckUserInterrupt(); double x=REAL(targets)[ix];
        for (int s=0;s<4;++s) {
            double ah=x/((s+1)*rh), bh=(1-x)/((s+1)*rh),
                   ab=x/((s+1)*rb), bb0=(1-x)/((s+1)*rb);
            double nh=lbeta(ah+1,bh+1), nb=lbeta(ab+1,bb0+1);
            if (!R_FINITE(nh) || !R_FINITE(nb)) error("beta plot kernel parameters not representable");
            for (R_xlen_t j=0;j<nv;++j) {
                out[ix+nx*j+nx*nv*s]=exp(ah*REAL(lv)[j]+bh*REAL(l1v)[j]-nh);
                out[ix+nx*j+nx*nv*(4+s)]=exp(ab*REAL(lv)[j]+bb0*REAL(l1v)[j]-nb);
            }
            int reflect=(x>.5); double xc=reflect?1-x:x;
            beta_rule(1+xc/((s+1)*rh),1+(1-xc)/((s+1)*rh),n,u,p,off,vec,work);
            if (reflect) for(int k=0;k<(n+1)/2;++k) {
                int z=n-1-k; double old=u[k],oldp=p[k];u[k]=1-u[z];p[k]=p[z];u[z]=1-old;p[z]=oldp;
            }
            for (int t=0;t<4;++t) for (int k=0;k<n;++k) {
                aa[t*n+k]=u[k]/((t+1)*rb);bb[t*n+k]=(1-u[k])/((t+1)*rb);
                normal[t*n+k]=lbeta(aa[t*n+k]+1,bb[t*n+k]+1);
            }
            for (int t=0;t<4;++t) for (R_xlen_t j=0;j<nv;++j) {
                long double sum=0;
                for (int k=0;k<n;++k) sum+=(long double)p[k]*exp(aa[t*n+k]*REAL(lv)[j]+bb[t*n+k]*REAL(l1v)[j]-normal[t*n+k]);
                out[ix+nx*j+nx*nv*(8+4*s+t)]=(double)sum;
            }
        }
    }
    SEXP dim=PROTECT(allocVector(INTSXP,3)); INTEGER(dim)[0]=(int)nx; INTEGER(dim)[1]=(int)nv; INTEGER(dim)[2]=24;
    setAttrib(ans,R_DimSymbol,dim); UNPROTECT(2);return ans;
}

SEXP C_np_plot_udens_beta_hull(SEXP center, SEXP lower, SEXP upper)
{
    if (TYPEOF(center)!=REALSXP || TYPEOF(lower)!=REALSXP || TYPEOF(upper)!=REALSXP ||
        XLENGTH(center)!=XLENGTH(lower) || XLENGTH(center)!=XLENGTH(upper) ||
        XLENGTH(center)>INT_MAX) error("invalid beta plot interval geometry");
    R_xlen_t n=XLENGTH(center); SEXP ans=PROTECT(allocMatrix(REALSXP,(int)n,2));
    for (R_xlen_t i=0;i<n;++i) {
        double c=REAL(center)[i], l=REAL(lower)[i], u=REAL(upper)[i];
        if (!R_FINITE(c)||!R_FINITE(l)||!R_FINITE(u)||l>u)
            error("invalid beta plot interval geometry");
        double d=fmax(c-l,u-c); if(d!=0) d=nextafter(d,INFINITY);
        double L=c-d,U=c+d,scale=fmax(fabs(c),fmax(fabs(L),fabs(U)));
        double spacing=nextafter(scale,INFINITY)-scale;
        if(!R_FINITE(L)||!R_FINITE(U)||L>l||U<u || fabs(L/2+U/2-c)>2*spacing)
            error("beta plot outward hull certificate failed");
        REAL(ans)[i]=L;REAL(ans)[n+i]=U;
    }
    UNPROTECT(1);return ans;
}

/* Exact weighted beta likelihood refits. Counts are supplied by the coordinator;
   this leaf never creates a draw or substitutes an influence update. */
SEXP C_np_plot_udens_beta_refits(SEXP data, SEXP counts)
{
    SEXP dim=getAttrib(counts,R_DimSymbol);
    if(TYPEOF(data)!=REALSXP || TYPEOF(counts)!=REALSXP || TYPEOF(dim)!=INTSXP ||
       XLENGTH(dim)!=2 || INTEGER(dim)[0]!=XLENGTH(data) || XLENGTH(data)<10)
        error("invalid beta plot refit dimensions");
    int n=INTEGER(dim)[0],B=INTEGER(dim)[1];
    if(B<1 || (double)n*B>33554432.0/sizeof(double))error("beta plot refit count batch budget");
    double *x=REAL(data),*lx=(double*)R_alloc(n,sizeof(double)),*l1=(double*)R_alloc(n,sizeof(double));
    for(int i=0;i<n;++i) {
        if(!(x[i]>0 && x[i]<1))error("beta plot refit requires interior observations");
        lx[i]=log(x[i]);l1[i]=log1p(-x[i]);
    }
    SEXP ans=PROTECT(allocMatrix(REALSXP,B,3));
    for(int b=0;b<B;++b) {
        R_CheckUserInterrupt();double mass=0,mean=0,t0=0,t1=0;
        for(int i=0;i<n;++i) {
            double c=REAL(counts)[i+(size_t)n*b];
            if(!R_FINITE(c)||c<0||c!=floor(c)||c>n)error("invalid beta plot multinomial count");
            double w=c/n;mass+=c;mean+=w*x[i];t0+=w*lx[i];t1+=w*l1[i];
        }
        if(mass!=n)error("beta plot resample size mismatch");
        double variance=0;
        for(int i=0;i<n;++i) {double d=x[i]-mean;variance+=(REAL(counts)[i+(size_t)n*b]/n)*d*d;}
        double s=mean*(1-mean)/variance-1;
        if(!(s>0)||!R_FINITE(s))error("beta plot refit moment start (draw %d)",b+1);
        double a=mean*s,c=(1-mean)*s;int done=0,iterations=0;
        for(int it=0;it<80;++it) {
            iterations=it+1;
            double z=trigamma(a+c),aa=trigamma(a)-z,bb=trigamma(c)-z;
            double g0=digamma(a)-digamma(a+c)-t0,g1=digamma(c)-digamma(a+c)-t1;
            double determinant=aa*bb-z*z;
            if(!(determinant>0)||!R_FINITE(determinant))error("beta plot refit information (draw %d)",b+1);
            double d0=(bb*g0+z*g1)/determinant,d1=(z*g0+aa*g1)/determinant;
            double dec=g0*d0+g1*d1;
            if(dec<1e-14 && fmax(fabs(g0),fabs(g1))<2e-8) {done=1;break;}
            double obj=lbeta(a,c)-(a-1)*t0-(c-1)*t1,step=1;int accepted=0;
            if(d0>0)step=fmin(step,.99*a/d0);
            if(d1>0)step=fmin(step,.99*c/d1);
            for(int ls=0;ls<40;++ls) {
                double ta=a-step*d0,tb=c-step*d1,val=lbeta(ta,tb)-(ta-1)*t0-(tb-1)*t1;
                if(R_FINITE(val)&&val<=obj-1e-4*step*dec+1e-14) {a=ta;c=tb;accepted=1;break;}
                step*=.5;
            }
            if(!accepted)error("beta plot refit line-search (draw %d)",b+1);
        }
        if(!done)error("beta plot refit stationarity (draw %d)",b+1);
        REAL(ans)[b]=a-1;REAL(ans)[B+b]=c-1;REAL(ans)[2*B+b]=iterations;
    }
    UNPROTECT(1);return ans;
}

/* Reuse exactly evaluated log weights/coordinates inside the same finite
   quadrature sum. No change to nodes, weights, loop order or normalization. */
static double beta_joint_rule_sum(const double *eta,int n,const double *u,const double *p,
                                  const double *v,const double *w)
{
    double ab[4],gamma=eta[4];if(!R_FINITE(gamma))error("invalid beta plot interaction");
    for(int j=0;j<4;++j) {ab[j]=exp(eta[j]);if(!R_FINITE(ab[j])||!(ab[j]>0))error("invalid beta plot joint shape");}
    double lp[2048],lw[2048],zu[2048],zv[2048];
    for(int i=0;i<n;++i) {lp[i]=log(p[i]);lw[i]=log(w[i]);zu[i]=2*u[i]-1;zv[i]=2*v[i]-1;}
    double peak=-INFINITY;
    for(int i=0;i<n;++i)for(int j=0;j<n;++j) {
        double logweight=lp[i]+lw[j]+gamma*zu[i]*zv[j];
        if(logweight>peak)peak=logweight;
    }
    long double mass=0;
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)
        mass+=exp(lp[i]+lw[j]+gamma*zu[i]*zv[j]-peak);
    double result=lbeta(ab[0],ab[1])+lbeta(ab[2],ab[3])+peak+log((double)mass);
    if(!R_FINITE(result))error("nonfinite beta plot joint normalizer");
    return result;
}

static double beta_joint_normalizer(const double *eta,int n,double *u,double *p,
                                   double *v,double *w,double *off,double *vec,double *work)
{
    double ab[4];
    for(int j=0;j<4;++j) {ab[j]=exp(eta[j]);if(!R_FINITE(ab[j])||!(ab[j]>0))error("invalid beta plot joint shape");}
    beta_rule(ab[0],ab[1],n,u,p,off,vec,work);beta_rule(ab[2],ab[3],n,v,w,off,vec,work);
    return beta_joint_rule_sum(eta,n,u,p,v,w);
}

SEXP C_np_plot_udens_beta_normalizer(SEXP eta, SEXP order)
{
    SEXP dim=getAttrib(eta,R_DimSymbol);
    if(TYPEOF(eta)!=REALSXP||TYPEOF(dim)!=INTSXP||XLENGTH(dim)!=2||INTEGER(dim)[1]!=5||
       TYPEOF(order)!=INTSXP||XLENGTH(order)!=1)error("invalid beta plot joint normalizer arguments");
    int B=INTEGER(dim)[0],n=INTEGER(order)[0];
    if(B<1||B>4096||n<2||n>2048)error("beta plot joint normalizer batch/order budget");
    double *u=(double*)R_alloc(n,sizeof(double)),*p=(double*)R_alloc(n,sizeof(double));
    double *v=(double*)R_alloc(n,sizeof(double)),*w=(double*)R_alloc(n,sizeof(double));
    double *off=(double*)R_alloc(n,sizeof(double)),*vec=(double*)R_alloc((size_t)n*n,sizeof(double));
    double *work=(double*)R_alloc(2*n,sizeof(double));SEXP ans=PROTECT(allocVector(REALSXP,B));
    for(int b=0;b<B;++b) {
        R_CheckUserInterrupt();double e[5];for(int j=0;j<5;++j)e[j]=REAL(eta)[j*B+b];
        REAL(ans)[b]=beta_joint_normalizer(e,n,u,p,v,w,off,vec,work);
    }
    UNPROTECT(1);return ans;
}

/* Same 16-node normalizer and five-point derivative as the retained sphere
   objective. This leaf only fuses its arithmetic and workspace preparation. */
SEXP C_np_plot_udens_beta_region_objective(SEXP eta,SEXP transform,SEXP target,SEXP direction,SEXP sign)
{
    if(TYPEOF(eta)!=REALSXP||XLENGTH(eta)!=5||TYPEOF(transform)!=REALSXP||XLENGTH(transform)!=25||
       TYPEOF(target)!=REALSXP||XLENGTH(target)!=5||TYPEOF(direction)!=REALSXP||XLENGTH(direction)!=5||
       TYPEOF(sign)!=REALSXP||XLENGTH(sign)!=1||!(REAL(sign)[0]==1||REAL(sign)[0]==-1))
        error("invalid beta plot region objective arguments");
    double e[5],th[5],gradient[5];
    for(int i=0;i<25;++i)if(!R_FINITE(REAL(transform)[i]))error("invalid beta plot region transform");
    for(int i=0;i<5;++i) {
        if(!R_FINITE(REAL(eta)[i])||!R_FINITE(REAL(target)[i])||!R_FINITE(REAL(direction)[i]))
            error("invalid beta plot region coordinates");
        e[i]=REAL(eta)[i];for(int j=0;j<5;++j)e[i]+=REAL(transform)[i+5*j]*REAL(direction)[j];
        th[i]=i<4?exp(e[i])-1:e[i];
    }
    double u[16],p[16],v[16],w[16],off[16],vec[256],work[32],step=.001;
    double normal=beta_joint_normalizer(e,16,u,p,v,w,off,vec,work),value=-normal;
    for(int j=0;j<5;++j) {
        double original=e[j],values[4],shifts[4]={2,1,-1,-2};
        for(int k=0;k<4;++k) {
            e[j]=original+step*shifts[k];double shifted[16],weight[16];
            if(j<2) {
                beta_rule(exp(e[0]),exp(e[1]),16,shifted,weight,off,vec,work);
                values[k]=beta_joint_rule_sum(e,16,shifted,weight,v,w);
            } else if(j<4) {
                beta_rule(exp(e[2]),exp(e[3]),16,shifted,weight,off,vec,work);
                values[k]=beta_joint_rule_sum(e,16,u,p,shifted,weight);
            } else values[k]=beta_joint_rule_sum(e,16,u,p,v,w);
        }
        e[j]=original;
        gradient[j]=(j<4?th[j]+1:1)*REAL(target)[j]-(-values[0]+8*values[1]-8*values[2]+values[3])/(12*step);
        value+=th[j]*REAL(target)[j];
    }
    SEXP ans=PROTECT(allocVector(VECSXP,4));SET_VECTOR_ELT(ans,0,ScalarReal(REAL(sign)[0]*value));
    for(int k=1;k<4;++k)SET_VECTOR_ELT(ans,k,allocVector(REALSXP,5));
    for(int i=0;i<5;++i) {
        double g=0;for(int j=0;j<5;++j)g+=REAL(transform)[j+5*i]*gradient[j];
        REAL(VECTOR_ELT(ans,1))[i]=REAL(sign)[0]*g;
        REAL(VECTOR_ELT(ans,2))[i]=e[i];REAL(VECTOR_ELT(ans,3))[i]=th[i];
    }
    UNPROTECT(1);return ans;
}

typedef struct {
    int n[2];const double *axis[2];double base[5],A;
    double mass[2][257],mom[2][257][9];
    int K;double probability[257],tail_bound;
} beta_series_workspace;

static void beta_series_state(beta_series_workspace *w,const double *theta,
                              double *A,double *mean,double *information)
{
    double r=fabs(theta[4]),maxT=0;int K=0;
    for(int j=0;j<4;++j)if(!R_FINITE(theta[j])||theta[j]<=-1)error("invalid beta plot refit shape");
    if(!R_FINITE(theta[4]))error("invalid beta plot refit interaction");
    for(int j=0;j<2;++j)for(int i=0;i<2*w->n[j];++i)maxT=fmax(maxT,fabs(w->axis[j][i]));
    double tolerance=1e-18/(1+maxT*maxT);
    if(r>0) {
        double z=4*r; if(z>256)error("beta plot positive-series term budget");K=(int)ceil(z);
        for(;;++K) {
            if(K>256)error("beta plot positive-series term budget");
            double bound=(K+1)*log(z)-lgammafn(K+2)-log1p(-z/(K+2));
            if(bound<log(tolerance))break;
        }
    }
    w->K=K;w->tail_bound=r==0?0:exp((K+1)*log(4*r)-lgammafn(K+2)-log1p(-4*r/(K+2)));
    for(int j=0;j<2;++j) {
        int n=w->n[j];const double *ax=w->axis[j];
        for(int k=0;k<=K;++k) {
            double peak=-INFINITY;
            for(int i=0;i<n;++i) {
                double ly=(j==1&&theta[4]<0)?ax[n+i]:ax[i];
                double z=ax[2*n+i]+ax[i]*(theta[2*j]-w->base[2*j])+ax[n+i]*(theta[2*j+1]-w->base[2*j+1])-2*r*exp(ly)+k*ly;
                peak=fmax(peak,z);
            }
            if(!R_FINITE(peak))error("beta plot positive-series mass");
            long double mass=0,moment[9]={0};
            for(int i=0;i<n;++i) {
                double lv=ax[i],l1=ax[n+i],ly=(j==1&&theta[4]<0)?l1:lv;
                double z=ax[2*n+i]+lv*(theta[2*j]-w->base[2*j])+l1*(theta[2*j+1]-w->base[2*j+1])-2*r*exp(ly)+k*ly;
                double p=exp(z-peak),u=2*exp(lv)-1;
                double T[9]={lv,l1,lv*lv,lv*l1,l1*l1,u,u*u,lv*u,l1*u};
                mass+=p;for(int c=0;c<9;++c)moment[c]+=(long double)p*T[c];
            }
            w->mass[j][k]=peak+log((double)mass);
            for(int c=0;c<9;++c)w->mom[j][k][c]=(double)(moment[c]/mass);
        }
    }
    double lc[257],peak=-INFINITY;long double mass=0;
    for(int k=0;k<=K;++k) {
        lc[k]=w->mass[0][k]+w->mass[1][k]+(r==0?0:k*log(4*r)-lgammafn(k+1));peak=fmax(peak,lc[k]);
    }
    for(int k=0;k<=K;++k)mass+=exp(lc[k]-peak);
    double lz=peak+log((double)mass);*A=w->A+r+lz;
    long double mm[5]={0},raw[25]={0};
    for(int k=0;k<=K;++k) {
        double p=exp(lc[k]-lz),*a=w->mom[0][k],*b=w->mom[1][k];
        w->probability[k]=p;
        double M[5]={a[0],a[1],b[0],b[1],a[5]*b[5]};
        for(int j=0;j<5;++j)mm[j]+=(long double)p*M[j];
        raw[0]+=p*a[2];raw[1]+=p*a[3];raw[5]+=p*a[3];raw[6]+=p*a[4];
        raw[12]+=p*b[2];raw[13]+=p*b[3];raw[17]+=p*b[3];raw[18]+=p*b[4];
        for(int i=0;i<2;++i) {
            for(int j=0;j<2;++j) {double v=p*a[i]*b[j];raw[i+5*(2+j)]+=v;raw[2+j+5*i]+=v;}
            double v=p*a[7+i]*b[5];raw[i+20]+=v;raw[4+5*i]+=v;
            v=p*b[7+i]*a[5];raw[2+i+20]+=v;raw[4+5*(2+i)]+=v;
        }
        raw[24]+=p*a[6]*b[6];
    }
    for(int i=0;i<5;++i)mean[i]=(double)mm[i];
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)information[i+5*j]=(double)(raw[i+5*j]-mm[i]*mm[j]);
}

static double beta_series_delta(beta_series_workspace *w,const double *old,
                                 const double *delta,const double *target)
{
    double size=0,shift=0;for(int k=0;k<5;++k) {size=fmax(size,fabs(delta[k]));shift+=target[k]*delta[k];}
    if(size>=1e-4) {
        double A,B,mean[5],I[25],trial[5];for(int k=0;k<5;++k)trial[k]=old[k]+delta[k];
        beta_series_state(w,old,&A,mean,I);beta_series_state(w,trial,&B,mean,I);return B-A-shift;
    }
    int n1=w->n[0],n2=w->n[1];const double *a=w->axis[0],*b=w->axis[1];
    double peak=-INFINITY,peakz=-INFINITY,maxz=0;
    for(int j=0;j<n2;++j)for(int i=0;i<n1;++i) {
        double T[5]={a[i],a[n1+i],b[j],b[n2+j],(2*exp(a[i])-1)*(2*exp(b[j])-1)};
        double lw=a[2*n1+i]+b[2*n2+j],z=-shift;
        for(int k=0;k<5;++k) {lw+=T[k]*(old[k]-w->base[k]);z+=T[k]*delta[k];}
        peak=fmax(peak,lw);peakz=fmax(peakz,lw+z);maxz=fmax(maxz,fabs(z));
    }
    long double mass=0,total=0;
    for(int j=0;j<n2;++j)for(int i=0;i<n1;++i) {
        double T[5]={a[i],a[n1+i],b[j],b[n2+j],(2*exp(a[i])-1)*(2*exp(b[j])-1)};
        double lw=a[2*n1+i]+b[2*n2+j],z=-shift;
        for(int k=0;k<5;++k) {lw+=T[k]*(old[k]-w->base[k]);z+=T[k]*delta[k];}
        double p=exp(lw-peak);mass+=p;
        total+=(maxz<.1)?(long double)p*expm1(z):exp(lw+z-peakz);
    }
    return maxz<.1?log1p((double)(total/mass)):peakz-peak+log((double)(total/mass));
}

SEXP C_np_plot_udens_beta_joint_refits(SEXP axes,SEXP base,SEXP normal,
                                      SEXP targets,SEXP starts)
{
    SEXP dm=getAttrib(targets,R_DimSymbol),ds=getAttrib(starts,R_DimSymbol);
    if(TYPEOF(axes)!=VECSXP||XLENGTH(axes)!=2||TYPEOF(base)!=REALSXP||XLENGTH(base)!=5||
       TYPEOF(normal)!=REALSXP||XLENGTH(normal)!=1||TYPEOF(targets)!=REALSXP||TYPEOF(starts)!=REALSXP||
       TYPEOF(dm)!=INTSXP||XLENGTH(dm)!=2||TYPEOF(ds)!=INTSXP||XLENGTH(ds)!=2||INTEGER(dm)[1]!=5||
       INTEGER(ds)[1]!=5||INTEGER(dm)[0]!=INTEGER(ds)[0])error("invalid beta plot joint refit arguments");
    int B=INTEGER(dm)[0];if(B<1||B>10000)error("beta plot joint refit batch budget");
    beta_series_workspace *w=(beta_series_workspace*)R_alloc(1,sizeof(beta_series_workspace));
    for(int j=0;j<5;++j) {w->base[j]=REAL(base)[j];if(!R_FINITE(w->base[j]))error("invalid beta plot joint base");}
    w->A=REAL(normal)[0];if(!R_FINITE(w->A))error("invalid beta plot joint base mass");
    for(int j=0;j<2;++j) {
        SEXP ax=VECTOR_ELT(axes,j),ad=getAttrib(ax,R_DimSymbol);
        if(TYPEOF(ax)!=REALSXP||TYPEOF(ad)!=INTSXP||XLENGTH(ad)!=2||INTEGER(ad)[1]!=3||INTEGER(ad)[0]<2||INTEGER(ad)[0]>8192)
            error("invalid beta plot joint refit axis");
        w->n[j]=INTEGER(ad)[0];w->axis[j]=REAL(ax);
        for(int i=0;i<3*w->n[j];++i)if(ISNAN(REAL(ax)[i])||REAL(ax)[i]>0||
            (i<2*w->n[j]&&!R_FINITE(REAL(ax)[i])))error("invalid beta plot joint refit log weight");
    }
    if((double)w->n[0]*w->n[1]>2560000)error("beta plot joint refit tensor budget");
    SEXP ans=PROTECT(allocMatrix(REALSXP,B,6));
    for(int b=0;b<B;++b) {
        R_CheckUserInterrupt();double theta[5],target[5];
        for(int k=0;k<5;++k) {theta[k]=REAL(starts)[b+B*k];target[k]=REAL(targets)[b+B*k];if(!R_FINITE(target[k]))error("invalid beta plot joint target");}
        int done=0,iterations=0;
        for(int it=0;it<50;++it) {
            double A,mean[5],I[25],di[5],dec=0,gradient=0;int piv[5],dim=5,nrhs=1,info;
            beta_series_state(w,theta,&A,mean,I);
            for(int k=0;k<5;++k) {di[k]=mean[k]-target[k];gradient=fmax(gradient,fabs(di[k]));}
            F77_CALL(dgesv)(&dim,&nrhs,I,&dim,piv,di,&dim,&info);
            if(info)error("beta plot joint information solve (draw %d)",b+1);
            for(int k=0;k<5;++k)dec+=(mean[k]-target[k])*di[k];
            iterations=it+1;
            if(dec<1e-14&&gradient<2e-8) {done=1;break;}
            double step=1;for(int k=0;k<4;++k)if(di[k]>0)step=fmin(step,.99*(theta[k]+1)/di[k]);int accepted=0;
            for(int ls=0;ls<35;++ls) {
                double trial[5],delta[5];for(int k=0;k<5;++k) {trial[k]=theta[k]-step*di[k];delta[k]=trial[k]-theta[k];}
                double val=beta_series_delta(w,theta,delta,target);
                if(val<=-1e-4*step*dec+1e-14) {for(int k=0;k<5;++k)theta[k]=trial[k];accepted=1;break;}
                step*=.5;
            }
            if(!accepted)error("beta plot joint likelihood line-search (draw %d)",b+1);
        }
        if(!done)error("beta plot joint likelihood stationarity (draw %d)",b+1);
        for(int k=0;k<5;++k)REAL(ans)[b+B*k]=theta[k];
        REAL(ans)[b+5*B]=iterations;
    }
    UNPROTECT(1);return ans;
}

SEXP C_np_plot_udens_beta_joint_block(SEXP parts,SEXP index,SEXP i1,SEXP i2)
{
    if(TYPEOF(parts)!=VECSXP||XLENGTH(parts)!=2||TYPEOF(index)!=VECSXP||XLENGTH(index)!=2||
       TYPEOF(i1)!=INTSXP||TYPEOF(i2)!=INTSXP||XLENGTH(i1)!=XLENGTH(i2))
        error("invalid beta plot joint block arguments");
    int nx[2],nv[2];const double *a[2];const int *ix[2];R_xlen_t m=0,t=XLENGTH(i1);
    for(int j=0;j<2;++j) {
        SEXP p=VECTOR_ELT(parts,j),d=getAttrib(p,R_DimSymbol),id=VECTOR_ELT(index,j);
        if(TYPEOF(p)!=REALSXP||TYPEOF(d)!=INTSXP||XLENGTH(d)!=3||INTEGER(d)[2]!=24||TYPEOF(id)!=INTSXP)
            error("invalid beta plot joint axis array");
        nx[j]=INTEGER(d)[0];nv[j]=INTEGER(d)[1];a[j]=REAL(p);ix[j]=INTEGER(id);
        if(nx[j]<1||nv[j]<1)error("empty beta plot joint axis");
        if(j==0)m=XLENGTH(id);else if(XLENGTH(id)!=m)error("beta plot joint target count mismatch");
        for(R_xlen_t k=0;k<m;++k)if(ix[j][k]<1||ix[j][k]>nx[j])error("beta plot joint target index");
    }
    if(m<1||t<1||m>INT_MAX||t>INT_MAX||(double)m*t>33554432.0/sizeof(double))
        error("beta plot joint block workspace budget");
    for(R_xlen_t k=0;k<t;++k)if(INTEGER(i1)[k]<1||INTEGER(i1)[k]>nv[0]||INTEGER(i2)[k]<1||INTEGER(i2)[k]>nv[1])
        error("beta plot joint donor index");
    double c[4]={4,-6,4,-1},coefficient[24];
    for(int k=0;k<4;++k)coefficient[k]=coefficient[k+4]=c[k];
    for(int s=0;s<4;++s)for(int t0=0;t0<4;++t0)coefficient[8+4*s+t0]=-c[s]*c[t0];
    SEXP ans=PROTECT(allocMatrix(REALSXP,(int)m,(int)t));double *out=REAL(ans);
    for(R_xlen_t col=0;col<t;++col) {
        if(col%32==0)R_CheckUserInterrupt();
        size_t p0=(size_t)nx[0]*(INTEGER(i1)[col]-1),p1=(size_t)nx[1]*(INTEGER(i2)[col]-1);
        for(R_xlen_t row=0;row<m;++row) {
            double value=0;
            for(int k=0;k<24;++k)value+=coefficient[k]*a[0][ix[0][row]-1+p0+(size_t)nx[0]*nv[0]*k]*a[1][ix[1][row]-1+p1+(size_t)nx[1]*nv[1]*k];
            if(!R_FINITE(value))error("nonfinite beta plot joint composition");
            out[row+m*col]=value;
        }
    }
    UNPROTECT(1);return ans;
}

/* Projection of accepted joint guides on cached composition axes. Each draw's
   positive-series weights are prepared once for the entire requested family. */
SEXP C_np_plot_udens_beta_series_means(SEXP parts,SEXP index,SEXP axes,
                                      SEXP base,SEXP normal,SEXP theta)
{
    SEXP td=getAttrib(theta,R_DimSymbol);
    if(TYPEOF(parts)!=VECSXP||XLENGTH(parts)!=2||TYPEOF(index)!=VECSXP||XLENGTH(index)!=2||
       TYPEOF(axes)!=VECSXP||XLENGTH(axes)!=2||TYPEOF(base)!=REALSXP||XLENGTH(base)!=5||
       TYPEOF(normal)!=REALSXP||XLENGTH(normal)!=1||TYPEOF(theta)!=REALSXP||TYPEOF(td)!=INTSXP||
       XLENGTH(td)!=2||INTEGER(td)[1]!=5)error("invalid beta plot series projection arguments");
    int B=INTEGER(td)[0],m=0,nx[2];const double *a[2];const int *ix[2];
    if(B<1||B>64)error("beta plot series projection batch budget");
    beta_series_workspace *w=(beta_series_workspace*)R_alloc(1,sizeof(beta_series_workspace));
    for(int k=0;k<5;++k) {w->base[k]=REAL(base)[k];if(!R_FINITE(w->base[k]))error("invalid beta plot series base");}
    w->A=REAL(normal)[0];if(!R_FINITE(w->A))error("invalid beta plot series mass");
    double bytes=0;
    for(int j=0;j<2;++j) {
        SEXP p=VECTOR_ELT(parts,j),pd=getAttrib(p,R_DimSymbol),id=VECTOR_ELT(index,j),ax=VECTOR_ELT(axes,j),ad=getAttrib(ax,R_DimSymbol);
        if(TYPEOF(p)!=REALSXP||TYPEOF(pd)!=INTSXP||XLENGTH(pd)!=3||INTEGER(pd)[2]!=24||TYPEOF(id)!=INTSXP||
           TYPEOF(ax)!=REALSXP||TYPEOF(ad)!=INTSXP||XLENGTH(ad)!=2||INTEGER(ad)[1]!=3||INTEGER(ad)[0]!=INTEGER(pd)[1])
            error("invalid beta plot series axes");
        nx[j]=INTEGER(pd)[0];w->n[j]=INTEGER(ad)[0];w->axis[j]=REAL(ax);a[j]=REAL(p);ix[j]=INTEGER(id);
        if(nx[j]<1||w->n[j]<2||w->n[j]>8192||XLENGTH(id)>10000)error("beta plot series geometry budget");
        if(j==0)m=(int)XLENGTH(id);else if(XLENGTH(id)!=m)error("beta plot series target mismatch");
        for(int k=0;k<m;++k)if(ix[j][k]<1||ix[j][k]>nx[j])error("beta plot series target index");
        for(R_xlen_t k=0;k<XLENGTH(p);++k)if(!R_FINITE(a[j][k]))error("nonfinite beta plot series composition");
        for(int k=0;k<3*w->n[j];++k)if(ISNAN(w->axis[j][k])||w->axis[j][k]>0||
            (k<2*w->n[j]&&!R_FINITE(w->axis[j][k])))error("invalid beta plot series log coordinates");
        bytes+=8.0*(257*w->n[j]+24.0*257*nx[j]+(double)XLENGTH(p));
    }
    if(m<1||bytes+8.0*B*m>128*1024*1024)error("beta plot series workspace budget");
    double *prob[2],*G[2];for(int j=0;j<2;++j) {
        prob[j]=(double*)R_alloc((size_t)257*w->n[j],sizeof(double));
        G[j]=(double*)R_alloc((size_t)24*257*nx[j],sizeof(double));
    }
    SEXP out=PROTECT(allocVector(VECSXP,4));SET_VECTOR_ELT(out,0,allocMatrix(REALSXP,B,m));
    for(int k=1;k<4;++k)SET_VECTOR_ELT(out,k,allocVector(REALSXP,B));
    double c[4]={4,-6,4,-1},coefficient[24],one=1,zero=0;char no='N';
    for(int k=0;k<4;++k)coefficient[k]=coefficient[k+4]=c[k];
    for(int s=0;s<4;++s)for(int t=0;t<4;++t)coefficient[8+4*s+t]=-c[s]*c[t];
    for(int b=0;b<B;++b) {
        R_CheckUserInterrupt();double th[5],A,mean[5],I[25];for(int k=0;k<5;++k)th[k]=REAL(theta)[b+B*k];
        beta_series_state(w,th,&A,mean,I);int terms=w->K+1;double r=fabs(th[4]);
        for(int j=0;j<2;++j) {
            int n=w->n[j];const double *ax=w->axis[j];
            for(int k=0;k<terms;++k)for(int i=0;i<n;++i) {
                double ly=(j==1&&th[4]<0)?ax[n+i]:ax[i];
                double lw=ax[2*n+i]+ax[i]*(th[2*j]-w->base[2*j])+ax[n+i]*(th[2*j+1]-w->base[2*j+1])-2*r*exp(ly)+k*ly;
                prob[j][i+(size_t)n*k]=exp(lw-w->mass[j][k]);
            }
            for(int k=0;k<24;++k)
                F77_CALL(dgemm)(&no,&no,&nx[j],&terms,&n,&one,a[j]+(size_t)nx[j]*n*k,&nx[j],prob[j],&n,
                                &zero,G[j]+(size_t)nx[j]*terms*k,&nx[j] FCONE FCONE);
        }
        for(int row=0;row<m;++row) {
            double result=0;
            for(int k=0;k<24;++k) {
                long double sum=0;
                for(int t=0;t<terms;++t)sum+=(long double)w->probability[t]*G[0][ix[0][row]-1+(size_t)nx[0]*(t+terms*k)]*
                    G[1][ix[1][row]-1+(size_t)nx[1]*(t+terms*k)];
                result+=coefficient[k]*(double)sum;
            }
            if(!R_FINITE(result))error("nonfinite beta plot series projection");
            REAL(VECTOR_ELT(out,0))[b+(size_t)B*row]=result;
        }
        REAL(VECTOR_ELT(out,1))[b]=A;REAL(VECTOR_ELT(out,2))[b]=terms;REAL(VECTOR_ELT(out,3))[b]=w->tail_bound;
    }
    UNPROTECT(1);return out;
}

/* Centered joint moments, with a diagonal/cross-only owner for SD inference.
   The two laws retain paired-donor weights. No product-guide approximation. */
SEXP C_np_plot_udens_beta_joint_moments(SEXP parts,SEXP index,SEXP axes,
                                       SEXP parameters,SEXP ratio,SEXP full)
{
    if(TYPEOF(parts)!=VECSXP||XLENGTH(parts)!=2||TYPEOF(index)!=VECSXP||XLENGTH(index)!=2||
       TYPEOF(axes)!=VECSXP||XLENGTH(axes)!=2||TYPEOF(parameters)!=REALSXP||XLENGTH(parameters)!=2||
       TYPEOF(ratio)!=REALSXP||TYPEOF(full)!=LGLSXP||XLENGTH(full)!=1||LOGICAL(full)[0]==NA_LOGICAL)
        error("invalid beta plot joint moment arguments");
    int nx[2],nv[2],m=0;const double *a[2],*ax[2];const int *ix[2];
    for(int j=0;j<2;++j) {
        SEXP p=VECTOR_ELT(parts,j),pd=getAttrib(p,R_DimSymbol),idx=VECTOR_ELT(index,j),v=VECTOR_ELT(axes,j),vd=getAttrib(v,R_DimSymbol);
        if(TYPEOF(p)!=REALSXP||TYPEOF(pd)!=INTSXP||XLENGTH(pd)!=3||INTEGER(pd)[2]!=24||
           TYPEOF(idx)!=INTSXP||TYPEOF(v)!=REALSXP||TYPEOF(vd)!=INTSXP||XLENGTH(vd)!=2||INTEGER(vd)[1]!=3||
           INTEGER(pd)[1]!=INTEGER(vd)[0]||XLENGTH(idx)>INT_MAX)
            error("invalid beta plot joint moment axes");
        nx[j]=INTEGER(pd)[0];nv[j]=INTEGER(pd)[1];a[j]=REAL(p);ax[j]=REAL(v);ix[j]=INTEGER(idx);
        if(nx[j]<1||nv[j]<2||nv[j]>1600)error("beta plot joint moment node cap");
        if(j==0)m=(int)XLENGTH(idx);else if(XLENGTH(idx)!=m)error("beta plot joint moment target mismatch");
        for(int k=0;k<m;++k)if(ix[j][k]<1||ix[j][k]>nx[j])error("beta plot joint moment target index");
        for(int k=0;k<3*nv[j];++k)if(ISNAN(ax[j][k])||ax[j][k]>0||(k<2*nv[j]&&!R_FINITE(ax[j][k])))
            error("invalid beta plot joint moment log axis");
    }
    SEXP rd=getAttrib(ratio,R_DimSymbol);int dense=LOGICAL(full)[0],d=m+5;
    if(m<1||m>10000||TYPEOF(rd)!=INTSXP||XLENGTH(rd)!=2||INTEGER(rd)[0]!=nv[0]||INTEGER(rd)[1]!=nv[1]||
       (double)nv[0]*nv[1]>2560000||!R_FINITE(REAL(parameters)[0])||!R_FINITE(REAL(parameters)[1]))
        error("invalid beta plot joint moment geometry");
    double bytes=8*(2.0*d*256+2.0*(2*d+5*d+(dense?(double)d*d:0))+
                   (double)nv[0]*nv[1]+24.0*(nx[0]*nv[0]+nx[1]*nv[1]));
    if(bytes>224*1024*1024)error("beta plot joint moment leaf workspace budget");
    for(R_xlen_t k=0;k<XLENGTH(ratio);++k)if(!R_FINITE(REAL(ratio)[k]))error("nonfinite beta plot ratio weights");
    SEXP ans=PROTECT(allocVector(VECSXP,2));double *mean[2],*variance[2],*cross[2],*cov[2]={NULL,NULL},logmass[2]={-INFINITY,-INFINITY};
    for(int law=0;law<2;++law) {
        SEXP v=allocVector(VECSXP,5);SET_VECTOR_ELT(ans,law,v);
        SET_VECTOR_ELT(v,0,allocVector(REALSXP,d));mean[law]=REAL(VECTOR_ELT(v,0));
        SET_VECTOR_ELT(v,1,allocVector(REALSXP,d));variance[law]=REAL(VECTOR_ELT(v,1));
        SET_VECTOR_ELT(v,2,allocMatrix(REALSXP,d,5));cross[law]=REAL(VECTOR_ELT(v,2));
        if(dense) {SET_VECTOR_ELT(v,3,allocMatrix(REALSXP,d,d));cov[law]=REAL(VECTOR_ELT(v,3));}
        SET_VECTOR_ELT(v,4,allocVector(REALSXP,1));
        for(int k=0;k<d;++k)mean[law][k]=variance[law][k]=0;
        for(int k=0;k<d*5;++k)cross[law][k]=0;
        if(dense)for(R_xlen_t k=0;k<(R_xlen_t)d*d;++k)cov[law][k]=0;
    }
    double *Z=(double*)R_alloc((size_t)d*256,sizeof(double)),*D=(double*)R_alloc((size_t)d*256,sizeof(double));
    double *bm=(double*)R_alloc(d,sizeof(double)),*delta=(double*)R_alloc(d,sizeof(double));
    double lw[2][256],prob[256],coefficient[24],cc[4]={4,-6,4,-1};
    for(int k=0;k<4;++k)coefficient[k]=coefficient[k+4]=cc[k];
    for(int s=0;s<4;++s)for(int t=0;t<4;++t)coefficient[8+4*s+t]=-cc[s]*cc[t];
    int total=nv[0]*nv[1],inc=1,five=5;double one=1,zero=0;
    char lower='L',normal='N',transpose='T';
    for(int start=0;start<total;start+=256) {
        R_CheckUserInterrupt();int count=total-start<256?total-start:256;
        for(int col=0;col<count;++col) {
            int flat=start+col,i=flat%nv[0],j=flat/nv[0];double interaction=(2*exp(ax[0][i])-1)*(2*exp(ax[1][j])-1);
            for(int row=0;row<m;++row) {
                double z=0;
                for(int k=0;k<24;++k)z+=coefficient[k]*a[0][ix[0][row]-1+(size_t)nx[0]*(i+nv[0]*k)]*
                    a[1][ix[1][row]-1+(size_t)nx[1]*(j+nv[1]*k)];
                if(!R_FINITE(z))error("nonfinite beta plot joint feature");
                Z[row+(size_t)d*col]=z;
            }
            Z[m+(size_t)d*col]=ax[0][i];Z[m+1+(size_t)d*col]=ax[0][nv[0]+i];
            Z[m+2+(size_t)d*col]=ax[1][j];Z[m+3+(size_t)d*col]=ax[1][nv[1]+j];Z[m+4+(size_t)d*col]=interaction;
            lw[0][col]=ax[0][2*nv[0]+i]+ax[1][2*nv[1]+j]+REAL(parameters)[0]*interaction-REAL(parameters)[1];
            lw[1][col]=lw[0][col]+REAL(ratio)[flat];
        }
        for(int law=0;law<2;++law) {
            double peak=-INFINITY;for(int k=0;k<count;++k)peak=fmax(peak,lw[law][k]);
            if(peak==-INFINITY)continue;
            long double mass=0;for(int k=0;k<count;++k) {prob[k]=exp(lw[law][k]-peak);mass+=prob[k];}
            double lb=peak+log((double)mass);
            for(int k=0;k<count;++k)prob[k]/=(double)mass;
            F77_CALL(dgemv)(&normal,&d,&count,&one,Z,&d,prob,&inc,&zero,bm,&inc FCONE);
            double maxlog=fmax(lb,logmass[law]),lm=maxlog+log(exp(lb-maxlog)+exp(logmass[law]-maxlog));
            double w=exp(lb-lm),v=exp(logmass[law]-lm);
            for(int row=0;row<d;++row)delta[row]=bm[row]-mean[law][row];
            for(int col=0;col<count;++col) {
                double p=sqrt(prob[col]);for(int row=0;row<d;++row)D[row+(size_t)d*col]=(Z[row+(size_t)d*col]-bm[row])*p;
            }
            if(dense) {
                for(int col=0;col<d;++col)for(int row=col;row<d;++row)
                    cov[law][row+(size_t)d*col]=v*cov[law][row+(size_t)d*col]+v*w*delta[row]*delta[col];
                F77_CALL(dsyrk)(&lower,&normal,&d,&count,&w,D,&d,&one,cov[law],&d FCONE FCONE);
            } else {
                for(int row=0;row<d;++row) {
                    long double sum=0;for(int col=0;col<count;++col)sum+=(long double)D[row+(size_t)d*col]*D[row+(size_t)d*col];
                    variance[law][row]=v*variance[law][row]+w*(double)sum+v*w*delta[row]*delta[row];
                }
                for(int col=0;col<5;++col)for(int row=0;row<d;++row)
                    cross[law][row+d*col]=v*cross[law][row+d*col]+v*w*delta[row]*delta[m+col];
                F77_CALL(dgemm)(&normal,&transpose,&d,&five,&count,&w,D,&d,D+m,&d,&one,cross[law],&d FCONE FCONE);
            }
            for(int row=0;row<d;++row)mean[law][row]=v*mean[law][row]+w*bm[row];
            logmass[law]=lm;
        }
    }
    for(int law=0;law<2;++law) {
        if(!R_FINITE(logmass[law]))error("beta plot joint moment mass is not finite");
        REAL(VECTOR_ELT(VECTOR_ELT(ans,law),4))[0]=logmass[law];
        if(dense) {
            for(int col=0;col<d;++col)for(int row=0;row<col;++row)cov[law][row+(size_t)d*col]=cov[law][col+(size_t)d*row];
            for(int row=0;row<d;++row)variance[law][row]=cov[law][row+(size_t)d*row];
            for(int col=0;col<5;++col)for(int row=0;row<d;++row)cross[law][row+d*col]=cov[law][row+(size_t)d*(m+col)];
        }
    }
    UNPROTECT(1);return ans;
}

/* A single exact residual covariance owner. It streams quadrature features and
   keeps one target covariance, without a target-by-node tensor or m-by-m map. */
SEXP C_np_plot_udens_beta_joint_covariance(SEXP parts,SEXP index,SEXP axes,
    SEXP parameters,SEXP ratio,SEXP projection,SEXP law,SEXP progress)
{
    SEXP kd=getAttrib(projection,R_DimSymbol),rd=getAttrib(ratio,R_DimSymbol);
    if(TYPEOF(parts)!=VECSXP||XLENGTH(parts)!=2||TYPEOF(index)!=VECSXP||XLENGTH(index)!=2||
       TYPEOF(axes)!=VECSXP||XLENGTH(axes)!=2||TYPEOF(parameters)!=REALSXP||XLENGTH(parameters)!=3||
       TYPEOF(ratio)!=REALSXP||TYPEOF(rd)!=INTSXP||XLENGTH(rd)!=2||TYPEOF(projection)!=REALSXP||
       TYPEOF(kd)!=INTSXP||XLENGTH(kd)!=2||INTEGER(kd)[1]!=5||TYPEOF(law)!=INTSXP||XLENGTH(law)!=1||
       INTEGER(law)[0]<0||INTEGER(law)[0]>1||(!isNull(progress)&&!isFunction(progress)))
        error("invalid beta plot covariance arguments");
    int m=INTEGER(kd)[0],nx[2],nv[2],selected=INTEGER(law)[0];const double *a[2],*ax[2];const int *ix[2];
    double bytes=8.0*((double)m*m+2.0*m*256+4*m);
    if(m<1||m>5000)error("beta plot covariance target budget");
    for(int j=0;j<3;++j)if(!R_FINITE(REAL(parameters)[j]))error("invalid beta plot covariance parameters");
    if(REAL(parameters)[2]<10)error("invalid beta plot covariance sample size");
    for(R_xlen_t k=0;k<XLENGTH(projection);++k)if(!R_FINITE(REAL(projection)[k]))error("invalid beta plot covariance projection");
    for(int j=0;j<2;++j) {
        SEXP p=VECTOR_ELT(parts,j),pd=getAttrib(p,R_DimSymbol),id=VECTOR_ELT(index,j),v=VECTOR_ELT(axes,j),vd=getAttrib(v,R_DimSymbol);
        if(TYPEOF(p)!=REALSXP||TYPEOF(pd)!=INTSXP||XLENGTH(pd)!=3||INTEGER(pd)[2]!=24||TYPEOF(id)!=INTSXP||XLENGTH(id)!=m||
           TYPEOF(v)!=REALSXP||TYPEOF(vd)!=INTSXP||XLENGTH(vd)!=2||INTEGER(vd)[1]!=3||INTEGER(vd)[0]!=INTEGER(pd)[1])
            error("invalid beta plot covariance axes");
        nx[j]=INTEGER(pd)[0];nv[j]=INTEGER(pd)[1];a[j]=REAL(p);ax[j]=REAL(v);ix[j]=INTEGER(id);
        if(nx[j]<1||nv[j]<2||nv[j]>1600||INTEGER(rd)[j]!=nv[j])error("invalid beta plot covariance geometry");
        for(int k=0;k<m;++k)if(ix[j][k]<1||ix[j][k]>nx[j])error("invalid beta plot covariance target index");
        for(R_xlen_t k=0;k<XLENGTH(p);++k)if(!R_FINITE(a[j][k]))error("nonfinite beta plot covariance axis");
        for(int k=0;k<3*nv[j];++k)if(ISNAN(ax[j][k])||ax[j][k]>0||(k<2*nv[j]&&!R_FINITE(ax[j][k])))
            error("invalid beta plot covariance log coordinates");
        bytes+=8.0*XLENGTH(p);
    }
    bytes+=8.0*XLENGTH(ratio);
    if(bytes>192*1024*1024||(double)nv[0]*nv[1]>2560000)error("beta plot covariance leaf workspace budget");
    for(R_xlen_t k=0;k<XLENGTH(ratio);++k)if(!R_FINITE(REAL(ratio)[k]))error("nonfinite beta plot covariance ratio");
    SEXP out=PROTECT(allocMatrix(REALSXP,m,m));double *cov=REAL(out);
    for(R_xlen_t k=0;k<(R_xlen_t)m*m;++k)cov[k]=0;
    double *Z=(double*)R_alloc((size_t)m*256,sizeof(double)),*D=(double*)R_alloc((size_t)m*256,sizeof(double));
    double *mean=(double*)R_alloc(m,sizeof(double)),*bm=(double*)R_alloc(m,sizeof(double)),*delta=(double*)R_alloc(m,sizeof(double));
    for(int k=0;k<m;++k)mean[k]=0;
    double lw[256],prob[256],cc[4]={4,-6,4,-1},coefficient[24],logmass=-INFINITY;
    for(int k=0;k<4;++k)coefficient[k]=coefficient[k+4]=cc[k];
    for(int s=0;s<4;++s)for(int t=0;t<4;++t)coefficient[8+4*s+t]=-cc[s]*cc[t];
    int total=nv[0]*nv[1],inc=1;double one=1,zero=0;char no='N',lower='L';
    for(int start=0;start<total;start+=256) {
        R_CheckUserInterrupt();int count=total-start<256?total-start:256;
        for(int col=0;col<count;++col) {
            int flat=start+col,i=flat%nv[0],j=flat/nv[0];
            double T[5]={ax[0][i],ax[0][nv[0]+i],ax[1][j],ax[1][nv[1]+j],(2*exp(ax[0][i])-1)*(2*exp(ax[1][j])-1)};
            for(int row=0;row<m;++row) {
                double value=0;
                for(int k=0;k<24;++k)value+=coefficient[k]*a[0][ix[0][row]-1+(size_t)nx[0]*(i+nv[0]*k)]*
                    a[1][ix[1][row]-1+(size_t)nx[1]*(j+nv[1]*k)];
                for(int k=0;k<5;++k)value-=REAL(projection)[row+(size_t)m*k]*T[k];
                if(!R_FINITE(value))error("nonfinite beta plot residual covariance feature");
                Z[row+(size_t)m*col]=value;
            }
            lw[col]=ax[0][2*nv[0]+i]+ax[1][2*nv[1]+j]+REAL(parameters)[0]*T[4]-REAL(parameters)[1]+(selected?REAL(ratio)[flat]:0);
        }
        double peak=-INFINITY;for(int col=0;col<count;++col)peak=fmax(peak,lw[col]);if(peak==-INFINITY)continue;
        long double mass=0;for(int col=0;col<count;++col) {prob[col]=exp(lw[col]-peak);mass+=prob[col];}
        double lb=peak+log((double)mass);for(int col=0;col<count;++col)prob[col]/=(double)mass;
        F77_CALL(dgemv)(&no,&m,&count,&one,Z,&m,prob,&inc,&zero,bm,&inc FCONE);
        double peakmass=fmax(lb,logmass),lm=peakmass+log(exp(lb-peakmass)+exp(logmass-peakmass)),w=exp(lb-lm),v=exp(logmass-lm);
        for(int row=0;row<m;++row)delta[row]=bm[row]-mean[row];
        for(int col=0;col<count;++col)for(int row=0;row<m;++row)D[row+(size_t)m*col]=(Z[row+(size_t)m*col]-bm[row])*sqrt(prob[col]);
        for(int col=0;col<m;++col)for(int row=col;row<m;++row)cov[row+(size_t)m*col]=v*cov[row+(size_t)m*col]+v*w*delta[row]*delta[col];
        F77_CALL(dsyrk)(&lower,&no,&m,&count,&w,D,&m,&one,cov,&m FCONE FCONE);
        for(int row=0;row<m;++row)mean[row]=v*mean[row]+w*bm[row];
        logmass=lm;
        if(!isNull(progress)&&(start%16384==0||start+count==total)) {
            SEXP call=PROTECT(lang4(progress,R_NilValue,R_NilValue,R_NilValue));
            SETCADR(call,mkString("exact covariance"));SETCADDR(call,ScalarInteger(start+count));
            SETCADDDR(call,ScalarInteger(total));
            eval(call,R_BaseEnv);UNPROTECT(1);
        }
    }
    if(!R_FINITE(logmass))error("beta plot residual covariance mass");
    for(int col=0;col<m;++col)for(int row=col;row<m;++row) {
        double value=cov[row+(size_t)m*col]/REAL(parameters)[2];
        cov[row+(size_t)m*col]=cov[col+(size_t)m*row]=value;
    }
    UNPROTECT(1);return out;
}

/* Complete eigensystem, in one owned matrix. No rank truncation or remainder
   envelope. DSYEV overwrites this copy, never the caller's covariance. */
SEXP C_np_plot_udens_beta_factor(SEXP covariance)
{
    SEXP dim=getAttrib(covariance,R_DimSymbol);
    if(TYPEOF(covariance)!=REALSXP||TYPEOF(dim)!=INTSXP||XLENGTH(dim)!=2||INTEGER(dim)[0]!=INTEGER(dim)[1])
        error("invalid beta plot factor covariance");
    int m=INTEGER(dim)[0];if(m<1||m>5000||16.0*m*m+8.0*64*m>192*1024*1024)
        error("beta plot factor workspace budget");
    const double *C=REAL(covariance);double *sd=(double*)R_alloc(m,sizeof(double));
    for(int i=0;i<m;++i) {
        if(!R_FINITE(C[i+(size_t)m*i])||C[i+(size_t)m*i]<=0)error("invalid beta plot factor variance");
        sd[i]=sqrt(C[i+(size_t)m*i]);
    }
    SEXP ans=PROTECT(allocMatrix(REALSXP,m,m));double *A=REAL(ans);
    for(int j=0;j<m;++j)for(int i=0;i<m;++i) {
        if(!R_FINITE(C[i+(size_t)m*j])||fabs(C[i+(size_t)m*j]-C[j+(size_t)m*i])>1e-12+1e-9*fabs(C[i+(size_t)m*j]))
            error("invalid beta plot factor symmetry");
        A[i+(size_t)m*j]=C[i+(size_t)m*j]/sd[i]/sd[j];
    }
    double *values=(double*)R_alloc(m,sizeof(double)),query;int lwork=-1,info;char vectors='V',lower='L';
    F77_CALL(dsyev)(&vectors,&lower,&m,A,&m,values,&query,&lwork,&info FCONE FCONE);
    if(info||!R_FINITE(query)||query<1||query>INT_MAX)error("beta plot eigensystem workspace query");
    lwork=(int)query;double *work=(double*)R_alloc(lwork,sizeof(double));
    R_CheckUserInterrupt();F77_CALL(dsyev)(&vectors,&lower,&m,A,&m,values,work,&lwork,&info FCONE FCONE);
    if(info)error("beta plot full eigensystem convergence");
    for(int j=0;j<m;++j) {
        if(!R_FINITE(values[j])||values[j]<-1e-10)error("beta plot correlation is not positive semidefinite");
        double scale=sqrt(fmax(values[j],0));for(int i=0;i<m;++i)A[i+(size_t)m*j]*=scale;
    }
    double *block=(double*)R_alloc((size_t)64*m,sizeof(double)),one=1,zero=0;char no='N',tr='T';
    for(int first=0;first<m;first+=64) {
        R_CheckUserInterrupt();int rows=m-first<64?m-first:64;
        F77_CALL(dgemm)(&no,&tr,&rows,&m,&m,&one,A+first,&m,A,&m,&zero,block,&rows FCONE FCONE);
        for(int j=0;j<m;++j)for(int i=0;i<rows;++i)
            if(!R_FINITE(block[i+(size_t)rows*j])||fabs(block[i+(size_t)rows*j]-C[first+i+(size_t)m*j]/sd[first+i]/sd[j])>=1e-9)
                error("beta plot full-rank factor reconstruction");
    }
    UNPROTECT(1);return ans;
}
