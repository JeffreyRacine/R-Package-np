/* Plot-only compact q2 computational leaves. Full covariance, no rank policy. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <limits.h>

/* The compact adapter uses a*b+a+b, embedded in the shared 24-slot format.
   Validate the entire encoding before decoding; other families never enter. */
static void compact_encoding(SEXP parts)
{
    for(int j=0;j<2;++j) {
        SEXP p=VECTOR_ELT(parts,j),dim=getAttrib(p,R_DimSymbol);
        R_xlen_t cells=(R_xlen_t)INTEGER(dim)[0]*INTEGER(dim)[1];
        const double *x=REAL(p);
        for(R_xlen_t i=0;i<cells;++i) {
            double a=x[i];if(!R_FINITE(a))error("nonfinite compact q2 feature");
            for(int k=0;k<24;++k) {
                double expected=k<4?a:k<8?(j?1:a):k%4==0?(j?-.25*a:1):0;
                if(x[i+cells*k]!=expected)
                    error("invalid compact q2 feature encoding");
            }
        }
    }
}

/* Conditional full-QR Gram: E Var(F|v) + Var(E(F|v)). Householder QR retains
   every basis row; its only dimension switch retains every original node when
   nodes <= columns. All residual and between-mean terms enter positive Grams. */
SEXP C_np_plot_udens_compact_q2_covariance(SEXP parts,SEXP index,SEXP axes,
    SEXP parameters,SEXP ratio,SEXP projection,SEXP law,SEXP progress)
{
    SEXP kd=getAttrib(projection,R_DimSymbol),rd=getAttrib(ratio,R_DimSymbol);
    if(TYPEOF(parts)!=VECSXP||XLENGTH(parts)!=2||TYPEOF(index)!=VECSXP||XLENGTH(index)!=2||
       TYPEOF(axes)!=VECSXP||XLENGTH(axes)!=2||TYPEOF(parameters)!=REALSXP||XLENGTH(parameters)!=3||
       TYPEOF(ratio)!=REALSXP||TYPEOF(rd)!=INTSXP||XLENGTH(rd)!=2||TYPEOF(projection)!=REALSXP||
       TYPEOF(kd)!=INTSXP||XLENGTH(kd)!=2||INTEGER(kd)[1]!=5||TYPEOF(law)!=INTSXP||XLENGTH(law)!=1||
       INTEGER(law)[0]<0||INTEGER(law)[0]>1||(!isNull(progress)&&!isFunction(progress)))
        error("invalid compact q2 covariance arguments");
    int m=INTEGER(kd)[0],nx[2],nv[2],selected=INTEGER(law)[0];
    const double *a[2],*ax[2],*K=REAL(projection);const int *ix[2];
    if(m<1||m>5000)error("compact q2 covariance target budget");
    for(int j=0;j<3;++j)if(!R_FINITE(REAL(parameters)[j]))error("invalid compact q2 covariance parameters");
    if(REAL(parameters)[2]<10)error("invalid compact q2 covariance sample size");
    for(R_xlen_t k=0;k<XLENGTH(projection);++k)if(!R_FINITE(K[k]))error("invalid compact q2 covariance projection");
    double bytes=8.0*((double)m*m+(double)m*256+8*m+XLENGTH(ratio));
    for(int j=0;j<2;++j) {
        SEXP p=VECTOR_ELT(parts,j),pd=getAttrib(p,R_DimSymbol),id=VECTOR_ELT(index,j),v=VECTOR_ELT(axes,j),vd=getAttrib(v,R_DimSymbol);
        if(TYPEOF(p)!=REALSXP||TYPEOF(pd)!=INTSXP||XLENGTH(pd)!=3||INTEGER(pd)[2]!=24||TYPEOF(id)!=INTSXP||XLENGTH(id)!=m||
           TYPEOF(v)!=REALSXP||TYPEOF(vd)!=INTSXP||XLENGTH(vd)!=2||INTEGER(vd)[1]!=3||INTEGER(vd)[0]!=INTEGER(pd)[1])
            error("invalid compact q2 covariance axes");
        nx[j]=INTEGER(pd)[0];nv[j]=INTEGER(pd)[1];a[j]=REAL(p);ax[j]=REAL(v);ix[j]=INTEGER(id);
        if(nx[j]<1||nv[j]<2||nv[j]>1600||INTEGER(rd)[j]!=nv[j])error("invalid compact q2 covariance geometry");
        for(int k=0;k<m;++k)if(ix[j][k]<1||ix[j][k]>nx[j])error("invalid compact q2 covariance target index");
        for(int k=0;k<3*nv[j];++k)if(ISNAN(ax[j][k])||ax[j][k]>0||(k<2*nv[j]&&!R_FINITE(ax[j][k])))
            error("invalid compact q2 covariance log coordinates");
        bytes+=8.0*XLENGTH(p);
    }
    for(R_xlen_t k=0;k<XLENGTH(ratio);++k)if(!R_FINITE(REAL(ratio)[k]))error("nonfinite compact q2 covariance ratio");
    compact_encoding(parts);
    int d=nx[0]+3,qr=nv[0]>d,rows=qr?d:nv[0],lwork=qr?64*d:1;
    bytes+=8.0*((double)d*nv[1]+(qr?(double)d*nv[0]:0)+lwork+d+4*nv[0]+4*nv[1]);
    if(bytes>192*1024*1024||(double)nv[0]*nv[1]>2560000)
        error("compact q2 covariance leaf workspace budget");
    double *cm=(double*)R_alloc((size_t)d*nv[1],sizeof(double));
    double *lc=(double*)R_alloc(nv[1],sizeof(double)),*pv=(double*)R_alloc(nv[1],sizeof(double));
    double *lw=(double*)R_alloc(nv[0],sizeof(double)),*prob=(double*)R_alloc(nv[0],sizeof(double));
    double *xu=(double*)R_alloc(nv[0],sizeof(double));
    long double *grand=(long double*)R_alloc(m,sizeof(long double));
    double gamma=REAL(parameters)[0],logz=REAL(parameters)[1];
    for(int i=0;i<nv[0];++i)xu[i]=2*exp(ax[0][i])-1;
    for(int t=0;t<m;++t)grand[t]=0;
    double peakall=-INFINITY;
    for(int j=0;j<nv[1];++j) {
        R_CheckUserInterrupt();double zv=2*exp(ax[1][j])-1,peak=-INFINITY;
        for(int i=0;i<nv[0];++i) {
            lw[i]=ax[0][2*nv[0]+i]+ax[1][2*nv[1]+j]+gamma*xu[i]*zv-logz+(selected?REAL(ratio)[i+(size_t)nv[0]*j]:0);
            peak=fmax(peak,lw[i]);
        }
        if(peak==-INFINITY){lc[j]=-INFINITY;for(int t=0;t<d;++t)cm[t+(size_t)d*j]=0;continue;}
        long double mass=0;for(int i=0;i<nv[0];++i){prob[i]=exp(lw[i]-peak);mass+=prob[i];}
        lc[j]=peak+log((double)mass);peakall=fmax(peakall,lc[j]);
        for(int i=0;i<nv[0];++i)prob[i]/=(double)mass;
        for(int t=0;t<d;++t) {
            long double value=0;
            for(int i=0;i<nv[0];++i) {
                double x=t<nx[0]?a[0][t+(size_t)nx[0]*i]:t==nx[0]?ax[0][i]:t==nx[0]+1?ax[0][nv[0]+i]:xu[i];
                value+=(long double)prob[i]*x;
            }
            cm[t+(size_t)d*j]=(double)value;
        }
    }
    if(!R_FINITE(peakall))error("compact q2 covariance mass is not finite");
    long double mass=0;for(int j=0;j<nv[1];++j){pv[j]=exp(lc[j]-peakall);mass+=pv[j];}
    for(int j=0;j<nv[1];++j) {
        pv[j]/=(double)mass;if(pv[j]==0)continue;
        double *muj=cm+(size_t)d*j,zv=2*exp(ax[1][j])-1;
        double T[5]={muj[nx[0]],muj[nx[0]+1],ax[1][j],ax[1][nv[1]+j],muj[nx[0]+2]*zv};
        for(int t=0;t<m;++t) {
            double b=a[1][ix[1][t]-1+(size_t)nx[1]*j];
            long double value=(long double)(b+1)*muj[ix[0][t]-1]+b;
            for(int k=0;k<5;++k)value-=(long double)K[t+(size_t)m*k]*T[k];
            grand[t]+=(long double)pv[j]*value;
        }
    }
    SEXP out=PROTECT(allocMatrix(REALSXP,m,m));double *cov=REAL(out);
    for(R_xlen_t k=0;k<(R_xlen_t)m*m;++k)cov[k]=0;
    double *Z=(double*)R_alloc((size_t)m*256,sizeof(double));int used=0,info;
    double *Q=qr?(double*)R_alloc((size_t)nv[0]*d,sizeof(double)):NULL;
    double *tau=qr?(double*)R_alloc(d,sizeof(double)):NULL;
    double *work=qr?(double*)R_alloc(lwork,sizeof(double)):NULL;
    double one=1;char lower='L',no='N';
    for(int j=0;j<nv[1];++j) {
        R_CheckUserInterrupt();if(pv[j]==0)continue;
        double *muj=cm+(size_t)d*j,zv=2*exp(ax[1][j])-1,peak=-INFINITY;
        for(int i=0;i<nv[0];++i) {
            lw[i]=ax[0][2*nv[0]+i]+ax[1][2*nv[1]+j]+gamma*xu[i]*zv-logz+(selected?REAL(ratio)[i+(size_t)nv[0]*j]:0);
            peak=fmax(peak,lw[i]);
        }
        long double sum=0;for(int i=0;i<nv[0];++i){prob[i]=exp(lw[i]-peak);sum+=prob[i];}
        for(int i=0;i<nv[0];++i)prob[i]/=(double)sum;
        if(qr) {
            for(int t=0;t<d;++t)for(int i=0;i<nv[0];++i) {
                double x=t<nx[0]?a[0][t+(size_t)nx[0]*i]:t==nx[0]?ax[0][i]:t==nx[0]+1?ax[0][nv[0]+i]:xu[i];
                Q[i+(size_t)nv[0]*t]=(x-muj[t])*sqrt(prob[i]);
            }
            F77_CALL(dgeqrf)(&nv[0],&d,Q,&nv[0],tau,work,&lwork,&info);
            if(info)error("compact q2 complete conditional QR failed (%d)",info);
        }
        double T[5]={muj[nx[0]],muj[nx[0]+1],ax[1][j],ax[1][nv[1]+j],muj[nx[0]+2]*zv};
        for(int r=0;r<=rows;++r) {
            for(int t=0;t<m;++t) {
                int ai=ix[0][t]-1;double b=a[1][ix[1][t]-1+(size_t)nx[1]*j];long double value;
                if(r==rows) {
                    value=(long double)(b+1)*muj[ai]+b;
                    for(int k=0;k<5;++k)value-=(long double)K[t+(size_t)m*k]*T[k];
                    value-=grand[t];
                } else if(qr) {
                    value=r<=ai?(long double)(b+1)*Q[r+(size_t)nv[0]*ai]:0;
                    for(int k=0;k<3;++k)if(r<=nx[0]+k)
                        value-=(long double)K[t+(size_t)m*(k==2?4:k)]*(k==2?zv:1)*Q[r+(size_t)nv[0]*(nx[0]+k)];
                } else {
                    value=(long double)(b+1)*(a[0][ai+(size_t)nx[0]*r]-muj[ai]);
                    value-=(long double)K[t]*(ax[0][r]-muj[nx[0]]);
                    value-=(long double)K[t+(size_t)m]*(ax[0][nv[0]+r]-muj[nx[0]+1]);
                    value-=(long double)K[t+(size_t)m*4]*zv*(xu[r]-muj[nx[0]+2]);
                    value*=sqrt(prob[r]);
                }
                Z[t+(size_t)m*used]=(double)(value*sqrt(pv[j]));
                if(!R_FINITE(Z[t+(size_t)m*used]))error("nonfinite compact q2 covariance feature");
            }
            ++used;
            if(used==256) {
                F77_CALL(dsyrk)(&lower,&no,&m,&used,&one,Z,&m,&one,cov,&m FCONE FCONE);used=0;
            }
        }
        if(!isNull(progress)&&(j%16==0||j+1==nv[1])) {
            SEXP call=PROTECT(lang4(progress,R_NilValue,R_NilValue,R_NilValue));
            SETCADR(call,mkString("exact covariance"));SETCADDR(call,ScalarInteger((j+1)*nv[0]));SETCADDDR(call,ScalarInteger(nv[0]*nv[1]));
            eval(call,R_BaseEnv);UNPROTECT(1);
        }
    }
    if(used)F77_CALL(dsyrk)(&lower,&no,&m,&used,&one,Z,&m,&one,cov,&m FCONE FCONE);
    for(int col=0;col<m;++col)for(int row=col;row<m;++row) {
        double value=cov[row+(size_t)m*col]/REAL(parameters)[2];
        cov[row+(size_t)m*col]=cov[col+(size_t)m*row]=value;
    }
#ifdef NP_COMPACT_Q2_TEST_WITNESS
    SEXP witness=PROTECT(allocVector(INTSXP,3));INTEGER(witness)[0]=qr;INTEGER(witness)[1]=rows;INTEGER(witness)[2]=nv[1];
    setAttrib(out,install("test.only.owner"),witness);UNPROTECT(1);
#endif
    UNPROTECT(1);return out;
}

SEXP C_np_plot_udens_compact_q2_moments(SEXP parts,SEXP index,SEXP axes,
                                       SEXP parameters,SEXP ratio,SEXP full)
{
    if(TYPEOF(parts)!=VECSXP||XLENGTH(parts)!=2||TYPEOF(index)!=VECSXP||XLENGTH(index)!=2||
       TYPEOF(axes)!=VECSXP||XLENGTH(axes)!=2||TYPEOF(parameters)!=REALSXP||XLENGTH(parameters)!=2||
       TYPEOF(ratio)!=REALSXP||TYPEOF(full)!=LGLSXP||XLENGTH(full)!=1||LOGICAL(full)[0]==NA_LOGICAL)
        error("invalid compact q2 joint moment arguments");
    int nx[2],nv[2],m=0;const double *a[2],*ax[2];const int *ix[2];
    for(int j=0;j<2;++j) {
        SEXP p=VECTOR_ELT(parts,j),pd=getAttrib(p,R_DimSymbol),idx=VECTOR_ELT(index,j),v=VECTOR_ELT(axes,j),vd=getAttrib(v,R_DimSymbol);
        if(TYPEOF(p)!=REALSXP||TYPEOF(pd)!=INTSXP||XLENGTH(pd)!=3||INTEGER(pd)[2]!=24||
           TYPEOF(idx)!=INTSXP||TYPEOF(v)!=REALSXP||TYPEOF(vd)!=INTSXP||XLENGTH(vd)!=2||INTEGER(vd)[1]!=3||
           INTEGER(pd)[1]!=INTEGER(vd)[0]||XLENGTH(idx)>INT_MAX)
            error("invalid compact q2 joint moment axes");
        nx[j]=INTEGER(pd)[0];nv[j]=INTEGER(pd)[1];a[j]=REAL(p);ax[j]=REAL(v);ix[j]=INTEGER(idx);
        if(nx[j]<1||nv[j]<2||nv[j]>1600)error("compact q2 joint moment node cap");
        if(j==0)m=(int)XLENGTH(idx);else if(XLENGTH(idx)!=m)error("compact q2 joint moment target mismatch");
        for(int k=0;k<m;++k)if(ix[j][k]<1||ix[j][k]>nx[j])error("compact q2 joint moment target index");
        for(int k=0;k<3*nv[j];++k)if(ISNAN(ax[j][k])||ax[j][k]>0||(k<2*nv[j]&&!R_FINITE(ax[j][k])))
            error("invalid compact q2 joint moment log axis");
    }
    compact_encoding(parts);
    SEXP rd=getAttrib(ratio,R_DimSymbol);int dense=LOGICAL(full)[0],d=m+5;
    if(m<1||m>10000||TYPEOF(rd)!=INTSXP||XLENGTH(rd)!=2||INTEGER(rd)[0]!=nv[0]||INTEGER(rd)[1]!=nv[1]||
       (double)nv[0]*nv[1]>2560000||!R_FINITE(REAL(parameters)[0])||!R_FINITE(REAL(parameters)[1]))
        error("invalid compact q2 joint moment geometry");
    double bytes=8*(2.0*d*256+2.0*(2*d+5*d+(dense?(double)d*d:0))+
                   (double)nv[0]*nv[1]+24.0*(nx[0]*nv[0]+nx[1]*nv[1]));
    if(bytes>224*1024*1024)error("compact q2 joint moment leaf workspace budget");
    for(R_xlen_t k=0;k<XLENGTH(ratio);++k)if(!R_FINITE(REAL(ratio)[k]))error("nonfinite compact q2 ratio weights");
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
    double lw[2][256],prob[256];
    int total=nv[0]*nv[1],inc=1,five=5;double one=1,zero=0;
    char lower='L',normal='N',transpose='T';
    for(int start=0;start<total;start+=256) {
        R_CheckUserInterrupt();int count=total-start<256?total-start:256;
        for(int col=0;col<count;++col) {
            int flat=start+col,i=flat%nv[0],j=flat/nv[0];double interaction=(2*exp(ax[0][i])-1)*(2*exp(ax[1][j])-1);
            for(int row=0;row<m;++row) {
                double av=a[0][ix[0][row]-1+(size_t)nx[0]*i];
                double bv=a[1][ix[1][row]-1+(size_t)nx[1]*j];
                double z=av*bv+av+bv;
                if(!R_FINITE(z))error("nonfinite compact q2 joint feature");
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
        if(!R_FINITE(logmass[law]))error("compact q2 joint moment mass is not finite");
        REAL(VECTOR_ELT(VECTOR_ELT(ans,law),4))[0]=logmass[law];
        if(dense) {
            for(int col=0;col<d;++col)for(int row=0;row<col;++row)cov[law][row+(size_t)d*col]=cov[law][col+(size_t)d*row];
            for(int row=0;row<d;++row)variance[law][row]=cov[law][row+(size_t)d*row];
            for(int col=0;col<5;++col)for(int row=0;row<d;++row)cross[law][row+d*col]=cov[law][row+(size_t)d*(m+col)];
        }
    }
    UNPROTECT(1);return ans;
}

