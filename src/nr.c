#include <float.h>
#include <math.h>
#include <stdio.h>

#include "headers.h"

#ifdef RCSID
static char rcsid[] = "$Id: nr.c,v 1.10 2006/11/02 16:56:49 tristen Exp $";
#endif

extern int int_VERBOSE;

/* The following routines (sort, powell, nrerror, vector, free_vector,
   linimn, brent, mnbrak, f1dim, erfun, ran3), are based on the
   routine(s) from the book Numerical Recipes in C (Cambridge
   University Press), Copyright (C) Copyright (C) 1987, 1988 by
   Numerical Recipes Software. Permission to include requested August
   1, 2006.  Use of these routines other than as an integral part of
   the np and npRmpi package requires an additional license from
   Numerical Recipes Software.  Further distribution in any form is
   prohibited.

   Note that some of these routines are heavily modified to admit
   restricted search, linmin with derivatives, and integer search
   using gradient information. */

typedef struct FCOMPLEX {double r,i;}
fcomplex;
typedef struct IMMENSE {unsigned long l,r;}
immense;
typedef struct GREAT {unsigned short l,c,r;}
great;

#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

#define GOLD 1.618034
#define GLIMIT 100.0
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define CGOLD 0.3819660
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349


/*
 * Heapsort algorithm based on Knuth's adaptation
 */

void sort(int n, double ra[])
{
    int l,j,ir,i;
    double rra;

    if(n == 0) return;

    l=(n >> 1)+1;
    ir=n;
    for (;;)
    {
        if (l > 1)
            rra=ra[--l];
        else
        {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1)
            {
                ra[1]=rra;
                return;
            }
        }
        i=l;
        j=l << 1;
        while (j <= ir)
        {
            if (j < ir && ra[j] < ra[j+1]) ++j;
            if (rra < ra[j])
            {
                ra[i]=ra[j];
                j += (i=j);
            }
            else j=ir+1;
        }
        ra[i]=rra;
    }
}


void  powell(int RESTRICT, int INTEGER, double *p_restrict, double *p, double **xi, int n, double ftol, double tol, double small, int itmax, int *iter,
double *fret, double (*func)(double *))
{
    int i,ibig,j;
    double t,fptt,fp,del;
    double *pt,*ptt,*xit;

    /* In powell(), want user to know instantly */

		spinner(4);

/* Only changes _non zero_ restrictions... not a general approach
   but valid for the current purpose */

    pt=vector(1,n);
    ptt=vector(1,n);
    xit=vector(1,n);

    *fret=(*func)(p);

    for (j=1;j<=n;j++)
    {
        if((RESTRICT == 1) && (p_restrict[j] > 0.0))
        {
            p[j]=p_restrict[j];                   /* Initialize parameters to restricted values */
        }
        pt[j]=p[j];                               /* pt contains initial values */
    }
    for (*iter=1;;(*iter)++)
    {
			spinner((*iter)-1);
/* Main body up to max number of iterations */
        fp=(*fret);
        ibig=0;
        del=(double)0.0;
        for (i=1;i<=n;i++)
        {
            for (j=1;j<=n;j++)
            {
                xit[j]=xi[j][i];
            }
            fptt=(*fret);
            if(RESTRICT == 1)
            {
                linmin(1,INTEGER,p_restrict,p,xit,n,tol,small,itmax,fret,func);
            }
            else
            {
                linmin(0,INTEGER,p,p,xit,n,tol,small,itmax,fret,func);
            }
            if (fabs(fptt-(*fret)) > del)
            {
                del=(double) fabs(fptt-(*fret));
                ibig=i;
            }
        }
/* Begin checking... */
        if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret)))
        {
            free_vector(xit,1,n);
            free_vector(ptt,1,n);
            free_vector(pt,1,n);
            return;
        }
        if (*iter == itmax)
        {
/*        nrerror("Too many iterations in routine POWELL");*/
            if(int_VERBOSE == 1)
            {
                printf("\n**Maximum number of iterations reached in routine POWELL\n");
            }
            pt=vector(1,n);
            ptt=vector(1,n);
            xit=vector(1,n);
            return;
        }
/* End checking and continue... */
        for (j=1;j<=n;j++)
        {
            if((RESTRICT == 1 ) && (p_restrict[j] > 0.0))
            {
                ptt[j]=pt[j]=p_restrict[j];
                xit[j]= (double) 0.0;
            }
            else
            {
                ptt[j]=(double) (2.0*p[j]-pt[j]);
                xit[j]=p[j]-pt[j];
                pt[j]=p[j];
            }
        }
        fptt=(*func)(ptt);
        if (fptt < fp)
        {
            t=(double)(2.0*(fp-2.0*(*fret)+fptt)*ipow((fp-(*fret)-del),2)-del*ipow((fp-fptt),2));
            if (t < 0.0)
            {
                if(RESTRICT == 1)
                {
                    linmin(1,INTEGER,p_restrict,p,xit,n,tol,small,itmax,fret,func);
                }
                else
                {
                    linmin(0,INTEGER,p,p,xit,n,tol,small,itmax,fret,func);
                }
                for (j=1;j<=n;j++)
                {
                    if((RESTRICT == 1 ) && (p_restrict[j] == 0.0))
                    {
/* Only change for unrestricted vars */
                        xi[j][ibig]=xi[j][n];
                        xi[j][n]=xit[j];
                    }
                    else
                    {
/* Going nowhere... */
                        xi[j][ibig]=(double)0.0;
                        xi[j][n]=(double)0.0;
                    }
                }
            }
        }
    }
}


void nrerror(char error_text[])
{
    printf("Numerical Recipes run-time error...\n");
    printf("%s\n",error_text);
    printf("...now exiting to system...\n");
    exit(1);
}


double *vector(int nl,int nh)
{
    double *v;

    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if (!v) nrerror("allocation failure in vector()");
    return(v-nl);
}


void free_vector(double *v, int nl, int nh)
{
    nh = nh;                                      /* Simply to avoid warning about unused variable */
    free((char*) (v+nl));
}




int ncom=0;                                       /* defining declarations */
double *pcom=0,*xicom=0,(*nrfunc)(double *);

/* 6/18/98 added DBL_MAX/2.0 for integer search - works like a charm */

void  linmin(int RESTRICT, int INTEGER, double *p_restrict, double *p, double *xi, int n, double tol, double small, int itmax, double *fret, double (*func)(double *))
{
    int j;
    double xx,xmin,fx,fb,fa,bx,ax;

/* Only changes _non zero_ restrictions... not a general approach
   but valid for the current purpose */

    if(RESTRICT == 1)
    {
        for(j=1;j<=n;j++)
        {
            if(p_restrict[j] > 0.0)
            {
                p[j]=p_restrict[j];               /* Initialize parameters to restricted values */
            }
        }
    }

    ncom=n;
    pcom=vector(1,n);
    xicom=vector(1,n);
    nrfunc=func;

    for (j=1;j<=n;j++)
    {
        pcom[j]=p[j];
        xicom[j]=xi[j];
    }

    if(INTEGER == 0)
    {

        ax=(double)0.0;
        xx=(double)1.0;
        bx=(double)2.0;

    }
    else
    {

        ax=(double)0.0;
        xx=(double)DBL_MAX/2.0;
        bx=(double)DBL_MAX;

    }

    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,small,f1dim);
    *fret=brent(ax,xx,bx,f1dim,tol,small,itmax,&xmin);
    for (j=1;j<=n;j++)
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_vector(xicom,1,n);
    free_vector(pcom,1,n);
}


double brent(double ax, double bx, double cx, double (*f)(double), double tol, double small, int itmax, double *xmin)
{
    int iter;
    double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    double d = 0.0;
    double e = 0.0;

    a=((ax < cx) ? ax : cx);
    b=((ax > cx) ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=(*f)(x);
    for (iter=1;iter<=itmax;iter++)
    {
        xm=(double) (0.5*(a+b));
        tol2=(double) (2.0*(tol1=(double) (tol*fabs(x)+small)));
        if (fabs(x-xm) <= (tol2-0.5*(b-a)))
        {
            *xmin=x;
            return fx;
        }
        if (fabs(e) > tol1)
        {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=(double)(2.0*(q-r));
            if (q > 0.0) p = -p;
            q=(double) fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                d=(double) (CGOLD*(e=(x >= xm ? a-x : b-x)));
            else
            {
                d=p/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=(double)SIGN(tol1,xm-x);
            }
        }
        else
        {
            d=(double) (CGOLD*(e=(x >= xm ? a-x : b-x)));
        }
        u=(double)(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        fu=(*f)(u);
        if (fu <= fx)
        {
            if (u >= x) a=x; else b=x;
            SHFT(v,w,x,u)
                SHFT(fv,fw,fx,fu)
        }
        else
        {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x)
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v=u;
                fv=fu;
            }
        }
    }
/*  nrerror("Too many iterations in BRENT");*/
    if(int_VERBOSE == 1)
    {
        printf("\n**Maximum number of iterations reached in routine BRENT\n");
    }
    *xmin=x;
    return fx;
}


#include <math.h>

void  mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
double *fc, double small, double (*func)(double))
{
    double ulim,u,r,q,fu,dum;

    *fa=(*func)(*ax);
    *fb=(*func)(*bx);
    if (*fb > *fa)
    {
        SHFT(dum,*ax,*bx,dum)
            SHFT(dum,*fb,*fa,dum)
    }
    *cx=(double)((*bx)+GOLD*(*bx-*ax));
    *fc=(*func)(*cx);
    while (*fb > *fc)
    {
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u= (double)((*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
            (2.0*SIGN(MAX(fabs(q-r),small),q-r)));
        ulim=(double)((*bx)+GLIMIT*(*cx-*bx));
        if ((*bx-u)*(u-*cx) > 0.0)
        {
            fu=(*func)(u);
            if (fu < *fc)
            {
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                return;
            }
            else if (fu > *fb)
            {
                *cx=u;
                *fc=fu;
                return;
            }
            u=(double) ((*cx)+GOLD*(*cx-*bx));
            fu=(*func)(u);
        }
        else if ((*cx-u)*(u-ulim) > 0.0)
        {
            fu=(*func)(u);
            if (fu < *fc)
            {
                SHFT(*bx,*cx,u,(double)(*cx+GOLD*(*cx-*bx)))
                    SHFT(*fb,*fc,fu,(*func)(u))
            }
        }
        else if ((u-ulim)*(ulim-*cx) >= 0.0)
        {
            u=ulim;
            fu=(*func)(u);
        }
        else
        {
            u= (double) ((*cx)+GOLD*(*cx-*bx));
            fu=(*func)(u);
        }
        SHFT(*ax,*bx,*cx,u)
            SHFT(*fa,*fb,*fc,fu)
    }
}


extern int ncom;                                  /* defined in LINMIN */
extern double *pcom,*xicom,(*nrfunc)(double *);

double f1dim(double x)
{
    int j;
    double f,*xt;

    xt=vector(1,ncom);
    for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    f=(*nrfunc)(xt);
    free_vector(xt,1,ncom);
    return f;
}
/* This is an approximation to the error function good to 1.2e-07.
     Compared to erfun() above it yields exact output for CDF
     estimation, and is substantially faster for computation. Now
     default replacement as of 10/6/03 */

double erfun(double x)
{
    double t,z,ans;

    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
        t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
        t*(-0.82215223+t*0.17087277)))))))));
    return(x >= 0.0 ? -(ans-1.0) : (ans-1.0));
}


/* This version uses dlinmin(), line search using derivatives */


#undef SIGN
#undef MOV3
#undef CGOLD

#undef GOLD
#undef GLIMIT
#undef MAX
#undef SHFT

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

#undef FREEALL


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

int iff = 0;

double ran3(int *idum)
{
    static int inext,inextp;
    static long ma[56];
    //static int iff=0;
    long mj,mk;
    int i,ii,k;

    if (*idum < 0 || iff == 0)
    {
        iff=1;
        mj=MSEED-(*idum < 0 ? -*idum : *idum);
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;i++)
        {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1;k<=4;k++)
            for (i=1;i<=55;i++)
        {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
        }
        inext=0;
        inextp=31;
        *idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;
    return(mj*FAC);
}


#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
