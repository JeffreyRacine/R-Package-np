/* This module contains the functions for the kernel function. */

/* Copyright (C) J. Racine, 1995-1997 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>

#include "headers.h"

#ifdef MPI2

#include "mpi.h"

extern  int my_rank;
extern  int source;
extern  int dest;
extern  int tag;
extern  int iNum_Processors;
extern  int iSeed_my_rank;
extern  MPI_Status status;
#endif

#ifdef RCSID
static char rcsid[] = "$Id: kernel.c,v 1.3 2006/11/02 16:56:49 tristen Exp $";
#endif

/*

The following kernel functions are supported.

To add additional kernel functions add the kernel
to kernel.c, add the necessary constants e.g.
{int K(z)^2 dz}, and modify kernelb.c to reflect the
convergence rate of the bandwidth.

0 = Second Order Gaussian kernel
1 = Fourth Order Gaussian kernel
2 = Sixth Order Gaussian kernel
3 = Eighth Order Gaussian kernel

4 = Second Order Epanechnikov
5 = Fourth Order Epanechnikov
6 = Sixth Order Epanechnikov
7 = Eighth Order Epanechnikov

8 = Rectangular kernel

The following techniques are currently supported.

The technique is input via disk files bandwidth_d.ini and bandwidth_r.ini
for density and regression.

0 = Fixed Bandwidth Estimator
1 = Generalized Nearest-Neighbor Bandwidth Estimator (Silverman, pg. 51)
2 = Adaptive Estimator Employing Nearest-Neighbor Bandwidth (Abramson)

*/

double kernel(int KERNEL, double z)
{

/* Evaluate the kernel function */

  double z_squared;
	double one_over_sqrt_two_pi = 0.39894228040143267794;
	double return_value = 0.0;

	switch(KERNEL)
	{

		case 0:

/* Second Order Gaussian (Standard Kernel) */

			return_value = one_over_sqrt_two_pi*exp(-0.5*ipow(z,2));

			break;

		case 1:

/* Fourth Order Gaussian */

      z_squared = ipow(z,2);

			return_value = one_over_sqrt_two_pi*(1.5-0.5*z_squared)*exp(-0.5*z_squared);

			break;

		case 2:

/* Sixth Order Gaussian */

      z_squared = ipow(z,2);

			return_value = one_over_sqrt_two_pi*(1.875-1.25*z_squared+0.125*ipow(z,4))*exp(-0.5*z_squared);

			break;

		case 3:

/* Eighth Order Gaussian */

      z_squared = ipow(z,2);

			return_value = one_over_sqrt_two_pi*(2.1875-2.1875*z_squared+0.4375*ipow(z,4)-0.02083333333*ipow(z,6))*exp(-0.5*z_squared);

			break;

		case 4:

/* Second Order Epanechnikov */
/* Note that return value is preset to 0 so no ifelse necessary */

      z_squared = ipow(z,2);

			if (z_squared < 5.0) return_value = (double)(0.33541019662496845446-0.067082039324993690892*z_squared);

			break;

		case 5:

/* Fourth Order Epanechnikov */

      z_squared = ipow(z,2);

			if (z_squared < 5.0) return_value = (double)(0.008385254916*(-15.0+7.0*z_squared)*(-5.0+z_squared));

			break;

		case 6:

/* Sixth Order Epanechnikov */

      z_squared = ipow(z,2);

			if (z_squared < 5.0) return_value = (double)(0.33541019662496845446*(2.734375+z_squared*(-3.28125+0.721875*z_squared))*(1.0-0.2*z_squared));

			break;

		case 7:

/* Eighth Order Epanechnikov */

      z_squared = ipow(z,2);

			if (z_squared < 5.0) return_value = (double)(0.33541019662496845446*(3.5888671875+z_squared*(-7.8955078125+z_squared*(4.1056640625-0.5865234375*z_squared)))
                                                   *(1.0-0.2*z_squared));

			break;

		case 8:

/* Rectangular Kernel */
/* Note that return value is preset to 0 so no ifelse necessary */

			if(fabs(z) < 1.0)	return_value = 0.5;


			break;

	}

	return(return_value);

}


double cdf_kernel(int KERNEL, double z)
{

/* Evaluate the cumulative kernel function */

/* Define z_squared used by all kernels */

	double z_squared;
	double sqrt_5 = 2.236067978;
	double return_value = 0.0;

	switch(KERNEL)
	{

		case 0:

/* Second Order Gaussian (Standard Kernel) */

			return_value = 0.5*erfun((double) (0.7071067810*z))+0.5;

			break;

		case 1:

/* Fourth Order Gaussian */

			return_value = 0.5*erfun((double)(0.7071067810*z))+0.1994711401*z*exp(-0.5*ipow(z,2))+0.5;

			break;

		case 2:

/* Sixth Order Gaussian */

      z_squared = ipow(z,2);

			return_value = 0.5*erfun((double)(0.7071067810*z))
				+0.3490744952*z*exp(-0.5*z_squared)
				-0.04986778504*exp(-0.5*z_squared)*ipow(z,3)+0.5;

			break;

		case 3:

/* Eighth Order Gaussian */

      z_squared = ipow(z,2);

			return_value = 0.5*erfun((double)(0.7071067810*z))+0.4737439578*z*exp(-0.5*z_squared)
				-0.1329807601*exp(-0.5*z_squared)*ipow(z,3)
				+0.008311297511*exp(-0.5*z_squared)*ipow(z,5) + 0.5;

			break;

		case 4:

/* Second Order Epanechnikov */

			if (z < -sqrt_5)
			{
				return_value = 0.0;
			}
			else if(z < sqrt_5)
			{
				return_value = (double) (0.3354101967*z-0.02236067978*ipow(z,3)+0.5);
			}
			else
			{
				return_value = 1.0;
			}

			break;

		case 5:

/* Fourth Order Epanechnikov */
			if (z < -sqrt_5)
			{
				return_value = 0.0;
			}
			else if(z < sqrt_5)
			{
				return_value = (double)(0.01173935688*ipow(z,5)
					-0.1397542486*ipow(z,3)+0.6288941188*z+0.5);
			}
			else
			{
				return_value = 1.0;
			}

			break;

		case 6:

/* Sixth Order Epanechnikov */
			if (z < -sqrt_5)
			{
				return_value = 0.0;
			}
			else if(z < sqrt_5)
			{
				return_value = (double)(-0.006917835307*ipow(z,7)
					+0.09244743547*ipow(z,5)
					-0.4279973864*ipow(z,3)+0.9171372566*z+0.5);
			}
			else
			{
				return_value = 1.0;
			}

			break;

		case 7:

/* Eighth Order Epanechnikov */
			if (z < -sqrt_5)
			{
				return_value = 0.0;
			}
			else if(z < sqrt_5)
			{
				return_value = (double)(0.004371687590*ipow(z,9)
					-0.06744889424*ipow(z,7)
					+0.3813456714*ipow(z,5)
					-0.9629941194*ipow(z,3)
					+1.203742649*z+0.5);
			}
			else
			{
				return_value = 1.0;
			}

			break;

		case 8:

/* Rectangular Kernel - bug Oct 10, 2007 - was using sqrt_5 in second call... */

			if(z < -1.0)
			{
				return_value = 0.0;
			}
			else if(z < 1.0)
			{
				return_value = 0.5 + 0.5*z;
			}
			else
			{
				return_value = 1.0;
			}

			break;

	}

	return(return_value);

}


/* This is the derivative of the kernel with respect to Z */

double kernel_deriv(int KERNEL, double z)
{

/* Evaluate the derivative of the kernel function */

/* Define z_squared used by all kernels */

	double z_squared = ipow(z,2);
	double return_value = 0.0;

	double one_over_sqrt_two_pi = 0.39894228040143267794;

	switch (KERNEL)
	{

		case 0:

/* Second Order Gaussian */

			return_value = -z*one_over_sqrt_two_pi*exp(-0.5*z_squared);

			break;

		case 1:

/* Fourth Order Gaussian */

			return_value = -one_over_sqrt_two_pi*(2.5*z-0.5*ipow(z,3))*exp(-0.5*z_squared);

			break;

		case 2:

/* Sixth Order Gaussian */

			return_value = -0.049867785050179084743*z*exp(-0.5*z_squared)*(35.0-14.0*z_squared+ipow(z,4));

			break;

		case 3:

/* Eighth Order Gaussian */

			return_value = -one_over_sqrt_two_pi*(6.5625*z-3.9375*ipow(z,3)+0.5625*ipow(z,4)*z-0.02083333333*ipow(z,6)*z)*exp(-0.5*z_squared);

			break;

		case 4:

/* Second Order Epanechnikov */

			if (z_squared < 5.0)
			{
				return_value = (double)(-0.13416407864998738178*z);
			}
			else
			{
				return_value = 0.0;
			}

			break;

		case 5:

/* Fourth Order Epanechnikov */
			if (z_squared < 5.0)
			{
				return_value = (double)(-.58696784409369479531*z*(1.0-0.2*z_squared)
					-0.13416407864998738178*(1.875-0.875*z_squared)*z);
			}
			else
			{
				return_value = 0.0;
			}

			break;

		case 6:

/* Sixth Order Epanechnikov */
			if (z_squared < 5.0)
			{
				return_value = -1.0022916396047925e-10*z*(z_squared*(2898847705.0*z_squared-18447212816.0)+25621128780.0);
			}
			else
			{
				return_value = 0.0;
			}

			break;

		case 7:

/* Eighth Order Epanechnikov */
			if (z_squared < 5.0)
			{
				return_value = (double)(-5.7779647152973081413*z
					+ 7.6269134241924467465 *ipow(z,3)
					- 2.8328535575571945058 *ipow(z,5)
					+ 0.31476150639524383398 *ipow(z,7));
			}
			else
			{
				return_value = 0.0;
			}

			break;

		case 8:

/* Rectangular kernel added 4/27/99 */
/* Derivative always zero except at +-1 where it is undefined. */

			return_value = 0.0;

			break;

	}

	return(return_value);

}


/* Completed 11/25/97 */

/* Revisited for two fucking days 5/19/99 - 5/20/99 for grades data set. */
/* No bugs at all - data has numerous ties!!! (Not continuous). */
/* However, beefed up by using 20 decimals for all kernels. */

double kernel_convol(int KERNEL, int BANDWIDTH, double z, double h1, double h2)
{

/* Evaluate the convolution kernel function */

/* Define z_squared used by all kernels */

	double z_squared = ipow(z,2);
	double return_value = 0.0;
	double s1;
	double s2;
	double s3;
	double s4;

	if((BANDWIDTH == 0)||(BANDWIDTH == 1))
	{

/* Convolution kernels when bandwidth is for evaluation point.
					 Fixed and Generalized NN */

		switch(KERNEL)
		{

			case 0:

/* Second Order Gaussian (Standard Kernel) */

				return_value = 0.28209479177387814348*exp(-0.25*z_squared);

				break;

			case 1:

/* Fourth Order Gaussian */

				return_value = 0.0044077311214668459918*exp(-0.25*z_squared)*(108.0-28.0*z_squared+ipow(z,4));

				break;

			case 2:

/* Sixth Order Gaussian */

				return_value = 0.00001721769969*exp(-0.25*z_squared)*(36240.0-19360.0*z_squared+2312.0*ipow(z,4)-88.0*ipow(z,6)+ipow(z,8));

				break;

			case 3:

/* Eighth Order Gaussian */

				return_value = 0.2989183974E-7*exp(-0.25*z_squared)*(25018560.0-20462400.0*z_squared+4202352.0*
					ipow(z,4)-331680.0*ipow(z,6)+11604.0*ipow(z,8)-180.0*ipow(z,10)+ipow(z,12));

				break;

			case 4:

/* Second Order Epanechnikov */

/* Error pointed out by Alicia Pérez-Alonso and colleagues May 22
   2009, previous version incorporated incorrect limits of
   integration. Confirmed May 27 and incorporated. New version
   requires additional test for negative and positive values, and
   reflects limits of -2*sqrt(5) and 2*sqrt(5) rather than -sqrt(5)
   and sqrt(5). Thanks ever so much Alicia! */

        return((z_squared < 20.0) ? 
               ((z < 0.0) ? 
                (5.579734404642339E-9*(26883*ipow(z,5)-2688300*ipow(z,3)-12022443*z_squared+48089773)) : 
                (-5.579734404642339E-9*(26883*ipow(z,5)-2688300*ipow(z,3)+12022443*z_squared-48089773)))
               : 0.0);
        
				break;

			case 5:

/* Fourth Order Epanechnikov */

        return((z_squared < 20.0) ? 
               ((z < 0.0) ? 
                (3.756009615384615e-9*(1456*ipow(z,9)-124800*ipow(z,7)+5491200*ipow(z,5)+15627432*ipow(z,4)-24960000*ipow(z,3)-111624513*z_squared+148832684))           :
                (-3.756009615384615e-9*(1456*ipow(z,9)-124800*ipow(z,7)+5491200*ipow(z,5)-15627432*ipow(z,4)-24960000*ipow(z,3)+111624513*z_squared-148832684)))
               : 0.0);

				break;

			case 6:

/* Sixth Order Epanechnikov */

        return((z_squared < 20.0) ? 
               ((z < 0.0) ? 
                (9.390024038461537E-11*(2079*ipow(z,13)-206388*ipow(z,11)+8867040*ipow(z,9)-255528000*ipow(z,7)-515705252*ipow(z,6)+1681680000*ipow(z,5)+4922641042*ipow(z,4)-3057600000*ipow(z,3)-13674002896*z_squared+9015826085)) :
                (-9.390024038461537E-11*(2079*ipow(z,13)-206388*ipow(z,11)+8867040*ipow(z,9)-255528000*ipow(z,7)+515705252*ipow(z,6)+1681680000*ipow(z,5)-4922641042*ipow(z,4)-3057600000*ipow(z,3)+13674002896*z_squared-9015826085)))
               : 0.0);
        
				break;

			case 7:

/* Eighth Order Epanechnikov */

        return((z_squared < 20.0) ? 
               ((z < 0.0) ? 
                (1.121969784007353E-13*(63063*ipow(z,17)-7351344*ipow(z,15)+373222080*ipow(z,13)-11040382080*ipow(z,11)+241727270400*ipow(z,9)+350679571413*ipow(z,8)-1900039680000*ipow(z,7)-4208154856956*ipow(z,6)+5757696000000*ipow(z,5)+16994471537707*ipow(z,4)-5757696000000*ipow(z,3)-25749199299557*z_squared+10097725215512)) :
                (-1.121969784007353E-13*(63063*ipow(z,17)-7351344*ipow(z,15)+373222080*ipow(z,13)-11040382080*ipow(z,11)+241727270400*ipow(z,9)-350679571413*ipow(z,8)-1900039680000*ipow(z,7)+4208154856956*ipow(z,6)+5757696000000*ipow(z,5)-16994471537707*ipow(z,4)-5757696000000*ipow(z,3)+25749199299557*z_squared-10097725215512)))
               : 0.0);
        
				break;

			case 8:
/* Rectangular kernel */
/* Need to verify 4/27/99 */

				if(z_squared < 1.0)
				{
					return_value = 0.5/(h1*h2);
				}
				else
				{
					return_value = 0.0;
				}

				break;

		}

	}
	else if(BANDWIDTH==2)
	{

/* Convolution kernels when bandwidth is for training point.
					 Adaptive */

/* h1 and h2 are the bandwidths for the evaluation and training points */

		switch(KERNEL)
		{

			case 0:

/* Second Order Gaussian (Standard Kernel) */

				return_value = 0.3989422803*exp(-0.5*z_squared/(h2*h2+h1*h1))/sqrt(h2*h2+h1*h1);

				break;

			case 1:

/* Fourth Order Gaussian */

				return_value = 0.9973557006E-1*exp(-0.5*z_squared/(h2*h2+h1*h1))*(-12.0*h1*h1*ipow(h2,4)
					*z_squared-12.0*ipow(h1,4)*z_squared*h2*h2+ipow(z,4)*h1*h1*h2*h2+27.0
					*ipow(h2,6)*h1*h1+42.0*ipow(h2,4)*ipow(h1,4)-2.0*ipow(h2,6)
					*z_squared+27.0*ipow(h1,6)*h2*h2-2.0*z_squared*ipow(h1,6)
					+6.0*ipow(h2,8)+6.0*ipow(h1,8))/sqrt(pow(h2*h2+h1*h1,9.0));

				break;

			case 2:

/* Sixth Order Gaussian */

				return_value = 0.6233473129E-2*exp(-0.5*z_squared/(h2*h2+h1*h1))*(-4.0*ipow(z,6)*ipow(h1,
					8)*h2*h2+ipow(z,8)*ipow(h1,4)*ipow(h2,4)-4.0*ipow(z,6)*h1*h1*ipow(h2,8)
					-40.0*ipow(z,6)*ipow(h1,4)*ipow(h2,6)+108.0*ipow(z,4)*ipow(h1,10)*h2*h2
					-40.0*ipow(z,6)*ipow(h1,6)*ipow(h2,4)-740.0*ipow(h2,12)*z_squared*h1*h1-3000.0*
					pow(h2,10.0)*z_squared*ipow(h1,4)-5860.0*ipow(h1,6)*ipow(h2,8)*z_squared+108.0*h1*h1*
					pow(h2,10.0)*ipow(z,4)-5860.0*ipow(h1,8)*ipow(h2,6)*z_squared+570.0*ipow(h1,4)*
					ipow(h2,8)*ipow(z,4)-3000.0*ipow(h1,10)*ipow(h2,4)*z_squared+940.0*ipow(h1,6)*
					ipow(h2,6)*ipow(z,4)-740.0*ipow(h1,12)*z_squared*h2*h2+570.0*ipow(h1,8)*ipow(z,4)
					*ipow(h2,4)+120.0*ipow(h2,16)+120.0*ipow(h1,16)+1020.0*ipow(h2,14)*h1*h1+
					3825.0*ipow(h2,12)*ipow(h1,4)+8040.0*ipow(h2,10)*ipow(h1,6)+10230.0*ipow(h2,
					8)*ipow(h1,8)-80.0*ipow(h2,14)*z_squared+8.0*ipow(h2,12)*ipow(z,4)+8040.0*ipow(
					h1,10)*ipow(h2,6)+3825.0*ipow(h1,12)*ipow(h2,4)+1020.0*ipow(h1,14)*h2*h2
					-80.0*z_squared*ipow(h1,14)+8.0*ipow(z,4)*ipow(h1,12))/sqrt(pow(h2*h2+h1*h1,17.0));

				break;

			case 3:

/* Eighth Order Gaussian */

				s1 = 0.1731520314E-3*exp(-0.5*z_squared/(h2*h2+h1*h1));
				s4 = 1008.0*ipow(z,4)*ipow(h1,20)+5040.0*ipow(h1,24)+5040.0*ipow(h2,
					24)-5040.0*ipow(h2,22)*z_squared-48.0*ipow(z,6)*ipow(h1,18)-5040.0*z_squared*ipow(h1,
					22)+2983050.0*ipow(h1,16)*ipow(h2,8)+362250.0*ipow(h2,20)*ipow(h1,4)+
					1008.0*ipow(h2,20)*ipow(z,4)+63000.0*ipow(h2,22)*h1*h1-5038110.0*ipow(h1,12
					)*ipow(h2,10)*z_squared-48.0*ipow(h2,18)*ipow(z,6)+2983050.0*ipow(h2,16)*ipow(h1,
					8)+362250.0*ipow(h1,20)*ipow(h2,4)+63000.0*ipow(h1,22)*h2*h2+5808600.0*
					ipow(h1,12)*ipow(h2,12)+4923765.0*ipow(h2,14)*ipow(h1,10)+1267875.0*ipow(h2,
					18)*ipow(h1,6)-105588.0*ipow(h1,8)*ipow(h2,10)*ipow(z,6)-49224.0*ipow(h1,
					6)*ipow(h2,12)*ipow(z,6)+1246266.0*ipow(h1,10)*ipow(h2,10)*ipow(z,4)
					-9876.0*ipow(h1,4)*ipow(h2,14)*ipow(z,6)+947520.0*ipow(h1,8)*ipow(h2,12)*
					ipow(z,4)-1104.0*h1*h1*ipow(h2,16)*ipow(z,6)+412335.0*ipow(h1,6)*ipow(h2,
					14)*ipow(z,4)-5038110.0*ipow(h1,10)*ipow(h2,12)*z_squared+102060.0*ipow(h2,16)*
					ipow(z,4)*ipow(h1,4)+15120.0*ipow(h2,18)*ipow(z,4)*h1*h1-3311280.0*ipow(h2,
					14)*z_squared*ipow(h1,8);
				s3 = s4-1420020.0*ipow(h2,16)*z_squared*ipow(h1,6)-391230.0*ipow(h2,18)*z_squared*
					ipow(h1,4)-65520.0*ipow(h2,20)*z_squared*h1*h1-49224.0*ipow(h1,12)*ipow(z,6)*ipow(
					h2,6)+102060.0*ipow(h1,16)*ipow(z,4)*ipow(h2,4)+412335.0*ipow(h1,14)*
					ipow(z,4)*ipow(h2,6)+947520.0*ipow(h1,12)*ipow(z,4)*ipow(h2,8)-65520.0*
					ipow(h1,20)*z_squared*h2*h2-391230.0*ipow(h1,18)*z_squared*ipow(h2,4)-1420020.0*ipow(h1,
					16)*z_squared*ipow(h2,6)-3311280.0*ipow(h1,14)*z_squared*ipow(h2,8)+1267875.0*ipow(h1,
					18)*ipow(h2,6)+4923765.0*ipow(h1,14)*ipow(h2,10)+ipow(z,12)*ipow(h1,6)*
					ipow(h2,6)-6.0*ipow(z,10)*ipow(h1,10)*ipow(h2,4)-84.0*ipow(z,10)*ipow(h1,
					8)*ipow(h2,6)-84.0*ipow(z,10)*ipow(h1,6)*ipow(h2,8)+24.0*ipow(z,8)*ipow(
					h1,14)*h2*h2+402.0*ipow(z,8)*ipow(h1,12)*ipow(h2,4)+2877.0*ipow(z,8)*
					ipow(h1,10)*ipow(h2,6)+4998.0*ipow(z,8)*ipow(h1,8)*ipow(h2,8)+2877.0*ipow(z
					,8)*ipow(h1,6)*ipow(h2,10)-1104.0*ipow(z,6)*ipow(h1,16)*h2*h2-9876.0*ipow(
					z,6)*ipow(h1,14)*ipow(h2,4)-105588.0*ipow(h1,10)*ipow(h2,8)*ipow(z,6)+
					15120.0*ipow(z,4)*ipow(h1,18)*h2*h2+402.0*ipow(z,8)*ipow(h1,4)*ipow(h2,12)
					-6.0*ipow(z,10)*ipow(h1,4)*ipow(h2,10)+24.0*ipow(z,8)*h1*h1*ipow(h2,14);
				s4 = 1/(sqrt(pow(h2*h2+h1*h1,25.0)));
				s2 = s3*s4;
				return_value = s1*s2;

				break;

			case 4:

/* Second Order Epanechnikov */
/* These are wrong and need to be corrected */

				if (z_squared < 5.0)
				{
					return_value = 0.3354101967E-1*(-3.0+5.0*h1*h1-z_squared+5.0*h2*h2-15.0*h1*h1*h2*h2
						+3.0*z_squared*h2*h2)/(h1*h1*h1)/(h2*h2*h2);
				}
				else
				{
					return_value = 0.0;
				}

				break;

			case 5:

/* Fourth Order Epanechnikov */
/* These are wrong and need to be corrected */
				if (z_squared < 5.0)
				{
					return_value = 0.1746928108E-3*(6125.0-6750.0*ipow(h2,4)*h1*h1*z_squared
						+7500.0*h2*h2*h1*h1*z_squared-11250.0*h1*h1+9450.0*z_squared
						-11250.0*h2*h2+4725.0*ipow(h1,4)+441.0*ipow(z,4)
						+4725.0*ipow(h2,4)-18900.0*z_squared*h2*h2-3150.0*h1*h1*z_squared
						+22500.0*h2*h2*h1*h1+10125.0*ipow(h1,4)*ipow(h2,4)
						+945.0*ipow(h2,4)*ipow(z,4)-11250.0*ipow(h1,4)*h2*h2
						-1050.0*h2*h2*ipow(z,4)-11250.0*ipow(h2,4)*h1*h1
						+9450.0*ipow(h2,4)*z_squared)/pow(h1,5.0)/pow(h2,5.0);

				}
				else
				{
					return_value = 0.0;
				}

				break;

			case 6:

/* Sixth Order Epanechnikov */
/* These are wrong and need to be corrected */
				if (z_squared < 5.0)
				{
					s2 = 0.4965438428137688E1*ipow(h1,4)*z_squared-0.1787557834129568E1*ipow(z,4)
						*h1*h1-0.4739736681404156E2*h2*h2*ipow(h1,4)-0.2681336751194351E2*h2*h2*ipow(z,
						4)-0.2843842008842494E2*h2*h2+0.6635631353965819E2*h2*h2*h1*h1
						-0.1042742069908914E3*h2*h2*z_squared-0.1489631528441306E2*ipow(h2,6)*z_squared+
						0.9479473362808313E1*ipow(h1,6)*h2*h2-0.2843842008842494E2*h1*h1+
						0.1931003833164656E2*ipow(h1,4)+0.3761695778892188E1*ipow(h1,6)*ipow(h2,6)
						-0.3546741734384063E1*ipow(h1,6)-0.4170968279635658E2*z_squared*h1*h1
						-0.8777290150748438E1*ipow(h1,6)*ipow(h2,4)+0.1092396454190291E2*ipow(z,4)+
						0.4468894585323919E2*z_squared+0.7448157642206531E2*z_squared*ipow(h2,4)+
						0.1931003833164656E2*ipow(h2,4)-0.4739736681404156E2*ipow(h2,4)*h1*h1;

					s1 = s2+0.1023783123183298E3*h2*h2*z_squared*h1*h1+0.1872679635754785*ipow(z,6
						)+0.3686461863314344E2*ipow(h2,4)*ipow(h1,4)+0.1895894672561663E1*ipow(h2,6)
						*ipow(z,4)*h1*h1-0.3546741734384063E1*ipow(h2,6)+0.1228820621104781E2*ipow(h2,
						4)*ipow(h1,4)*z_squared-0.4423754235977213E1*ipow(h2,4)*ipow(z,4)*h1*h1+
						0.4634409199595175*ipow(h2,4)*ipow(z,6)-0.1327126270793164E2*h2*h2*ipow(h1,4
						)*z_squared-0.5266374090449063E1*ipow(h2,6)*ipow(h1,4)*z_squared+0.1895894672561663E2*
						ipow(h2,6)*z_squared*h1*h1+0.477765457485539E1*h2*h2*ipow(z,4)*h1*h1
						-0.7962757624758983E2*ipow(h2,4)*z_squared*h1*h1-0.8777290150748438E1*ipow(h2,6)*
						ipow(h1,4)-0.4965438428137688E1*ipow(h2,6)*ipow(z,4)-0.1986175371255075*ipow(
						h2,6)*ipow(z,6)-0.5005161935562789*h2*h2*ipow(z,6)+0.9479473362808313E1*
						ipow(h2,6)*h1*h1+0.2085484139817829E2*ipow(h2,4)*ipow(z,4)+
						0.1260457447142644E2;

					s2 = 1/pow(h1,7.0)/pow(h2,7.0);

					return_value = s1*s2;

				}
				else
				{
					return_value = 0.0;
				}

				break;

			case 7:

/* Eighth Order Epanechnikov */
/* These are wrong and need to be corrected */
				if (z_squared < 5.0)
				{

					s3 = 0.622090439547381E2*ipow(h2,6)*ipow(h1,6)*z_squared+0.5824692032190657E3
						*ipow(z,4)+0.1470881826310772E2*ipow(h1,8)-0.1444138520377849E3*ipow(h1,6)+
						0.4032479099208916E3*ipow(h1,4)+0.550698155770753E2*ipow(z,6)
						-0.4326914081055916E3*h1*h1-0.4326914081055916E3*h2*h2+0.4032479099208916E3*
						ipow(h2,4)-0.1444138520377849E3*ipow(h2,6)+0.4236139659775023E3*h2*h2*ipow(h1,
						6)+0.9884325872808387E1*ipow(h2,8)*ipow(z,6)+0.1129094147778497E4*ipow(h2,
						4)*ipow(h1,4)-0.5769218774741222E1*ipow(z,6)*h1*h1+0.1630913769013384E4*
						ipow(h2,4)*ipow(z,4)-0.3540202429954841E3*ipow(z,4)*h1*h1+
						0.1470881826310772E2*ipow(h2,8)-0.1143757708139256E4*ipow(h2,4)*h1*h1;

					s2 = s3+0.2668767985658264E4*ipow(h2,4)*z_squared-0.1615381256927542E3*h2*h2*
						ipow(z,6)+0.4236139659775023E3*ipow(h2,6)*h1*h1-0.9884325872808387E3*z_squared*ipow(
						h2,6)-0.3530116383145853E2*ipow(h1,6)*z_squared-0.149777795113474E4*z_squared*h1*h1
						-0.1143757708139256E4*h2*h2*ipow(h1,4)-0.1652094467312259E4*h2*h2*ipow(z,4)+
						0.5718788540696281E3*z_squared*ipow(h1,4)-0.4538721064044668E2*h2*h2*ipow(h1,8)
						-0.1483513399219171E1*h2*h2*ipow(z,8)-0.4399068108227909E3*ipow(h2,4)*ipow(h1,
						6)+0.1677511305270909E3*ipow(h2,4)*ipow(z,6)-0.4399068108227909E3*ipow(h2,
						6)*ipow(h1,4)-0.6354209489662535E3*ipow(h2,6)*ipow(z,4)+
						0.1059034914943756E3*ipow(h2,8)*z_squared-0.4538721064044668E2*ipow(h2,8)*h1*h1+
						0.6480108745285219E1*ipow(h2,8)*ipow(h1,8)+0.2118069829887512*ipow(h2,8)*
						ipow(z,8);
					s3 = s2-0.2592043498114088E2*ipow(h2,6)*ipow(h1,8)-0.8472279319550046*
						ipow(h2,6)*ipow(z,8)-0.2795852175451515E4*h2*h2*z_squared+0.1198222360907792E4*h2*
						h2*h1*h1-0.2592043498114088E2*ipow(h2,8)*ipow(h1,6)+0.2329876812876263E2*
						ipow(z,4)*ipow(h1,4)+0.1677511305270909E1*ipow(h2,4)*ipow(z,8)+
						0.1866271318642143E3*ipow(h2,6)*ipow(h1,6)-0.7116714628422039E2*ipow(h2,6)*
						ipow(z,6)+0.5132246126265893E2*ipow(h2,8)*ipow(h1,4)+0.741324440460629E2*
						ipow(h2,8)*ipow(z,4)+0.4807682312284352*ipow(z,8)+0.1590777236077458E3
						-0.1555226098868453E2*ipow(h2,8)*ipow(h1,6)*z_squared+0.4575030832557025E3*ipow(h2,
						6)*ipow(z,4)*h1*h1-0.190626284689876E3*ipow(h2,8)*z_squared*h1*h1
						-0.7189334165446753E2*h2*h2*ipow(z,4)*ipow(h1,4)+0.1780216079063006E2*h2*h2*
						ipow(z,6)*h1*h1;
					s1 = s3+0.1742030970858252E4*ipow(h2,4)*z_squared*ipow(h1,4)
						-0.1078400124817013E4*ipow(h2,4)*ipow(z,4)*h1*h1+0.1026449225253179E2*ipow(h2,
						8)*ipow(z,4)*ipow(h1,4)-0.2541683795865014E1*ipow(h2,8)*ipow(z,6)*h1*h1
						-0.4105796901012715E2*ipow(h2,6)*ipow(z,4)*ipow(h1,4)+0.108929305537072E3*h2
						*h2*ipow(h1,6)*z_squared+0.5132246126265893E2*ipow(h2,4)*ipow(h1,8)+
						0.1026449225253179E3*ipow(h2,8)*z_squared*ipow(h1,4)+0.1016673518346006E2*ipow(h2,
						6)*ipow(z,6)*h1*h1-0.6354209489662535E2*ipow(h2,8)*ipow(z,4)*h1*h1+
						0.4248242915945809E4*h2*h2*z_squared*h1*h1+0.103845937945342E4*h2*h2*ipow(z,4)*h1*h1
						-0.4193778263177273E4*ipow(h2,4)*z_squared*h1*h1-0.1677511305270909E4*h2*h2*z_squared*ipow(
						h1,4)-0.1231739070303814E3*ipow(h2,4)*ipow(h1,6)*z_squared-0.7390434421822886E3*
						ipow(h2,6)*z_squared*ipow(h1,4)+0.8129477864005175E2*ipow(h2,4)*ipow(z,4)*ipow(h1,
						4)-0.2013013566325091E2*ipow(h2,4)*ipow(z,6)*h1*h1+0.163393958305608E4*ipow(
						h2,6)*z_squared*h1*h1+0.1009613285579714E4*z_squared;

					s2 = 1/pow(h1,9.0)/pow(h2,9.0);
					return_value = s1*s2;

				}
				else
				{
					return_value = 0.0;
				}

				break;

			case 8:
/* Rectangular kernel */
/* Really need to verify!!! This one is almost certainly incorrect 4/27/99 */

				if(z_squared < 1.0)
				{
					return_value = 0.5/(h1*h2);
				}
				else
				{
					return_value = 0.0;
				}

				break;

		}

	}

	return(return_value);

}


int initialize_kernel_regression_asymptotic_constants(
int KERNEL,
int num_reg_continuous,
double *INT_KERNEL_P,
double *K_INT_KERNEL_P,
double *INT_KERNEL_PM_HALF,
double *DIFF_KER_PPM)
{

/* Initialize constants for various kernels required for asymptotic standard errors */

	switch(KERNEL)
	{

		case 0:

/* Second-order Gaussian */

			*INT_KERNEL_P = 0.28209479177387814348;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.21969564473386119853;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.06239914704001694495);

			break;

		case 1:

/* Fourth-order Gaussian */

			*INT_KERNEL_P = 0.47603496111841936711;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.27805230036629307938;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.19798266075212628773);

			break;

		case 2:

/* Sixth-order Gaussian */

			*INT_KERNEL_P = 0.62396943688265038571;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.25618196366213489976;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.36778747322051548595);

			break;

		case 3:

/* Eighth-order Gaussian */

			*INT_KERNEL_P = 0.74785078617543927990;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.19644083574560137818;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.55140995042983790172);

			break;

		case 4:

/* Second-order Epanechnikov */

			*INT_KERNEL_P = 0.26832815729997476357;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.20250390621232470438;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.06582425108765005919);

			break;

		case 5:

/* Fourth-order Epanechnikov */

			*INT_KERNEL_P = 0.55901699437494742410;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.25635637709255874475;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.30266061728238867935);

			break;

		case 6:

/* Sixth-order Epanechnikov */

			*INT_KERNEL_P = 0.84658823667359826246;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.27428761935713012265;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.57230061731646813981);

			break;

		case 7:

/* Eighth-order Epanechnikov */

			*INT_KERNEL_P = 1.1329342579014329689;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.15585854498586945817;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.97707571291556351073);

			break;

		case 8:

			*INT_KERNEL_P = 0.5;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_reg_continuous);
			*INT_KERNEL_PM_HALF = 0.25;
			*DIFF_KER_PPM = (2.0 * (*K_INT_KERNEL_P/ *INT_KERNEL_P) * 0.25);

			break;

	}

	return(0);

}


int initialize_kernel_density_asymptotic_constants(
int KERNEL,
int num_var_continuous,
double *INT_KERNEL_P,
double *K_INT_KERNEL_P)
{

/* Initialize constants for various kernels required for asymptotic standard errors */

	switch(KERNEL)
	{

		case 0:

/* Second-order Gaussian */

			*INT_KERNEL_P = 0.28209479177387814348;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

		case 1:

/* Fourth-order Gaussian */

			*INT_KERNEL_P = 0.47603496111841936711;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

		case 2:

/* Sixth-order Gaussian */

			*INT_KERNEL_P = 0.62396943688265038571;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

		case 3:

/* Eighth-order Gaussian */

			*INT_KERNEL_P = 0.74785078617543927990;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

		case 4:

/* Second-order Epanechnikov */

			*INT_KERNEL_P = 0.26832815729997476357;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

		case 5:

/* Fourth-order Epanechnikov */

			*INT_KERNEL_P = 0.55901699437494742410;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

		case 6:

/* Sixth-order Epanechnikov */

			*INT_KERNEL_P = 0.84658823667359826246;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

		case 7:

/* Eighth-order Epanechnikov */

			*INT_KERNEL_P = 1.1329342579014329689;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

		case 8:

			*INT_KERNEL_P = 0.5;
			*K_INT_KERNEL_P = ipow(*INT_KERNEL_P, num_var_continuous);

			break;

	}

	return(0);

}


double kernel_ordered(int KERNEL, double x, double y, double lambda)
{

	double return_value = 0.0;

	switch(KERNEL)
	{

		case 0:

/* Geometric weights in Ahmad & Cerrito 1994 Journal of Statistical Planning and Inference */
/* Taken from Wang & Van Ryzin (1981) - verified sums to 1 (via their claim that it does so) */

			if (x == y)
			{
				return_value = 1.0-lambda;
			}
			else
			{
				return_value = ipow(lambda,(int)fabs(x-y))*(1.0-lambda)*0.5;
			}

			break;

		case 1:

/* Qi's Kernel (Habbemma?) */

			if (x == y)
			{
				return_value = 1.0;
			}
			else
			{
				return_value = ipow(lambda,(int)fabs(x-y));
			}

			break;

	}

	return(return_value);

}


/* Cumulative version of Aitchenson and Aitken's beautiful categorical
kernel. */

double cdf_kernel_ordered(int KERNEL, double x, double y, double lambda, int c, double *categorical_vals)
{

/* X is eval, Y is train */

	double l;
	double return_value = 0.0;

/* Had used unordered kernel as template. However, this baby did not integrate to one */
/* Now going from -max to max in steps of 1 - Ahmad & Cerrito claim that this must */
/* integrate to onc from -infty to infty - using sample analog */

	for(l = categorical_vals[0]-fabs(categorical_vals[0]-categorical_vals[c-1]); l <= categorical_vals[c-1]; l++)
	{
		if(l <= x)
		{
			return_value += kernel_ordered(KERNEL, l, y, lambda);
		}
	}

	return(return_value);

}


/* My implementation of a categorical convolution kernel - thanks Qi! */

double kernel_ordered_convolution(int KERNEL, double x, double y, double lambda, int c, double *c_vals)
{

	int i;
	double *pc_vals = &c_vals[0];

	double return_value = 0.0;

	for(i=0; i < c; i++)
	{
		return_value += kernel_ordered(KERNEL, x, *pc_vals, lambda)*kernel_ordered(KERNEL, y, *pc_vals, lambda);
		pc_vals++;
	}

	return(return_value);

}


/* Aitchenson and Aitken's beautiful categorical kernel */

double kernel_unordered(int KERNEL, double x, double y, double lambda, int c)
{

	double return_value = 0.0;

	switch(KERNEL)
	{

		case 0:

			if (x == y)
			{
				return_value = 1.0 - lambda;
			}
			else
			{
				return_value = lambda/((double) c - 1.0);
			}

			break;
        case 1:

          if (x == y)
            {
              return_value = 1.0;
            }
          else
            {
              return_value = lambda;
            }
          
          break;


	}

	return(return_value);

}


/* 
   This version of the kernel is suitable for some estimators constructed as
   ratios of kernel sums, such as regressions and conditional densities.
   The kernels return 1 outside of their support.
 */

double kernel_unordered_ratio(int KERNEL, double x, double y, double lambda, int c){

	double return_value = 1.0;
  const double dc = (double) c;

	switch(KERNEL)	{

  case 0: /* Aitchenson and Aitken */
    if(x != y){
      if(x <= c){
				return_value = lambda/(dc - 1.0);
      }
    } else {
      return_value = 1.0 - lambda;
    }
    
    break;
  case 1: /* Li and Racine */
    if((x <= c) && (x != y))
      return_value = lambda;
          
    break;
	}

	return(return_value);

}


/* Cumulative version of Aitchenson and Aitken's beautiful categorical
kernel. */

double cdf_kernel_unordered(int KERNEL, double x, double y, double lambda, int c, double *categorical_vals)
{

/* X is eval, Y is train */

	int l;
	double return_value = 0.0;

	for(l = 0; l < c; l++)
	{
		if(categorical_vals[l] <= x)
		{
			return_value += kernel_unordered(KERNEL, categorical_vals[l], y, lambda, c);
		}
	}

	return(return_value);

}


/* My implementation of a categorical convolution kernel - thanks Qi! */

double kernel_unordered_convolution(int KERNEL, double x, double y, double lambda, int c, double *c_vals)
{

	int i;
	double *pc_vals = &c_vals[0];

	double return_value = 0.0;

	for(i=0; i < c; i++)
	{
		return_value += kernel_unordered(KERNEL, x, *pc_vals, lambda, c)*kernel_unordered(KERNEL, y, *pc_vals, lambda, c);
		pc_vals++;
	}

	return(return_value);

}
