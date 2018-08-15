/*--------------------------------------------------------------------
 * $Id: regress.c 2623 2011-12-23 10:52:38Z robert.buras $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numeric.h"


/***********************************************************************************/
/* Function: gauss                                                        @31_30i@ */
/* Description:                                                                    */
/*  Calculate a Gauss function for given average and standard deviation.           */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double x:      x value where the Gauss function is to be evaluated.            */
/*  double mu:     Average.                                                        */
/*  double sigma:  Standard deviation.                                             */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

double gauss (double x, double mu, double sigma)
{
  double result=0.39894228 / sigma * exp (-(x-mu)*(x-mu) / 2.0 / sigma / sigma);
  return result;
}


/**************************************************************/
/* calculate distorted gausscurve;                            */
/* a distorted gausscurve is an invention of Bernhard Mayer;  */
/* normal gausscurve:  exp (a2*x^2 + a3*x + a4), a2<0;        */ 
/* distorted gausscurve: added a term x^3 to allow asymmetry; */
/**************************************************************/

double gauss_distorted (double x, double a0, double a1, double a2, double a3)
{
  double dummy=0;
  dummy = exp (a3*x*x*x + a2*x*x + a1*x + a0);
  
  return dummy;
}





/***********************************************************************************/
/* Function: regression                                                   @31_30i@ */
/* Description:                                                                    */
/*  Calculate coefficients for y = a + b*x by linear regression.                   */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:            Vector (0..n-1) of x-values.                             */
/*  double *y:            Vector (0..n-1) of y-values.                             */
/*  int n:                Number of data points.                                   */
/*  double *a:            Coefficient a.                                           */ 
/*  double *b:            Coefficient b.                                           */ 
/*  double *sigma_a:      Standard deviation of a.                                 */
/*  double *sigma_b:      Standard deviation of b.                                 */
/*  double *correlation:  Correlation coefficient.                                 */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int regression (double *x, double *y, int n, double *a, double *b,
		   double *sigma_a, double *sigma_b, double *correlation)
{
  int i=0;
  double delta_x=0, delta_y=0;
  double sum_x=0, sum_y=0, sum_xx=0, sum_yy=0, sum_xy=0;
  double s_2=0;
  
  for (i=0; i<n; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_xx += (x[i]*x[i]);
    sum_yy += (y[i]*y[i]);
    sum_xy += (x[i]*y[i]);
  }
  
  delta_x = n*sum_xx - sum_x*sum_x;
  delta_y = n*sum_yy - sum_y*sum_y;
  
  *a = (sum_xx*sum_y - sum_x*sum_xy) / delta_x;
  *b = (n*sum_xy - sum_x*sum_y) / delta_x;
  
  *correlation = n* sum_xy - sum_x*sum_y;
  *correlation /= sqrt (delta_x);
  *correlation /= sqrt (delta_y);
  
  /* calculate errors of coefficients */
  for (i=0; i<n; i++)
    s_2 += ( (y[i]-*a-*b*x[i]) * (y[i]-*a-*b*x[i]) );
  
  s_2 /= (n-2);
  
  *sigma_a = sqrt (sum_xx * s_2 / delta_x);
  *sigma_b = sqrt (n * s_2 / delta_x);
  
  
  return 0;
}



/***********************************************************************************/
/* Function: weight_regression                                            @31_30i@ */
/* Description:                                                                    */
/*  Calculate coefficients for y = a + b*x by @emph{weighted} linear regression.   */
/*  Each data point is weighted with (1 / sigma[i]**2)                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:            Vector (0..n-1) of x-values.                             */
/*  double *y:            Vector (0..n-1) of y-values.                             */
/*  double *sigma:        Vector (0..n-1) of weighting coefficients.               */
/*  int n:                Number of data points.                                   */
/*  double *a:            Coefficient a.                                           */ 
/*  double *b:            Coefficient b.                                           */ 
/*  double *sigma_a:      Standard deviation of a.                                 */
/*  double *sigma_b:      Standard deviation of b.                                 */
/*  double *correlation:  Correlation coefficient.                                 */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

double weight_regression (double *x, double *y, double *sigma, int n,
			  double *a, double *b,
			  double *sigma_a, double *sigma_b)
{
  int i=0;
  double delta_x=0, delta_y=0;
  double correlation=0;
  double sum_x=0, sum_y=0, sum_xx=0, sum_yy=0, sum_xy=0, sum_1=0;
  
  
  /* every term is weighted with factor 1/sigma[i]^2 */
  for (i=0; i<n; i++) {
    sum_x += x[i] / (sigma[i] * sigma[i]);
    sum_y += y[i] / (sigma[i] * sigma[i]);
    sum_xx += (x[i]*x[i]) / (sigma[i] * sigma[i]);
    sum_yy += (y[i]*y[i]) / (sigma[i] * sigma[i]);
    sum_xy += (x[i]*y[i]) / (sigma[i] * sigma[i]);
    sum_1 += 1 / (sigma[i] * sigma[i]);
  }
  
  
  delta_x = sum_1*sum_xx - sum_x*sum_x;
  delta_y = sum_1*sum_yy - sum_y*sum_y;
  
  *a = (sum_xx*sum_y - sum_x*sum_xy) / delta_x;
  *b = (sum_1*sum_xy - sum_x*sum_y) / delta_x;
  
  correlation = sum_1* sum_xy - sum_x*sum_y;
  correlation /= sqrt (delta_x);
  correlation /= sqrt (delta_y);
  
  /* calculate errors of coefficients */
  /* see Bevington, pg. 116 !!!! */
  
  *sigma_a = sqrt (sum_xx / delta_x);
  *sigma_b = sqrt (sum_1 / delta_x);
  
  
  return (correlation);
}


/***************************************/
/* fit a gaussian curve to data points */
/***************************************/

int gaussfit (double *x, double *y, long number,
	      double *mu, double *sigma, double *area)
{
  double *A[3]    = {NULL,NULL,NULL};
  double b[3]     = {0,0,0};
  double *res = NULL;
  double *gamma   = (double *) calloc (number, sizeof(double));
  double *lny     = (double *) calloc (number, sizeof(double));
  long i=0;
  
  
  for (i=0; i<3; i++)
    A[i] = (double *) calloc (3, sizeof(double)); 

  for (i=0; i<number; i++)  
    lny[i] = log (y[i]);
  
  for (i=0; i<number; i++)
    gamma[i] = y[i] * y[i];
  
  for (i=0; i<number; i++) {
    A[0][0] += (gamma[i] * x[i] * x[i] * x[i] * x[i]);
    A[0][1] += (gamma[i] * x[i] * x[i] * x[i]);
    A[0][2] += (gamma[i] * x[i] * x[i]);
    A[1][2] += (gamma[i] * x[i]);
    A[2][2] += gamma[i];
    b[0] += (gamma[i] * lny[i] * x[i] * x[i]);
    b[1] += (gamma[i] * lny[i] * x[i]);
    b[2] += (gamma[i] * lny[i]);
  }
  
  A[1][0] = A[0][1];    /* use symmetries */
  A[1][1] = A[0][2];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];
  

  free(gamma);
  free(lny);


  if (solve_gauss (A, b, 3, &res) != 0)  {
    for (i=0; i<3; i++)
      free(A[i]); 
    return FIT_NOT_POSSIBLE;
  }

  for (i=0; i<3; i++)
    free(A[i]); 


  if (res[0] >= 0)
    return FIT_NOT_POSSIBLE;

  
  
  *mu = -res[1] / 2.0 / res[0];
  *sigma = sqrt ( -1.0 / 2.0 / res[0] );
  *area = *sigma * 2.506628275 * exp (res[2] - res[1]*res[1] / 4.0 / res[0]);
  
  free(res);

  return 0;
}


/**************************************************************/
/* fit a curve exp (a1x^3 + a2x^2 + a3x +a4) to data points;  */
/* that's a kind of distorted gaussian curve                  */
/**************************************************************/

int gaussfit_distorted (double *x, double *y, long number,
			double *a0, double *a1, double *a2, double *a3)
{
  double *A[4]   = {NULL,NULL,NULL,NULL};
  double b[4]    = {0,0,0,0};
  double *res = NULL;

  double *gamma = (double *) calloc (number, sizeof(double));
  double *lny   = (double *) calloc (number, sizeof(double));

  long i=0;
  
  for (i=0; i<4; i++)
    A[i] = (double *) calloc (4, sizeof(double)); 


  
  for (i=0; i<number; i++)  
    lny[i] = log (y[i]);
  
  for (i=0; i<number; i++)
    gamma[i] = y[i] * y[i];
  
  for (i=0; i<number; i++) {
    A[0][0]  += (gamma[i] * x[i] * x[i] * x[i] * x[i] * x[i] * x[i]);
    A[0][1]  += (gamma[i] * x[i] * x[i] * x[i] * x[i] * x[i]);
    A[0][2]  += (gamma[i] * x[i] * x[i] * x[i] * x[i]);
    A[0][3]  += (gamma[i] * x[i] * x[i] * x[i]);
    A[1][3]  += (gamma[i] * x[i] * x[i]);
    A[2][3]  += (gamma[i] * x[i]);
    A[3][3]  += (gamma[i]);
    b[0]  += (gamma[i] * lny[i] * x[i] * x[i] * x[i]);
    b[1]  += (gamma[i] * lny[i] * x[i] * x[i]);
    b[2]  += (gamma[i] * lny[i] * x[i]);
    b[3]  += (gamma[i] * lny[i]);
  }
  
  A[1][0]  = A[0][1];    /* use symmetries */
  A[1][1]  = A[0][2];
  A[1][2]  = A[0][3];
  A[2][0]  = A[0][2];
  A[2][1]  = A[1][2];
  A[2][2]  = A[1][3];
  A[3][0]  = A[0][3];
  A[3][1]  = A[1][3];
  A[3][2]  = A[2][3];
  
  free(gamma);
  free(lny);
  

  if (solve_gauss (A, b, 4, &res) != 0)  {
    for (i=0; i<4; i++)
      free (A[i]);
    return FIT_NOT_POSSIBLE;
  }

  for (i=0; i<4; i++)
    free (A[i]);

  *a3 = res[0];
  *a2 = res[1];
  *a1 = res[2];
  *a0 = res[3];
  
  free(res);

  return 0;
}


/*************************************************/
/* fit a curve  y = a * exp (b/x) to data points */
/*              ?????????????????                */      
/*************************************************/

int exponentialfit (double *x, double *y, 
		    long number, double *a, double *b)
{
  double *gamma = (double *) calloc (number, sizeof(double));

  double a1=0, a2=0, a3=0, a4=0, a5=0, a6=0;
  double denominator=0;
  long i=0;
  
  double c=0;
  
  
  for (i=0; i<number; i++) 
    gamma[i] = y[i] * y[i];
  
  
  for (i=0; i<number; i++)  {
    if (y[i] != 0)  {
      c = gamma[i] * log (y[i]);
      
      a1 += (gamma[i] * x[i] * x[i]);
      a2 += (gamma[i] * x[i]);
      a3 += (c * x[i]);
      a4  = a2;
      a5 += gamma[i];
      a6 += c;
    }
  }
  
  free(gamma);


  if ( (denominator = a1 * a5 - a2 * a4) == 0)
    return FIT_NOT_POSSIBLE;
  
  *a = (a3*a5 - a2*a6) / denominator;
  *b = (a1*a6 - a3*a4) / denominator;
  
  return 0;
}



/************************************************/
/* fit a curve  y = exp (ax + b) to data points */
/************************************************/

int boltzmannfit (double *x, double *y, long number, double *a, double *b)
{
  double *gamma = (double *) calloc (number, sizeof(double));

  int i=0;
  double a1=0, a2=0, a3=0, a4=0, a5=0, a6=0;
  double denominator=0;
  double alpha=0, beta=0;
  
  
  for (i=0; i<number; i++)
    gamma[i] = y[i] * y[i];
  
  for (i=0; i<number; i++)  {
    if (x[i] != 0 && y[i] != 0)  {
      a1 += ( gamma[i] * log ( y[i] ));
      a2 += ( gamma[i] / x[i] );
      a3 += ( gamma[i] * log (y[i]) / x[i]);
      a4 += ( gamma[i] * y[i] / x[i]);
      a5 += ( gamma[i] / x[i] / x[i]);
      a6 += gamma[i];
    }
  }
  
  free(gamma);

  denominator = a6 * a5 - a2 * a2;
  
  if (denominator == 0)
    return FIT_NOT_POSSIBLE;

  alpha = (a1 * a5 - a2 * a3) / denominator;
  beta  = (a6 * a3 - a1 * a2) / denominator;
  
  *a = exp (alpha);
  *b = beta;
  
  return 0;
}



/******************************************************/
/* fit a parabola  a0 + a1*x + a2*x^2  to data points */ 
/******************************************************/

int parabolafit (double *x, double *y, long number, 
		 double *a0, double *a1, double *a2)
{
  double *A[3]   = {NULL,NULL,NULL};
  double b[3]    = {0,0,0};
  double *result = NULL;
  long i=0;
  
  for (i=0; i<3; i++)
    A[i] = (double *) calloc (3, sizeof(double));


  A[0][0] = number;
  
  for (i=0; i<number; i++)  {
    A[0][1] += x[i];
    A[0][2] += (x[i] * x[i]);
    A[1][2] += (x[i] * x[i] * x[i]);
    A[2][2] += (x[i] * x[i] * x[i] * x[i]);
    b[0] += y[i];
    b[1] += (y[i]*x[i]);
    b[2] += (y[i]*x[i]*x[i]);
  }
  
  A[1][0] = A[0][1];
  A[1][1] = A[0][2];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];
  
  if (solve_gauss (A, b, 3, &result) != 0)  {
    for (i=0; i<3; i++)
      free(A[i]);
    return FIT_NOT_POSSIBLE;
  }

  for (i=0; i<3; i++)
    free(A[i]);

  *a0 = result[0];
  *a1 = result[1];
  *a2 = result[2];
  
  free(result);

  return 0;
}



/***********************************************************/
/* fit a cubic  a0 + a1*x + a2*x^2 + a3*x^3 to data points */ 
/***********************************************************/

int cubicfit (double *x, double *y, long number, 
	      double *a0, double *a1, double *a2, double *a3)
{
  double *A[4]   = {NULL,NULL,NULL,NULL};
  double b[4]    = {0,0,0,0};
  double *result = NULL;
  long i=0;
  
  for (i=0; i<4; i++)
    A[i] = (double *) calloc (4, sizeof(double));
    


  A[0][0] = number;
  
  for (i=0; i<number; i++)  {
    A[0][1] += x[i];
    A[0][2] += (x[i] * x[i]);
    A[0][3] += (x[i] * x[i] * x[i]);
    A[1][3] += (x[i] * x[i] * x[i] * x[i]);
    A[2][3] += (x[i] * x[i] * x[i] * x[i] * x[i]);
    A[3][3] += (x[i] * x[i] * x[i] * x[i] * x[i] * x[i]);
    
    b[0] += y[i];
    b[1] += (y[i]*x[i]);
    b[2] += (y[i]*x[i]*x[i]);
    b[3] += (y[i]*x[i]*x[i]*x[i]);
  }
  
  A[1][0] = A[0][1];
  A[1][1] = A[0][2];
  A[1][2] = A[0][3];
  A[2][0] = A[0][2];
  A[2][1] = A[0][3];
  A[2][2] = A[1][3];
  A[3][0] = A[0][3];
  A[3][1] = A[1][3];
  A[3][2] = A[2][3];
  
  
  
  
  if (solve_gauss (A, b, 4, &result) != 0)  {
    for (i=0; i<4; i++)
      free(A[i]);
    return FIT_NOT_POSSIBLE;
  }
  
  for (i=0; i<4; i++)
    free(A[i]);

  *a0 = result[0];
  *a1 = result[1];
  *a2 = result[2];
  *a3 = result[3];
  
  free(result);

  return 0;
}


/*********************************************/
/* fit a hyperbola  1/(a+b*x) to data points */ 
/*********************************************/

int hyperbolafit (double *x, double *y, long number, 
		  double *a, double *b)
{
  long i=0;
  double *x_trans = (double *) calloc(number,sizeof(double)); 
  double *y_trans = (double *) calloc(number,sizeof(double)); 
  double *sigma   = (double *) calloc(number,sizeof(double)); 
  long counter=0;
  double sigma_a=0;
  double sigma_b=0;
  
  /* take only not-zero values of y */
  for (i=0; i<number; i++)  {
    if (y[i] != 0)  {
      x_trans[counter] = x[i];
      y_trans[counter] = 1/y[i];
      sigma[counter] = 1/y[i]/y[i];
      counter++;
    }
  }
  
  weight_regression (x_trans, y_trans, sigma, counter,
		     a, b, &sigma_a, &sigma_b);
  
  free(x_trans);
  free(y_trans);
  free(sigma);




  return 0;
}


/**************************************************************/
/* fit an inverse parabola  1/(a0+a1*x+a2*x^2) to data points */ 
/**************************************************************/

int inv_parabolafit (double *x, double *y, long number, 
		     double *a0, double *a1, double *a2)
{
  long i=0;
  double *x_trans = (double *) calloc (number, sizeof(double));  
  double *y_trans = (double *) calloc (number, sizeof(double));  
  double *weight  = (double *) calloc (number, sizeof(double));  
  long counter=0;
  double *A[3]     = {NULL,NULL,NULL};
  double b[3]      = {0,0,0};
  double *result   = NULL;
  
  for (i=0; i<3; i++)
    A[i] = (double *) calloc (3, sizeof(double));


  /* take only not-zero values of y */
  for (i=0; i<number; i++)  {
    if (y[i] != 0)  {
      x_trans[counter] = x[i];
      y_trans[counter] = 1/y[i];
      weight[counter] = y[i]*y[i]*y[i]*y[i];
      counter++;
    }
  }
  
  
  for (i=0; i<counter; i++)  {
    A[0][0] += (weight[i]);
    A[0][1] += (weight[i] * x_trans[i]);
    A[0][2] += (weight[i] * x_trans[i]*x_trans[i]);
    A[1][2] += (weight[i] * x_trans[i]*x_trans[i]*x_trans[i]);
    A[2][2] += (weight[i] * x_trans[i]*x_trans[i]*x_trans[i]*x_trans[i]);
    b[0] += (weight[i] * y_trans[i]);
    b[1] += (weight[i] * y_trans[i]*x_trans[i]);
    b[2] += (weight[i] * y_trans[i]*x_trans[i]*x_trans[i]);
  }
  
  free(x_trans);
  free(y_trans);
  free(weight);

  A[1][0] = A[0][1];
  A[1][1] = A[0][2];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];
  
  if (solve_gauss (A, b, 3, &result) != 0)  {
    for (i=0; i<3; i++)
      free (A[i]);
    return FIT_NOT_POSSIBLE;
  }

  for (i=0; i<3; i++)
    free (A[i]);

  *a0 = result[0];
  *a1 = result[1];
  *a2 = result[2];
  
  free(result);

  return 0;
  
}
