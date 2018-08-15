/*--------------------------------------------------------------------
 * $Id: spl.c 2623 2011-12-23 10:52:38Z robert.buras $
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
/* Function: spline                                                       @31_30i@ */
/* Description:                                                                    */
/*  Interpolate between given data points using natural cubic splines.             */
/*  The input data (x,y) are interpolated to an equidistant grid (start, step).    */
/*  Memory for result vectors will be allocated automatically.                     */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:      x[0..number-1], x-values of the input data.                    */
/*  double *y:      y[0..number-1], y-values of the input data.                    */
/*  int number:     Number of input data points (x[i], y[i]).                      */
/*  double start:   Start  point for the output grid.                              */
/*  double step:    Step of the output grid; data are interpolated to              */
/*                  start, start+step, start+2*step, ...                           */
/*  int *newnumber: Number of output data points.                                  */
/*  double *x_new:  x_new[0..newnumber], x-values of the interpolated data points. */
/*  double *y_new:  y_new[0..newnumber], y-values of the interpolated data points. */
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

int spline (double *x, double *y, int number, double start, double step,
	    int *newnumber, double **new_x, double **new_y)
{
  int i=0, j=0;
  int status=0;
  double xx=0, xnew=0;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double poly1=0, poly2=0, poly3=0;

  *newnumber=0;

  if (number<2)              /* if too few data points */
    return TOO_FEW_DATA_POINTS;


  /* check for first value */
  if (start >= x[0])
    xx = start;
  else  
    xx = (ceil (x[0]/step)) * step;

  *newnumber = (int) ((x[number-1] - xx) / step + 1.0 + 1E-8);
  
  if (*newnumber<=0)
    return NO_SPLINED_VALUES;

  
  /* calculate coefficients */
  if ((status = spline_coeffc (x, y, number, &a0, &a1, &a2, &a3)) != 0)
    return status;


  /* allocate memory for new arrays */
  *new_x = (double *) calloc (*newnumber, sizeof(double));
  *new_y = (double *) calloc (*newnumber, sizeof(double));

  
  /* calculate new values */
  i=0;
  for (j=0; j<*newnumber; j++)  {
    xnew = xx + (double) j * step;
    if (xnew > x[i+1])
      while (xnew > x[i+1] && i < number-2)
	i++;

    (*new_x)[j] = xnew; 
    poly1 = xnew-x[i];
    poly2 = poly1*poly1;
    poly3 = poly2*poly1;
    (*new_y)[j] = a3[i] * poly3 + a2[i] * poly2 + a1[i] * poly1 + a0[i];
  }


  /* free memory of coefficients */
  free(a0);
  free(a1);
  free(a2);
  free(a3);
  
  return 0;
}





/***********************************************************************************/
/* Function: spline_coeffc                                                @31_30i@ */
/* Description:                                                                    */
/*  Calculate coefficients for natural cubic spline. These coefficients may        */
/*  be used as input to calc_splined_value() in order to interpolate the data      */
/*  points at arbitrary x values. Memory for the vectors a0, a1, a2, and a3        */
/*  will be allocated automatically.                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:      x[0..number-1], x-values of the input data.                    */
/*  double *y:      y[0..number-1], y-values of the input data.                    */
/*  int number:     Number of input data points (x[i], y[i]).                      */
/*  double **a0:    Pointer to vector of 0th order spline coefficients.            */
/*  double **a1:    Pointer to vector of 1st order spline coefficients.            */
/*  double **a2:    Pointer to vector of 2nd order spline coefficients.            */
/*  double **a3:    Pointer to vector of 3rd order spline coefficients.            */
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

int spline_coeffc (double *x, double *y, int number, 
		   double **a0, double **a1, double **a2, double **a3)
{
  int i=0, status=0;
  double *h=NULL;
  double *n=NULL, *res=NULL;
  double **m=NULL;

  if (number<2) {
    fprintf(stderr,"Error, too few data points during spline_coeffc (in src_c/spl.c) \n");
    return TOO_FEW_DATA_POINTS;
  }

  if (number==2)  {     /* linear interpolation */

    if (x[1] <= x[0]) {
      fprintf(stderr,"x not ascending %d %f %f\n", i, x[i], x[i+1]);
      return X_NOT_ASCENDING;
    }   

    *a0 = (double *) calloc (number, sizeof(double));
    *a1 = (double *) calloc (number, sizeof(double));
    *a2 = (double *) calloc (number, sizeof(double));
    *a3 = (double *) calloc (number, sizeof(double));
    
    (*a0)[0] = y[0];
    (*a1)[0] = (y[1]-y[0])/(x[1]-x[0]);

    return 0;
  }


  h = (double *) calloc (number-1, sizeof(double));
  
  for (i=0; i<number-1; i++)  {
    h[i] = x[i+1] - x[i];
    if ( h[i] <= 0 )  {
      free (h);
      fprintf(stderr,"x not ascending %d %f %f\n", i, x[i], x[i+1]);
      return X_NOT_ASCENDING;
    }
  }
  
  /* create linear equation system  mx = n */
  n   = (double *) calloc (number-2, sizeof(double));
  
  m   = (double **) calloc (number-2, sizeof(double *)); 
  for (i=0; i<number-2; i++)  
    m[i] = (double *) calloc (3, sizeof(double));

  for (i=1; i<number-2; i++)
    m[i][0] = h[i];
  
  for (i=0; i<number-2; i++)
    m[i][1] = 2.0*(h[i]+h[i+1]);
  
  for (i=0; i<number-3; i++)
    m[i][2] = h[i+1];
  
  for (i=0; i<number-2; i++)
    n[i] = 3.0 / h[i+1] * (y[i+2]-y[i+1]) - 3.0 / h[i] * (y[i+1]-y[i]);
  
  status = solve_three_ms (m, n, number-2, &res);
  
  /* free memory */
  for (i=0; i<number-2; i++)  
    free (m[i]);
  free(m);
  free(n);
  
  /* return if error solving equation system */
  if (status != 0)  {
    free(h);
    fprintf(stderr,"Error solving equation system during spline_coeffc (in src_c/spl.c) \n");
    return SPLINE_NOT_POSSIBLE ;
  }
  
  *a0 = (double *) calloc (number, sizeof(double)); /* spline coefficients */
  *a1 = (double *) calloc (number, sizeof(double));
  *a2 = (double *) calloc (number, sizeof(double));
  *a3 = (double *) calloc (number, sizeof(double));


  (*a2)[0] = 0;

  for (i=1; i<number-1; i++)
    (*a2)[i] = res[i-1];

  free(res);
  
  for (i=0; i<number-1; i++)  {
    (*a3)[i] = ((*a2)[i+1] - (*a2)[i]) / 3.0 / h[i];
    (*a1)[i] = (y[i+1] - y[i]) / h[i]  - h[i] / 3.0 * ((*a2)[i+1] + 2.0*(*a2)[i]);
    (*a0)[i] = y[i];
  }
  
  /* free memory */
  free(h);
  
  return 0;
}



/***********************************************************************************/
/* Function: fspline_coeffc                                               @31_30i@ */
/* Description:                                                                    */
/*  Similar to spline_coeffc, only all data are of type float instead of           */
/*  double. Calculate coefficients for natural cubic spline. These coefficients    */
/*  may be used as input to calc_splined_value() in order to interpolate the data  */
/*  points at arbitrary x values. Memory for the vectors a0, a1, a2, and a3        */
/*  will be allocated automatically.                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*  float *x:      x[0..number-1], x-values of the input data.                     */
/*  float *y:      y[0..number-1], y-values of the input data.                     */
/*  int number:    Number of input data points (x[i], y[i]).                       */
/*  float **a0:    Pointer to vector of 0th order spline coefficients.             */
/*  float **a1:    Pointer to vector of 1st order spline coefficients.             */
/*  float **a2:    Pointer to vector of 2nd order spline coefficients.             */
/*  float **a3:    Pointer to vector of 3rd order spline coefficients.             */
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

int fspline_coeffc (float *x, float *y, int number, 
		    float **a0, float **a1, float **a2, float **a3)
{
  int i=0, status=0;
  float *h=NULL;
  float *n=NULL, *res=NULL;
  float **m=NULL;

  if (number<2)
    return TOO_FEW_DATA_POINTS;


  if (number==2)  {     /* linear interpolation */

    if (x[1] <= x[0]) {
      fprintf(stderr,"x not ascending %d %f %f\n", i, x[i], x[i+1]);
      return X_NOT_ASCENDING;
    }   

    *a0 = (float *) calloc (number, sizeof(float));
    *a1 = (float *) calloc (number, sizeof(float));
    *a2 = (float *) calloc (number, sizeof(float));
    *a3 = (float *) calloc (number, sizeof(float));
    
    (*a0)[0] = y[0];
    (*a1)[0] = (y[1]-y[0])/(x[1]-x[0]);

    return 0;
  }


  h = (float *) calloc (number-1, sizeof(float));
  
  for (i=0; i<number-1; i++)  {
    h[i] = x[i+1] - x[i];
    if ( h[i] <= 0 )  {
      free (h);
      fprintf(stderr,"x not ascending %d %f %f\n", i, x[i], x[i+1]);
      return X_NOT_ASCENDING;
    }
  }
  
  /* create linear equation system  mx = n */
  n   = (float *) calloc (number-2, sizeof(float));
  
  m   = (float **) calloc (number-2, sizeof(float *)); 
  for (i=0; i<number-2; i++)  
    m[i] = (float *) calloc (3, sizeof(float));

  for (i=1; i<number-2; i++)
    m[i][0] = h[i];
  
  for (i=0; i<number-2; i++)
    m[i][1] = 2.0*(h[i]+h[i+1]);
  
  for (i=0; i<number-3; i++)
    m[i][2] = h[i+1];
  
  for (i=0; i<number-2; i++)
    n[i] = 3.0 / h[i+1] * (y[i+2]-y[i+1]) - 3.0 / h[i] * (y[i+1]-y[i]);
  
  status = fsolve_three_ms (m, n, number-2, &res);
  
  /* free memory */
  for (i=0; i<number-2; i++)  
    free (m[i]);
  free(m);
  free(n);
  
  /* return if error solving equation system */
  if (status != 0)  {
    free(h);
    return SPLINE_NOT_POSSIBLE ;
  }
  
  *a0 = (float *) calloc (number, sizeof(float)); /* spline coefficients */
  *a1 = (float *) calloc (number, sizeof(float));
  *a2 = (float *) calloc (number, sizeof(float));
  *a3 = (float *) calloc (number, sizeof(float));


  for (i=1; i<number-1; i++)
    (*a2)[i] = res[i-1];

  free(res);
  

  (*a1)[0] = (y[1] - y[0]) / h[0];
  (*a0)[0] = y[0];
  
  for (i=1; i<number-1; i++)  {
    (*a3)[i] = ((*a2)[i+1] - (*a2)[i]) / 3.0 / h[i];
    (*a1)[i] = (y[i+1] - y[i]) / h[i]  - h[i] / 3.0 * ((*a2)[i+1] + 2.0*(*a2)[i]);
    (*a0)[i] = y[i];
  }
  
  /* free memory */
  free(h);
  
  return 0;
}




/***********************************************************************************/
/* Function: appspl                                                       @31_30i@ */
/* Description:                                                                    */
/*  @strong{Approximate} data points using natural cubic splines.                  */
/*  The input data (x,y) are interpolated to an equidistant grid (start, step).    */
/*  Memory for result vectors will be allocated automatically.                     */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:      x[0..number-1], x-values of the input data.                    */
/*  double *y:      y[0..number-1], y-values of the input data.                    */
/*  double *w:      y[0..number-1], weighting coefficients. Increasing w will      */
/*                  increase the degree of approximation; extreme cases are        */
/*                  w=0 (interpolating spline) and                                 */ 
/*                  w=infinity (linear regression).                                */
/*  int number:     Number of input data points (x[i], y[i]).                      */
/*  double start:   Start  point for the output grid.                              */
/*  double step:    Step of the output grid; data are interpolated to              */
/*                  start, start+step, start+2*step, ...                           */
/*  int *newnumber: Number of output data points.                                  */
/*  double *x_new:  x_new[0..newnumber], x-values of the interpolated data points. */
/*  double *y_new:  y_new[0..newnumber], y-values of the interpolated data points. */
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

int appspl (double *x, double *y, double *w, int number, 
	    double start, double step, 
	    int *newnumber, double **new_x, double **new_y)
     
{
  int i=0, j=0;
  int status=0;
  double xx=0, xnew=0;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double poly1=0, poly2=0, poly3=0;

  *newnumber=0;



  /* check for first value */
  if (start >= x[0])
    xx = start;
  else  
    xx = (ceil (x[0]/step)) * step;

  *newnumber = (int) ((x[number-1] - xx) / step + 1.0 + 1E-8);
  
  if (*newnumber<=0)
    return NO_SPLINED_VALUES;

  
  /* calculate coefficients */
  if ((status = appspl_coeffc (x, y, w, number, &a0, &a1, &a2, &a3)) != 0)
    return status;


  /* allocate memory for new arrays */
  *new_x = (double *) calloc (*newnumber, sizeof(double));
  *new_y = (double *) calloc (*newnumber, sizeof(double));

  
  /* calculate new values */
  i=0;
  for (j=0; j<*newnumber; j++)  {
    xnew = xx + (double) j * step;
    if (xnew > x[i+1])
      while (xnew > x[i+1] && i < number-2)
	i++;

    (*new_x)[j] = xnew; 
    poly1 = xnew-x[i];
    poly2 = poly1*poly1;
    poly3 = poly2*poly1;
    (*new_y)[j] = a3[i] * poly3 + a2[i] * poly2 + a1[i] * poly1 + a0[i];
  }


  /* free memory of coefficients */
  free(a0);
  free(a1);
  free(a2);
  free(a3);
  
  return 0;


}





/***********************************************************************************/
/* Function: appspl_coeffc                                                @31_30i@ */
/* Description:                                                                    */
/*  Calculate coefficients for approximating cubic spline. These coefficients may  */
/*  be used as input to calc_splined_value() in order to interpolate the data      */
/*  points at arbitrary x values. Memory for the vectors a0, a1, a2, and a3        */
/*  will be allocated automatically.                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:      x[0..number-1], x-values of the input data.                    */
/*  double *y:      y[0..number-1], y-values of the input data.                    */
/*  double *w:      w[0..number-1], weighting coefficients. Increasing w will      */
/*                  increase the degree of approximation; extreme cases are        */
/*                  w=0 (interpolating spline) and                                 */ 
/*                  w=infinity (linear regression).                                */
/*  int number:     Number of input data points (x[i], y[i]).                      */
/*  double **a0:    Pointer to vector of 0th order spline coefficients.            */
/*  double **a1:    Pointer to vector of 1st order spline coefficients.            */
/*  double **a2:    Pointer to vector of 2nd order spline coefficients.            */
/*  double **a3:    Pointer to vector of 3rd order spline coefficients.            */
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

int appspl_coeffc (double *x, double *y, double *w, int number, 
		   double **a0, double **a1, double **a2, double **a3)
{
  int n = number-1;
  int i=0;
  double help1=0, help2=0;
  double *a=NULL, *b=NULL, *c=NULL, *d=NULL;
  double *h1=NULL, *h2=NULL, *H=NULL;    
  double *res=NULL;   
  double **mat=NULL;  



  if (n<5)              /* if too few data points */
    return TOO_FEW_DATA_POINTS;
  
  for (i=0; i<n; i++)   /* data sorted by x-value ? */
    if (x[i] >= x[i+1])
      return DATA_NOT_SORTED;
  
  for (i=0; i<=n; i++)  /* weighting factors positive ? */ 
    if (w[i] < 0) 
      return NEGATIVE_WEIGHTING_FACTORS;
  



  a      = (double *) calloc (n+1, sizeof(double)); 
  b      = (double *) calloc (n+1, sizeof(double)); 
  c      = (double *) calloc (n+1, sizeof(double)); 
  d      = (double *) calloc (n+1, sizeof(double));
  h1     = (double *) calloc (n+1, sizeof(double));
  h2     = (double *) calloc (n+1, sizeof(double));
  H      = (double *) calloc (n+1, sizeof(double));

  /* matrix mat[n+1][n+1] */ 
  mat   = (double **) calloc(n+1, sizeof(double *));
  for (i=0; i<=n; i++)
    mat[i] = (double *) calloc(5, sizeof(double));


  for (i=1; i<n; i++)
    d[i] = h2[i-1] + h2[i];
  
  /* calculate coefficients */
  for (i=0; i<n; i++)  {
    h1[i] = x[i+1] - x[i];
    h2[i] = 1.0 / h1[i];
    b[i]  = 6.0 / w[i];
  }
  b[n] = 6.0 / w[n];
  
  for (i=0; i<n-1; i++)
    H[i] = h2[i] + h2[i+1];
  
  for (i=1; i<n; i++)
    d[i] = h2[i-1] + h2[i];
  
  
  help1 = (y[1]-y[0]) * h2[0];
  
  for (i=0; i<n-1; i++)  {
    help2 = (y[i+2] - y[i+1]) * h2[i+1];
    c[i] = (help2-help1) * 3.0;
    help1 = help2;
  }
  
  for (i=0; i<n-1; i++)
    mat[i][2] = 2 * (h1[i]+h1[i+1]) + b[i]*h2[i]*h2[i] + 
      b[i+1]*H[i]*H[i] + b[i+2]*h2[i+1]*h2[i+1]; 
  
  for (i=0; i<n-2; i++)
    mat[i][3] = h1[i+1] - b[i+1]*h2[i+1]*H[i] - b[i+2]*h2[i+1]*H[i+1];
  
  for (i=0; i<n-3; i++)
    mat[i][4] = b[i+2]*h2[i+1]*h2[i+2];
  
  for (i=1; i<n-1; i++) 
    mat[i][1] = h1[i] - b[i]*h2[i]*H[i-1] - b[i+1]*h2[i]*H[i];
  
  for (i=2; i<n-1; i++)
    mat[i][0] = b[i]*h2[i-1]*h2[i];
  
  
  /* free memory */
  free(H);      


  if (solve_five_ms (mat, c, n-1, &res) != 0)  {
    free(a);      
    free(b);     
    free(c);      
    free(d);      
    free(h1);     
    free(h2);     
    
    for (i=0; i<5; i++)
      free (mat[i]);
    free(mat);

    return SPLINE_NOT_POSSIBLE ;
  }

  for (i=0; i<=n; i++)
    free (mat[i]);
  free(mat);


  for (i=0; i<n; i++) 
    c[i] = res[i-1];
  
  c[0] = 0.0;
  c[n] = 0.0;
  
  a[0] = y[0] + b[0] / 3.0 * h2[0] * (c[0] - c[1]);
  for (i=1; i<n; i++)  
    a[i] = y[i] - b[i] / 3.0 * (c[i-1] * h2[i-1] - d[i] * c[i] + 
				c[i+1] * h2[i]);
  a[n] = y[n] - b[n] / 3.0 * h2[n-1] * (c[n-1] - c[n]);
  
  
  for (i=0; i<n; i++)  {
    b[i] = h2[i] * (a[i+1] - a[i]) - h1[i] / 3.0 * (c[i+1] +2.0 * c[i]);
    d[i] = h2[i] / 3.0 * (c[i+1] - c[i]);
  }
  
  /* free memory */
  free(h1);     
  free(h2);     
  free(res);    

  
  *a0=a;
  *a1=b;
  *a2=c;
  *a3=d;

  
  return 0;
}




/***********************************************************************************/
/* Function: calc_splined_value                                           @31_30i@ */
/* Description:                                                                    */
/*  Interpolate/approximate data to x_new. Basis are the spline coefficients       */
/*  which have been determined either by spline_coeffc or by appspl_coeffc,        */
/*  or even by linear_coeffc.                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double xnew:      x-value to be interpolated.                                  */
/*  double *ynew:     Pointer to the calculated y-value.                           */
/*  double *x:        x[0..number-1], original x-values.                           */
/*  int number:       Number of original x-values.                                 */
/*  double **a0:      Pointer to vector of 0th order spline coefficients.          */
/*  double **a1:      Pointer to vector of 1st order spline coefficients.          */
/*  double **a2:      Pointer to vector of 2nd order spline coefficients.          */
/*  double **a3:      Pointer to vector of 3rd order spline coefficients.          */
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

int calc_splined_value (double xnew, double *ynew, 
			double *x, int number, 
			double *a0, double *a1, double *a2, double *a3)
{
  int i=0;
  double poly1=0, poly2=0, poly3=0;

  *ynew=0;

  if (number==0)
    return NO_INTERPOLATION;

  /* no return value, but programm can still continue */
  if (xnew < x[0] || xnew > x[number-1])
    return NO_EXTRAPOLATION;

  i=locate (x, number, xnew);

  poly1 = xnew-x[i];
  poly2 = poly1*poly1;
  poly3 = poly2*poly1;
  *ynew = a3[i] * poly3 + a2[i] * poly2 + a1[i] * poly1 + a0[i];

  return 0;
}
  

/***********************************************************************************/
/* Function: fcalc_splined_value                                          @31_30i@ */
/* Description:                                                                    */
/*  Similar to calc_splined_value, only all data are of type float instead of      */
/*  double. Interpolate/approximate data to x_new. Basis are the spline            */
/*  coefficients which have been determined by fspline_coeffc.                     */
/*                                                                                 */
/* Parameters:                                                                     */
/*  float xnew:      x-value to be interpolated.                                   */
/*  float *ynew:     Pointer to the calculated y-value.                            */
/*  float *x:        x[0..number-1], original x-values.                            */
/*  int number:      Number of original x-values.                                  */
/*  float **a0:      Pointer to vector of 0th order spline coefficients.           */
/*  float **a1:      Pointer to vector of 1st order spline coefficients.           */
/*  float **a2:      Pointer to vector of 2nd order spline coefficients.           */
/*  float **a3:      Pointer to vector of 3rd order spline coefficients.           */
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

int fcalc_splined_value (float xnew, float *ynew, 
			 float *x, int number, 
			 float *a0, float *a1, float *a2, float *a3)
{
  int i=0;
  float poly1=0, poly2=0, poly3=0;

  *ynew=0;

  if (xnew < x[0] || xnew > x[number-1])
    return NO_EXTRAPOLATION;

  while (xnew>x[i])
    i++;
  
  if (i>0)
    i--;

  poly1 = xnew-x[i];
  poly2 = poly1*poly1;
  poly3 = poly2*poly1;
  *ynew = a3[i] * poly3 + a2[i] * poly2 + a1[i] * poly1 + a0[i];

  return 0;
}
  



/***********************************************************************************/
/* Function: linear_eqd                                                   @31_30i@ */
/* Description:                                                                    */
/*  Interpolate linearely between given data points (x[i], y[i]).                  */
/*  The input data (x,y) are interpolated to an equidistant grid (start, step).    */
/*  Memory for result vectors will be allocated automatically. The syntax of       */ 
/*  linear_eqd is identical to spline.                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:      x[0..number-1], x-values of the input data.                    */
/*  double *y:      y[0..number-1], y-values of the input data.                    */
/*  int number:     Number of input data points (x[i], y[i]).                      */
/*  double start:   Start  point for the output grid.                              */
/*  double step:    Step of the output grid; data are interpolated to              */
/*                  start, start+step, start+2*step, ...                           */
/*  int *newnumber: Number of output data points.                                  */
/*  double *x_new:  x_new[0..newnumber], x-values of the interpolated data points. */
/*  double *y_new:  y_new[0..newnumber], y-values of the interpolated data points. */
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

int linear_eqd (double *x, double *y, int number, double start, double step,
		int *newnumber, double **new_x, double **new_y)
{
  int i=0, j=0;
  int status=0;
  double xx=0, xnew=0;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double poly1=0, poly2=0, poly3=0;

  *newnumber=0;

  if (number<2)              /* if too few data points */
    return TOO_FEW_DATA_POINTS;


  /* check for first value */
  if (start >= x[0])
    xx = start;
  else  
    xx = (ceil (x[0]/step)) * step;

  *newnumber = (int) ((x[number-1] - xx) / step + 1.0 + 1E-8);
  
  if (*newnumber<=0)
    return NO_SPLINED_VALUES;

  
  /* calculate coefficients */
  if ((status = linear_coeffc (x, y, number, &a0, &a1, &a2, &a3)) != 0)
    return status;


  /* allocate memory for new arrays */
  *new_x = (double *) calloc (*newnumber, sizeof(double));
  *new_y = (double *) calloc (*newnumber, sizeof(double));

  
  /* calculate new values */
  i=0;
  for (j=0; j<*newnumber; j++)  {
    xnew = xx + (double) j * step;
    if (xnew > x[i+1])
      while (xnew > x[i+1] && i < number-2)
	i++;

    (*new_x)[j] = xnew; 
    poly1 = xnew-x[i];
    poly2 = poly1*poly1;
    poly3 = poly2*poly1;
    (*new_y)[j] = a3[i] * poly3 + a2[i] * poly2 + a1[i] * poly1 + a0[i];
  }


  /* free memory of coefficients */
  free(a0);
  free(a1);
  free(a2);
  free(a3);
  
  return 0;
}



