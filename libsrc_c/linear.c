/*--------------------------------------------------------------------
 * $Id: linear.c 2623 2011-12-23 10:52:38Z robert.buras $
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
/* Function: linear_coeffc                                                @31_30i@ */
/* Description:                                                                    */
/*  Calculate coefficients for linear interpolation;                               */
/*  memory for coefficients will be allocated automatically!                       */
/*  These function has been created for compatibility with the spline              */
/*  interpolation functions; for this reason four coefficients are calculated,     */
/*  but a2[] and a3[] are set to zero. The interpolation may be done with          */
/*  calc_spline_values().                                                          */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:    x values of the data points, i=0...number-1                      */
/*  double *y:    y values of the data points, i=0...number-1                      */
/*  int number:   number of datapoints                                             */
/*  double **a0:  array of coefficients, i=0...number-1                            */ 
/*  double **a1:  array of coefficients, i=0...number-1                            */ 
/*  double **a2:  array of coefficients, i=0...number-1, set to zero               */ 
/*  double **a3:  array of coefficients, i=0...number-1, set to zero               */ 
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

int linear_coeffc (double *x, double *y, int number, 
		   double **a0, double **a1, double **a2, double **a3)
{
  int i=0;


  /* check if x ascending */
  for (i=0; i<number-1; i++)  {
    if (x[i]>=x[i+1]) {
      fprintf(stderr,"x not ascending %d %f %f\n", i, x[i], x[i+1]);
      return X_NOT_ASCENDING;
    }
  }

  
  /* allocate memory for coefficients */
  *a0 = (double *) calloc (number, sizeof(double));
  *a1 = (double *) calloc (number, sizeof(double));
  *a2 = (double *) calloc (number, sizeof(double));
  *a3 = (double *) calloc (number, sizeof(double));


  /* calculate coefficients */
  for (i=0; i<number-1; i++)  {
    (*a0)[i] = y[i];
    (*a1)[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]);
  }

  if (number==1) {
    (*a0)[0] = y[0];
    (*a1)[0] = 0;
  }

  return 0;   /* if o.k. */
}


/***********************************************************************************/
/* Function: linear_coeffc_fast                                           @31_30i@ */
/* Description:                                                                    */
/*  Calculate coefficients for linear interpolation;                               */
/*  memory for coefficients will be allocated automatically!                       */
/*  In contrast to linear_coeffc(), only a0 and a1 are used. The interpolation     */
/*  may be done with calc_linear_values().                                         */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:    x values of the data points, i=0...number-1                      */
/*  double *y:    y values of the data points, i=0...number-1                      */
/*  int number:   number of datapoints                                             */
/*  double **a0:  array of coefficients, i=0...number-1                            */ 
/*  double **a1:  array of coefficients, i=0...number-1                            */ 
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

int linear_coeffc_fast (double *x, double *y, int number, 
			double **a0, double **a1)
{
  int i=0;

  /* check if x ascending */
  for (i=0; i<number-1; i++)  {
    if (x[i]>=x[i+1]) {
      fprintf(stderr,"x not ascending %d %f %f\n", i, x[i], x[i+1]);
      return X_NOT_ASCENDING;
    }
  }
  
  /* allocate memory for coefficients */
  *a0 = (double *) calloc (number, sizeof(double));
  *a1 = (double *) calloc (number, sizeof(double));

  /* calculate coefficients */
  for (i=0; i<number-1; i++)  {
    (*a0)[i] = y[i];
    (*a1)[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]);
  }

  if (number==1) {
    (*a0)[0] = y[0];
    (*a1)[0] = 0;
  }

  return 0;   /* if o.k. */
}




/***********************************************************************************/
/* Function: slinear_coeffc                                               @31_30i@ */
/* Description:                                                                    */
/*  Like linear_coeffc(), but sort data before calculating coefficients. Please    */
/*  note that the fields x and y themselves are sorted by slinear_coeffc()  which  */
/*  is an important pre-requisite when they are passed to calc_splined_value()     */
/*  later.                                                                         */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x:    x values of the data points, i=0...number-1                      */
/*  double *y:    y values of the data points, i=0...number-1                      */
/*  int number:   number of datapoints                                             */
/*  double **a0:  array of coefficients, i=0...number-1                            */ 
/*  double **a1:  array of coefficients, i=0...number-1                            */ 
/*  double **a2:  array of coefficients, i=0...number-1, set to zero               */ 
/*  double **a3:  array of coefficients, i=0...number-1, set to zero               */ 
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

int slinear_coeffc (double *x, double *y, int number, 
		    double **a0, double **a1, double **a2, double **a3)
{
  int i=0, status=0;
  double **value=NULL;

  value = calloc (number, sizeof(double *));

  for (i=0; i<number; i++) {
    value[i] = calloc (2, sizeof(double));
    value[i][0] = x[i];
    value[i][1] = y[i];
  }

  status = ASCII_sortarray (value, number, 2, 0, 0);
  if (status!=0) {
    fprintf (stderr, "Error %d sorting data\n", status);
    return status;
  }

  for (i=0; i<number; i++) {
    x[i] = value[i][0];
    y[i] = value[i][1];
    free(value[i]);
  }

  free(value);
  
  /* allocate memory for coefficients */
  *a0 = (double *) calloc (number, sizeof(double));
  *a1 = (double *) calloc (number, sizeof(double));
  *a2 = (double *) calloc (number, sizeof(double));
  *a3 = (double *) calloc (number, sizeof(double));


  /* calculate coefficients */
  for (i=0; i<number-1; i++)  {
    (*a0)[i] = y[i];
    (*a1)[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]);
  }

  if (number==1) {
    (*a0)[0] = y[0];
    (*a1)[0] = 0;
  }

  return 0;   /* if o.k. */
}




/***********************************************************************************/
/* Function: calc_linear_value                                            @31_30i@ */
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
/*                                                                                 */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:  Bernhard Mayer                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int calc_linear_value (double xnew, double *ynew, 
		       double *x, int number, 
		       double *a0, double *a1)
{
  int i=0;

  *ynew=0;

  if (number==0)
    return NO_INTERPOLATION;

  /* no return value, but programm can still continue */
  if (xnew < x[0] || xnew > x[number-1])
    return NO_EXTRAPOLATION;

  i=locate (x, number, xnew);

  *ynew = a0[i] + a1[i] * (xnew-x[i]);

  return 0;
}
