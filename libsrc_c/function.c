/*--------------------------------------------------------------------
 * $Id: function.c 2623 2011-12-23 10:52:38Z robert.buras $
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



/* local definitions */
#define DOUBLE_RELATIVE_ERROR 1e-10


/***********************************************************************************/
/* Function: double_equal                                                 @31_30i@ */
/* Description:                                                                    */
/*  Compare two float values; returns 0, if the relative difference is             */
/*  bigger than DOUBLE_RELATIVE_ERROR, which is 1E-10 here. The intention of       */
/*  this function is to avoid roundoff errors when comparing two floats.           */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double a:      First float to be compared.                                     */
/*  double b:      Second float to be compared.                                    */
/*                                                                                 */
/* Return value:                                                                   */
/*  1  if 'equal', 0 if not equal.                                                 */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int double_equal (double a, double b)  
{
  double diff=0, temp=0;

  /* return 1, if a==b */
  if ((diff=fabs(a-b)) == 0.0)
    return 1;


  /* if one of both == 0, the relative difference cannot be calculated */
  if (a==0 || b==0) 
    return 0;

  a = fabs(a);
  b = fabs(b);

  /* sort a,b */
  if (a>b)  {
    temp=a;
    a=b;
    b=temp;
  }

  if ( diff/a < DOUBLE_RELATIVE_ERROR)
    return 1;

  return 0;
}




/***********************************************************************************/
/* Function: sort_long                                                    @31_30i@ */
/* Description:                                                                    */
/*  Sort two long integers in ascending order.                                     */
/*                                                                                 */
/* Parameters:                                                                     */
/*  long *x1:      Pointer to first integer.                                       */
/*  long *x2:      Pointer to second integer.                                      */
/*                                                                                 */
/* Return value:                                                                   */
/*  None.                                                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

void sort_long (long *x1, long *x2)
{
  long dummy=0;
  
  if (*x1>*x2)  {
    dummy = *x2;
    *x2 = *x1;
    *x1 = dummy;
  }
}


/***********************************************************************************/
/* Function: sort_double                                                  @31_30i@ */
/* Description:                                                                    */
/*  Sort two doubles in ascending order.                                           */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x1:      Pointer to first double.                                      */
/*  double *x2:      Pointer to second double.                                     */
/*                                                                                 */
/* Return value:                                                                   */
/*  None.                                                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

void sort_double (double *x1, double *x2)
{
  double dummy=0;
  
  if (*x1>*x2)  {
    dummy = *x2;
    *x2 = *x1;
    *x1 = dummy;
  }
}



/***********************************************************************************/
/* Function: fak                                                          @31_30i@ */
/* Description:                                                                    */
/*  Calculate the faculty n! of an integer number n.                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*  long n:   Function input.                                                      */
/*                                                                                 */
/* Return value:                                                                   */
/*  Result n!                                                                      */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

double fak (long n)
{
  double dummy=0;
  
  if (n<=0) 
    return 1;
  else  {
    dummy = (double) n * fak(n-1);
    return dummy;
  }
}


/***********************/
/* calculates n over m */
/***********************/

long over (long n, long m)  
{
  long dummy = (long) ( fak (n) / fak (m) / fak (n-m) + 0.5);
  return dummy;
}




/************************************************************/ 
/* average a data set of n points with a gaussian function  */
/* of base width "width"                                    */
/************************************************************/ 

void average (long width, long *y, long n)
{
  long i=0, j=0;
  long counter=0;
  double sum=0;
  double *weight = (double *) calloc (width+1, sizeof(double));
  double *new_y  = (double *) calloc (n      , sizeof(double));
  
  for (i=0; i<n; i++) 
    new_y [i] = 0;
  
  for (i=0; i<=width; i++)  
    weight [i] = (double) over (width, i);
  
  
  for (i=0; i<n; i++)  {
    sum=0;
    for (j=0; j<=width; j++)  {
      counter = i-width/2+j;
      if (counter >= 0 && counter < n)  {  
	new_y[i] += ( (double) y[counter] * weight[j]);
	sum += weight[j];     
      }
    }  
    new_y[i] /= sum;  
  }
  
  for (i=0; i<n; i++)  
    y[i] = (long) new_y[i];

  free(weight);
  free(new_y);
}


/***********************************/
/* calculate mean of n data points */
/***********************************/

double mean (double *x, int n)
{
  int i=0;
  double tmp=0;
  
  for (i=0; i<n; i++)
    tmp += x[i];
  
  tmp /= n;
  return tmp;
}


/********************************************/
/* calculate WEIGHTED mean of n data points */
/********************************************/

double weight_mean (double *x, double *sigma, int n)
{
  int i=0;
  double tmp=0;
  double numerator=0;
  
  for (i=0; i<n; i++)   {
    tmp += x[i] / (sigma[i]*sigma[i]);
    numerator += 1/(sigma[i]*sigma[i]);
  }
  
  tmp /= numerator;
  return tmp;
}


/*************************************************/
/* calculate standard deviation of n data points */
/*************************************************/

double standard_deviation (double *x, int n)
{
  int i=0;
  double tmp=0;
  
  double mu = mean (x, n);
  
  for (i=0; i<n; i++)
    tmp += ( (x[i] - mu) * (x[i] - mu) );
  
  tmp /= (float) (n-1);
  tmp = sqrt(tmp);
  
  return tmp;
}


/**********************************************************/
/* calculate weighted standard deviation of n data points */
/**********************************************************/

double weight_standard_deviation (double *x, double *sigma, int n)
{
  int i=0;
  double tmp=0;
  double numerator=0;
  
  double mu = weight_mean (x, sigma, n);
  
  for (i=0; i<n; i++)  {
    tmp += ( (x[i] - mu) * (x[i] - mu) / (sigma[i]*sigma[i]) );
    numerator += (1 / (sigma[i] * sigma[i]) );
  }
  
  tmp *= ( (float) n / (float) (n-1));
  tmp /= numerator;
  
  tmp = sqrt(tmp);
  
  return tmp;
}
