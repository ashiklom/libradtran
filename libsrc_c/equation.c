/*--------------------------------------------------------------------
 * $Id: equation.c 2623 2011-12-23 10:52:38Z robert.buras $
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
#include <float.h>
#include "numeric.h"


/***********************************************************************************/
/* Function: solve_gauss                                                  @31_30i@ */
/* Description:                                                                    */
/*  Solve a system of n linear equations, A*x = b, using the Gauss algorithm       */
/*  The pivot element is determined using 'relative column maximum strategy'.      */
/*  For a description of the algorithm see H.R.Schwarz: "Numerische Mathematik",   */
/*  pg. 21. Memory for the result vector res is allocated automatically.           */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double **A:            Matrix[n x n]  (see above).                             */
/*  double *b:             Vector[n] (see above).                                  */
/*  double n:              Number of equations.                                    */
/*  double **res:          Pointer to the result vector[n]; if no unique solution  */ 
/*                         exists, *res will be set to NULL.                       */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if no unique solution.                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int solve_gauss (double **A, double *b, int n, double **res)
{
  int i=0, j=0, k=0, p=0;
  double div=0.0;
  double *tmp1=NULL;
  double tmp2=0.0;
  double sum=0.0;
  double max = (0.0 - DBL_MAX);
  

  for (k=0; k<n; k++)  {
    /* get pivot element (line p) */
    max = (0.0 - DBL_MAX);
    p = 0;
    
    for (i=k; i<n; i++)  {
      
      sum = 0.0;
      for (j=k; j<n; j++)
	sum += fabs(A[i][j]);
      if (sum == 0.0)
	return GAUSS_SINGULAR;
      
      if ( (tmp2 = fabs (A[i][k]) / sum) > max )  {
	max = tmp2;
	p = i;
      }
    }
    
    /* exchange lines k and p */
    tmp1 = A[k];
    A[k] = A[p];
    A[p] = tmp1;

    tmp2 = b[k];
    b[k] = b[p];
    b[p] = tmp2;
    
    
    if ( (div = A[k][k]) == 0)   /* no definite solution */
      return GAUSS_SINGULAR;
    
    for (i=k; i<n; i++)
      A[k][i] /= div;
    
    b[k] /= div;
    
    for (i=k+1; i<n; i++)  {
      div = A[i][k];
      for (j=k; j<n; j++)
	A[i][j] = A[i][j] - A[k][j] * div;
      b[i] = b[i] - b[k] * div;
    }
  }
  

  /* allocate memory for result vector */
  *res = (double *) calloc (n, sizeof(double));

  for (i=n-1; i>=0; i--)  {
    (*res)[i] += b[i];
    for (k=i+1; k<n; k++)
      (*res)[i] -= A[i][k] * (*res)[k];
  }

  return 0;
}





/***********************************************************************************/
/* Function: solve_five                                                   @31_30i@ */
/* Description:                                                                    */
/*  Solve a system of n linear equations, A*x = b, where A is a five-diagonal      */
/*  matrix; for details see Engeln-Muellges, pg.95ff.                              */
/*  Memory for the result vector is allocated automatically.                       */
/*  If possible, solve_five_ms() should be preferred to solve_five(), because      */
/*  much less memory is required by the latter.                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double **A:            Matrix[n x n]  (see above).                             */
/*  double *b:             Vector[n] (see above).                                  */
/*  double n:              Number of equations.                                    */
/*  double **res:          Pointer to the result vector[n]; if no unique solution  */ 
/*                         exists, *res will be set to NULL.                       */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if no unique solution.                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int solve_five (double **A, double *b, int n, double **res)
{
  int i=0;

  double *a = (double *) calloc (n+1, sizeof(double));

  double *c=NULL, *d=NULL, *e=NULL, *f=NULL, *g=NULL;
  double *r=NULL;
  double *alpha=NULL, *beta=NULL, *gamma=NULL, *delta=NULL, *epsilon=NULL;

  /* check diagonal elements */
  for (i=0; i<n; i++)
    if (A[i][i] == 0)
      return GAUSS_SINGULAR;


  c = (double *) calloc (n+1, sizeof(double));
  d = (double *) calloc (n+1, sizeof(double));
  e = (double *) calloc (n+1, sizeof(double));
  f = (double *) calloc (n+1, sizeof(double));
  g = (double *) calloc (n+1, sizeof(double));




  /* ATTENTION: arrays a, c, d, e, f, g range from 1 to n !! */

  for (i=2; i<=n; i++)  
    c[i] = A[i-1][i-2];

  for (i=1; i<=n; i++)  
    d[i] = A[i-1][i-1];

  for (i=1; i<=n-1; i++)  
    e[i] = A[i-1][i];

  for (i=1; i<=n-2; i++)  
    f[i] = A[i-1][i+1];

  for (i=3; i<=n; i++)  
    g[i] = A[i-1][i-3];

    
  for (i=1; i<=n; i++) 
    a[i] = b[i-1];

  alpha   = (double *) calloc (n+1, sizeof(double));
  beta    = (double *) calloc (n+1, sizeof(double));
  gamma   = (double *) calloc (n+1, sizeof(double));
  delta   = (double *) calloc (n+1, sizeof(double));
  epsilon = (double *) calloc (n+1, sizeof(double));


  alpha[1] = d[1];
  gamma[1] = e[1]/alpha[1];
  delta[1] = f[1]/alpha[1];
  beta[2]  = c[2];
  alpha[2] = d[2] - beta[2]*gamma[1];
  gamma[2] = (e[2]-beta[2]*delta[1])/alpha[2];
  delta[2] = f[2]/alpha[2];
  
  for (i=3; i<=n-2; i++)  {
    beta[i]  = c[i] - g[i]*gamma[i-2];
    alpha[i] = d[i] - g[i]*delta[i-2] - beta[i]*gamma[i-1];
    gamma[i] = (e[i] - beta[i]*delta[i-1])/alpha[i];
    delta[i] = f[i]/alpha[i];
  }

  beta[n-1]  = c[n-1] - g[n-1]*gamma[n-3];
  alpha[n-1] = d[n-1] - g[n-1]*delta[n-3] - beta[n-1]*gamma[n-2];
  gamma[n-1] = (e[n-1] - beta[n-1]*delta[n-2])/alpha[n-1];
  beta[n]    = c[n] - g[n]*gamma[n-2];
  alpha[n]   = d[n] - g[n]*delta[n-2] - beta[n]*gamma[n-1];

  for (i=3; i<=n; i++)
    epsilon[i] = g[i];

  r = (double *) calloc (n+1, sizeof(double));
  
  r[1] = a[1]/alpha[1];
  r[2] = (a[2] - beta[2]*r[1])/alpha[2];

  for (i=3; i<=n; i++)  
    r[i] = (a[i] - epsilon[i]*r[i-2] - beta[i]*r[i-1])/alpha[i];

  free(c);
  free(d);
  free(e);
  free(f);
  free(g);

  *res = (double *) calloc(n, sizeof(double));
  
  (*res)[n-1] = r[n];
  (*res)[n-2] = r[n-1] - gamma[n-1] * (*res)[n-1];
  
  for (i=n-2; i>=1; i--)
    (*res)[i-1] = r[i] - gamma[i]*(*res)[i] - delta[i]*(*res)[i+1];


  free(r);

  free(alpha);
  free(beta);
  free(gamma);
  free(delta);
  free(epsilon);

  return 0;
}




/***********************************************************************************/
/* Function: solve_five_ms                                                @31_30i@ */
/* Description:                                                                    */
/*  Same function as solve_five(), but much more memory efficient on the cost of   */
/*  some user-friendlyness. In contrast to solve_five(), A is defined as a         */
/*  n*5-matrix, rather than a n*n-matrix.                                          */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double **A:            Matrix[n x 5]  (see above).                             */
/*  double *b:             Vector[n] (see above).                                  */
/*  double n:              Number of equations.                                    */
/*  double **res:          Pointer to the result vector[n]; if no unique solution  */ 
/*                         exists, *res will be set to NULL.                       */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if no unique solution.                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int solve_five_ms (double **A, double *b, int n, double **res)
{
  int i=0;

  double *a=NULL;
  double *c=NULL, *d=NULL, *e=NULL, *f=NULL, *g=NULL;
  double *r = NULL;
  double *alpha=NULL, *beta=NULL, *gamma=NULL, *delta=NULL, *epsilon=NULL;


  /* check diagonal elements */
  for (i=0; i<n; i++)
    if (A[i][2] == 0)
      return GAUSS_SINGULAR;


  a = (double *) calloc (n+1, sizeof(double));
  c = (double *) calloc (n+1, sizeof(double));
  d = (double *) calloc (n+1, sizeof(double));
  e = (double *) calloc (n+1, sizeof(double));
  f = (double *) calloc (n+1, sizeof(double));
  g = (double *) calloc (n+1, sizeof(double));


  /* ATTENTION: arrays a, c, d, e, f, g range from 1 to n !! */

  for (i=2; i<=n; i++)  
    c[i] = A[i-1][1];

  for (i=1; i<=n; i++)  
    d[i] = A[i-1][2];

  for (i=1; i<=n-1; i++)  
    e[i] = A[i-1][3];

  for (i=1; i<=n-2; i++)  
    f[i] = A[i-1][4];

  for (i=3; i<=n; i++)  
    g[i] = A[i-1][0];

    
  for (i=1; i<=n; i++) 
    a[i] = b[i-1];


  alpha   = (double *) calloc (n+1, sizeof(double));
  beta    = (double *) calloc (n+1, sizeof(double));
  gamma   = (double *) calloc (n+1, sizeof(double));
  delta   = (double *) calloc (n+1, sizeof(double));
  epsilon = (double *) calloc (n+1, sizeof(double));


  alpha[1] = d[1];
  gamma[1] = e[1]/alpha[1];
  delta[1] = f[1]/alpha[1];
  beta[2]  = c[2];
  alpha[2] = d[2] - beta[2]*gamma[1];
  gamma[2] = (e[2]-beta[2]*delta[1])/alpha[2];
  delta[2] = f[2]/alpha[2];
  
  for (i=3; i<=n-2; i++)  {
    beta[i]  = c[i] - g[i]*gamma[i-2];
    alpha[i] = d[i] - g[i]*delta[i-2] - beta[i]*gamma[i-1];
    gamma[i] = (e[i] - beta[i]*delta[i-1])/alpha[i];
    delta[i] = f[i]/alpha[i];
  }

  beta[n-1]  = c[n-1] - g[n-1]*gamma[n-3];
  alpha[n-1] = d[n-1] - g[n-1]*delta[n-3] - beta[n-1]*gamma[n-2];
  gamma[n-1] = (e[n-1] - beta[n-1]*delta[n-2])/alpha[n-1];
  beta[n]    = c[n] - g[n]*gamma[n-2];
  alpha[n]   = d[n] - g[n]*delta[n-2] - beta[n]*gamma[n-1];

  for (i=3; i<=n; i++)
    epsilon[i] = g[i];

  r = (double *) calloc (n+1, sizeof(double));

  r[1] = a[1]/alpha[1];
  r[2] = (a[2] - beta[2]*r[1])/alpha[2];

  for (i=3; i<=n; i++)  
    r[i] = (a[i] - epsilon[i]*r[i-2] - beta[i]*r[i-1])/alpha[i];
    
  free(a);
  free(c);
  free(d);
  free(e);
  free(f);
  free(g);
    
  /* allocate memory for result vector */
  *res = (double *) calloc(n, sizeof(double));

  (*res)[n-1] = r[n];
  (*res)[n-2] = r[n-1] - gamma[n-1] * (*res)[n-1];
  
  for (i=n-2; i>=1; i--)
    (*res)[i-1] = r[i] - gamma[i] * (*res)[i] - delta[i] * (*res)[i+1];

  free(r);
  
  free(alpha);
  free(beta);
  free(gamma);
  free(delta);
  free(epsilon);

  return 0;
}





/***********************************************************************************/
/* Function: solve_three                                                  @31_30i@ */
/* Description:                                                                    */
/*  Solve a system of n linear equations, A*x = b, where A is a three-diagonal     */
/*  matrix; for details see Engeln-Muellges, pg.95ff.                              */
/*  Memory for the result vector is allocated automatically.                       */
/*  If possible, solve_three_ms() should be preferred to solve_three(), because    */
/*  much less memory is required by the latter.                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double **A:            Matrix[n x n]  (see above).                             */
/*  double *b:             Vector[n] (see above).                                  */
/*  double n:              Number of equations.                                    */
/*  double **res:          Pointer to the result vector[n]; if no unique solution  */ 
/*                         exists, *res will be set to NULL.                       */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if no unique solution.                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int solve_three (double **A, double *b, int n, double **res)
{
  int i=0;

  double *a = (double *) calloc (n+1, sizeof(double));

  double *c=NULL, *d=NULL, *e=NULL;
  double *r=NULL;
  double *alpha=NULL, *gamma=NULL;


  /* check diagonal elements */
  for (i=0; i<n; i++)
    if (A[i][i] == 0)
      return GAUSS_SINGULAR;


  c = (double *) calloc (n+1, sizeof(double));
  d = (double *) calloc (n+1, sizeof(double));
  e = (double *) calloc (n+1, sizeof(double));


  /* ATTENTION: arrays a, c, d, e range from 1 to n !! */

  for (i=1; i<=n-1; i++)  
    c[i] = A[i-1][i];

  for (i=1; i<=n; i++)  
    d[i] = A[i-1][i-1];

  for (i=2; i<=n; i++)  
    e[i] = A[i-1][i-2];


  for (i=1; i<=n; i++) 
    a[i] = b[i-1];

  alpha = (double *) calloc (n+1, sizeof(double));
  gamma = (double *) calloc (n+1, sizeof(double));

  alpha[1] = d[1];
  gamma[1] = c[1]/alpha[1];
  
  for (i=2; i<=n-1; i++)  {
    if ( (alpha[i] = d[i] - e[i]*gamma[i-1]) == 0)  {
      free(a);
      free(c);
      free(d);
      free(e);
      free(alpha);
      free(gamma);
      return GAUSS_SINGULAR;
    }
    gamma[i] = c[i]/alpha[i];
  }

  alpha[n] = d[n] - e[n]*gamma[n-1];

  r = (double *) calloc (n+1, sizeof(double));

  r[1] = a[1]/d[1];

  for (i=2; i<=n; i++)  
    r[i] = (a[i] - e[i]*r[i-1])/alpha[i];

  free(a);
  free(c);
  free(d);
  free(e);

  *res = (double *) calloc(n, sizeof(double));

  (*res)[n-1] = r[n];
  
  for (i=n-1; i>=1; i--)
    (*res)[i-1] = r[i] - gamma[i] * (*res)[i];


  free(r);

  free(alpha);
  free(gamma);

  return 0;
}




/***********************************************************************************/
/* Function: solve_three_ms                                               @31_30i@ */
/* Description:                                                                    */
/*  Same function as solve_three(), but much more memory efficient on the cost of  */
/*  some user-friendlyness. In contrast to solve_three(), A is defined as a        */
/*  n*3-matrix, rather than a n*n-matrix.                                          */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double **A:            Matrix[n x 3]  (see above).                             */
/*  double *b:             Vector[n] (see above).                                  */
/*  double n:              Number of equations.                                    */
/*  double **res:          Pointer to the result vector[n]; if no unique solution  */ 
/*                         exists, *res will be set to NULL.                       */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if no unique solution.                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int solve_three_ms (double **A, double *b, int n, double **res)
{
  int i=0;

  double *r=NULL;
  double *alpha=NULL, *gamma=NULL;

  /* check diagonal elements */
  for (i=0; i<n; i++)
    if (A[i][1] == 0)
      return GAUSS_SINGULAR;

  alpha   = (double *) calloc (n+1, sizeof(double));
  gamma   = (double *) calloc (n+1, sizeof(double));

  alpha[1] = A[0][1];
  gamma[1] = A[0][2]/alpha[1];
  
  for (i=2; i<=n-1; i++)  {
    if ( (alpha[i] = A[i-1][1] - A[i-1][0]*gamma[i-1]) == 0)  {
      free(alpha);
      free(gamma);
      return GAUSS_SINGULAR;
    }
    gamma[i] = A[i-1][2]/alpha[i];
  }

  alpha[n] = A[n-1][1] - A[n-1][0]*gamma[n-1];


  r = (double *) calloc (n+1, sizeof(double));


  r[1] = b[0]/A[0][1];

  for (i=2; i<=n; i++)  
    r[i] = (b[i-1] - A[i-1][0]*r[i-1])/alpha[i];

  /* allocate memory for result vector */
  *res = (double *) calloc (n, sizeof(double));

  (*res)[n-1] = r[n];
  
  for (i=n-1; i>=1; i--)
    (*res)[i-1] = r[i] - gamma[i] * (*res)[i];

  free(r);
  
  free(alpha);
  free(gamma);

  return 0;
}



/***********************************************************************************/
/* Function: fsolve_three_ms                                              @31_30i@ */
/* Description:                                                                    */
/*  Similar to solve_three_ms, only all data are of type float instead of          */
/*  float. Same function as solve_three(), but much more memory efficient on       */
/*  the cost of some user-friendlyness. In contrast to solve_three(), A is defined */
/*  as a n*3-matrix, rather than a n*n-matrix.                                     */
/*                                                                                 */
/* Parameters:                                                                     */
/*  float **A:            Matrix[n x 3]  (see above).                              */
/*  float *b:             Vector[n] (see above).                                   */
/*  float n:              Number of equations.                                     */
/*  float **res:          Pointer to the result vector[n]; if no unique solution   */ 
/*                         exists, *res will be set to NULL.                       */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if no unique solution.                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

int fsolve_three_ms (float **A, float *b, int n, float **res)
{
  int i=0;

  float *r=NULL;
  float *alpha=NULL, *gamma=NULL;

  /* check diagonal elements */
  for (i=0; i<n; i++)
    if (A[i][1] == 0)
      return GAUSS_SINGULAR;

  alpha   = (float *) calloc (n+1, sizeof(float));
  gamma   = (float *) calloc (n+1, sizeof(float));

  alpha[1] = A[0][1];
  gamma[1] = A[0][2]/alpha[1];
  
  for (i=2; i<=n-1; i++)  {
    if ( (alpha[i] = A[i-1][1] - A[i-1][0]*gamma[i-1]) == 0)  {
      free(alpha);
      free(gamma);
      return GAUSS_SINGULAR;
    }
    gamma[i] = A[i-1][2]/alpha[i];
  }

  alpha[n] = A[n-1][1] - A[n-1][0]*gamma[n-1];


  r = (float *) calloc (n+1, sizeof(float));


  r[1] = b[0]/A[0][1];

  for (i=2; i<=n; i++)  
    r[i] = (b[i-1] - A[i-1][0]*r[i-1])/alpha[i];

  /* allocate memory for result vector */
  *res = (float *) calloc (n, sizeof(float));

  (*res)[n-1] = r[n];
  
  for (i=n-1; i>=1; i--)
    (*res)[i-1] = r[i] - gamma[i] * (*res)[i];

  free(r);
  
  free(alpha);
  free(gamma);

  return 0;
}

