/*--------------------------------------------------------------------
 * $Id: cnv.c 2623 2011-12-23 10:52:38Z robert.buras $
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


/* @31c@ */
/* @c31@ */

/****************************************************************************/
/* In order to use the functions provided by the numeric library,  @31_10c@ */
/* #include <numeric.h> in your source code and link with libRadtran_c.a.   */
/*                                                                          */
/* @strong{Example:}                                                        */
/* Example for a source file:                                               */
/* @example                                                                 */
/*                                                                          */
/*   ...                                                                    */
/*   #include "../src_c/numeric.h"                                          */
/*   ...                                                                    */
/*                                                                          */
/* @end example                                                             */
/*                                                                          */
/* Linking of the executable, using the GNU compiler gcc:                   */
/* @example                                                                 */
/*                                                                          */
/*   gcc -o test test.c -lRadtran_c -L../lib                                */
/*                                                                          */
/* @end example                                                   @c31_10@  */
/****************************************************************************/


/****************************************************************************/
/* The numeric library provides various numeric functions. The     @31_20c@ */
/* source code is split up into the following source files:                 */
/* @table @asis                                                             */
/* @item cnv.c                                                              */
/*    Data convolution.                                                     */
/* @item equation.c                                                         */
/*    Solve systems of linear equations.                                    */
/* @item function.c                                                         */
/*    Miscellaneous functions.                                              */
/* @item integrat.c                                                         */
/*    Numerical integration.                                                */
/* @item linear.c                                                           */
/*    Linear interpolation.                                                 */
/* @item regress.c                                                          */
/*    Linear regression and related things.                                 */
/* @item spl.c                                                              */
/*    Spline interpolation and approximation.                               */
/* @end table                                                               */
/*                                                                 @c31_20@ */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "numeric.h"

/***********************************************************************************/
/* Function: convolute                                                    @31_30i@ */
/* Description:                                                                    */
/*  Convolute a dataset (x_spec[i], y_spec[i]) with another data set               */
/*  (x_conv[i], y_conv[i]). The data sets must obey the following principles:      */
/*  (1) Both datasets must be defined in equidistant steps; (2) the step width     */
/*  must be the same for both datasets; and (3) 0 must be a point of the grid of   */
/*  the convolution function, x_conv[]. The results is stored in                   */
/*  (x_spec_conv[i], y_spec_conv[i]), i=0...spec_conv_num, the memory of which     */
/*  is allocated automatically.                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x_spec:        x values of the data points, i=0...spec_num-1           */
/*  double *y_spec:        y values of the data points, i=0...spec_num-1           */
/*  int spec_num:          number of data points                                   */
/*  double *x_conv:        x values of the convolution function, i=0...conv_num-1  */
/*  double *y_conv:        y values of the convolution function, i=0...conv_num-1  */
/*  int conv_num:          number of convolution function data points              */
/*  double **x_spec_conv:  x values of the convoluted spectrum                     */
/*  double **y_spec_conv:  y values of the convoluted spectrum                     */
/*  int spec_conv_num:     number of result data points                            */
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

int convolute (double *x_spec, double *y_spec, int spec_num,
	       double *x_conv, double *y_conv, int conv_num,
	       double **x_spec_conv, double **y_spec_conv, int *spec_conv_num)
{

  double spec_delta=0, conv_delta=0;
  double sum=0;
  int i=0, index=0, spec_index=0;
  int mid=0;

  
  /* check for equidistant values */

  if (spec_num > 1)
    spec_delta = x_spec[1] - x_spec[0];

  for (i=2; i<spec_num; i++)  
    if (!double_equal(x_spec[i]-x_spec[i-1], spec_delta)) 
      return SPEC_NOT_EQUIDISTANT;

  if (conv_num > 1)
    conv_delta = x_conv[1] - x_conv[0];

  for (i=2; i<conv_num; i++)  
    if (!double_equal(x_conv[i]-x_conv[i-1], conv_delta)) 
      return CONV_NOT_EQUIDISTANT;

  if (!double_equal(conv_delta, spec_delta))  
    return (SPEC_CONV_DIFFERENT);




  /* look for center wavelength of convolution function */
  mid=0;
  while (x_conv[mid++] != 0.0)
    if (mid == conv_num) 
      return CONV_NOT_CENTERED;

  mid--;


  /* number of values */
  *spec_conv_num = spec_num;

  /* allocate memory for convoluted function */
  *x_spec_conv = (double *) calloc (*spec_conv_num, sizeof(double));
  *y_spec_conv = (double *) calloc (*spec_conv_num, sizeof(double));



  /* do convolution */

  for (i=0; i<spec_num; i++)  {
    sum=0.0;
    for (index=0; index<conv_num; index++)  {
      spec_index = i - mid + index;
      if (spec_index >= 0 && spec_index < spec_num)  {
	(*y_spec_conv)[i] += y_conv[index] * y_spec[spec_index];
	sum+=y_conv[index];
      }
    }

    if (sum != 0.0)
      (*y_spec_conv)[i] /= sum;
    else
      (*y_spec_conv)[i] = 0;
      
    (*x_spec_conv)[i] = x_spec[i];
  }


  return 0;
}



/***********************************************************************************/
/* Function: int_convolute                                                @31_30i@ */
/* Description:                                                                    */
/*  Convolute a dataset (x_spc[i], y_spc[i]) with another data set                 */
/*  (x_conv[i], y_conv[i]). In contrast to the conv() function provided by         */
/*  this library, the dataset to be convoluted and the convolution function        */
/*  may be defined on different grids. While x_spc[] may be an arbitrary grid,     */
/*  x_conv[] must obey the following principles: (1) The spacing must be           */
/*  equidistant; and (2) 0 must be part of the grid. During the convolution        */
/*  process, the function to be convoluted is interpolated to the grid defined by  */
/*  x_conv[]. It is therefore necessary to specify a fine enough grid x_conv[],    */
/*  even if it only describes, e.g., a triangle. The output is stored in           */
/*  y_spec_conv[], which is evaluated at the original data points x_spc[].         */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *x_spc:         x values of the data points, i=0...spc_num-1            */
/*  double *y_spc:         y values of the data points, i=0...spc_num-1            */
/*  int spc_num:           number of data points                                   */
/*  double *x_conv:        x values of the convolution function, i=0...conv_num-1  */
/*  double *y_conv:        y values of the convolution function, i=0...conv_num-1  */
/*  int conv_num:          number of convolution function data points              */
/*  double **y_spec_conv:  convoluted spectrum on the original grid x_spc,         */
/*                         i=0...spc_num-1                                         */
/*  int std:               this is a special feature for covolution of MYSTIC      */
/*                         standard deviations - for those we need to square       */
/*                         the coefficients, but AFTER normalization of the        */
/*                         not squared convolution function: that is, convolute    */
/*                         with (y_conf[index] / sum)**2                           */
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

int int_convolute (double *x_spc,  double *y_spc,  int spc_num,
		   double *x_conv, double *y_conv, int conv_num,
		   double **y_spec_conv, int std)
{

  double interval=0, ratio=0;
  double *x_spec=NULL, *y_spec=NULL; 
  double *x_spec_tmp=NULL, *y_spec_tmp=NULL;
  int spec_tmp_num=0, spec_num=0;

  int i=0, index=0, spec_index=0, status=0;
  int mid=0;

  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;

  double sum=0, stepwidth=0;

  /* check for equidistant values */
  if (conv_num > 1)
    stepwidth = x_conv[1] - x_conv[0];

  for (i=2; i<conv_num; i++)  
    if (!double_equal(x_conv[i]-x_conv[i-1], stepwidth)) 
      return CONV_NOT_EQUIDISTANT;

  /* interpolate spectrum to stepwidth of convolution function */

  /* calculate LINEAR interpolation coefficients */
  status = linear_coeffc (x_spc, y_spc, spc_num, &a0,  &a1,  &a2,  &a3);
  if (status!=0)  {
    fprintf (stderr, " ... int_convolute(): error calculating interpolation coefficients\n");
    return status;
  }
  
  /* calculate number of interpolated values */
  interval = x_spc[spc_num-1] - x_spc[0];
  ratio    = interval/stepwidth;
  
  spec_num = (int) ceil(ratio) + 1;

  /* allocate memory for interpolated spectrum */
  x_spec  = calloc (spec_num, sizeof(double));
  y_spec  = calloc (spec_num, sizeof(double));

  
  /* calculate interpolated values */
  for (i=0; i<spec_num; i++)  {
    x_spec[i] = x_spc[0] + (double) i * stepwidth;

    status = calc_splined_value (x_spec[i], &(y_spec[i]), 
				 x_spc, spc_num, 
				 a0, a1, a2, a3);

    if (status!=0)  {
      if (i==spec_num-1)     /* if last data point */
	y_spec[i] = y_spec[i-1];
      else  {
	fprintf (stderr, " ... int_convolute(): error interpolating spectrum at %g\n", x_spec[i]);
	return status;
      }
    }
  }
   
  /* free interpolation coefficients */
  free(a0);
  free(a1);
  free(a2);
  free(a3);
  


  /* look for center wavelength of convolution function */
  mid=0;
  while (x_conv[mid++] != 0)
    if (mid == conv_num) 
      return CONV_NOT_CENTERED;

  mid--;


  /* number of values */
  spec_tmp_num = spec_num;

  /* allocate memory for convoluted function */
  x_spec_tmp = (double *) calloc (spec_tmp_num, sizeof(double));
  y_spec_tmp = (double *) calloc (spec_tmp_num, sizeof(double));



  /* do convolution */

  for (i=0; i<spec_num; i++)  {
    sum=0.0;
    for (index=0; index<conv_num; index++)  {
      spec_index = i - mid + index;
      if (spec_index >= 0 && spec_index < spec_num)  {

	if (!std)    /* "normal" convolution */
	  y_spec_tmp[i] += y_conv[index] * y_spec[spec_index];
	else         /* convolution of standard deviation */
	  y_spec_tmp[i] += y_conv[index] * y_conv[index] * y_spec[spec_index];
	
	sum+=y_conv[index];
      }
    }

    if (sum != 0.0) {
      if (!std)    /* "normal" convolution */
	y_spec_tmp[i] /= sum;
      else         /* convolution of standard deviation */
	y_spec_tmp[i] /= (sum*sum);
    }	
    else
      y_spec_tmp[i] = 0;
      
    x_spec_tmp[i] = x_spec[i];
  }


  /* now again interpolate convoluted spectra to original steps */

  /* calculate LINEAR interpolation coefficients */
  status = linear_coeffc (x_spec_tmp, y_spec_tmp, spec_tmp_num, &a0,  &a1,  &a2,  &a3);
  if (status!=0)  {
    fprintf (stderr, " ... int_convolute(): error calculating interpolation coefficients\n");
    return status;
  }
  

  /* allocate memory for convoluted function */
  *y_spec_conv = (double *) calloc (spc_num, sizeof(double));

  for (i=0; i<spc_num; i++)  {

    status = calc_splined_value (x_spc[i], &((*y_spec_conv)[i]), 
				 x_spec_tmp, spec_tmp_num, 
				 a0, a1, a2, a3);
    
    if (status!=0)  {
      fprintf (stderr, " ... int_convolute(): error interpolating spectrum\n");
      return status;
    }
  }
    
  free(a0);
  free(a1);
  free(a2);
  free(a3);
  free(x_spec);
  free(y_spec);
  free(x_spec_tmp);
  free(y_spec_tmp);

  
  return 0;
}
