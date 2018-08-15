/*--------------------------------------------------------------------
 * $Id: wavelength.c 2623 2011-12-23 10:52:38Z robert.buras $
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

/* @69c@ */
/* @c69@ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numeric.h"
#include "wavelength.h"


/*********************************************************************************/
/* In order to use the functions provided by the misc library,          @69_10c@ */
/* #include <misc.h> in your source code and link with libRadtran_c.a.  @c69_10@ */
/*********************************************************************************/





/***********************************************************************************/
/* Function: air_refraction                                               @69_30i@ */
/* Description:                                                                    */
/*  Calculate index of refraction of air for `standard air' according to the       */
/*  1997/98 `CRC handbook of Chemistry and Physics'. `Standard air' refers to      */
/*  a temperature of 15 deg C and a pressure of 1013.25 mbar. The index of         */
/*  refraction n is defined as a function of vacuum wavelength, but due to the     */
/*  slow variation of n, the errors are negligibly small when using the air        */
/*  wavelength as input for air_refraction(). The formula is valid between         */
/*  200 and 2000 nm.                                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double lambda:  Wavelength in nm.                                              */
/*                                                                                 */
/* Return value:                                                                   */
/*  Index of refraction (double)                                                   */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i69_30@ */
/***********************************************************************************/

double air_refraction (double lambda)
{	
  double sigma2 = 1.0 / (lambda/1000.0) / (lambda/1000.0);
  double n = 1.0 + 1.0E-8 * (8342.13 + 2406030.0/(130.0-sigma2) + 15997.0/(38.9-sigma2));
  
  return n;
}





/***********************************************************************************/
/* Function: vac2air                                                      @69_30i@ */
/* Description:                                                                    */
/*  Shift a spectrum from vacuum to air or vice versa. The shifted data are        */
/*  re-gridded to the original wavelength grid. The index of refraction is         */
/*  calculated with air_refraction().                                              */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double *lambda:    Original wavelength grid, i = 0..rows                       */
/*  double *irradince: Original irradiance data, i = 0..rows                       */
/*  int rows:          Number of data pairs.                                       */
/*  char *reverse:     If 0, convert from vacuum to air, else vice versa.          */
/*  char *linear:      If 0, use spline interpolation, else linear.                */
/*  double **irradiance\_shifted:  Shifted irradiance, i = 0..rows, referring      */
/*                                 to the original wavelength grid; memory is      */
/*                                 allocated automatically.                        */ 
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i69_30@ */
/***********************************************************************************/

int vac2air (double *lambda, double *irradiance, int rows,
	     char reverse, char linear,
	     double **irradiance_shifted)
{
  int j=0, status=0;
  double *d0=NULL, *d1=NULL, *d2=NULL, *d3=NULL;
  double indr=0, di_shifted=0;
  double lambda_shifted=0;

  
  /* allocate memory for shifted spectrum */
  *irradiance_shifted = calloc (rows, sizeof(double));


  
  /* calculate spline coefficients for irradiance */
  if (!linear)   /* spline interpolation */
    status = spline_coeffc (lambda, irradiance, rows, &d0, &d1, &d2, &d3);
  else           /* linear interpolation */
    status = linear_coeffc (lambda, irradiance, rows, &d0, &d1, &d2, &d3);
  
  if (status!=0)
    return status;
  
  
  for (j=0; j<rows; j++)  {
    
    /* calculate shifted wavelength */
    indr = air_refraction (lambda[j]);
    
    
    if (!reverse)    /* vacuum to air */
      lambda_shifted = lambda[j] / indr;
    else             /* air to vacuum */
      lambda_shifted = lambda[j] * indr;
    
    
    /* interpolate irradiance */
    status = calc_splined_value (lambda[j] - (lambda_shifted - lambda[j]), &di_shifted, 
				 lambda, rows, d0, d1, d2, d3);
    
    
    /* store value in a shifted array */
    if (status != 0) {
      di_shifted = NAN;
      fprintf (stderr, "Warning: cannot interpolate at %g nm, inserting %g\n", 
	       lambda[j], di_shifted);
    }

    (*irradiance_shifted)[j] = di_shifted;
  }
  
  /* free memory of coefficients */
  free(d0); free(d1); free(d2); free(d3);

  /* if o.k. */
  return 0;
}
