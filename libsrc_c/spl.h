/*--------------------------------------------------------------------
 * $Id: spl.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __spl_h
#define __spl_h

#if defined (__cplusplus)
extern "C" {
#endif


#include "numeric.h"


#define NO_SPLINED_VALUES           -1
#define X_NOT_ASCENDING             -2
#define SPLINE_NOT_POSSIBLE         -3
#define TOO_FEW_DATA_POINTS         -4
#define DATA_NOT_SORTED             -5
#define NEGATIVE_WEIGHTING_FACTORS  -6
#define NO_EXTRAPOLATION            -7
#define NO_INTERPOLATION            -8



/* prototypes */
int spline (double *x, double *y, int number, double start, double step,
	    int *newnumber, double **new_x, double **new_y);
int spline_coeffc (double *x, double *y, int number, 
		   double **a0, double **a1, double **a2, double **a3);
int fspline_coeffc (float *x, float *y, int number, 
		    float **a0, float **a1, float **a2, float **a3);
int appspl (double *x, double *y, double *w, int number, 
	    double start, double step, 
	    int *newnumber, double **new_x, double **new_y);
int appspl_coeffc (double *x, double *y, double *w, int number, 
		   double **a0, double **a1, double **a2, double **a3);
int calc_splined_value (double xnew, double *ynew, 
			double *x, int number, 
			double *a0, double *a1, double *a2, double *a3);
int fcalc_splined_value (float xnew, float *ynew, 
			 float *x, int number, 
			 float *a0, float *a1, float *a2, float *a3);
int linear_eqd (double *x, double *y, int number, double start, double step,
		int *newnumber, double **new_x, double **new_y);

int locate (double *xx, int n, double x);
int flocate (float *xx, int n, float x);

#if defined (__cplusplus)
}
#endif

#endif
