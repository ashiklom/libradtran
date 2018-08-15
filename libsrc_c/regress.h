/*--------------------------------------------------------------------
 * $Id: regress.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __regress_h
#define __regress_h

#if defined (__cplusplus)
extern "C" {
#endif

#include "numeric.h"


#define FIT_NOT_POSSIBLE -30


/* prototypes */

double gauss             (double x, double mu, double sigma);
double gauss_distorted   (double x, double a0, double a1, double a2, double a3);

int regression        (double *x, double *y, int n, 
		       double *alpha, double *beta,
		       double *sigma_a, double *sigma_b,
		       double *correlation);
double weight_regression (double *x, double *y, double *sigma, int n,
                          double *alpha, double *beta,
                          double *sigma_a, double *sigma_b);
int gaussfit             (double *x, double *y, long number,
			  double *mu, double *sigma, double *area);
int gaussfit_distorted   (double *x, double *y, long number,
			  double *a0, double *a1, double *a2, double *a3);
int boltzmannfit         (double *x, double *y, long number, 
			  double *a, double *b);
int exponentialfit       (double *x, double *y, long number, 
			  double *a, double *b);
int parabolafit          (double *x, double *y, long number, 
			  double *a0, double *a1, double *a2);
int cubicfit             (double *x, double *y, long number, 
			  double *a0, double *a1, double *a2, double *a3);
int hyperbolafit         (double *x, double *y, long number, 
			  double *a, double *b);
int inv_parabolafit      (double *x, double *y, long number, 
			  double *a0, double *a1, double *a2);


#if defined (__cplusplus)
}
#endif

#endif

