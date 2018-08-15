/*--------------------------------------------------------------------
 * $Id: function.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __function_h
#define __function_h

#if defined (__cplusplus)
extern "C" {
#endif

#include "numeric.h"


/* prototypes */

int double_equal (double a, double b);
void sort_long (long *x1, long *x2);
void sort_double (double *x1, double *x2);
double fak (long n);
long over (long n, long m);
void average (long width, long *y, long n);
double mean (double *x, int n);
double weight_mean (double *x, double *sigma, int n);
double standard_deviation (double *x, int n);
double weight_standard_deviation (double *x, double *sigma, int n);


#if defined (__cplusplus)
}
#endif

#endif


