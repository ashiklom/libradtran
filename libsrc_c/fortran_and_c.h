/*--------------------------------------------------------------------
 * $Id: fortran_and_c.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __fortran_and_c_h
#define __fortran_and_c_h

#if defined (__cplusplus)
extern "C" {
#endif

#include <stdlib.h>
#include "ascii.h"

/* prototypes */
  
float* c2fortran_2D_float_ary(int dim1, int dim2, float** ary_2D );
double *c2fortran_2D_double_ary (int dim1, int dim2, double** ary_2D);
float *c2fortran_3D_float_ary (int dim1, int dim2, int dim3, float*** ary_3D); 
double *c2fortran_3D_double_ary (int dim1, int dim2, int dim3, double*** ary_3D); 
float** fortran2c_2D_float_ary(int dim1, int dim2, float* ary_1D );
float*** fortran2c_3D_float_ary(int dim1, int dim2, int dim3, float* ary_1D );
float**** fortran2c_4D_float_ary(int dim1, int dim2, int dim3, int dim4, float* ary_1D );

void fortran2c_2D_float_ary_noalloc (int dim1, int dim2,           float *ary_1D, float **ary_2D);
void fortran2c_3D_float_ary_noalloc (int dim1, int dim2, int dim3, float *ary_1D, float ***ary_3D);
void fortran2c_2D_double_ary_noalloc (int dim1, int dim2, double* ary_1D, double** ary_2D);
void fortran2c_4D_double_ary_noalloc (int dim1, int dim2, int dim3, int dim4, double* ary_1D,
				      double ****ary_4D);


double** fortran2c_2D_double_ary(int dim1, int dim2, double* ary_1D );
double**** fortran2c_4D_double_ary(int dim1, int dim2, int dim3, int dim4, double* ary_1D );

#if defined (__cplusplus)
}
#endif

#endif
