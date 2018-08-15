/*--------------------------------------------------------------------
 * $Id: fortran_and_c.c 2623 2011-12-23 10:52:38Z robert.buras $
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


/* @32c@ */
/* @c32@ */

/****************************************************************************/
/* In order to use the functions provided by the fortran_and_c    @32_10c@  */
/* library, #include <fortran_and_c.h> in your source code and link         */
/* with libRadtran_c.a.                                                     */
/*                                                                          */
/* @strong{Example:}                                                        */
/* Example for a source file:                                               */
/* @example                                                                 */
/*                                                                          */
/*   ...                                                                    */
/*   #include <fortran_and_c.h>                                             */
/*   ...                                                                    */
/*                                                                          */
/* @end example                                                             */
/*                                                                          */
/* Linking of the executable, using the GNU compiler gcc:                   */
/* @example                                                                 */
/*                                                                          */
/*   gcc -o test test.c -lRadtran_c -L../lib                                */
/*                                                                          */
/* @end example                                                   @c32_10@  */
/****************************************************************************/


/****************************************************************************/
/* The fortran_and_c library provides functions converting        @32_20c@  */
/* from multidimensional C arrays to one-dimensional fortran arrays         */
/* that can be input to fortran routines, and functions for converting      */
/* one-dimensional fortran compatible arrays to multidimensional C arrays.  */
/*                                                                          */
/* These functions are useful when calling fortran functions and            */
/* subroutines from C.                                                      */
/*                                                                 @c32_20@ */
/****************************************************************************/

#include "fortran_and_c.h"


/***********************************************************************************/
/* Function: c2fortran_2D_float_ary                                       @32_30i@ */
/* Description:                                                                    */
/* Convert 2D C float array to 1D column float array and map it such that it       */
/* can be passed to a fortran routine.                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  float **ary_2D:    pointer to the 2D C array                                   */
/*                                                                                 */
/* Return value:                                                                   */
/*  Pointer to the 1D fortran compatible array. Memory allocation is done          */
/*  automatically and can be freed with a simple free().                           */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

float *c2fortran_2D_float_ary (int dim1, int dim2, float** ary_2D)
{
  int i1=0, i2=0, j=0;
  float *tmp = NULL;

  if (dim1*dim2 > 0) {
    tmp = (float *) calloc (dim1*dim2, sizeof(float));
    
    for (i1=0; i1<dim1; i1++)
      for (i2=0;i2<dim2;i2++)
	tmp[j++] = ary_2D[i1][i2];
  }

  return tmp;
}


double *c2fortran_2D_double_ary (int dim1, int dim2, double** ary_2D)
{
  int i1=0, i2=0, j=0;
  double *tmp = NULL;

  if (dim1*dim2 > 0) {
    tmp = (double *) calloc (dim1*dim2, sizeof(double));
    
    for (i1=0; i1<dim1; i1++)
      for (i2=0;i2<dim2;i2++)
	tmp[j++] = ary_2D[i1][i2];
  }

  return tmp;
}

float *c2fortran_3D_float_ary (int dim1, int dim2, int dim3, float*** ary_3D)
{
  int i1=0, i2=0, i3=0, j=0;
  float *tmp = NULL;

  if (dim1*dim2*dim3 > 0) {
    tmp = (float *) calloc (dim1*dim2*dim3, sizeof(float));
    
    for (i1=0; i1<dim1; i1++)
      for (i2=0; i2<dim2; i2++)
        for (i3=0; i3<dim3; i3++)
          tmp[j++] = ary_3D[i1][i2][i3];
  }

  return tmp;
}

double *c2fortran_3D_double_ary (int dim1, int dim2, int dim3, double*** ary_3D)
{
  int i1=0, i2=0, i3=0, j=0;
  double *tmp = NULL;

  if (dim1*dim2*dim3 > 0) {
    tmp = (double *) calloc (dim1*dim2*dim3, sizeof(double));
    
    for (i1=0; i1<dim1; i1++)
      for (i2=0; i2<dim2; i2++)
        for (i3=0; i3<dim3; i3++)
          tmp[j++] = ary_3D[i1][i2][i3];
  }

  return tmp;
}


/***********************************************************************************/
/* Function: fortran2c_2D_float_ary                                       @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible float array to 2D C float array. Space for the    */
/* returned ary is automatically allocated.                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  float *ary_1D:     pointer to the 1D fortran compatible array                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  Pointer to the 2D C array. Memory allocation is done automatically.            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

float **fortran2c_2D_float_ary (int dim1, int dim2, float* ary_1D)
{
  int i1=0, i2=0, j=0, status=0;
  float **tmp=NULL;

  status = ASCII_calloc_float (&tmp, dim1, dim2);
  if (status==0) {
    for (i1=0; i1<dim1; i1++)
      for (i2=0; i2<dim2; i2++)
	tmp[i1][i2] = ary_1D[j++];
  }

  return tmp;
}

/***********************************************************************************/
/* Function: fortran2c_2D_float_ary_noalloc                               @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible float array to a pre-allocated 2D C float array.  */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  float *ary_1D:     pointer to the 1D fortran compatible array                  */
/*  float *ary_2D:     pointer to the 2D C compatible array                        */  
/*                                                                                 */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

void fortran2c_2D_float_ary_noalloc (int dim1, int dim2, float* ary_1D, float** ary_2D)
{
  int i1=0, i2=0, j=0;
  
  for (i1=0; i1<dim1; i1++)
    for (i2=0; i2<dim2; i2++)
      ary_2D[i1][i2] = ary_1D[j++];
}

/***********************************************************************************/
/* Function: fortran2c_2D_double_ary_noalloc                              @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible double array to 2D C double array.                */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  double *ary_1D:    pointer to the 1D fortran compatible array                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  Pointer to the 2D C array. Memory allocation is done automatically.            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

void fortran2c_2D_double_ary_noalloc (int dim1, int dim2, double* ary_1D, double** ary_2D)
{
  int i1=0, i2=0, j=0;
  
  for (i1=0; i1<dim1; i1++)
    for (i2=0; i2<dim2; i2++) 
      ary_2D[i1][i2] = ary_1D[j++];
}


/***********************************************************************************/
/* Function: fortran2c_2D_double_ary                                      @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible double array to 2D C double array. Space for the  */
/* returned ary is automatically allocated.                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  double *ary_1D:    pointer to the 1D fortran compatible array                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  Pointer to the 2D C array. Memory allocation is done automatically.            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

double **fortran2c_2D_double_ary (int dim1, int dim2, double* ary_1D )
{
  int i1=0, i2=0, j=0, status=0;
  double **tmp=NULL;

  status = ASCII_calloc_double (&tmp, dim1, dim2);
  if (status==0) {
    for (i1=0; i1<dim1; i1++)
      for (i2=0; i2<dim2; i2++)
	tmp[i1][i2] = ary_1D[j++];
  }

  return tmp;
}

/***********************************************************************************/
/* Function: fortran2c_3D_float_ary_noalloc                               @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible float array to a pre-allocated 3D C float array.  */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  float *ary_1D:     pointer to the 1D fortran compatible array                  */
/*                                                                                 */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

void fortran2c_3D_float_ary_noalloc (int dim1, int dim2, int dim3, float *ary_1D, float ***ary_3D)
{
  int i1=0, i2=0, i3=0, j=0;

  for (i1=0; i1<dim1; i1++)
    for (i2=0; i2<dim2; i2++)
      for (i3=0; i3<dim3; i3++)
	ary_3D[i1][i2][i3] = ary_1D[j++];
}

/***********************************************************************************/
/* Function: fortran2c_3D_float_ary                                       @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible float array to 3D C float array. Space for the    */
/* returned ary is automatically allocated.                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  float *ary_1D:     pointer to the 1D fortran compatible array                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  Pointer to the 3D C array. Memory allocation is done automatically.            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

float ***fortran2c_3D_float_ary (int dim1, int dim2, int dim3, float* ary_1D)
{
  int i1=0, i2=0, i3=0, j=0, status=0;
  float ***tmp=NULL;

  status = ASCII_calloc_float_3D (&tmp, dim1, dim2, dim3);
  if (status==0) {
    for (i1=0;i1<dim1;i1++)
      for (i2=0;i2<dim2;i2++)
	for (i3=0;i3<dim3;i3++)
	  tmp[i1][i2][i3] = ary_1D[j++];
  }

  return tmp;
}

/***********************************************************************************/
/* Function: fortran2c_4D_float_ary                                       @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible float array to 4D C float array. Space for the    */
/* returned ary is automatically allocated.                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  float *ary_1D:     pointer to the 1D fortran compatible array                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  Pointer to the 4D C array. Memory allocation is done automatically.            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

float ****fortran2c_4D_float_ary (int dim1, int dim2, int dim3, int dim4, float* ary_1D)
{
  int i1=0, i2=0, i3=0, i4=0, j=0, status=0;
  float ****tmp=NULL;

  status = ASCII_calloc_float_4D (&tmp, dim1, dim2, dim3, dim4);
  if (status==0) {
    for (i1=0; i1<dim1; i1++) 
      for (i2=0; i2<dim2; i2++) 
	for (i3=0; i3<dim3; i3++) 
	  for (i4=0; i4<dim4; i4++) 
	    tmp[i1][i2][i3][i4] = ary_1D[j++];
  }

  return tmp;
}

/***********************************************************************************/
/* Function: fortran2c_4D_double_ary_noalloc                              @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible double array to 4D C double array.                */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  double *ary_1D:    pointer to the 1D fortran compatible array                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  Pointer to the 4D C array. Memory allocation is done automatically.            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

void fortran2c_4D_double_ary_noalloc (int dim1, int dim2, int dim3, int dim4, double* ary_1D,
				      double ****ary_4D)
{
  
  int i1=0, i2=0, i3=0, i4=0, j=0;

  for (i1=0;i1<dim1;i1++)
    for (i2=0;i2<dim2;i2++)
      for (i3=0;i3<dim3;i3++)
	for (i4=0;i4<dim4;i4++) 
	  ary_4D[i1][i2][i3][i4] = ary_1D[j++];
}


/***********************************************************************************/
/* Function: fortran2c_4D_double_ary                                      @32_30i@ */
/* Description:                                                                    */
/* Convert 1D fortran compatible double array to 4D C double array. Space for the  */
/* returned ary is automatically allocated.                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/*  int dim1:          first dimension of C array                                  */
/*  int dim2:          second dimension of C array                                 */
/*  double *ary_1D:    pointer to the 1D fortran compatible array                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  Pointer to the 4D C array. Memory allocation is done automatically.            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                            */
/*                                                                        @i32_30@ */
/***********************************************************************************/

double ****fortran2c_4D_double_ary (int dim1, int dim2, int dim3, int dim4, double* ary_1D)
{
  
  int i1=0, i2=0, i3=0, i4=0, j=0, status=0;
  double ****tmp=NULL;

  status = ASCII_calloc_double_4D (&tmp, dim1, dim2, dim3, dim4);
  if (status==0) {

    /* ??? I do not understand why the mapping between fortran 1-D ary and */
    /* ??? a multidimensional ary should be like this, I would believe the */
    /* ??? order of the loops should be i1, i2, i3, i4 and not as follows. */
    
    for (i1=0;i1<dim1;i1++)
      for (i3=0;i3<dim3;i3++)
	for (i2=0;i2<dim2;i2++)
	  for (i4=0;i4<dim4;i4++)
	    tmp[i1][i2][i3][i4] = ary_1D[j++];
  }

  return tmp;
}

