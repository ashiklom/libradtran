/*--------------------------------------------------------------------
 * $Id: cnv.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __cnv_h
#define __cnv_h

#if defined (__cplusplus)
extern "C" {
#endif

#include "numeric.h"


#define SPEC_NOT_EQUIDISTANT   (-10)
#define CONV_NOT_EQUIDISTANT   (-11)
#define CONV_NOT_CENTERED      (-12)
#define SPEC_CONV_DIFFERENT    (-13)


/* prototypes */

int convolute (double *x_spec, double *y_spec, int spec_num,
	       double *x_conv, double *y_conv, int conv_num,
	       double **x_spec_conv, double **y_spec_conv, int *spec_conv_num);

int int_convolute (double *x_spec, double *y_spec, int spec_num,
		   double *x_conv, double *y_conv, int conv_num,
		   double **y_spec_conv, int std);

#if defined (__cplusplus)
}
#endif

#endif





