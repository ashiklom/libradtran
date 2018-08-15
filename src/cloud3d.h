/*--------------------------------------------------------------------
 * $Id: cloud3d.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __cloud3d_h
#define __cloud3d_h

#include <stdio.h>
#include "uvspec.h"

profile *calloc_profile (int n);
profile3D *calloc_profile3D (int Nz, int Nx, int Ny, int *threed);

int setup_caoth3D (input_struct input, output_struct *output);
int setup_atmosphere3D (input_struct input, output_struct *output);

int convert_caoth_mystic ( input_struct        input,
			   caoth_inp_struct    input_caoth,
			   output_struct      *output,
			   caoth_out_struct   *output_caoth,
			   caoth3d_out_struct *caoth3d,
			   int                 iv );

int free_caoth_mystic ( int                 properties,
			caoth3d_out_struct *caoth3d );

#endif



