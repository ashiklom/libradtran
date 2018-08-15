/*--------------------------------------------------------------------
 * $Id: molecular3d.h 3308 2017-09-18 13:01:41Z Claudia.Emde $
 *
 * This file is part of libRadtran.
 * Copyright (c) 1997-2017 by Arve Kylling, Bernhard Mayer,
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

#ifndef __molecular3d_h
#define __molecular3d_h

#include <stdio.h>
#include "uvspec.h"

int setup_molecular3d(input_struct input, output_struct *output);

int setup_optprop_molecular3d(input_struct input, output_struct *output);

int optical_properties_molecular3d (input_struct   input,
				    output_struct *output,
				    caoth3d_out_struct *caoth3d,
				    int            iv,
				    int            iq); 



#endif
