/*--------------------------------------------------------------------
 * $Id: molecular.h 2627 2011-12-29 08:41:28Z svn-kylling $
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

#ifndef __molecular_h
#define __molecular_h

#include <stdio.h>
#include "uvspec.h"

int setup_crs (input_struct input, output_struct *output);
int setup_raman (input_struct input, output_struct *output);
int setup_rayleigh (input_struct input, output_struct *output);
int closest_above(float lambda, float* lambda_raw, int n_crs); 
int closest_below(float lambda, float* lambda_raw, int n_crs);
int closest_above_double(double lambda, double* lambda_raw, int n_crs); 
int closest_below_double(double lambda, double* lambda_raw, int n_crs);

#endif



