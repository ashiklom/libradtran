/*--------------------------------------------------------------------
 * $Id: atmosphere.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __atmosphere_h
#define __atmosphere_h

#include <stdio.h>
#include "uvspec.h"

int setup_altitude    (input_struct input, output_struct *output);
int check_hydrostatic_equation ( float *z, float *p, float *T, int nlev, char *profile_name, int quiet, int verbose );
int setup_atmosphere  (input_struct input, output_struct *output);
int setup_zout_levels (input_struct input, output_struct *output, float **z_all, int *nlev_all);
int setup_gases       (input_struct input, output_struct *output);
int setup_temperature (input_struct input, output_struct *output);

int column (float *dens, float *zd, int nlev, float altitude, int interpol_method, float *dens_air, 
            float **dens_avg, int allocate, float *column);

double vapor_pressure (double t);
double vapor_pressure_over_ice (double t);

char* gas_number2string(int gas_number);

int read_molecular_absorption_lambda (char* filename, int quiet,
				      float **lambda, int *nlambda, int *monochromatic);

#define M_H2O 18.015
#define N_MOL_CKDFU 6.0235e23
#define DU2particle 2.687e16   /* DU to column [particles per cm^2] */

#endif



