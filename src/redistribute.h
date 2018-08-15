/*--------------------------------------------------------------------
 * $Id: redistribute.h 3279 2017-07-07 20:35:28Z Claudia.Emde $
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

#ifndef __redistribute_h
#define __redistribute_h

#include <stdio.h>

#define REDISTRIBUTE_FLOAT  1
#define REDISTRIBUTE_DOUBLE 2

int setup_redistribute (input_struct input, output_struct *output);

int regrid (double *xin, double  *yin, int nin, 
	    double *xout, double *yout, int nout,
	    int scale);

int regrid_float (float *xin, float  *yin,  int nin,
		  float *xout, float *yout, int nout,
		  int scale);

int sort_vec (float *xin, int nin, float *xout, int nout, int **i_loc);

int redistribute_optprop (optprop_struct *optprop, int nlambda, 
			  int nlambda_lower, int nlambda_upper,
			  float *zd, int nlyr, int nphamat,
			  float *zd_common, int nlyr_common, int alloc);

int redistribute_1D (void **data, 
		     float *z_old, int nlyr_old, 
		     float *z_new, int nlyr_new,
		     int scale, int type, int alloc);

int redistribute_2D (void ***data, int nlambda, 
		     int nlambda_lower, int nlambda_upper,
		     float *z_old, int nlyr_old, 
		     float *z_new, int nlyr_new,
		     int scale, int type, int alloc);

int redistribute_molecular (atm_optprop_struct *optprop, int nlambda,
			    int nlambda_lower, int nlambda_upper,
			    int *nq,
                            int nipa,
			    float *zd, int nlyr, int Nx, int Ny,
			    float *zd_common, int nlyr_common);


void set_common_z_grid (float *zd, int nlev, float **zd_common, int *nlev_common);
int add_file2grid (float **zd, int *nlev, char *filename);
int add_sigma_nc_file2grid (float **zd, int *nlev, char *filename, 
                            float latitude, float longitude, struct tm UTC, int time_interpolate,
                            float *z_atm, float *press, int nlev_atm,
                            int verbose, int quiet);

#endif



