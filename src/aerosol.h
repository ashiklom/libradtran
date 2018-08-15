/*--------------------------------------------------------------------
 * $Id: aerosol.h 2859 2013-02-25 17:04:24Z robert.buras $
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

#ifndef __aerosol_h
#define __aerosol_h

#include <stdio.h>
#include "uvspec.h"


#define AEROSOL_HAZE_NN            7
#define AEROSOL_HAZE_RURAL         1
#define AEROSOL_HAZE_MARITIME      4
#define AEROSOL_HAZE_URBAN         5
#define AEROSOL_HAZE_TROPOSPHERIC  6

#define AEROSOL_VULCAN_NN          5
#define AEROSOL_VULCAN_BACKGROUND  1
#define AEROSOL_VULCAN_MODERATE    2
#define AEROSOL_VULCAN_HIGH        3
#define AEROSOL_VULCAN_EXTREME     4

#define AEROSOL_SEASON_NN           3
#define AEROSOL_SEASON_SPRINGSUMMER 1
#define AEROSOL_VULCAN_FALLWINTER   2


int setup_aerosol (input_struct input, output_struct *output);

int read_optprop_files (char *filename,
			float *lambda_r, int nlambda_r, 
			int nlambda_rte_lower, int nlambda_rte_upper,
			float *zd, int nlyr, 
			float **dtau, float **ssa, float **gg,
			float ****mom, int **nmom, int verbose);

int calloc_aer_out (aer_out_struct *out, int nlambda, int nlev, int nphamat);
int cp_aer_out (aer_out_struct *source, aer_out_struct *target, int nlambda, int alloc_moment, int quiet);
int free_aer_out (aer_out_struct *out);

#endif



