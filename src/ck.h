/*--------------------------------------------------------------------
 * $Id: ck.h 3276 2017-07-04 14:16:36Z Claudia.Emde $
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

#ifndef __ck_h
#define __ck_h

#include <stdio.h>
#include "uvspec.h"

/* Kato et al. [1999]           */
/* components: CO2, O3, O2, H2O */
#define KATO_COMPONENTS 4

#define KATO_CO2  1
#define KATO_O3   2
#define KATO_O2   3
#define KATO_H2O  4

#define KATO_MAXINT         20
#define KATO_MAXINT_CO2     18
#define KATO_MAXINT_O3      14
#define KATO_MAXINT_O2       6
#define KATO_MAXINT_H2O     11


/* Fu and Liou [1992/93]               */
/* components: CO2, O3, H2O, CH4, N2O; */
/* band overlap handled within ck.c,   */
/* therefore only one component needs  */
/* to be considered here               */    

#define FU_COMPONENTS  1
#define FU_WGHT        1
#define FU_MAXINT     14


/* AVHRR parameterization by David Kratz             */
/* components: CO2, O3, O2, H2O, CH4, N2O, F11, F12; */
/* band overlap handled within ck.c,                 */
/* therefore only one component needs                */
/* to be considered here                             */    

#define AVHRR_KRATZ_COMPONENTS   1
#define AVHRR_KRATZ_WGHT         1
#define AVHRR_KRATZ_MAXINT      20



/* LOWTRAN/SBDART parameterization */

#define LOWTRAN_COMPONENTS   1
#define LOWTRAN_WGHT         1
#define LOWTRAN_MAXINT       3



/* Generic ck                          */
/* band overlap handled internally     */
/* therefore only one component needs  */
/* to be considered here               */    

#define GENERIC_COMPONENTS  1
#define GENERIC_WGHT        1



int kato_readtables (int scheme, ck_struct *crs, char *path, char *photon_filename);
int fu_readtables (ck_struct *crs, char *path, char *photon_filename);
int avhrr_kratz_readtables (ck_struct *ck, char *path, char *photon_filename);
int ck_generic_readtables (ck_struct *ck, char *filename, char *photon_filename);
int sbdart_readtables (ck_struct *ck, int nwvl, char *path, char *photon_filename, int quiet);

int sbdart_profile (float ***temper, float ***press, float *z, 
		    float ***h2o, float ***o3, 
		    float xn2, float xo2, float xco2,
		    float xch4, float xn2o, float xno2,
		    int ck_o4abs, int ck_n2abs, int ck_coabs, int ck_so2abs, int ck_nh3abs, int ck_noabs, int ck_hno3abs,
		    float *dt, float *ssa, int nlev, 
		    float sza, float lambda, int iv, ck_profile *ck, int spectral_is);

int crs_ck (ck_struct *ck, 
	    float *lambda, int n_lambda, 
	    float ***temp, float ***temp_avg, float ***press, float *z, 
            float ****dens, float ****dens_avg, int nlev, int Nx, int Ny,
	    float *mixing_ratio, int h2ocont, 
	    int ck_o4abs, int ck_n2abs, int ck_coabs, int ck_so2abs, 
	    int ck_nh3abs, int ck_noabs, int ck_hno3abs,
	    float *sza,
	    char *path, char *filename, int first, int verbose,
            int nipa,
	    crs_ck_out_struct *crs_ck,
            input_struct input, output_struct *output);

double fu_rayleigh (double pressure, double temperature, int iv);

#endif

