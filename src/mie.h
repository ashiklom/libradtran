/*--------------------------------------------------------------------
 * $Id: mie.h 3194 2015-11-25 18:19:46Z svn-kylling $
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

#ifndef __mie_h
#define __mie_h

#include "miecalc.h"

#define OUTPUT_USER_NONE         0
#define OUTPUT_USER_WAVE         1
#define OUTPUT_USER_WAVENUMBER   2
#define OUTPUT_USER_R_EFF        3
#define OUTPUT_USER_REFRAC_REAL  4
#define OUTPUT_USER_REFRAC_IMAG  5
#define OUTPUT_USER_QEXT         6
#define OUTPUT_USER_QSCA         7
#define OUTPUT_USER_OMEGA        8
#define OUTPUT_USER_GG           9
#define OUTPUT_USER_SPIKE       10
#define OUTPUT_USER_PMOM        11
#define OUTPUT_USER_SFORW       12
#define OUTPUT_USER_SBACK       13
#define OUTPUT_USER_QBACK       14

#define OUTPUT_FORMAT_NORMAL     0
#define OUTPUT_FORMAT_CLOUDPROP  1
#define OUTPUT_FORMAT_AERPROP    2
#define OUTPUT_FORMAT_NETCDF     3

#define CK_CRS          0
#define CK_KATO         1
#define CK_KATO2        2
#define CK_KATO2_96     3
#define CK_FU           4
#define CK_AVHRR_KRATZ  5
#define CK_GENERIC      6
#define CK_LOWTRAN      7
#define CK_KATO2ANDWANDJI  12

#define DIST_NONE       0
#define DIST_GAMMA      1
#define DIST_LOGNORMAL  2
#define DIST_FILE       3
#define DIST_AER        4

#define RHO_ICE   0.917
#define RHO_WATER 1.000

#define TYPE_NONE 0
#define TYPE_INSO 1
#define TYPE_WASO 2
#define TYPE_SOOT 3
#define TYPE_SSAM 4
#define TYPE_SSCM 5
#define TYPE_MINM 6
#define TYPE_MIAM 7
#define TYPE_MICM 8
#define TYPE_MITR 9
#define TYPE_SUSO 10


/* DROPLET RADIUS DISTRIBUTION STRUCTURE */

typedef struct {
  double *radius;
  double *number_dens;
  int   n_dens;
} dist_struct;


/* Aerosol size distributions for an aerosol type. It contains 
 size distributions for each relative humidity */
typedef struct {
  int type; 
  int n_hum;
  double *rel_hum; 
  double *r_min; 
  double *r_max;
  double *rmod;
  double *rho;
  double sigma;
dist_struct *size_distr; 
} aerdist_struct;

/* INPUT STRUCTURE */

typedef struct {
  int         verbose; 
  int         program;
  mie_inp_struct  mie;
  mie_wl_inp_struct wl;
  int         lambda_i_start;
  int         lambda_i_end;
  double      temperature; 
  int         medium; 
  int         aerosol_type;
  mie_complex crefin;         /* User-defined refractive index */
  char        *ri_filename;
  double      rho_user;       /* User-defined mass density */ 
  double      r_eff_min;
  double      r_eff_max;
  double      r_eff_step; 
  int         r_eff_log;
  int         distribution;
  float       alpha;
  char        *sd_filename;
  int         n_output_user;
  int         *output_user;
  int         output_format;
  char        *basename;
  char        *data_files_path;
  int         ck_scheme;
  char        *filename_ck_generic;
  int         nstokes;  /*Number of Stokes parameters.*/
  int         nthetamax;
  int         nmom_netcdf;
  double      n_r_max;
  double      dx_max;
  double      accuracy;
} input_struct;


/* OUTPUT STRUCTURE */

typedef struct {
  mie_complex       ref; 
  dist_struct       dist;
  aerdist_struct    aerdistr;
  mie_out_struct    mie;
  mie_wl_out_struct wl;
  double            beta;
  double            omega; 
  double            g;
  double            back;
  int               nmommax;
  float**           pmom;
  int*              ntheta; 
  float**           phase; 
  float**           theta;
  int               n_r_eff;
  double            *r_eff;
  double            rho_medium;
  FILE*             ptr_output;
  int**             warnings_nmom;
  int**             warnings_ntheta;
  char*             outputfile; 
} output_struct;


typedef struct {
  int ncid;
  int id_reff; 
  int id_nreff;
  int id_lam; 
  int id_nlam; 
  int id_theta; 
  int id_nthetamax; 
  int id_ntheta; 
  int id_phase; 
  int id_nphamat;
  int id_pmom; 
  int id_nmom; 
  int id_nmommax; 
  int id_ext; 
  int id_ssa; 
  int id_gg; 
  int id_refre; 
  int id_refim; 
  int id_rho;
  int id_nrho;
} netcdf_ids;
    

#endif
