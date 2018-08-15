/*--------------------------------------------------------------------
 * $Id: atmosphere.c 3305 2017-09-15 13:44:02Z Claudia.Emde $
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "uvspec.h"
#include "ckdfu.h"
#include "cloud.h"
#include "ascii.h"

#include "fortran_and_c.h"
#include "f77-uscore.h"
#define EPSILON 1E-6

#ifndef MAX_LENGTH_OF_LINE
#define MAX_LENGTH_OF_LINE     65536
#endif

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

#define QSORT_CAST (int (*)(const void *, const void *))

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif


/************************************/
/* prototypes of internal functions */
/************************************/

static int read_atmosphere (char *filename, atm_out_struct *out, 
                            /* int interpol_method_press, int interpol_method_temper, int *interpol_method_gas, */
                            int cdf, int quiet, int verbose);

static int extrapolate_atmosphere (float  altitude, atm_out_struct *out,
                                   int quiet, int verbose);

/* static int read_density (char *filename, */
/* 			 atm_out_struct *atm_out, int dens_id, int *interpol_method, */
/* 			 int quiet, int verbose);  */

static int read_radiosonde_file (input_struct input, output_struct *output);

static int read_ECMWF_atmosphere (float lat, float lon, struct tm UTC, int time_interpolate, char *filename, int ECMWF_new_format, 
                                  float ***atm_data, int *rows, float altitude, float **z_model_layer, int *nz_model_layer, 
                                  float *z_cpt_sea, float *press_cpt, float *temper_cpt, int verbose, int quiet);

static int read_ECHAM_atmosphere (float lat, float lon, struct tm UTC, int
                                  time_interpolate, char *filename, 
                                  float ***atm_data, int *rows, float altitude,
                                  float **z_model_layer, int *nz_model_layer,
                                  int verbose, int quiet);

static int find_cold_point_tropopause ( float *press, float *temper, float *z, int nlev, 
                                        float *z_CPT_above_sea, float *p_CPT, float *T_CPT, int quiet, int verbose );

static int read_and_combine_density (char *filename, atm_out_struct *out,
				     int dens_id, int *interpol_method_gas, int interpol_method_temper, 
                                     int unit_profile, float *mol_mass,
				     float *dens_air_all, float *temper_all, float *z_all, int nlev_all,
                                     float *dens_air_atm, float *z_atm, int nlev_atm,
				     int solver, int quiet, int verbose);

static int scale_pressure (float ****pressure, float *****dens, int nlev, float scale_factor, 
                           int *well_mixed_gas, int ialt, int verbose);

static int scale_density_profile (float **dens, float **dens_avg, int alloc_avg, int i,
				  float input_column, float *output_column, 
				  float *scale_factor, int interpol_method_gas, 
				  input_struct input, output_struct *output);


static int convert_profile_units    (float **ptr_dens, float *z_dens, int n, int dens_id,   
                                     float *dens_air_atm, float *temper_atm, float *z_atm, int n_atm,
                                     int interpol_method_air, int interpol_method_temper,
                                     int unit_profile, float *mol_mass, int quiet);

static int calculate_or_read_refind (char *filename,
                                     float ***refind, float *z_all, int nlev_all,
                                     int interpol_method_refind, 
                                     float *t, float *p,
                                     int nlambda, float *lambda,
                                     int quiet);
 
static int combine_profiles (float **dens, float *z, int nlev, int sza,
                             float *dens_all, float *dens_air_all, float *z_all, int nlev_all, 
                             int interpol_method_dens, int interpol_method_air,
                             float ***dens_comb, float **z_comb, int *nlev_comb, int quiet);

static int read_molecular_absorption (char *filename, 
				      float **zd, int *nlyr,
				      int quiet,
				      double **dtau_mono, 
				      float ***dtau_spec, float **lambda, int *nlambda,
				      int *monochromatic, int wl_start_index, int wl_end_index);

static int read_molecular_absorption_z (char* filename, int quiet,
					float **zd, int *nlyr);

static char* interpol_number2string(int method_number);

static int compare_float_inv (void *ap, void *bp);

static int read_dtheta_dx_from_ECMWF_file (input_struct input, output_struct *output);

static int read_dtheta_dy_from_ECMWF_file (input_struct input, output_struct *output);

static int calculate_dtheta_dz            (input_struct input, output_struct *output);

static int  read_ECMWF_ozone_climatology (char *data_directory, float lat, struct tm UTC,
                                          float *dens_air_atm, int nlev_atm, float *press_atm, float **ozone_atm, 
                                          int verbose, int quiet );

static int scale_profile_with_mixing_ratio(input_struct input, output_struct *output, int i_mx, int i_mol, 
					   int ialt, double default_ppm, char *mol_name);


/**************************************************************/
/* Setup model atmosphere.                                    */
/**************************************************************/

int setup_altitude (input_struct input, output_struct *output)
{
  int status=0, cdf=0, ix=0, iy=0, lc=0;
  double  fact = 0;
  float min_z=NOT_DEFINED_FLOAT, max_z=NOT_DEFINED_FLOAT;

  char  filename[FILENAME_MAX];

  /* copy altitude input to output structure*/
  output->alt.source      = input.alt.source;
  output->alt.altitude    = input.alt.altitude;

  /* where to get the atmosphere from */
  if (input.ck_scheme == CK_FILE) {
    strcpy (filename, input.ck_scheme_filename);
    cdf=1;

    if (!input.quiet)
      if (strlen(input.filename[FN_ATMOSPHERE]) > 0)
	fprintf (stderr, " ... ignoring %s, reading atmosphere from %s\n",
		 input.filename[FN_ATMOSPHERE], input.ck_scheme_filename);
  }
  else {
    strcpy (filename, input.filename[FN_ATMOSPHERE]);
    cdf=0;
  }

  /* read atmosphere file */
  status = read_atmosphere (filename, &(output->atm), 
			    /* input.atm.interpol_method_press,  */
			    /* input.atm.interpol_method_temper,  */
			    /* input.atm.interpol_method_gas, */
			    cdf, input.quiet, input.verbose);
  
  if (status!=0) {
    fprintf (stderr, "Error %d reading atmosphere from %s\n",
	     status, filename);
    return status;
  }
  
  output->atm.Nxatm=1; output->atm.Nyatm=1;
  //output->atm.Nzatm=output->atm.nlev-1;
  
  /* 3DAbs moved to molecular3d.c setup_optprop_molecular3d else{ */
  /*   /\* save number of layers from 1D absorption calculation *\/ */
  /*   output->atm.Nzatm=output->atm.nlev-1; */
    
  /*   status = read_atmosphere_3d (input.atmosphere3d_filename, &(output->atm), */
  /* 				 input.atm.mol_mass,  */
  /* 				 input.quiet, input.verbose); */
  /*   if (status!=0) { */
  /*     fprintf (stderr, "Error %d reading 3D atmosphere from %s\n", */
  /* 	       status, input.atmosphere3d_filename); */
  /*     return status; */
  /*   } */
    
  /* } */ 

  /* fprintf (stderr, "setup_altitude 3DAbs Nzatm nlyr %d %d \n",   output->atm.Nzatm, output->atm.nlyr); */
  
  /* Density of O4 as used in the calculations of optical depth is given as below. ak20110404 */
  for (ix=0; ix<output->atm.Nxatm; ix++) {
    for (iy=0; iy<output->atm.Nyatm; iy++) {
      for (lc=0; lc<output->atm.nlev; lc++) {
	fact = output->atm.microphys.dens[MOL_O2][ix][iy][lc]* 1e-23; /* O4 cross section is scaled by 1e46 due */
	/* to float precision limitations         */
	/* account for that here. ak 20110404     */
	output->atm.microphys.dens[MOL_O4][ix][iy][lc] = fact*fact;
      }
    }
  }
  
  /* altitude options only available for 1D atmosphere */
  switch (input.alt.source) {
  case ALT_NOT_DEFINED:
    output->alt.altitude = output->atm.zd[output->atm.nlev-1];
    output->alt.source = ALT_FROM_ATM;
    break;
  case ALT_DIRECT_INPUT:
    output->alt.altitude = input.alt.altitude; /* copy altitude data */
    output->alt.source   = input.alt.source;   /* copy source of altitude data */
    if (output->molecular3d){
      fprintf(stderr, "Error, altitude can not be specified for 3D atmosphere");
      return -1;
    }
    break;
  case ALT_FROM_MAP:
      
    /* read surface elevation from map */
    status = get_number_from_netCDF_map (input.latitude, input.longitude, input.UTC, input.atm.time_interpolate,
					 input.filename[FN_ALTITUDE_MAP], &(output->alt.altitude), TYPE_FLOAT, input.alt.netCDF_alt_name,
					 input.verbose, input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d reading altitude map %s (line %d, function %s in %s)\n",
	       status, input.filename[FN_ALTITUDE_MAP], __LINE__, __func__, __FILE__);
      return status;
    }
      
    output->alt.altitude = output->alt.altitude / (input.alt.scale_factor);
    if (input.verbose)
      fprintf (stderr, "     set altitude to %8.3f km\n", output->alt.altitude);
      
    output->alt.source = input.alt.source;    /* copy source of altitude data */
      
    if (output->molecular3d){
      fprintf(stderr, "Error, altitude_map can not be specified for 3D molecular atmosphere.");
      return -1;
    }
    break;
  default:
    fprintf (stderr, "Error determining surface altitude! (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -1;
  }
    
  

  /* if altitude lower than atmosphere_file */
  /* extrapolate the profiles with constant mixing ratios and constant temperature gradient */
  if (output->alt.altitude < output->atm.zd[output->atm.nlev-1] ) 
    status = extrapolate_atmosphere(output->alt.altitude, &(output->atm),
				    input.quiet, input.verbose);

  if (status!=0) {
    fprintf (stderr, "Error %d returned by extrapolate_atmosphere() (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /* obsolete scince extrapolate atmosphere */
  /*   if (output->alt.altitude < output->atm.zd[output->atm.nlev-1] || output->alt.altitude > output->atm.zd[0]) { */
  /*     fprintf (stderr, "Error -5, altitude = %7.3f km is outside\n", output->alt.altitude); */
  /*     fprintf (stderr, "          of the range of the atmosphere file = [%8.3f km,%8.3f km]\n",  */
  /*                       output->atm.zd[output->atm.nlev-1], output->atm.zd[0]); */
  /*     return -5; */
  /*   } */

  /* read radiosonde                                                     */
  /* (needs surface elevation in order to integrate from there)          */
  /* (must be before setup_clouds, as top of atmosphere is needed there) */
  if (input.atm.rs_source != RS_NO_DATA){
    status = read_radiosonde_file(input, output);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by read_radiosonde_file() (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
      return status;
    }
  }
  

  /*********************************************************/
  /* reduce atmosphere to forced (by user) atmosphere grid */
  /*********************************************************/

  if ( input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER ) {

    /* check if  z_atm_forced_sea are inside the range of the atmosphere */
    if (output->atm.zd[output->atm.nlev-1] < output->atm.zd[0] ) {
      min_z=output->atm.zd[output->atm.nlev-1];
      max_z=output->atm.zd[0];
    }
    else {
      min_z=output->atm.zd[0];
      max_z=output->atm.zd[output->atm.nlev-1];
    }
    /* check only upper and lower boundary */
    for ( lc=0; lc<input.atm.nz_atm_forced_sea; lc=lc+input.atm.nz_atm_forced_sea-1 ) {
      if ( input.atm.z_atm_forced_sea[lc] < min_z || max_z < input.atm.z_atm_forced_sea[lc] ) {
	fprintf (stderr, "Error, atmosphere_zgrid level %f outside the range of the atmosphere file %f %f\n",
		 input.atm.z_atm_forced_sea[lc], min_z, max_z );
	return -1;
      }
    }

    /* new surface altitude */
    output->alt.altitude = input.atm.z_atm_forced_sea[input.atm.nz_atm_forced_sea-1];
    output->alt.source   = ALT_FROM_ATM;
    if (input.verbose)
      fprintf (stderr, " ... forced new altitude = %f\n", output->alt.altitude);

    
    /* interpolate the atmospheric profiles to the forced z-grid */
    status = interpolate_atmosphere (output->atm.zd,
				     &(output->atm.microphys.press),
				     &(output->atm.microphys.temper),
				     &(output->atm.microphys.dens),
				     &(output->atm.microphys.temper_avg),
				     &(output->atm.microphys.dens_avg),
				     output->atm.nlev,
				     input.atm.z_atm_forced_sea, 
				     input.atm.nz_atm_forced_sea,
				     input.atm.interpol_method_press,
				     input.atm.interpol_method_temper,
				     input.atm.interpol_method_gas,
				     input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by interpolate_atmosphere()\n", status);
      return status;
    }

    /* redefine output->atm.zd */
    free(output->atm.zd);
    /* this is important in order to allocate output->atm.zd in set_common_z_grid again */
    output->atm.nlev=0; 
    /* copy z_forced_sea to the grid, where to interpolate */
    set_common_z_grid (input.atm.z_atm_forced_sea, input.atm.nz_atm_forced_sea, &(output->atm.zd), &(output->atm.nlev));

    /* check hydrostatic equation of the interpolated atmosphere */
    status = check_hydrostatic_equation ( output->atm.zd, output->atm.microphys.press[0][0], output->atm.microphys.temper[0][0], output->atm.nlev,
					  "the atmosphere forced to 'atmosphere_zgrid'-levels", input.quiet, input.verbose );
    if (status!=0) {
      fprintf (stderr, "Error %d returned by check_hydrostatic_equation()\n", status);
      return status;
    }
  }
  

  return status;
}


/**************************************************************/
/* Setup model atmosphere.                                    */
/**************************************************************/

int setup_atmosphere (input_struct input, output_struct *output)
{
  int i=0, status=0, lc=0, iv=0, iz=0, ialt=0, ipa=0, iscale=0, is=0, i_gas=0, isp=0;
  float   scale_factor=0;
  float  *zuser     = NULL;
  float  *zd_user   = NULL;

  int     nlyr_user = -999;
  int     nuser     = -999;

  float   kappa = (C_P_DRY_STD-C_V_DRY_STD) / C_P_DRY_STD;  /* == R / cp */;

  char   *method    = NULL;

  float  *z_all     = NULL;
  int     nlev_all  = -999;

  float  *tmp_lambda=NULL, *tmp_dtauf=NULL, *tmp_dtauf_r=NULL;
  float **tmp_dtau=NULL;
  int     tmp_nlambda=0;

  float   tmpdens=0, idealdiff=0;
  int     idealgas=0, firstnonideal=0;
  
  float   mmr_h2o     = 0.0;
  float   c_p_dry_air = 0.0;
  float   c_p_wv      = 0.0;
  float   precip_water = NOT_DEFINED_FLOAT;

  int     monochromatic=0;

  float  *dens_air_atm=NULL;
  float  *z_atm=NULL;
  char   *gas=NULL;
  float   press_h2o = -999;
  int     io=0;
  

  
  /* allocate memory for layer average densities */
  output->atm.microphys.dens_avg = calloc(MOL_NN, sizeof(float *));
  if (output->atm.microphys.dens_avg == NULL) {
    fprintf (stderr, "Error allocating memory for dens_avg\n");
    return -1;
  }
  
  if (input.verbose) {
    /* all method names MUST have the same numer of characters */
    method = (char *) calloc (strlen("spline             ")+1, sizeof (char));
    fprintf (stderr, " ... Interpolation methods\n");
    /* write the interpolate into the verbose file */ 
    method = interpol_number2string(input.atm.interpol_method_press);
    fprintf(stderr, "     z_interpolate pressure          %s interpolation\n", method);
    method = interpol_number2string(input.atm.interpol_method_temper);
    fprintf(stderr, "     z_interpolate temperature       %s interpolation\n", method);
    method = interpol_number2string(input.atm.interpol_method_refind);
    fprintf(stderr, "     z_interpolate refraction index  %s interpolation\n", method);
    for (i=0; i<MOL_NN; i++){
      method = interpol_number2string(input.atm.interpol_method_gas[i]);
      fprintf(stderr, "     z_interpolate %s density      %s interpolation\n", gas_number2string(i), method);
    }
    free(method);
  }
  
  /*******************************************/
  /* read additional atmospheric input files */
  /*******************************************/
  
  /* read ECMWF ozone climatology */
  if ( input.atm.ECMWF_ozone_climatology ) {
    status = read_ECMWF_ozone_climatology (input.filename[FN_PATH], input.latitude, 
					   input.UTC,
					   output->atm.microphys.dens[MOL_AIR][0][0],
					   output->atm.nlev,
					   output->atm.microphys.press[0][0], 
					   &(output->atm.microphys.dens[MOL_O3][0][0]),
					   input.verbose, input.quiet );

    if (status != 0) {
      fprintf (stderr, "Error %d returned by 'read_ECMWF_ozone_climatology' (line %d, function %s in %s) \n", 
	       status, __LINE__, __func__, __FILE__);
      return status;
    }
  }
  
  /* copy the z-grid of the atmosphere file */
  z_atm = (float *) calloc (output->atm.nlev, sizeof(float));
  for (lc=0; lc<output->atm.nlev; lc++)
    z_atm[lc] = output->atm.zd[lc];
  
  /* copy the air number density of the atmosphere file for linmix interpolations*/
  if ( (dens_air_atm =    calloc (output->atm.nlev, sizeof(float))) == NULL ) {
    fprintf (stderr, "Error allocating memory for dens_air_atm\n");
    return -44;
  }
  
  for (lc=0; lc<output->atm.nlev; lc++)
    dens_air_atm[lc] = output->atm.microphys.dens[MOL_AIR][0][0][lc];
  
  /********************************************************************************/
  /* combine all levels from other additional input files to a common grid: z_all */ 
  /********************************************************************************/
  if (input.atm.nz_atm_forced_sea == NOT_DEFINED_INTEGER) {
    
    /* initialisation of the common z-grid: z_all */
    nlev_all = 0; 
    set_common_z_grid (output->atm.zd, output->atm.nlev, 
		       &(z_all), &(nlev_all));
    
    /* dens                  */
    /* refind                */
    /* rayscat               */
    /* molabs                */
    /* zout                  */
    /* altitude 2nd argument */
    /* altitude 1st argument */
    
    
    /* densities (including rh-file) */
    for (i=0; i<MOL_NN; i++)
      if (strlen(input.atm.filename[i]) > 0) {
	status = add_file2grid (&(z_all), &(nlev_all),
				input.atm.filename[i]);
	if (status!=0) {
	  fprintf (stderr, "Error %d returned by add_file2grid()\n", status);
	  return status;
	}
      }
    
    /* refractive index */
    if (strlen(input.filename[FN_REFIND]) > 0) {
      status = add_file2grid (&(z_all), &(nlev_all), 
			      input.filename[FN_REFIND]);
      if (status!=0) {
	fprintf (stderr, "Error %d returned by add_file2grid()\n", status);
	return status;
      }
    }
  }
  else {
    /* user wants a forced z-grid          */
    /* so we do not include further levels */
    nlev_all = 0; 
    set_common_z_grid (input.atm.z_atm_forced_sea, input.atm.nz_atm_forced_sea, &(z_all), &(nlev_all));
  }
  
  /* rayleigh scattering cross section profile, user-defined */
  if (strlen(input.filename[FN_MOL_TAU_SCA]) > 0) {
    status = read_molecular_absorption_z (input.filename[FN_MOL_TAU_SCA], input.quiet,
					  &zd_user, &nlyr_user);
    if (status!=0) {
      fprintf (stderr, "Error reading z levels from user-defined Rayleigh scattering file %s\n",
	       input.filename[FN_MOL_TAU_SCA]);
      fprintf (stderr, "Please check format and retry!\n");
      return -1;
    }

    /* add new levels to the common atmospheric profile */
    set_common_z_grid (zd_user, nlyr_user+1, &z_all, &nlev_all);

    /* free memory */
    free (zd_user);
  }
  
  
  /* molecular absorption cross section profile, user-defined  */
  if (strlen(input.filename[FN_MOL_TAU_ABS]) > 0) {
    status = read_molecular_absorption_z (input.filename[FN_MOL_TAU_ABS], input.quiet,
					  &zd_user, &nlyr_user);
    if (status!=0) {
      fprintf (stderr, "Error reading z levels from user-defined molecular absorption file %s\n", 
	       input.filename[FN_MOL_TAU_ABS]);
      fprintf (stderr, "Please check format and retry!\n");
      return -1;
    }
    
    /* add new levels to the common atmospheric profile */
    set_common_z_grid (zd_user, nlyr_user+1, &z_all, &nlev_all);

    /* free memory */
    free (zd_user);
  }
  
  
  /* SOLVER_POLRADTRAN requires zout levels as atmosphere levels */
  if (input.rte.solver == SOLVER_POLRADTRAN && input.atm.zout_interpolate == NO_ZOUT_INTERPOLATE) {
    input.atm.zout_interpolate = ZOUT_INTERPOLATE;
    if (!input.quiet)
      fprintf (stderr, " ... zout_interpolate is turned on by default for rte_solver POLRADTRAN\n");
  }


  /* add user-defined (ground) altitude   */
  /* as level to the common z-grid        */
  zuser = (float *) calloc (1, sizeof(float));
  zuser[0] = output->alt.altitude;
  set_common_z_grid (zuser, 1, &(z_all), &nlev_all);
  free (zuser);
        
  /* add user-defined equidistant altitude grid */
  if (input.atm.zout_interpolate == ZOUT_INTERPOLATE) {
    if (input.alt.altitude_dz > 0) {
      nuser = (int) ((z_all[0] - z_all[nlev_all-1]) / 
		     input.alt.altitude_dz          + 0.5) + 1;
      
      zuser = (float *) calloc (nuser, sizeof(float));
      
      for (i=0; i<nuser; i++) 
	zuser[i] = z_all[nlev_all-1] + input.alt.altitude_dz * (float) i;
    
      while (zuser[nuser-1] > z_all[0])
	nuser--;

      /* add user-defined grid to atmospheric grid */
      set_common_z_grid (zuser, nuser, &z_all, &nlev_all);

      free (zuser);
    }
  }
  /* else if (input.atm.zout_interpolate == NO_ZOUT_INTERPOLATE) */
  /* altitude second argument heights are added to the z-grid    */
  /* before setup_redistribution                                 */

  /* add zout levels to z_all */
  status = setup_zout_levels( input, output, &z_all, &nlev_all);
  if (status!=0) {
    fprintf (stderr, "Error %d returned by setup_zout_levels() (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }
  

  /************************************************************************/
  /* READING ADDITIONAL INPUT FILES                                       */ 
  /* INTERPOLATION TO THE COMMON Z-GRID: z_all                            */
  /* if input only covers a part of the vertical range of atmosphere_file */
  /* information of the atmosphere_file is assumed in the regions,        */
  /* which are not covered by the additional input file                   */
  /************************************************************************/

  /* first we need pressure, temperature and air density on the new grid z_all for interpolation */
  
  /* Pressure */ 
  status = interpolate_profile (z_atm, &(output->atm.microphys.press[0][0]),
				output->atm.nlev, 
				z_all, nlev_all, 
				input.atm.interpol_method_press, input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating pressure (setup_atmosphere, in atmosphere.c)\n", status);
    return status;
  }
  
  /* Temperature */
  status = interpolate_profile (z_atm, &(output->atm.microphys.temper[0][0]),
				output->atm.nlev, 
				z_all, nlev_all,
				input.atm.interpol_method_temper,
				input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating temperature\n", status);
    return status;
  }
  
  /* replace interpolated cold point temperature, which would smear the extrem of cold temperature */
  /* with double linear interpolated CP-temperature calculated in setup_zout                       */
  if ( output->atm.microphys.temper_cpt != NOT_DEFINED_FLOAT ) {
    switch (input.atm.zout_source) {
    case OUTLEVEL_ZOUT_ABOVE_SUR:
    case OUTLEVEL_ZOUT_ABOVE_SEA:
      for (lc=0; lc<nlev_all; lc++) {
	if ( z_all[lc] == output->atm.microphys.z_cpt_sea ) {
	  output->atm.microphys.temper[lc][0][0] = output->atm.microphys.temper_cpt;
	  break;
	}
      }
      break;
    case OUTLEVEL_PRESS:
      for (lc=0; lc<nlev_all; lc++) {
	if ( fabs(output->atm.microphys.press[lc][0][0] - output->atm.microphys.press_cpt) < 1.e-6*output->atm.microphys.press_cpt ) {
	  output->atm.microphys.temper[lc][0][0] = output->atm.microphys.temper_cpt;
	  break;
	}
      }
      break;
    case OUTLEVEL_ATM_LEVELS:
    case OUTLEVEL_ALL_LEVELS:
    case OUTLEVEL_MODEL_LEVELS:
    case OUTLEVEL_MODEL_LAYERS:
    case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
      break;
    default:
      fprintf (stderr, "Error determining output levels zout.\n");
      fprintf (stderr, "       unknown option %d for zout_source! (line %d, function '%s' in '%s')\n", input.atm.zout_source, __LINE__, __func__, __FILE__);
      return -1;
    }
  }
  
  
  /* Air */
  /* interpolation of air dens is not 100% consistent with interpolation of p and T and ideal gas law */
  /* but otherwise all informations of air density given by the user will be thrown away              */
  status = interpolate_profile (z_atm, &(output->atm.microphys.dens[MOL_AIR][0][0]), 
				output->atm.nlev,
				z_all, nlev_all, 
				input.atm.interpol_method_gas[MOL_AIR], input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating dens_air\n", status);
    return status;
  }
  
  /* correct the cold point tropopause */
  if (output->atm.microphys.z_cpt_sea != NOT_DEFINED_FLOAT)
    for (lc=0; lc<nlev_all; lc++)
      if ( z_all[lc] - output->atm.microphys.z_cpt_sea == 0.0 ) {

	output->atm.microphys.dens[MOL_AIR][lc][0][0] =
	  output->atm.microphys.press[lc][0][0] /
	  (BOLTZMANN * output->atm.microphys.temper[lc][0][0]) * 100.0 / 1e6;
	if ( input.verbose )
	  fprintf (stderr, " ... correct air density at the cold point tropopause in %f km\n", output->atm.microphys.z_cpt_sea);
	break;
      }
  
  /*   /\* determine air number density from pressure and temperature *\/ */
  /*   free(output->atm.microphys.dens[MOL_AIR]); */
  /*   output->atm.microphys.dens[MOL_AIR] = calloc (nlev_all, sizeof(float)); */
  /*   for (lc=0; lc<nlev_all; lc++) */
  /*     output->atm.microphys.dens[MOL_AIR][lc] =  output->atm.microphys.press[lc] /  */

  /*                                                (BOLTZMANN * output->atm.microphys.temper[lc]) * 100.0 / 1e6; */
  
  
  /* calculate layer average temperature  */ 
  /* this is here, as average_dens() needs a suitable dens_air */
  ASCII_calloc_float_3D( &output->atm.microphys.temper_avg, 1, 1, nlev_all);
  status = average_dens(output->atm.microphys.temper[0][0],
			output->atm.microphys.dens[MOL_AIR][0][0], 
			z_all, nlev_all,
			input.atm.interpol_method_temper, 
			&(output->atm.microphys.temper_avg[0][0]), NO);
  
  if (status!=0) {
    fprintf (stderr, "Error %d calculating average temperature (line %d, function '%s' in '%s') \n", 
	     status, __LINE__, __func__, __FILE__ );
    return status;
  }
  
  /* now read and combine additional data */

  /* pressure    */
  /* temperature */
  /* dens (dens_files might also be dens_tab files) */
  /* refind      */

  /*****************************************/
  /* density profiles (including rh-file)  */
  /*****************************************/
  
  for (i=0; i<MOL_NN; i++) {
    
    if (i != MOL_AIR ) {
      /* call it in any case as densities are also interpolated onto the z_all grid */ 
      /* it also calculates the otherwise empty rh profile */
      /* even if it doesn't read rh from the file          */
      status = read_and_combine_density (input.atm.filename[i], &output->atm, 
					 i, input.atm.interpol_method_gas,
					 input.atm.interpol_method_temper, 
					 input.atm.unit_profile[i], input.atm.mol_mass,
					 output->atm.microphys.dens[MOL_AIR][0][0],
					 output->atm.microphys.temper[0][0],
					 z_all, nlev_all,
					 dens_air_atm, z_atm, output->atm.nlev,
					 input.rte.solver, input.quiet, input.verbose);
      if (status !=0) {
	fprintf (stderr, "Error %d reading density profile (line %d, function '%s' in '%s')\n", status, __LINE__, __func__, __FILE__);
	return status;
      }
    }
  }
  free(dens_air_atm);
  free(z_atm);
  
  /* allocate memory for scale factors */
  output->scale_factor = calloc(MOL_NN, sizeof(float));
  if (output->scale_factor == NULL) {
    fprintf (stderr, "Error allocating memory for 'scale_factor' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -1;
  }
  
  /* allocate memory for vertical columns */  
  output->column = calloc(MOL_NN, sizeof(float));
  if (output->column == NULL) {
    fprintf (stderr, "Error allocating memory for 'column' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -1;
  }
    
  if (output->atm.microphys.nsza_denstab > 0) {
    /* allocate memory for scale factors */
    
    status  = ASCII_calloc_float (&output->scale_factor_denstab,
				  MOL_NN, output->atm.microphys.nsza_denstab);
    if (status != 0) {
      fprintf (stderr, "Error allocating memory for 'scale_factor_denstab' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return status;
    }
    
    /* allocate memory for vertical columns */  
    status  = ASCII_calloc_float (&output->column_denstab, MOL_NN, output->atm.microphys.nsza_denstab);
    if (status != 0) {
      fprintf (stderr, "Error allocating memory for 'column_denstab' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return status;
    }
  }
  
  
  /* the following 20 lines may have become obsolete after deleting saturation */
  
  /* calculate h2o column */
  if ((ASCII_calloc_float_3D(&output->atm.microphys.dens_avg[MOL_H2O], 1, 1, output->atm.nlev))
      !=0) {
    fprintf (stderr, "Error: Allocation of memory for dens_avg[MOL_H2O] failed (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -1;
  }
  
  status = column (output->atm.microphys.dens[MOL_H2O][0][0],
		   output->atm.zd, output->atm.nlev, output->alt.altitude,
		   input.atm.interpol_method_gas[MOL_H2O], 
		   output->atm.microphys.dens[MOL_AIR][0][0],
		   &(output->atm.microphys.dens_avg[MOL_H2O][0][0]),
		   NO, &(output->column[MOL_H2O]));
  if (status != 0) {
    fprintf (stderr, "Error %d returned by function 'column' (line %d, function '%s' in '%s')\n", status, __LINE__, __func__, __FILE__);
    return status;
  }
  
  free(output->atm.microphys.dens_avg[MOL_H2O][0][0]);
  /* this is alloacted again, when dens columns are scaled */
  
  /*********************/
  /* refraction index  */
  /*********************/

  /* calculate_or_read_refind() needs to be called in any case */
  /* because it calculates the otherwise empty refind profile  */
  /* even if it doesn't read refind from the file              */ 
  
  status = calculate_or_read_refind (input.filename[FN_REFIND],
				     &(output->atm.microphys.refind),
				     z_all, nlev_all,
				     input.atm.interpol_method_refind, 
				     output->atm.microphys.temper[0][0],
				     output->atm.microphys.press[0][0],
				     output->wl.nlambda_r, output->wl.lambda_r,
				     input.quiet);
  
  if (status!=0) {
    fprintf (stderr, "Error %d returned by function 'calculate_or_read_refind' (line %d, function '%s' in '%s')\n", status, __LINE__, __func__, __FILE__);
    return status;
  }
  
  /***********************************************/
  /*   copy altitude grid to final destination   */
  /***********************************************/

  free (output->atm.zd);
  output->atm.zd = calloc (nlev_all, sizeof(float));
  
  for (lc=0; lc<nlev_all; lc++)
    output->atm.zd[lc] = z_all[lc];
  
  output->atm.nlev = nlev_all;
  output->atm.nlyr = nlev_all-1; 
  
  free (z_all); 
  
  /*****************************************************************************************/
  /*                               Scaling profiles                                        */
  /*****************************************************************************************/ 

  if (input.verbose)
    fprintf (stderr, "*** Scaling profiles\n");
  
  /* Test if surface altitude is really in the zd-grid */
  for (ialt=0; ialt<output->atm.nlev; ialt++) 
    if (output->alt.altitude == output->atm.zd[ialt]) 
      break;
  
  if (ialt == output->atm.nlev) {
    fprintf (stderr, "Error, user-defined altitude %f not found in altitude grid (line %d, function '%s' in '%s')\n", 
	     output->alt.altitude, __LINE__, __func__, __FILE__);
    return -1;
  }
  
  /* allocate memory for layer average densities */
  for (i=0; i<MOL_NN; i++) {
    status = 
      ASCII_calloc_float_3D(&output->atm.microphys.dens_avg[i], 1, 1, output->atm.nlyr);
    if (status!=0) {
      fprintf (stderr, "Error allocating memory for dens_avg[%3d] (line %d, function '%s' in '%s')\n", i, __LINE__, __func__, __FILE__);
      return -1;
    }
  }
  
  /*****************************************************************************************/
  /*                      Scaling profiles with COLUMN VALUEs                              */
  /*                and calculating average densities in the layers                        */
  /*****************************************************************************************/ 

  /* BrO, HCHO, OCLO are not included in normal calculations */
  /* but might be added by dens- or denstab files            */
  for (i=0; i<MOL_NN; i++) {
    
    gas = gas_number2string(i);
    
    switch (i) {
      /************************************/
      /* Scale pressure (and air density) */
      /************************************/
    case(MOL_AIR):
      if (input.pressure >= 0.0) {
	
	scale_factor = input.pressure / output->atm.microphys.press[0][0][ialt];
	scale_pressure (&(output->atm.microphys.press),
			&(output->atm.microphys.dens),output->atm.nlev, 
			scale_factor, input.atm.well_mixed_gas, ialt, input.verbose);
	
      }

      /* call to column () is essential because it also includes calculation of dens_avg[MOL_AIR] */
      status = column (output->atm.microphys.dens[i][0][0], output->atm.zd, 
		       output->atm.nlev, output->alt.altitude, 
		       input.atm.interpol_method_gas[i], 
		       output->atm.microphys.dens[MOL_AIR][0][0],
		       &(output->atm.microphys.dens_avg[i][0][0]),
		       NO, &(output->column[i]));  /* in atmosphere.c */

      if (status!=0) {
	fprintf (stderr, "\nError %d returned by column()\n", status);
	fprintf (stderr, "      during calculation of %s column (line %d, function '%s' in '%s')\n", gas, __LINE__, __func__, __FILE__);
	return status;
      }
      break;
      
      /********************/
      /* Scale gas column */
      /********************/
    case(MOL_O3):
    case(MOL_O2):
      /* H2O is a special case, see below */
    case(MOL_CO2):
    case(MOL_NO2):
    case(MOL_BRO):
    case(MOL_OCLO):
    case(MOL_HCHO):
    case(MOL_O4):
    case(MOL_SO2):
    case(MOL_CH4):
    case(MOL_N2O):
    case(MOL_CO):
    case(MOL_N2):
      
      status = 
	ASCII_calloc_float_3D(&output->atm.microphys.dens_avg[i], 1, 1, output->atm.nlyr);
      if (status!=0) {
	fprintf (stderr, "Error allocating memory for dens_avg[%3d] (line %d, function '%s' in '%s')\n", i, __LINE__, __func__, __FILE__);
	return -1;
      }
      
      status = scale_density_profile (&(output->atm.microphys.dens[i][0][0]), 
				      &(output->atm.microphys.dens_avg[i][0][0]), NO, i,
				      input.atm.column[i], &(output->column[i]),
				      &(output->scale_factor[i]), input.atm.interpol_method_gas[i],
				      input, output);   /* in atmosphere.c */
      
      if (status!=0) {
	fprintf (stderr, "Error %d returned by scale_density_profile()\n", status);
	fprintf (stderr, "      while scaling dens (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
	return status;
      }
      
      /* air mass factor profiles */
      if (i == output->atm.microphys.denstab_id) {
	
	for (is=0; is<output->atm.microphys.nsza_denstab; is++) {
	  
	  status = scale_density_profile (&(output->atm.microphys.denstab[i][is]), 
					  &(output->atm.microphys.denstab_avg[i][is]), 
					  NO, i,
					  input.atm.column[i], 
					  &(output->column_denstab[i][is]),
					  &(output->scale_factor_denstab[i][is]), 
					  input.atm.interpol_method_gas[i],
					  input, output);   /* in atmosphere.c */
	  
	  if (status!=0) {
	    fprintf (stderr, "Error %d returned by scale_density_profile()\n", status);
	    fprintf (stderr, "      while scaling denstab (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
	    return status;
	  }
	  
	  for (lc=0; lc<output->atm.nlev; lc++)
	    output->atm.microphys.denstab_amf[is][lc] *= output->scale_factor_denstab[i][is];
	  
	} /* for (is=0; is<output->atm.microphys.nsza_denstab; is++) { */
      }
      
      break;
      
    case MOL_H2O:
      /****************************/
      /* Scale precipitable water */
      /****************************/
      if (strcmp(" -1 ",gas)==0) {
	fprintf (stderr, "Error, unknown gas id number %d (line %d, function '%s' in '%s')\n", i, __LINE__, __func__, __FILE__);
	return -1;
      }
      
      status = column (output->atm.microphys.dens[i][0][0],
		       output->atm.zd, output->atm.nlev, output->alt.altitude,
		       input.atm.interpol_method_gas[i],
		       output->atm.microphys.dens[MOL_AIR][0][0],
		       &(output->atm.microphys.dens_avg[i][0][0]),YES,&(output->column[i]));  /* in atmosphere.c */
      if (status!=0) {
	fprintf (stderr, "\nError %d returned by column()\n", status);
	fprintf (stderr, "      during calculation of %s column (line %d, function '%s' in '%s')\n", gas, __LINE__, __func__, __FILE__);
	return status;
      }

      output->precip_water = output->column[i] * DU2particle * (M_H2O / N_MOL_CKDFU * 10.);
      /*     DU         *  DU->cm-2   *       cm-2->kg/m2         */
      
      /* alternative input option dens_column, convert to kg/m2 */
      if (input.atm.column[MOL_H2O] >= 0.0) {
	switch (input.atm.unit_column[MOL_H2O]) {
	case MOL_UNIT_DU:
	  precip_water = input.atm.column[MOL_H2O] * DU2particle * (M_H2O / N_MOL_CKDFU * 10.);
	  /*         DU             *  DU->cm-2   *       cm-2->kg/m2         */
	  break;
	case MOL_UNIT_CM_2:
	  precip_water = input.atm.column[MOL_H2O] * (M_H2O / N_MOL_CKDFU * 10.);
	  /*        cm-2            *       cm-2->kg/m2         */ 
	  break;
	case MOL_UNIT_MM:
	  precip_water = input.atm.column[MOL_H2O];
	  /*   kg/m2              */
	  break;
	default:
	  fprintf (stderr, "\nError, unknown unit for %s column (line %d, function '%s' in '%s')\n", gas, __LINE__, __func__, __FILE__);
	  return -1;
	}
      }
      
      if (precip_water >= 0.0) {
	
	/* verbose output */
	if (!input.quiet)
	  fprintf (stderr, " ... scaling precipitable water from %10.6f kg/m2 to ", 
		   output->precip_water);
	
	/* determine scale factor and scale profile */
	/* scaling profile for right column value */
	/* problem: changing level properties, but you want sum of layer properties */
	/* and a variable integration/interpolation method has to be considered */
	
	/* doing this with an iteration, until column value is right, 20 iterations should be more than enough */
	for (iscale=0;iscale<20;iscale++) { 
	  
	  /* determine scale factor */
	  if (precip_water == 0) {
	    output->scale_factor[i] = 0;
	  }
	  else if (output->precip_water == 0.0) {
	    fprintf (stderr, "\nError, can't scale a H2O profile which is == 0, column = 0 kg/m2\n");
	    return -1;
	  }
	  else
	    /* (also valid for no clouds at all) */
	    output->scale_factor[i] = precip_water/output->precip_water;
	  
	  /* scale profile */
	  for (lc = 0; lc<output->atm.nlev; ++lc) {
	    
	    /* /\* for additional verbose output *\/ */
	    /* old_h2o = output->atm.microphys.dens[MOL_H2O][lc]; */
	    
	    output->atm.microphys.dens[i][0][0][lc] *= output->scale_factor[i];
	    
	    /* additional verbose output */
	    /* fprintf (stderr, "%5d  %6.2f  %.5e ", lc, output->atm.zd[lc], old_h2o); */
	    /* fprintf (stderr, " %.5e  %8.4f\n", output->atm.microphys.dens[MOL_H2O][lc],  */
	    /*                                   output->atm.microphys.dens[i][lc]/old_h2o); */
	  } /* end of: for (lc = 0; lc<output->atm.nlev; ++lc) { ... */
	  
	    /* recalculate the h2o column AND dens_avg */
	  status = column (output->atm.microphys.dens[i][0][0],
			   output->atm.zd, output->atm.nlev, output->alt.altitude, 
			   input.atm.interpol_method_gas[i], 
			   output->atm.microphys.dens[MOL_AIR][0][0],
			   &(output->atm.microphys.dens_avg[i][0][0]),
			   NO,&(output->column[i]));
	  
	  if (status!=0) {
	    fprintf (stderr, "Error %d calculating %s column!\n", status, gas);
	    return status;
	  }
	  output->precip_water = output->column[i] * 
	    DU2particle * (M_H2O / N_MOL_CKDFU * 10.);   
	  /*     DU         *  DU->cm-2   *       cm-2->kg/m2         */
	  
	  if (fabs(output->precip_water-precip_water) < 0.00001)
	    break;
	} /* end of: for (iscale=0;iscale<10;iscale++) { ... */
	
	  /* recalculate the h2o column */
	status = column (output->atm.microphys.dens[i][0][0], output->atm.zd,
			 output->atm.nlev, output->alt.altitude, 
			 input.atm.interpol_method_gas[i], 
			 output->atm.microphys.dens[MOL_AIR][0][0],
			 &(output->atm.microphys.dens_avg[i][0][0]),
			 NO,&(output->column[i]));
	if (status!=0) {
	  fprintf (stderr, "Error %d returned by column()\n", status);
	  fprintf (stderr, "      during calculation of %s column (line %d, function '%s' in '%s')\n", gas, __LINE__, __func__, __FILE__);
	  return status;
	}
	output->precip_water = output->column[i] * DU2particle * (M_H2O / N_MOL_CKDFU * 10.);   
	/*     DU         *  DU->cm-2   *       cm-2->kg/m2         */
	
	if (fabs(output->precip_water-precip_water) >= 0.0001) {
	  fprintf (stderr, "Error, iteration for scaling water column (in setup_atmosphere) to %10.6f kg/m2 failed.\n",precip_water);
	  fprintf (stderr, "Reached end of iteration: %10.6f kg/m2\n",output->precip_water);
	  return -1;
	}
	
	/* verbose output */
	if (!input.quiet)
	  fprintf (stderr, "%10.6f kg/m2\n", output->precip_water);
	
      } /* end of: if (input.precip_water >= 0.0) {... */
      
      break;      
    default:
      fprintf (stderr, "Error, unknown gas_number %d when scaling gas profile (line %d, function '%s' in '%s')\n",
	       i, __LINE__, __func__, __FILE__);
      return -1;
    }
    free(gas);
  }
  

  /*****************************************************************************************/
  /*                      Scaling profiles with MIXING RATIOs                              */
  /*                         at the user-defined altitude                                  */
  /*****************************************************************************************/ 
  
  for (ialt=0; ialt<output->atm.nlev; ialt++)
    if (output->alt.altitude == output->atm.zd[ialt])
      break;
  
  if (ialt == output->atm.nlev) {
    fprintf (stderr, "Error, user-defined altitude %f not found in altitude grid (line %d, function '%s' in '%s')\n", 
	     output->alt.altitude, __LINE__, __func__, __FILE__);
    return -1;
  }
  
  output->mixing_ratio = calloc (MX_NN, sizeof(float));
  
  status =scale_profile_with_mixing_ratio(input,output,MX_O2 ,MOL_O2 ,ialt,-1.0    ,"O2" );
  status+=scale_profile_with_mixing_ratio(input,output,MX_H2O,MOL_H2O,ialt,-1.0    ,"H2O");
  status+=scale_profile_with_mixing_ratio(input,output,MX_CO2,MOL_CO2,ialt,-1.0    ,"CO2");
  status+=scale_profile_with_mixing_ratio(input,output,MX_NO2,MOL_NO2,ialt,-1.0    ,"NO2");
  
  if(output->wl.use_reptran){
    status+=scale_profile_with_mixing_ratio(input,output,MX_CH4,MOL_CH4,ialt,-1.0    ,"CH4");
    status+=scale_profile_with_mixing_ratio(input,output,MX_N2O,MOL_N2O,ialt,-1.0    ,"N2O");
  }
  else{
    status+=scale_profile_with_mixing_ratio(input,output,MX_CH4,-1     ,ialt,1.6     ,"CH4");
    status+=scale_profile_with_mixing_ratio(input,output,MX_N2O,-1     ,ialt,0.28    ,"N2O");
  }
  
  status+=scale_profile_with_mixing_ratio(input,output,MX_F11,-1     ,ialt,0.268E-3,"F11");
  status+=scale_profile_with_mixing_ratio(input,output,MX_F12,-1     ,ialt,0.503E-3,"F12");
  status+=scale_profile_with_mixing_ratio(input,output,MX_F22,-1     ,ialt,0.105E-3,"F22");
  
  if (status){
    fprintf (stderr, "Error during scale_profile_with_mixing_ratio()!\n");
    return -1;
  }
  
  /* some warnings */
  if(!input.quiet && output->wl.use_reptran) {
    if (input.mixing_ratio[MX_F11]>=0)
      fprintf (stderr, "Warning: 'mixing_ratio f11' is ignored!\n");
    if (input.mixing_ratio[MX_F12]>=0)
      fprintf (stderr, "Warning: 'mixing_ratio f12' is ignored!\n");
    if (input.mixing_ratio[MX_F22]>=0)
      fprintf (stderr, "Warning: 'mixing_ratio f22' is ignored!\n");
  }
  
  /********************************************/
  /* Read Rayleigh cross section if specified */
  /********************************************/
  if (strlen(input.filename[FN_MOL_TAU_SCA]) > 0) {
    status = read_molecular_absorption (input.filename[FN_MOL_TAU_SCA],
					&zd_user, &nlyr_user,
					input.quiet,
					&(output->atm.optprop.tau_rayleigh_user), 
					&tmp_dtau, &tmp_lambda, &tmp_nlambda,
					&monochromatic, output->wl.nlambda_rte_lower,
					output->wl.nlambda_rte_upper);

    if (status!=0) {
      fprintf (stderr, "Error %d returned by read_molecular_absorption (\"%s\")\n", 
	       status, input.filename[FN_MOL_TAU_SCA]);
      return status;
    }
    /* FIXCE should be no longer necessary here ... */
    /*if(input.rte.mc.spectral_is){ */
    /*       /\* MYSTIC spectral calculation should be done for all wavelengths in the {molecular,rayleigh}_tau_file. *\/ */
    /*       output->mc.alis.nlambda_abs=tmp_nlambda; */
    /*       output->mc.alis.lambda = calloc(output->wl.nlambda_r, sizeof(float)); */
    
    /*       for (iv=0; iv<output->wl.nlambda_r; iv++){ */
    /*         output->mc.alis.lambda[iv]=tmp_lambda[iv]; */
    /*       } */
    /*     } */
    
    if (monochromatic) {
      output->atm.rayleigh = RAYLEIGH_FILE_MONO;

      if (!input.quiet)
	fprintf (stderr, " ... read monochromatic molecular absorption from %s\n", 
		 input.filename[FN_MOL_TAU_SCA]);
    }
    else {
      output->atm.rayleigh = RAYLEIGH_FILE_SPEC;
      
      if (!input.quiet)
	fprintf (stderr, " ... read Rayleigh scattering cross section from %s\n", 
		 input.filename[FN_MOL_TAU_SCA]);
      
      /*********************************************************************/
      /* interpolate optical thickness to internal wavelength grid, linear */
      /*********************************************************************/

      /* allocate memory for optical thickness on internal wavelength grid */
      output->atm.optprop.tau_rayleigh_user_r = calloc (output->wl.nlambda_r, sizeof(double *));
      for (iv=0; iv<output->wl.nlambda_r; iv++)
	output->atm.optprop.tau_rayleigh_user_r[iv] = calloc (nlyr_user, sizeof(double));
      
      switch (output->wl.type) {
      case WLGRID_CK:
      case WLGRID_USER:
      case WLGRID_BANDS:
      case WLGRID_UVSPEC:
      case WLGRID_MOLABS:
	
	/* memory for temporary arrays */
	tmp_dtauf   = calloc (tmp_nlambda,          sizeof (float));
	tmp_dtauf_r = calloc (output->wl.nlambda_r, sizeof (float));
	
	/* loop over layers */
	for (lc=0; lc<nlyr_user; lc++) {
	  
	  /* interpolate tmp_dtauf to internal wavelength grid */

	  /* copy double array to temporary float array */
	  /* ??? should improve that by providing a ??? */
	  /* ??? an arb_wvn for arrays of double    ??? */
	  for (iv=0; iv<tmp_nlambda; iv++)
	    tmp_dtauf[iv] = tmp_dtau[lc][iv];
	  
	  status = arb_wvn (tmp_nlambda, tmp_lambda, tmp_dtauf, 
			    output->wl.nlambda_r, output->wl.lambda_r, 
			    tmp_dtauf_r,
			    1, 0);
          if (status!=0) {
	    fprintf (stderr, "Error %d during interpolation of dtau (line %d, function '%s' in '%s')\n", 
		     status, __LINE__, __func__, __FILE__);
	    return status;
	  }
	  
	  /* copy temporary float array to final destination */
	  for (iv=0; iv<output->wl.nlambda_r; iv++)
	    output->atm.optprop.tau_rayleigh_user_r[iv][lc] = tmp_dtauf_r[iv];
	  
	  if (status != 0) {
	    fprintf (stderr, "Error %d interpolating molecular absorption grid to internal\n", status);
	    fprintf (stderr, "wavelength grid (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
	    return status;
	  }
	}

	free(tmp_dtauf); free(tmp_dtauf_r);
	break;
	  
      default:
	fprintf (stderr, "Error, unknown wavelength grid type %d (line %d, function '%s' in '%s')\n", 
		 output->wl.type, __LINE__, __func__, __FILE__);
	return -1;
      }

      free(tmp_lambda); 
      for (lc=0; lc<nlyr_user; lc++)
	free (tmp_dtau[lc]);
      free (tmp_dtau);
    }


    /* redistribute tau_rayleigh_user and tau_rayleigh_user_r to common atmospheric grid */
    status = 0;
    if (monochromatic)
      status = redistribute_1D ((void *) (&(output->atm.optprop.tau_rayleigh_user)),
				zd_user, nlyr_user, 
				output->atm.zd, output->atm.nlev-1, 
				1, REDISTRIBUTE_DOUBLE, 1);
    else 
      status = redistribute_2D ((void *) (&(output->atm.optprop.tau_rayleigh_user_r)), 
				output->wl.nlambda_r, 
				output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
				zd_user, nlyr_user, 
				output->atm.zd, output->atm.nlev-1,
				1, REDISTRIBUTE_DOUBLE, 1);
    
    if (status!=0) {
      fprintf (stderr, "Error interpolating user-defined molecular absorption coefficients\n");
      fprintf (stderr, "to common atmospheric vertical grid (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return status;
    }
  }			

  /* Read molecular absorption cross section if specified */ 
  if (strlen(input.filename[FN_MOL_TAU_ABS]) > 0) {
    status = read_molecular_absorption (input.filename[FN_MOL_TAU_ABS], 
					&zd_user, &nlyr_user,
					input.quiet,
					&(output->atm.optprop.tau_molabs_user), 
					&tmp_dtau, &tmp_lambda, &tmp_nlambda,
					&monochromatic, output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper);

    if (status!=0) {
      fprintf (stderr, "Error %d returned by read_molecular_absorption (\"%s\")\n", 
	       status, input.filename[FN_MOL_TAU_ABS]);
      return status;
    }
    
    if (monochromatic) {
      output->atm.molabs = MOLABS_FILE_MONO;
      if (!input.quiet)
        fprintf (stderr, " ... read monochromatic molecular absorption from %s\n", 
                 input.filename[FN_MOL_TAU_ABS]);
    }
    else {
      output->atm.molabs = MOLABS_FILE_SPEC;
        
      if (!input.quiet)
        fprintf (stderr, " ... read spectral molecular absorption from %s\n", 
                 input.filename[FN_MOL_TAU_ABS]);
        
      /*********************************************************************/
      /* interpolate optical thickness to internal wavelength grid, linear */
      /*********************************************************************/

      /* allocate memory for optical thickness on internal wavelength grid */
      output->atm.optprop.tau_molabs_user_r = calloc (output->wl.nlambda_r, sizeof(double *));
      for (iv=0; iv<output->wl.nlambda_r; iv++)
        output->atm.optprop.tau_molabs_user_r[iv] = calloc (nlyr_user, sizeof(double));
        
      switch (output->wl.type) {
      case WLGRID_CK:
      case WLGRID_USER:
      case WLGRID_BANDS:
      case WLGRID_UVSPEC:
          
        /* memory for temporary arrays */
        tmp_dtauf   = calloc (tmp_nlambda,          sizeof (float));
        tmp_dtauf_r = calloc (output->wl.nlambda_r, sizeof (float));
          
        /* loop over layers */
        for (lc=0; lc<nlyr_user; lc++) {
            
          /* interpolate tmp_dtauf to internal wavelength grid */
            
          /* copy double array to temporary float array */
          /* ??? should improve that by providing a ??? */
          /* ??? an arb_wvn for arrays of double    ??? */
          for (iv=0; iv<tmp_nlambda; iv++)
            tmp_dtauf[iv] = tmp_dtau[lc][iv];
            
          status = arb_wvn (tmp_nlambda, tmp_lambda, tmp_dtauf, 
                            output->wl.nlambda_r, output->wl.lambda_r, 
                            tmp_dtauf_r,
                            1, 0);
            
          /* copy temporary float array to final destination */
          for (iv=0; iv<output->wl.nlambda_r; iv++)
            output->atm.optprop.tau_molabs_user_r[iv][lc] = tmp_dtauf_r[iv];
            
            
          if (status != 0) {
            fprintf (stderr, "Error %d interpolating molecular absorption grid to internal\n", status);
            fprintf (stderr, "wavelength grid\n");
            return status;
          }
        }

        free(tmp_dtauf); free(tmp_dtauf_r);
        break;
	  
      case WLGRID_MOLABS:
        /* no need to interpolate, just copy */
        /* copy temporary float array to final destination */
        for (lc=0; lc<nlyr_user; lc++)
          for (iv=0; iv<output->wl.nlambda_r; iv++)
            output->atm.optprop.tau_molabs_user_r[iv][lc] = tmp_dtau[lc][iv];
          
        break;
	
      default:
        fprintf (stderr, "Error, unknown wavelength grid type %d\n", 
                 output->wl.type);
        return -1;
      }
        
      free(tmp_lambda); 
      for (lc=0; lc<nlyr_user; lc++)
        free (tmp_dtau[lc]);
      free (tmp_dtau);
    }
      

    /* redistribute tau_molabs_user and tau_molabs_user_r to common atmospheric grid */
    status = 0;
    if (monochromatic) 
      status = redistribute_1D ((void *) (&(output->atm.optprop.tau_molabs_user)),
                                zd_user, nlyr_user, 
                                output->atm.zd, output->atm.nlev-1, 
                                1, REDISTRIBUTE_DOUBLE, 1);
    else
      status = redistribute_2D ((void *) (&(output->atm.optprop.tau_molabs_user_r)), 
                                output->wl.nlambda_r, 
                                output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
                                zd_user, nlyr_user, 
                                output->atm.zd, output->atm.nlev-1,
                                1, REDISTRIBUTE_DOUBLE, 1);
      
    if (status!=0) {
      fprintf (stderr, "Error interpolating user-defined molecular absorption coefficients\n");
      fprintf (stderr, "to common atmospheric vertical grid.\n");
      return status;
    }
      
  }
  
  /**************************************************************************************/
  /* When using representative wavelengths, molecular absorption lookup tables are used */
  /**************************************************************************************/
  if (input.ck_scheme == CK_REPTRAN || input.ck_scheme == CK_REPTRAN_CHANNEL)
    if (output->atm.molabs != MOLABS_NONE)
      output->atm.molabs = MOLABS_LOOKUP;

  /******************/
  /* Verbose output */
  /******************/

  if (input.verbose) {
 
    fprintf (stderr, " --------------------------------------------------------------------------------------------------------------------------");
    if (output->column[MOL_BRO]  !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_OCLO] !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_HCHO] !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_O4]   !=0)
      fprintf (stderr, "--------------");
    fprintf (stderr, "\n");
    fprintf (stderr, "  lc |  z[km]  |  Pressure  | Temp.  |    Air      |   Ozone     |     O2      | Water vap.  |    CO2      |    NO2      | ");
    if (output->column[MOL_BRO]  !=0)
      fprintf (stderr, "    BRO      |");
    if (output->column[MOL_OCLO] !=0)
      fprintf (stderr, "    OCLO     |");
    if (output->column[MOL_HCHO] !=0)
      fprintf (stderr, "    HCHO     |");
    if (output->column[MOL_O4]   !=0)
      fprintf (stderr, "     O4      |");
    fprintf (stderr, "\n");
    fprintf (stderr, "     |         |   [hPa]    |  [K]   |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    | ");
    if (output->column[MOL_BRO]  !=0)
      fprintf (stderr, "   [cm-3]    |");
    if (output->column[MOL_OCLO] !=0)
      fprintf (stderr, "   [cm-3]    |");
    if (output->column[MOL_HCHO] !=0)
      fprintf (stderr, "   [cm-3]    |");
    if (output->column[MOL_O4]   !=0)
      fprintf (stderr, "[1.0e+46cm-6]|");
    fprintf (stderr, "\n");
    fprintf (stderr, " --------------------------------------------------------------------------------------------------------------------------");
    if (output->column[MOL_BRO]  !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_OCLO] !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_HCHO] !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_O4]   !=0)
      fprintf (stderr, "--------------");
    fprintf (stderr, "\n");


    for (lc=0; lc<output->atm.nlev; lc++) {
      if (output->atm.zd[lc] >= output->alt.altitude) {
        fprintf (stderr, "%5d  %8.4f  %10.5f   %6.2f   %.5e   %.5e   %.5e   %.5e   %.5e   %.5e",
	         lc, output->atm.zd[lc],
	         output->atm.microphys.press[0][0][lc], output->atm.microphys.temper[0][0][lc],output->atm.microphys.dens[MOL_AIR][0][0][lc],
	         output->atm.microphys.dens[MOL_O3][0][0][lc],output->atm.microphys.dens[MOL_O2][0][0][lc], output->atm.microphys.dens[MOL_H2O][0][0][lc],
                 output->atm.microphys.dens[MOL_CO2][0][0][lc],output->atm.microphys.dens[MOL_NO2][0][0][lc]);
        if (output->column[MOL_BRO]  !=0)
          fprintf (stderr, "   %.5e",output->atm.microphys.dens[MOL_BRO][0][0][lc]);
        if (output->column[MOL_OCLO] !=0)
          fprintf (stderr, "   %.5e",output->atm.microphys.dens[MOL_OCLO][0][0][lc]);
        if (output->column[MOL_HCHO] !=0)
          fprintf (stderr, "   %.5e",output->atm.microphys.dens[MOL_HCHO][0][0][lc]);
        if (output->column[MOL_O4]   !=0)
          fprintf (stderr, "   %.5e",output->atm.microphys.dens[MOL_O4][0][0][lc]);
        fprintf (stderr, "\n");
      }
    }

    fprintf (stderr, " -------------------------------------------------------------------------------------------------------------------------");
    if (output->column[MOL_BRO]  !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_OCLO] !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_HCHO] !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_O4]   !=0)
      fprintf (stderr, "--------------");
    fprintf (stderr, "\n");
    fprintf (stderr, " %4s | %6.2f | %10.5f | %6.2f | %7.1ecm-2 |  %7.3f DU | %7.1ecm-2 | %6.3f kg/m2 | %7.1ecm-2 | %7.5f DU  |",
	     "sum", 0.0/0.0, 0.0/0.0, 0.0/0.0, output->column[MOL_AIR]*DU2particle, output->column[MOL_O3], 
	     output->column[MOL_O2]*DU2particle, output->precip_water, output->column[MOL_CO2]*DU2particle, output->column[MOL_NO2]);
    if (output->column[MOL_BRO]  !=0)
      fprintf (stderr, " %7.5f DU |",output->column[MOL_BRO]);
    if (output->column[MOL_OCLO] !=0)
      fprintf (stderr, " %7.5f DU |",output->column[MOL_OCLO]);
    if (output->column[MOL_HCHO] !=0)
      fprintf (stderr, " %7.5f DU |",output->column[MOL_HCHO]);
    if (output->column[MOL_O4] !=0)
      fprintf (stderr, " %7.1ecm-5|",output->column[MOL_O4]*DU2particle*1.0e+46);
    fprintf (stderr, "\n");
    fprintf (stderr, " -------------------------------------------------------------------------------------------------------------------------");
    if (output->column[MOL_BRO]  !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_OCLO] !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_HCHO] !=0)
      fprintf (stderr, "--------------");
    if (output->column[MOL_O4]   !=0)
      fprintf (stderr, "--------------");
    fprintf (stderr, "\n");

  }

  /* read wind data, if given */
  if (strlen(input.filename[FN_ECMWF_WIND_MAP]) > 0) {
    status = read_wind_from_ECMWF_file ( input.latitude, input.longitude, input.UTC, input.atm.time_interpolate, 
                                         input.filename[FN_ECMWF_WIND_MAP], output->alt.altitude,
                                         &(output->wind.u), &(output->wind.v), &(output->wind.w), &(output->wind.z), &(output->wind.nlev), 
                                         input.verbose, input.quiet);                                         
    if (status!=0) {
#if HAVE_LIBNETCDF
      fprintf (stderr, "Error '%s' reading '%s' from %s (line %d, function '%s' in '%s')\n", 
	       nc_strerror(status), "u(z) and v(z)", input.filename[FN_ECMWF_WIND_MAP], __LINE__, __func__, __FILE__ );
#else
      fprintf (stderr, "Error %d reading '%s' from %s (line %d, function '%s' in '%s')\n", 
	       status,              "u(z) and v(z)", input.filename[FN_ECMWF_WIND_MAP], __LINE__, __func__, __FILE__ );
#endif
      return status;
    }
  }


  /*****************************************************************/
  /* calculate special VERTICAL FIELD AT "zout levels" ZOUT-LEVLES */
  /* for heating rate calculations and for output user             */
  /*****************************************************************/
  if (input.verbose) {
    fprintf (stderr, " ... calculate output properties on zout levels\n"); 
  }


  /* allocation */
  /* gas densities */
  if ((output->atm.microphys.dens_zout = calloc (MOL_NN, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for dens_zout (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -30;
  }
  for (i_gas=0;i_gas<MOL_NN;i_gas++) {
    /* gas density */
    if ((output->atm.microphys.dens_zout[i_gas] = calloc (output->atm.nlev, sizeof (float))) == NULL) {
      fprintf (stderr,"Error allocating memory for dens_zout[%d] (line %d, function '%s' in '%s')\n", i_gas, __LINE__, __func__, __FILE__ );
      return -31;
    }
  }

  /* pressure */
  switch (input.atm.zout_source) {
  case OUTLEVEL_ZOUT_ABOVE_SUR:
  case OUTLEVEL_ZOUT_ABOVE_SEA:
  case OUTLEVEL_ATM_LEVELS:
  case OUTLEVEL_ALL_LEVELS:
  case OUTLEVEL_MODEL_LEVELS:
  case OUTLEVEL_MODEL_LAYERS:
  case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
    if ((output->atm.microphys.press_zout = calloc (output->atm.nlev, sizeof (float))) == NULL) {
      fprintf (stderr,"Error allocating memory for press_zout (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return -34;
    }
    break;
  case OUTLEVEL_PRESS:
    if ((output->atm.microphys.press_zout = calloc (output->atm.nzout, sizeof (float))) == NULL) {
      fprintf (stderr,"Error allocating memory for press_zout (line %d, function '%s' in '%s')\n",  __LINE__, __func__, __FILE__ );
      return -34;
    }
    break;
  default:
    fprintf (stderr, "Error, unknown output level type %d! (line %d, function '%s' in '%s')\n", input.atm.zout_source, __LINE__, __func__, __FILE__ );
    fprintf (stderr, "This is a program bug, please report this bug to the programmers and include the input file. Thanx.\n");
    return -1;
  }

  /* temperature */
  if ((output->atm.microphys.temper_zout      = calloc (output->atm.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for temper_zout (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -33;
  }
  /* specific heat capacity at constant pressure */
  if ((output->atm.microphys.c_p              = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for c_p (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -34;
  }
  /* potential temperature theta */
  if ((output->atm.microphys.theta            = calloc (output->atm.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for theta (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -35;
  }
  /* potential temperature theta */
  if ((output->atm.microphys.theta_zout       = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for theta_zout (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -35;
  }
  /* dewpoint temperature */
  if ((output->atm.microphys.temper_d_zout    = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for temper_d_zout (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -36;
  }
  /* equivalent potential temperature theta_e */
  if ((output->atm.microphys.theta_e_zout = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for theta_equiv_zout (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -37;
  }
  /* "cloud liquid water" content CLWC */
  for (isp=0; isp<input.n_caoth; isp++) {
    if ((output->caoth[isp].microphys.lwc_zout = calloc (output->atm.nzout, sizeof (float))) == NULL) {
      fprintf (stderr,"Error allocating memory for mass concentration (LWC) for %s (line %d, function '%s' in '%s')\n", input.caoth[isp].fullname, __LINE__, __func__, __FILE__);
      return -38;
    }
    /* caoth effective radius effr */
    if ((output->caoth[isp].microphys.effr_zout = calloc (output->atm.nzout, sizeof (float))) == NULL) {
      fprintf (stderr,"Error allocating memory for effective radius for %s (line %d, function '%s' in '%s')\n", input.caoth[isp].fullname, __LINE__, __func__, __FILE__);
      return -39;
    }
  }
  /* cloud fraction */
  if ((output->cf.cf_zout = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for cloud fraction (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
    return -42;
  }

  /* wind data */
  if (strlen(input.filename[FN_ECMWF_WIND_MAP]) > 0) {
    /* wind component u_zout */
    if ((output->wind.u_zout = calloc (output->atm.nzout, sizeof (float))) == NULL) {
      fprintf (stderr,"Error allocating memory for 'wind.u_zout'(line %d, function '%s' in '%s') \n", __LINE__, __func__, __FILE__);
      return -43;
    }
    /* wind component v_zout */
    if ((output->wind.v_zout = calloc (output->atm.nzout, sizeof (float))) == NULL) {
      fprintf (stderr,"Error allocating memory for 'wind.v_zout' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return -44;
    }
    /* wind component w_zout */
    if ((output->wind.w_zout = calloc (output->atm.nzout, sizeof (float))) == NULL) {
      fprintf (stderr,"Error allocating memory for 'wind.w_zout' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return -45;
    }
  }

  /* DEBUGGING LINES */
  /* for (lc=0;lc<output->atm.nlev;lc++) */
  /*   fprintf(stderr," DEBUG !!! zd[%3d] = %11.6f\n",lc,output->atm.zd[lc]); */
  /* for (lc=0;lc<output->atm.nzout;lc++) */
  /*   fprintf(stderr," DEBUG !!! zout_sea[%3d] = %11.6f\n",lc,output->atm.zout_sea[lc]); */

  /* calculation of number density of trave gases (or dummy value for MYSTIC 2D surface ) */
  if (output->atm.zout_sur[0] != ZOUT_MYSTIC && input.rte.solver!=SOLVER_SSLIDAR) { /* all cases except MYSTIC 2D surface */

    /* AIR DENSITY (must be before other gases, because result */
    /*  is needed for linmix interpolation of the other gases) */

    /* copy original values */
    for (iz=0; iz<output->atm.nlev; iz++)
      output->atm.microphys.dens_zout[MOL_AIR][iz] = output->atm.microphys.dens[MOL_AIR][0][0][iz];

    /* interpolation */
    status = interpolate_profile (output->atm.zd, &(output->atm.microphys.dens_zout[MOL_AIR]), output->atm.nlev,
                                  output->atm.zout_sea, output->atm.nzout, 
                                  input.atm.interpol_method_gas[MOL_AIR],input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error interpolating dens_zout[MOL_AIR] (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
      return status;
    }

    /* GAS DENSITY */
    for (i_gas=0;i_gas<MOL_NN;i_gas++) {

      /* other than air */
      if (i_gas != MOL_AIR ) {

        /* copy atmosphere values */
        for (iz=0; iz<output->atm.nlev; iz++)
          output->atm.microphys.dens_zout[i_gas][iz] = output->atm.microphys.dens[i_gas][0][0][iz];

        /* interpolation */
        status = interpolate_density (output->atm.zd, &(output->atm.microphys.dens_zout[i_gas]), output->atm.nlev, 
				      output->atm.zout_sea, output->atm.nzout, input.atm.interpol_method_gas[i_gas],
				      output->atm.microphys.dens[MOL_AIR][0][0], output->atm.microphys.dens_zout[MOL_AIR], input.quiet);
        if (status!=0) {
          fprintf (stderr, "Error %d interpolate_density for dens_zout[%d] (line %d, function '%s' in '%s')\n",
		   status, i_gas, __LINE__, __func__, __FILE__);
          return status;
        }
      }
    }

    /* pressure */
    switch (input.atm.zout_source) {
    case OUTLEVEL_ZOUT_ABOVE_SUR:
    case OUTLEVEL_ZOUT_ABOVE_SEA:
    case OUTLEVEL_ATM_LEVELS:
    case OUTLEVEL_ALL_LEVELS:
    case OUTLEVEL_MODEL_LEVELS:
    case OUTLEVEL_MODEL_LAYERS:
    case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
      /* interpolate PRESSURE to zout levels */
      /* copy values from atmosphere grid */
      for (iz=0; iz<output->atm.nlev; iz++)
        output->atm.microphys.press_zout[iz] = output->atm.microphys.press[0][0][iz];
        
      /* interpolation, size is adjusted automagically in interpolate_profile */
      status = interpolate_profile (output->atm.zd, &output->atm.microphys.press_zout,
				    output->atm.nlev,
                                    output->atm.zout_sea, output->atm.nzout, 
                                    input.atm.interpol_method_press, input.quiet);
      if (status!=0) {
        fprintf (stderr, "Error interpolating press_zout\n");
        return status;
      }
      break;
    case OUTLEVEL_PRESS:
      /* copy values from the input (just in order to avoid rounding inaccuracies) */
      for (iz=0; iz<output->atm.nzout; iz++) {
        output->atm.microphys.press_zout[iz] = output->atm.press_zout_org[output->atm.zout_index[iz]];
      }
      break;
    default:
      fprintf (stderr, "Error, unknown output level type %d (line %d, function '%s' in '%s').\n", input.atm.zout_source, __LINE__, __func__, __FILE__ );
      fprintf (stderr, "This is a program bug, please report this bug to the programmers and include the input file. Thanx.\n");
      return -1;
    }


    /* TEMPERATURE */
    /* copy original values, size is adjusted automagically in interpolate_profile */
    for (iz=0; iz<output->atm.nlev; iz++)
      output->atm.microphys.temper_zout[iz] = output->atm.microphys.temper[0][0][iz];
        
    status = interpolate_profile (output->atm.zd, &output->atm.microphys.temper_zout,
				  output->atm.nlev,
                                  output->atm.zout_sea, output->atm.nzout, 
                                  input.atm.interpol_method_temper,input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error interpolating temper_zout\n");
      return status;
    }

    /* replace interpolated cold point temperature, which would smear the extrem of cold temperature */
    /* with double linear interpolated CP-temperature calculated in setup_zout                       */
    if ( output->atm.microphys.temper_cpt != NOT_DEFINED_FLOAT ) {
      switch (input.atm.zout_source) {
      case OUTLEVEL_ZOUT_ABOVE_SUR:
      case OUTLEVEL_ZOUT_ABOVE_SEA:
        for (iz=0; iz<output->atm.nzout; iz++) {
          if ( output->atm.zout_sea[iz] == output->atm.microphys.z_cpt_sea ) {
            output->atm.microphys.temper_zout[iz] = output->atm.microphys.temper_cpt;
            break;
          }
        }
        break;
      case OUTLEVEL_PRESS:
        for (iz=0; iz<output->atm.nzout; iz++) {
          if ( fabs(output->atm.microphys.press_zout[iz] - output->atm.microphys.press_cpt) < 1.e-6*output->atm.microphys.press_cpt ) {
            output->atm.microphys.temper_zout[iz] = output->atm.microphys.temper_cpt;
            break;
          }
        }
        break;
      case OUTLEVEL_ATM_LEVELS:
      case OUTLEVEL_ALL_LEVELS:
      case OUTLEVEL_MODEL_LEVELS:
      case OUTLEVEL_MODEL_LAYERS:
      case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
        break;
      default:
        fprintf (stderr, "Error determining output levels zout.\n");
        fprintf (stderr, "       unknown option %d for zout_source! (line %d, function '%s' in '%s')\n", 
                 input.atm.zout_source, __LINE__, __func__, __FILE__);
        return -1;
      }
    }


    /* specific heating capacity of moist air */
    /* dependent on h2o-dens and T            */
    for (iz=0;iz<output->atm.nzout;iz++) {

      mmr_h2o = (output->atm.microphys.dens_zout[MOL_H2O][iz] * input.atm.mol_mass[MOL_H2O]) /
	(output->atm.microphys.dens_zout[MOL_AIR][iz] * input.atm.mol_mass[MOL_AIR]);

      /* temperature dependent heat capacity of dry air */
      c_p_dry_air = specific_heat_capacity_dry_air (output->atm.microphys.temper_zout[iz], input.quiet);
      if (c_p_dry_air < 0.0) {
        fprintf (stderr, "Error returned by specific_heat_capacity_dry_air()\n");
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -1;
      }

      /* temperature dependent heat capacity of water vapour */
      c_p_wv = specific_heat_capacity_water_vapour(output->atm.microphys.temper_zout[iz], input.quiet);
      if (c_p_wv < 0.0) {
        fprintf (stderr, "Error returned by specific_heat_capacity_water_vapour()\n");
        return -1;
      }

      /* mass weighted mean of specific heat capacity (moist air), e.g. Etling page 34 */
      output->atm.microphys.c_p[iz] = ( 1 - mmr_h2o ) * c_p_dry_air  +  mmr_h2o * c_p_wv;
    }

    /* potential temperature THETA from p and T */
    /* on atmosphere levels */
    for (lc=0;lc<output->atm.nlev;lc++)
      output->atm.microphys.theta[lc] = output->atm.microphys.temper[0][0][lc] * 
	pow((1000./output->atm.microphys.press[0][0][lc]),kappa);
    /* on zout levels */
    for (iz=0;iz<output->atm.nzout;iz++)
      output->atm.microphys.theta_zout[iz] = output->atm.microphys.temper_zout[iz] *
	pow((1000./output->atm.microphys.press_zout[iz]),kappa);

    /* calculate DEWPOINT TEMPERATURE from dens_h2o */
    for (iz=0;iz<output->atm.nzout;iz++) {

      press_h2o = output->atm.microphys.dens_zout[MOL_H2O][iz] * BOLTZMANN * output->atm.microphys.temper_zout[iz] * 1E+04; 
      /* 10E4 == J/cm-3 -> hPa */

      output->atm.microphys.temper_d_zout[iz] = dewpoint(press_h2o);
    }

    /* calculate EQUIVALENT POTENTIAL TEMPERATURE from  T, T_D, and p    */
    /* approximately theta_e = (T + L_v R/c_p h2o_mix) (p_0/p(z))^kappa) */
    /* but L_v is temperature dependent                                  */
    for (lc=0;lc<output->atm.nzout;lc++) {
      output->atm.microphys.theta_e_zout[lc] = EPT (output->atm.microphys.temper_zout[lc], 
                                                    output->atm.microphys.temper_d_zout[lc],
						    output->atm.microphys.press_zout[lc],
                                                    output->atm.microphys.dens_zout[MOL_AIR][lc],
                                                    output->atm.microphys.dens_zout[MOL_H2O][lc]);
    }

    /* Cloud LWC and reff r_eff */
    /* pay attention to different order */
    /* zout_sea in ascending  order 0,1,2,3 ... */
    /* wc.zd    in descending order 6,4,3,2 ... */
    for (isp=0; isp<input.n_caoth; isp++) {
      for (iz=0; iz<output->atm.nzout; iz++) {
	if ( output->atm.zout_sea[iz] < output->caoth[isp].zd[output->caoth[isp].nlev-1] )
	  /* zout below caoths */
	  output->caoth[isp].microphys.lwc_zout[iz] = 0.0;
	else if ( output->atm.zout_sea[iz] >= output->caoth[isp].zd[0] )
	  /* zout above caoths */
	  output->caoth[isp].microphys.lwc_zout[iz] = 0.0;
	else {
	  /* zout in caoth range */
	  for (lc=0; lc<output->caoth[isp].nlev-1; lc++) { 
	    if ( output->caoth[isp].zd[lc+1] <= output->atm.zout_sea[iz] &&
		 output->atm.zout_sea[iz] < output->caoth[isp].zd[lc] ) {
	      output->caoth[isp].microphys.lwc_zout[iz]
		= output->caoth[isp].microphys.lwc_layer[lc];
	      output->caoth[isp].microphys.effr_zout[iz]
		= output->caoth[isp].microphys.effr[lc+1]; /* +1 == conversion to layer */
	      break;
	    }
	  }
	}
      }
    }

    /* Cloud fraction */  
    if (output->cf.nlev>0) { 
      for (iz=0; iz<output->atm.nzout; iz++) {
        if      (output->atm.zout_sea[iz] < output->cf.zd[output->cf.nlev-1]) {
          output->cf.cf_zout[iz] = 0.0;
        }
        else if (output->atm.zout_sea[iz] >= output->cf.zd[0]) {
          output->cf.cf_zout[iz] = 0.0;
        }
        else {
          for (lc=0; lc<output->cf.nlev; lc++) {
            if (output->cf.zd[lc+1] <= output->atm.zout_sea[iz] && output->atm.zout_sea[iz] < output->cf.zd[lc] ) {
              output->cf.cf_zout[iz]  = output->cf.cf[lc];
              break;
            }
          }
        }
      }
    }
    else
      for (iz=0; iz<output->atm.nzout; iz++)
        output->cf.cf_zout[iz] = 0.0/0.0;

    /* wind data */
    if (strlen(input.filename[FN_ECMWF_WIND_MAP]) > 0) {

      /* interpolation of wind.u */
      status = arb_wvn (output->wind.nlev, output->wind.z, output->wind.u, 
                        output->atm.nzout, output->atm.zout_sea, output->wind.u_zout, INTERP_METHOD_LINEAR, 1);

      if (status!=0) {
        fprintf (stderr, "Error %d interpolating u_zout\n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }

      /* interpolation of wind.v */
      status = arb_wvn (output->wind.nlev, output->wind.z, output->wind.v, 
                        output->atm.nzout, output->atm.zout_sea, output->wind.v_zout, INTERP_METHOD_LINEAR, 1);

      if (status!=0) {
        fprintf (stderr, "Error %d interpolating v_zout\n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }

      /* interpolation of wind.w */
      status = arb_wvn (output->wind.nlev, output->wind.z, output->wind.w, 
                        output->atm.nzout, output->atm.zout_sea, output->wind.w_zout, INTERP_METHOD_LINEAR, 1);

      if (status!=0) {
        fprintf (stderr, "Error %d interpolating w_zout\n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }

    }

    /* read and interpolate derivatives of theta */
    for (io=0; io<input.n_output_user; io++) {
      if ( input.output_user[io] == OUTPUT_USER_DTDX || input.output_user[io] == OUTPUT_USER_HEAT_AD_X || input.output_user[io] == OUTPUT_USER_HEAT_AD ) { 
        status = read_dtheta_dx_from_ECMWF_file (input, output);
        if (status!=0) {
          fprintf (stderr, "Error %d returned by read_dtheta_dx_from_netCDF()\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
        break;
      }
    }
    for (io=0; io<input.n_output_user; io++) {
      if ( input.output_user[io] == OUTPUT_USER_DTDY || input.output_user[io] == OUTPUT_USER_HEAT_AD_Y || input.output_user[io] == OUTPUT_USER_HEAT_AD ) { 
        status = read_dtheta_dy_from_ECMWF_file (input, output);
        if (status!=0) {
          fprintf (stderr, "Error %d returned by read_dtheta_dy_from_netCDF()\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
        break;
      }
    }
    for (io=0; io<input.n_output_user; io++) {
      if ( input.output_user[io] == OUTPUT_USER_DTDZ || input.output_user[io] == OUTPUT_USER_HEAT_AD_Z || input.output_user[io] == OUTPUT_USER_HEAT_AD ) { 
        status = calculate_dtheta_dz (input, output);
        if (status!=0) {
          fprintf (stderr, "Error %d returned by calculate_dtheta_dz()\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
        break;
      }
    }
  }
  else { /* for MYSTIC 2D surface set dummy value */

    for (i_gas=0;i_gas<MOL_NN;i_gas++)
      output->atm.microphys.dens_zout[i_gas][0] = -999.0;
    output->atm.microphys.press_zout[0]         = -999.0;
    output->atm.microphys.temper_zout[0]        = -999.0;
    output->atm.microphys.theta_zout[0]         = -999.0;
  }

  /* -------------------------------------------------------------------------------------------*/
  /*      every z-level added to zd_common will be taken care of in setup_redistribute()        */
  /* -------------------------------------------------------------------------------------------*/

  /*************************************************************/
  /*  Initialise atm.zd_common grid, all properties will be at */
  /*  this grid latest after setup_redistribute().             */
  /*  output->atm.microphys.xxx is defined on output->atm.zd   */ 
  /*  until setup_redistribute()                               */ 
  /*************************************************************/

  /* Set common vertical resolution */
  output->atm.nlev_common = 0;
  set_common_z_grid (output->atm.zd, output->atm.nlev,
		     &output->atm.zd_common, &output->atm.nlev_common);

  /* if zout_interpolate is not specified, data for                  */
  /* altitude (2nd argument) levels will be added in redistribute.c  */
  /* and zout-Levels will not be added anywhere to the z-grid        */
  /* solvers do that for us                                          */

  /* add caoth levels */
  /* (no interpolation, but redistribution on these levels */
  /* add the caoth levels to the common profile */
  for (isp=0; isp<input.n_caoth; isp++) {
    if (input.caoth[isp].ipa) {
      for (ipa=0; ipa<output->nipa; ipa++)
	set_common_z_grid (output->caoth_ipa[isp][ipa].zd, output->caoth_ipa[isp][ipa].nlev, 
			   &output->atm.zd_common, &output->atm.nlev_common);
    }
    else {
      set_common_z_grid (output->caoth[isp].zd, output->caoth[isp].nlev, 
			 &output->atm.zd_common, &output->atm.nlev_common);
    }
  }

  /* this is just a check for stupid programming errors, in usual case this is unnessesary */
  if (input.verbose) 
    fprintf (stderr, " ... check ideal gas law for pressure, temperature and number density\n"); 

  /* Finally, check if profiles of pressure, temperature and air density are consistent */
  /* with ideal gas law over 6 orders of magnitude in pressure; check only first six    */
  /* orders of magnitude because for very small pressures roundoff errors due to too    */
  /* few digits cause problems                                                          */
  firstnonideal=1;
  for (lc=0; lc<output->atm.nlev; lc++) {

    idealgas=1;


    if (output->atm.microphys.press[0][0][lc]/output->atm.microphys.press[0][0][output->atm.nlev-1]>1e-6) {
      tmpdens = output->atm.microphys.press[0][0][lc]*1E-04 / (BOLTZMANN*output->atm.microphys.temper[0][0][lc]);
      
      if ((tmpdens==0 && output->atm.microphys.dens[MOL_AIR][0][0][lc]!=0) || 
	  (tmpdens!=0 && output->atm.microphys.dens[MOL_AIR][0][0][lc]==0))
	idealgas=0;
      
      if (output->atm.microphys.dens[MOL_AIR][0][0][lc]>0)
	if ((idealdiff = fabs ((tmpdens - output->atm.microphys.dens[MOL_AIR][0][0][lc]) / output->atm.microphys.dens[MOL_AIR][0][0][lc])) > 0.01)
	  idealgas=0;
      
      if (!idealgas && firstnonideal) {
	firstnonideal=0;
	fprintf (stderr, "\n");
	fprintf (stderr, "ERROR, pressure, temperature, and air density not\n");
	fprintf (stderr, "                  consistent with ideal gas law at the following levels\n");
	fprintf (stderr, "                  (assuming a Boltzmann constant %.6e):\n", BOLTZMANN);
	fprintf (stderr, "    z [km]  density [cm-3]  p/kT [cm-3]  difference  press[hPa]  T[K]\n");
      }
      
      if (!idealgas) {
	fprintf (stderr, "   %7.2f  %.4e      %.4e    %6.1f%%   %10.5f  %8.3f\n", 
		 output->atm.zd[lc], output->atm.microphys.dens[MOL_AIR][0][0][lc], tmpdens, idealdiff*100.0, 
                 output->atm.microphys.press[0][0][lc], output->atm.microphys.temper[0][0][lc]);
      }
    }
  }
  
  if (!firstnonideal) {
    fprintf (stderr, "*** Please check pressure, temperature, and air density in\n");
    fprintf (stderr, "*** %s for consistency!\n", input.filename[FN_ATMOSPHERE]);
    fprintf (stderr, "\n");
    return -1;
  }

  for (isp=0; isp<input.n_caoth; isp++)
    free(output->caoth[isp].z3D);

  return 0;
}
 

/***********************************************************************************/
/* Function: scale_profile_with_mixing_ratio                                       */
/*                                                                                 */
/* Description:                                                                    */
/* Initializes output->mixing_ratio and scales the profile with an user specified  */
/* mixing ratio at an user-defined altitude for a single species.                  */
/*                                                                                 */
/* Parameters:                                                                     */
/* i_mx: MX_XX value of species                                                    */
/* i_mol: MOL_XX value of species (if species has no MOL_XX: i_mol=-1)             */
/* ialt: index of user-specified altitude level                                    */
/* default_ppm: default value for output->mixing_ratio which is used if no         */
/*                 mixing_ratio was specified by the user;                         */
/*                 default_ppm=-1: there is no default value                       */
/* mol_name: String with species name                                              */                
/***********************************************************************************/
int scale_profile_with_mixing_ratio(input_struct input, output_struct *output, int i_mx, int i_mol, 
				    int ialt, double default_ppm, char *mol_name){
  
  int status, lc;
  double scale_factor;

  if (input.mixing_ratio[i_mx]<0) {

    /* If no mixing_ratio not specified, use the default_ppm or the mixing ratio at the user-defined altitude */
    if(default_ppm>=0)
      output->mixing_ratio[i_mx] = default_ppm;
    else if(i_mol>0) 
      output->mixing_ratio[i_mx] = output->atm.microphys.dens[i_mol][0][0][ialt] / output->atm.microphys.dens[MOL_AIR][0][0][ialt] * 1.0e6;
    else{
      fprintf(stderr, "Error: internal problem with scale_profile_with_mixing_ratio() for %s!\n", mol_name);
      return -1;
    }

  }
  else {

    output->mixing_ratio[i_mx] = input.mixing_ratio[i_mx];

    /* do the actual scaling */
    if(i_mol>0){

      if(input.mixing_ratio[i_mx]>0 && output->atm.microphys.dens[i_mol][0][0][ialt]==0){
        fprintf(stderr, "Error, unable to scale profile to %s-mixing ratio at user-defined altitude because the density is zero at this altitude!\n", mol_name);
        return -1;
      }
    
      scale_factor = (input.mixing_ratio[i_mx] / 1.0e6) / 
        (output->atm.microphys.dens[i_mol][0][0][ialt] / 
         output->atm.microphys.dens[MOL_AIR][0][0][ialt]);
    
      if (!input.quiet)
        fprintf(stderr, " ... scaling %s-mixing ratio from %10.5e ppm to %10.5e ppm\n", mol_name,
		(output->atm.microphys.dens[i_mol][0][0][ialt]/output->atm.microphys.dens[MOL_AIR][0][0][ialt])*1e6,
		input.mixing_ratio[i_mx]); 

      for (lc=0; lc<output->atm.nlev; lc++)
        output->atm.microphys.dens[i_mol][0][0][lc] *= scale_factor;

      /* calculate new column and new dens_avg values */
      status = column (output->atm.microphys.dens[i_mol][0][0], output->atm.zd, output->atm.nlev, output->alt.altitude, 
		       input.atm.interpol_method_gas[i_mol], output->atm.microphys.dens[MOL_AIR][0][0],
		       &(output->atm.microphys.dens_avg[i_mol][0][0]), NO, &(output->column[i_mol]));  /* in atmosphere.c */
      if (status!=0) {
        fprintf (stderr, "\nError %d returned by column()\n", status);
        fprintf (stderr, "      during calculation of %s column (line %d, function '%s' in '%s')\n", "O2", __LINE__, __func__, __FILE__);
        return status;
      }

      if(i_mol==MOL_H2O) 
        output->precip_water = output->column[MOL_H2O] * DU2particle * (M_H2O / N_MOL_CKDFU * 10.); 

    }
    else{

      if (!input.quiet)
        fprintf(stderr, " ... using %s-mixing ratio %10.5e ppm\n", mol_name,output->mixing_ratio[i_mx]);
    
    }
 
  }

  return 0;

}

/***********************************************************************************/
/* Function: setup_zout_levels                                                     */
/*                                                                                 */
/* Description:                                                                    */
/*  * create list of output levels: zout_sea, and zout_sur                         */
/*  * add all zout_sea-levels to z_all and nlev_all                                */
/*                                                                                 */
/*  In detail:                                                                     */
/*  Calculate the output levels 'zout' from different input options:               */
/*  zout, zout_sea, pressure_out, atm_levels, or all_levels                        */
/*  Replace place holder ZOUT_TOA (top of atmosphere), and ZOUT_SURFACE.           */
/*  Throw away output levels, which are below the surface or above TOA.            */
/*  Remember original zout levels in (nessesary for empty levels in netCDF output) */
/*  output->atm.zout_sea_org, output->atm.zout_sur_org, and output->atm.nzout_org. */
/*  If user wants local heating rates, than add some additional levels above and   */
/*  below the zout levels, for best possible calculation of absorption             */
/*  coefficients on levels (instead of layers).                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/*     input_struct input                                                          */
/*     output_struct *output                                                       */
/*                                                                                 */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    atmosphere.c                                                          */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*    Jan 2008   U. Hamann     Additional atm. levels for local heating            */
/*    Feb 2008   U. Hamann     Additional zout level for cold point tropopause     */
/*                                                                                 */
/***********************************************************************************/

int setup_zout_levels(input_struct input, output_struct *output, 
                      float **z_all, int *nlev_all)
{

  int status =0;

  int lc=0, i=0, iz=0;
  float  *tmp_zout  = NULL;

  int nuser = -999;
  float *zuser = NULL;

  float  *log_p_atm  = NULL;
  float  *log_p_zout = NULL;

  float z_below = NOT_DEFINED_FLOAT, z_above = NOT_DEFINED_FLOAT;

  int status_CPT = NOT_DEFINED_INTEGER;
  int skipp_CPT = 0;

  int shift = 0;

  if (input.verbose)
    fprintf (stderr, " ... setup_zout_levels\n");

  /*********************************************************/
  /* determine the number of zout levels                   */
  /* allocating output-levels fields                       */
  /* (relative to sea level and relative to surface level) */
  /*********************************************************/
  if (input.rte.solver==SOLVER_SSLIDAR) {
    /* special case! */
    output->atm.nzout_org = input.sslidar_nranges;
    tmp_zout=calloc (output->atm.nzout_org, sizeof(float));
    for (iz=0;iz<output->atm.nzout_org;iz++)
      tmp_zout[iz] = input.sslidar[SSLIDAR_RANGE] * ( (double) iz + 0.5 ) + output->alt.altitude;
  }
  else {
    switch (input.atm.zout_source) {
    case OUTLEVEL_ZOUT_ABOVE_SUR:
    case OUTLEVEL_ZOUT_ABOVE_SEA:
    case OUTLEVEL_PRESS:
      output->atm.nzout_org = input.atm.nzout;
      break;
    case OUTLEVEL_ATM_LEVELS:
      output->atm.nzout_org = output->atm.nlev;
      break;
    case OUTLEVEL_ALL_LEVELS:
      output->atm.nzout_org = 0;
      set_common_z_grid ( (*z_all), (*nlev_all), &tmp_zout, &output->atm.nzout_org );

      if (input.alt.altitude_dz > 0) {
        /* add altitude levels for 2nd argument of altitude to tmp_zout */
        nuser = (int) (((*z_all)[0] - (*z_all)[(*nlev_all)-1]) / 
	               input.alt.altitude_dz          + 0.5) + 1;
        zuser = (float *) calloc (nuser, sizeof(float));
        for (i=0; i<nuser; i++) 
          zuser[i] = (*z_all)[(*nlev_all)-1] + input.alt.altitude_dz * (float) i;
        while (zuser[nuser-1] > (*z_all)[0])
          nuser--;

        /* add user-defined grid to atmospheric grid */
        set_common_z_grid (zuser, nuser, &tmp_zout, &output->atm.nzout_org);
        free (zuser);
      } 
      break;
    case OUTLEVEL_MODEL_LEVELS:
      if ( input.atm.nzout == NOT_DEFINED_INTEGER )
        output->atm.nzout_org = output->atm.nz_model_level;
      else 
        output->atm.nzout_org = input.atm.nzout;
      break;
    case OUTLEVEL_MODEL_LAYERS:
      if ( input.atm.nzout == NOT_DEFINED_INTEGER )
        output->atm.nzout_org = output->atm.nz_model_layer;
      else 
        output->atm.nzout_org = input.atm.nzout;
      break;
    case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
      /* first get all level, add layers, and reduce later, if wanted */
      output->atm.nzout_org = output->atm.nz_model_level;
      break;
    default:
      fprintf (stderr, "Error determining output levels zout.\n");
      fprintf (stderr, "       unknown option %d for zout_source! (line %d, function '%s' in '%s')\n", input.atm.zout_source, __LINE__, __func__, __FILE__);
      return -1;
    }
  }

  /* allocate memory for zout_org levels */    
  if ((output->atm.zout_sur_org  = calloc (output->atm.nzout_org, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for zout_sur_org (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -1;
  }
  if ((output->atm.zout_sea_org = calloc (output->atm.nzout_org, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for zout_sea_org (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -1;
  }

  /***********************************************************************************/
  /* copy the values from input.atm.zout_* / press_out ... -> output->atm.zout_*_org */
  /* here we replace the marker for surface, top of atmosphere                       */
  /* and cold point tropopause with real altitude or pressure values                 */
  /***********************************************************************************/
  switch (input.atm.zout_source) {
  case OUTLEVEL_ZOUT_ABOVE_SUR:
    for (iz=0;iz<output->atm.nzout_org;iz++) {
      if (input.atm.zout_sur[iz] == ZOUT_MYSTIC) {   /* marker for 2D-SURFACE IN MYSTIC */
        output->atm.zout_sea_org[iz - skipp_CPT] = ZOUT_MYSTIC;
        output->atm.zout_sur_org[iz - skipp_CPT] = ZOUT_MYSTIC;
      }
      else if (input.atm.zout_sur[iz] == ZOUT_SURFACE) { /* marker for 2D-SURFACE IN MYSTIC */
        output->atm.zout_sea_org[iz - skipp_CPT] = output->alt.altitude;
        output->atm.zout_sur_org[iz - skipp_CPT] = 0;
      }
      else if (input.atm.zout_sur[iz] == ZOUT_TOA) {   /* marker for TOP_OF_ATMOSPHERE */
        output->atm.zout_sea_org[iz - skipp_CPT] = output->atm.zd[0];
        output->atm.zout_sur_org[iz - skipp_CPT] = output->atm.zd[0] - output->alt.altitude;
      }
      else if (input.atm.zout_sur[iz] == ZOUT_CPT) {   /* marker for COLD_POINT_TROPOPAUSE */
        /* if z_cpt_sea is DEFINED, than CPT was calculated with more accurate data before */
        if ( output->atm.microphys.z_cpt_sea == NOT_DEFINED_FLOAT ) {
          status_CPT = find_cold_point_tropopause ( output->atm.microphys.press[0][0],
						    output->atm.microphys.temper[0][0],
						    output->atm.zd, output->atm.nlev, 
                                                    &(output->atm.microphys.z_cpt_sea), 
						    &(output->atm.microphys.press_cpt), 
                                                    &(output->atm.microphys.temper_cpt),
						    input.quiet, input.verbose );
        }
        if ( status_CPT == 0 ) {
          output->atm.zout_sea_org[iz - skipp_CPT] = output->atm.microphys.z_cpt_sea;
          output->atm.zout_sur_org[iz - skipp_CPT] = output->atm.microphys.z_cpt_sea - output->alt.altitude;
        }
        else {
          /* we did not find a CPT, but we like to continue anyway          */
          skipp_CPT = skipp_CPT + 1;
        }
      }
      else { /* normal case */
        output->atm.zout_sea_org[iz - skipp_CPT] = input.atm.zout_sur[iz] + output->alt.altitude;
        output->atm.zout_sur_org[iz - skipp_CPT] = input.atm.zout_sur[iz]; 
      }
      output->atm.nzout_org = output->atm.nzout_org - skipp_CPT;
    } 
    break;
  case OUTLEVEL_ZOUT_ABOVE_SEA: /* BCA this is almost identical to ABOVE_SUR */
    for (iz=0;iz<output->atm.nzout_org;iz++) {
      if      (input.atm.zout_sea[iz] == ZOUT_MYSTIC) {    /* marker for MYSTIC */
        output->atm.zout_sea_org[iz - skipp_CPT] = ZOUT_MYSTIC;  
        output->atm.zout_sur_org[iz - skipp_CPT] = ZOUT_MYSTIC;
      }
      else if (input.atm.zout_sea[iz] == ZOUT_SURFACE) { /* marker for SURFACE_ELEVATION */
        output->atm.zout_sea_org[iz - skipp_CPT] = output->alt.altitude;
        output->atm.zout_sur_org[iz - skipp_CPT] = 0;
      }
      else if (input.atm.zout_sea[iz] == ZOUT_TOA) {     /* marker for TOP_OF_ATMOSPHERE */
        output->atm.zout_sea_org[iz - skipp_CPT] = output->atm.zd[0];
        output->atm.zout_sur_org[iz - skipp_CPT] = output->atm.zd[0] - output->alt.altitude;
      }
      else if (input.atm.zout_sea[iz] == ZOUT_CPT) {   /* marker for COLD_POINT_TROPOPAUSE */
        /* if z_cpt_sea is DEFINED, than CPT was calculated with more accurate data before */
        if ( output->atm.microphys.z_cpt_sea == NOT_DEFINED_FLOAT ) {
          status_CPT = find_cold_point_tropopause ( output->atm.microphys.press[0][0],
						    output->atm.microphys.temper[0][0],
						    output->atm.zd, output->atm.nlev, 
                                                    &(output->atm.microphys.z_cpt_sea),
						    &(output->atm.microphys.press_cpt), 
                                                    &(output->atm.microphys.temper_cpt),
						    input.quiet, input.verbose);
        }

        if ( output->atm.microphys.z_cpt_sea != NOT_DEFINED_FLOAT ) {
          output->atm.zout_sea_org[iz - skipp_CPT] = output->atm.microphys.z_cpt_sea;
          output->atm.zout_sur_org[iz - skipp_CPT] = output->atm.microphys.z_cpt_sea - output->alt.altitude;
        }
        else {
          /* we did not find an CPT, but we like to continue anyway          */
          skipp_CPT = skipp_CPT + 1;
        }
      }
      else { /* normal case */
        output->atm.zout_sea_org[iz - skipp_CPT] = input.atm.zout_sea[iz];
        output->atm.zout_sur_org[iz - skipp_CPT] = input.atm.zout_sea[iz] - output->alt.altitude;
      }
    }
    output->atm.nzout_org = output->atm.nzout_org - skipp_CPT;
    break;
  case OUTLEVEL_PRESS:
    /* log pressure for atmosphere grid */
    log_p_atm = calloc (output->atm.nlev, sizeof (float));

    if ( ( log_p_atm = calloc (output->atm.nlev, sizeof (float)) )==NULL) {
      fprintf (stderr, "Error allocating memory for 'log_p_atm' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }
    for ( lc=0; lc<output->atm.nlev; lc++ ) {
      log_p_atm[lc] = log(output->atm.microphys.press[0][0][lc]);
    }

    /* log p for zout grid */
    if ( ( log_p_zout = calloc (output->atm.nzout_org, sizeof (float)) )==NULL) {
      fprintf (stderr, "Error allocating memory for 'log_p_atm' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }
    if ( ( output->atm.press_zout_org = calloc (output->atm.nzout_org, sizeof (float)) )==NULL) {
      fprintf (stderr, "Error allocating memory for 'log_p_atm' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    /* replace marker, where possible */
    for (iz=0;iz<output->atm.nzout_org;iz++) {
      if      (input.atm.press_zout[iz] == ZOUT_MYSTIC)    /* marker for MYSTIC */
        output->atm.press_zout_org[iz - skipp_CPT] = ZOUT_MYSTIC;
      else if (input.atm.press_zout[iz] == ZOUT_SURFACE) { /* marker for SURFACE_ELEVATION */
        output->atm.press_zout_org[iz - skipp_CPT] = output->atm.microphys.press[output->atm.nlev-1][0][0];
      }
      else if (input.atm.press_zout[iz] == ZOUT_TOA) {     /* marker for TOP_OF_ATMOSPHERE */
        output->atm.press_zout_org[iz - skipp_CPT] = output->atm.microphys.press[0][0][0];
      }
      else if (input.atm.press_zout[iz] == ZOUT_CPT) {     /* marker for COLD_POINT_TROPOPAUSE */
        /* if z_cpt_sea is DEFINED, than CPT was calculated with more accurate data before */
        if ( output->atm.microphys.z_cpt_sea == NOT_DEFINED_FLOAT ) {
          status_CPT = find_cold_point_tropopause ( output->atm.microphys.press[0][0],
						    output->atm.microphys.temper[0][0],
						    output->atm.zd, output->atm.nlev, 
                                                    &(output->atm.microphys.z_cpt_sea), 
						    &(output->atm.microphys.press_cpt), 
                                                    &(output->atm.microphys.temper_cpt), 
						    input.quiet, input.verbose);
        }
        if ( output->atm.microphys.z_cpt_sea != NOT_DEFINED_FLOAT ) {
          output->atm.press_zout_org[iz - skipp_CPT] = output->atm.microphys.press_cpt;
        }
        else {
          /* we did not find an CPT, but we like to continue anyway          */
          /* therefor we write some values below the ground into the profile */
          skipp_CPT = skipp_CPT + 1;
        }
      }
      else { /* normal case, just copy to output structure */
        output->atm.press_zout_org[iz - skipp_CPT] = input.atm.press_zout[iz];
      }
    }

    output->atm.nzout_org = output->atm.nzout_org - skipp_CPT;

    /* take logarithm for interpolation of z */
    for (iz=0;iz<output->atm.nzout_org;iz++)
      log_p_zout[iz] = log(output->atm.press_zout_org[iz]);

    /* interpolate z-values to the press_zout grid, assuming linear variation of z with log_p */
    status = arb_wvn (output->atm.nlev,  log_p_atm,  output->atm.zd, 
                      output->atm.nzout_org, log_p_zout, output->atm.zout_sea_org, 
                      INTERP_METHOD_LINEAR, 0);
    if (status!=0) {
      fprintf (stderr, "Error %d while determining zout-levels from press_zout (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
      return status;
    }

    free(log_p_zout);
    free(log_p_atm);

    /* correct numberical inaccuracy of interpolation (outputlevel epsilon below altitude) */
    for (iz=0;iz<output->atm.nzout_org;iz++)
      if ( fabs(output->atm.zout_sea_org[iz]-output->alt.altitude) < 1.0e-6 )
        output->atm.zout_sea_org[iz] = output->alt.altitude;

    /* create zout_sur_org from zout_sea_org */
    for(iz=0;iz<output->atm.nzout_org;iz++) {
      output->atm.zout_sur_org[iz] = output->atm.zout_sea_org[iz] - output->alt.altitude;
    }
    break;

  case OUTLEVEL_ATM_LEVELS:
    /* if user wants zout-levels == atm levels -> copy atm to zout */
    if (input.verbose)
      fprintf (stderr, "     copy atmosphere levels to zout levels\n");
    for (lc=0; lc < output->atm.nzout_org; lc++) {
      output->atm.zout_sea_org[lc] = output->atm.zd[lc];
      output->atm.zout_sur_org[lc] = output->atm.zout_sea_org[lc] - output->alt.altitude;
    }
    break;
  case OUTLEVEL_ALL_LEVELS:
    /* if user wants all levels -> copy all to zout */
    if (input.verbose)
      fprintf (stderr, "     copy all levels to zout levels\n");
    for (lc=0; lc < output->atm.nzout_org; lc++) {
      output->atm.zout_sea_org[lc] = tmp_zout[lc];
      output->atm.zout_sur_org[lc] = output->atm.zout_sea_org[lc] - output->alt.altitude;
    }
    free(tmp_zout);
    break;
  case OUTLEVEL_MODEL_LEVELS:
    /* if user wants zout-levels -> copy them to zout */
    if (input.verbose)
      fprintf (stderr, "     copy atmosphere levels from model (ECMWF/ECHAM) to zout levels\n");

    /* copy levels to zout-grid */
    if ( output->atm.nzout_org > output->atm.nz_model_level ) {
      fprintf (stderr, "Error, user specified more zout-levels (%d) with third argument of the input option 'zout'\n", output->atm.nzout_org); 
      fprintf (stderr, "       than there are levels of the atmosphere (%d)!\n", output->atm.nz_model_level);
      fprintf (stderr, "       (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }
    shift = output->atm.nz_model_level - output->atm.nzout_org;
    for (lc=0; lc < output->atm.nzout_org; lc++) {
      output->atm.zout_sea_org[lc] = output->atm.z_model_level[lc + shift] ;
      output->atm.zout_sur_org[lc] = output->atm.zout_sea_org[lc] - output->alt.altitude ;
    }
    break;
  case OUTLEVEL_MODEL_LAYERS:
    /* if user wants zout-layer -> copy them to zout */
    if (input.verbose)
      fprintf (stderr, "     copy atmosphere layers from model (ECMWF/ECHAM) to zout levels\n");

    /* copy levels to zout-grid */
    if ( output->atm.nzout_org > output->atm.nz_model_layer ) {
      fprintf (stderr, "Error, user specified more zout-levels (%d) with third argument of the input option 'zout'\n", output->atm.nzout_org); 
      fprintf (stderr, "       than there are layers of the atmosphere (%d)!\n", output->atm.nz_model_layer);
      fprintf (stderr, "       (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }
    shift = output->atm.nz_model_layer - output->atm.nzout_org;
    for (lc=0; lc < output->atm.nzout_org; lc++) {
      output->atm.zout_sea_org[lc] = output->atm.z_model_layer[lc + shift] ;
      output->atm.zout_sur_org[lc] = output->atm.zout_sea_org[lc] - output->alt.altitude ;
    }
    break;
  case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:

    /* copy all levels */
    for (lc=0; lc < output->atm.nz_model_level; lc++)
      output->atm.zout_sea_org[lc] = output->atm.z_model_level[lc];

    /* add layer to z-grid */
    set_common_z_grid ( output->atm.z_model_layer, output->atm.nz_model_layer, 
                        &(output->atm.zout_sea_org), &(output->atm.nzout_org));

    /* adjust combined zout-grid to user defined size (input.atm.zout) */
    if ( input.atm.nzout != NOT_DEFINED_INTEGER ) {
      if ( input.atm.nzout > output->atm.nzout_org ) {
        fprintf (stderr, "Error, user specified more zout-levels (%d) with third argument of the input option 'zout'\n", input.atm.nzout); 
        fprintf (stderr, "       are more than there are levels and layers of the atmosphere (%d)!\n", output->atm.nzout_org);
        fprintf (stderr, "       (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -1;
      }
      shift = output->atm.nzout_org - input.atm.nzout;
      for (lc=0; lc < input.atm.nzout; lc++)
        output->atm.zout_sea_org[lc] = output->atm.zout_sea_org[lc + shift] ;

      output->atm.nzout_org = input.atm.nzout;
      if ((output->atm.zout_sea_org = realloc (output->atm.zout_sea_org, output->atm.nzout_org * sizeof(float))) == NULL) {
        fprintf (stderr, "Error reallocating memory for zout_sur_org (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
        return -1;
      }
    }
   
    /* reallocate zout_SUR for more levels, as nzout_org increased during set_common_z_grid  */
    if ((output->atm.zout_sur_org = realloc (output->atm.zout_sur_org, output->atm.nzout_org * sizeof(float))) == NULL) {
      fprintf (stderr,"Error reallocating memory for zout_sur_org (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
      return -1;
    }

    for (lc=0; lc < output->atm.nzout_org; lc++)
      output->atm.zout_sur_org[lc] = output->atm.zout_sea_org[lc] - output->alt.altitude ;
    
    break;
  default:
    fprintf (stderr, "Error determining output levels zout.\n");
    fprintf (stderr, "       unknown option %d for zout_source! (line %d, function '%s' in '%s')\n", 
	     input.atm.zout_source, __LINE__, __func__, __FILE__);
    return -1;
  }


  /*******************************************************************/
  /* at this point, the zout_*_org are done                          */
  /*******************************************************************/

  /* DEBUGGING LINES */
  /* for (iz=0;iz<output->atm.nzout_org;iz++) */
  /*   fprintf (stderr, " DEBUG!!! zout_org[%3d] = %11.6f %11.6f\n",  */
  /*     iz, output->atm.zout_sea_org[iz], output->atm.zout_sur_org[iz]); */

  /*******************************************************************/
  /* copy zout_sea_org -> tmp_zout                                   */
  /* skip values below ground and take care of surface and toa flag  */
  /*******************************************************************/

  /* copy zout_sea_org to tmp_zout */
  tmp_zout = calloc (output->atm.nzout_org, sizeof (float));
  for (iz=0;iz<output->atm.nzout_org;iz++) {
    tmp_zout[iz] = output->atm.zout_sea_org[iz];
  }

  /* sort tmp_zout in ascending order */
  qsort ( tmp_zout, output->atm.nzout_org, sizeof(float), QSORT_CAST compare_float_inv);

  /* now check, if there are values below the ground or double occuring values */
  output->atm.nzout=0;
  for (iz=0; iz<output->atm.nzout_org; iz++) {
    if ( output->alt.altitude <= tmp_zout[iz] || tmp_zout[iz] == ZOUT_MYSTIC ) {
      if ( iz == 0 ) {
        output->atm.nzout++;
      }
      else if (tmp_zout[iz] != tmp_zout[iz-1]) {
        output->atm.nzout++;
      }
    }
  }

  /* override check for sslidar */
  if (input.rte.solver==SOLVER_SSLIDAR)
    output->atm.nzout=output->atm.nzout_org;

  output->atm.nzout_user = output->atm.nzout;

  /* Warning */
  if ( output->atm.nzout != output->atm.nzout_org )
    if (!input.quiet)
      fprintf (stderr, "*** Warning, some zout levels (%3d of %3d) are not inside the atmosphere\n",
	       output->atm.nzout_org - output->atm.nzout, output->atm.nzout_org );
  if (input.verbose) {
    fprintf (stderr, "     only following output levels are considered:\n"); 
    fprintf (stderr, "  lc    zout[km]\n");
    fprintf (stderr, "----------------\n");
    for (iz=0; iz<output->atm.nzout_org; iz++) {
      if ( output->alt.altitude <= tmp_zout[iz] || tmp_zout[iz] == ZOUT_MYSTIC )
        fprintf (stderr, " %3d  %9.5f\n", iz, tmp_zout[iz]);
    }
  }

  /* Error */
  if ( output->atm.nzout == 0 ) {
    fprintf (stderr,  "Error, no zout levels of %d levels is inside the atmosphere height range\n", output->atm.nzout_org );
    for (iz=0; iz<output->atm.nzout_org; iz++)
      fprintf (stderr, "     level = %3d, zout = %f\n", iz, output->atm.zout_sea_org[iz] ); 
    return -2;
  }

  /* allocate memory for final zout grids */
  if ((output->atm.zout_sea    = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for zout_sea\n");
    return -1;
  }
  if ((output->atm.zout_sur    = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for zout_sur\n");
    return -1;
  }

  /* sort tmp_zout into final grid */
  /* avoid double occuring values  */
  lc=0;
  for (iz=0; iz<output->atm.nzout_org; iz++)
    if ( output->alt.altitude <= tmp_zout[iz] || tmp_zout[iz] == ZOUT_MYSTIC ) {
      if ( tmp_zout[iz] == ZOUT_MYSTIC ) {
        output->atm.zout_sea[lc]=ZOUT_MYSTIC;
        output->atm.zout_sur[lc]=ZOUT_MYSTIC;
        lc++;
      }
      else {
        if ( iz == 0 ) {
          output->atm.zout_sea[lc]=tmp_zout[iz];
          output->atm.zout_sur[lc]=output->atm.zout_sea[lc] - output->alt.altitude;
          /* fprintf (stderr, " output->atm.zout_sea[%d] = %f\n", lc, output->atm.zout_sea[lc] ); */
          lc++;
        }
        else if ( tmp_zout[iz-1]  != tmp_zout[iz] ) {
          output->atm.zout_sea[lc] = tmp_zout[iz];
          output->atm.zout_sur[lc] = output->atm.zout_sea[lc] - output->alt.altitude;
          /* fprintf (stderr, " output->atm.zout_sea[%d] = %f\n", lc, output->atm.zout_sea[lc] ); */
          lc++;
        }
      }
    }

  free(tmp_zout);

  /* searching for index relationship */
  output->atm.zout_index  = calloc (output->atm.nzout, sizeof (int));
  for (lc=0;lc<output->atm.nzout;lc++) {
    output->atm.zout_index[lc] = -999;
    for (iz=0; iz<output->atm.nzout_org; iz++) {
      if ( fabs(output->atm.zout_sur[lc] - output->atm.zout_sur_org[iz]) < 0.000001 ) {
        output->atm.zout_index[lc] = iz; 
        break;
      }
    }
  }

  /* restore marker for netCDF writing */
  switch (input.atm.zout_source) {
  case OUTLEVEL_ZOUT_ABOVE_SUR:
    for (iz=0;iz<output->atm.nzout_org;iz++) {
      if (input.atm.zout_sur[iz] == ZOUT_TOA) {          /* marker for TOP_OF_ATMOSPHERE */
        output->atm.zout_sea_org[iz] = ZOUT_TOA;
        output->atm.zout_sur_org[iz] = ZOUT_TOA;
      }
      else if (input.atm.zout_sur[iz] == ZOUT_SURFACE) { /* marker for SURFACE */
        output->atm.zout_sea_org[iz] = ZOUT_SURFACE;
        output->atm.zout_sur_org[iz] = ZOUT_SURFACE;
      }
      else if (input.atm.zout_sur[iz] == ZOUT_CPT) {     /* marker for COLD_POINT_TROPOPAUSE */
        output->atm.zout_sea_org[iz] = ZOUT_CPT;
        output->atm.zout_sur_org[iz] = ZOUT_CPT;
      }
    }
    break;
  case OUTLEVEL_ZOUT_ABOVE_SEA:
    for (iz=0;iz<output->atm.nzout_org;iz++) {
      if (input.atm.zout_sea[iz] == ZOUT_TOA) {            /* marker for TOP_OF_ATMOSPHERE */
        output->atm.zout_sea_org[iz] = ZOUT_TOA;
        output->atm.zout_sur_org[iz] = ZOUT_TOA;
      }
      else if (input.atm.zout_sea[iz] == ZOUT_SURFACE) {   /* marker for TOP_OF_ATMOSPHERE */
        output->atm.zout_sea_org[iz] = ZOUT_SURFACE;
        output->atm.zout_sur_org[iz] = ZOUT_SURFACE;
      }
      else if (input.atm.zout_sea[iz] == ZOUT_CPT) {       /* marker for COLD_POINT_TROPOPAUSE */
        output->atm.zout_sea_org[iz] = ZOUT_CPT;
        output->atm.zout_sur_org[iz] = ZOUT_CPT;
      }
    }
    break;
  case OUTLEVEL_PRESS:
    for (iz=0;iz<output->atm.nzout_org;iz++) {
      if (input.atm.press_zout[iz] == ZOUT_TOA) {          /* marker for TOP_OF_ATMOSPHERE */
        output->atm.zout_sea_org[iz] = ZOUT_TOA;
        output->atm.zout_sur_org[iz] = ZOUT_TOA;
      }
      else if (input.atm.press_zout[iz] == ZOUT_SURFACE) { /* marker for TOP_OF_ATMOSPHERE */
        output->atm.zout_sea_org[iz] = ZOUT_SURFACE;
        output->atm.zout_sur_org[iz] = ZOUT_SURFACE;
      }
      else if (input.atm.press_zout[iz] == ZOUT_CPT) {     /* marker for COLD_POINT_TROPOPAUSE */
        output->atm.zout_sea_org[iz] = ZOUT_CPT;
        output->atm.zout_sur_org[iz] = ZOUT_CPT;
      }
    }
    break;
  case OUTLEVEL_ATM_LEVELS:
  case OUTLEVEL_ALL_LEVELS:
  case OUTLEVEL_MODEL_LEVELS:
  case OUTLEVEL_MODEL_LAYERS:
  case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
    /* nothing to do here */
    break;
  default:
    fprintf (stderr, "Error determining output levels zout.\n");
    fprintf (stderr, "       unknown option %d for zout_source! (line %d, function '%s' in '%s')\n", 
	     input.atm.zout_source, __LINE__, __func__, __FILE__);
    return -1;
  }


  /* security check for values below surface and above top of atmosphere */
  if (input.rte.solver!=SOLVER_SSLIDAR) {
    for (iz=0;iz<output->atm.nzout;iz++) { 
      /* check if values are below surface */
      if (output->atm.zout_sea[iz] < output->alt.altitude && output->atm.zout_sea[iz] != ZOUT_MYSTIC) {
	fprintf (stderr, "Error zout[%d] = %f km is below the ground = %f km\n",
		 iz, output->atm.zout_sea[iz],output->alt.altitude);
	fprintf (stderr, "      This is a program bug. Please contact the programmers and include the input file. \n");
	fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
	return -2;
      }
      /* check if values are above atmosphere */
      if (output->atm.zout_sea[iz] > output->atm.zd[0]) {
	fprintf (stderr, "Error zout[%d] = %f km is above top of atmosphere = %f km\n", 
		 iz, output->atm.zout_sea[iz],output->atm.zd[0]);
	return -3;
      }
    } /* for (iz=0; iz<output->atm.nzout; iz++) { ... */
  }

  /* DEBUGGING LINES */
  /* for (iz=output->atm.nzout-1; iz>=0; iz--) */
  /*   fprintf (stderr, " DEBUG!!! zout_sea[%3d] = %11.6f, zout_sur =%11.6f, zout_index = %3d\n",  */
  /*     iz, output->atm.zout_sea[iz], output->atm.zout_sur[iz], output->atm.zout_index[iz]); */

  /* add new levels to the common atmospheric profile */
  /* before interpolation is done */
  if (input.atm.zout_interpolate == ZOUT_INTERPOLATE) {
    for (iz=0; iz<output->atm.nzout; iz++) {
      if (output->atm.zout_sur[iz] != ZOUT_MYSTIC) {
        set_common_z_grid (&(output->atm.zout_sea[iz]), 1, z_all, nlev_all);
      }
    }
  }
  /* else if (input.atm.zout_interpolate == NO_ZOUT_INTERPOLATE) */
  /* zout level are not added to any z-grid */
  /* when input.atm.zout_interpolate == NO_ZOUT_INTERPOLATE */
  /* solver get the zout levels separat and evaluate there by themselfes*/


  /* if simulation result is a local heating rate,                */
  /* it is nessesary to  have small layers around the zout levels */
  /* add 10m above and 10m below zout-levels also levels          */
  if (input.heating == HEAT_LOCAL) {
    for (iz=0; iz<output->atm.nzout; iz++) {
      z_below=output->atm.zout_sea[iz]-0.1;
      if (z_below > output->alt.altitude) {
        set_common_z_grid (&(z_below), 1, z_all, nlev_all);
        /* fprintf (stderr, " add internal interpolation level (below zout): %7.3f\n", z_below); */
      }
      z_above=output->atm.zout_sea[iz]+0.1;
      if (z_above < output->atm.zd[0]) {
        set_common_z_grid (&(z_above), 1, z_all, nlev_all);
        /* fprintf (stderr, " add internal interpolation level (above zout): %7.3f\n", z_above); */
      }
    }
  }

  /* DEBUGGING LINES */
  /* for (lc=0; lc<(*nlev_all); lc++) */
  /*   fprintf (stderr, " z[%d] = %f\n", lc, (*z_all)[lc]); */

  return status;
}

/***********************************************************************************/
/* Function: find_cold_point_tropopause                                            */
/*                                                                                 */
/* Description:                                                                    */
/*  this function searches for the cold point tropopause (CPT)                     */
/*  (here level between 40 and 500hPa with 2 warmer levels above)                  */
/*  then it determines the height and temperature of the CPT                       */
/*                                                                                 */
/* Parameters:                                                                     */
/*     input_struct input                                                          */
/*     output_struct *output                                                       */
/*                                                                                 */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    atmosphere.c                                                          */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Feb 2008   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

static int find_cold_point_tropopause ( float *press, float *temper, float *z, int nlev, 
                                        float *z_CPT_sea, float *p_CPT, float *T_CPT, int quiet, int verbose ) 
{

  int i_CPT_2=NOT_DEFINED_INTEGER,i_CPT_3=NOT_DEFINED_INTEGER;
  int i_CPT = NOT_DEFINED_INTEGER, di_CPT = NOT_DEFINED_INTEGER;
  int i1,i2,i3,i4;
  int lc=0;
  int status=0;

  /* search minimum temperature between 40 and 500hPa */    
  for (lc = nlev-1; lc>2; lc--) {
    /* fprintf (stderr, " check: T[%3d]=%10.5f, p=%10.5f, z=%10.5f\n", lc, temper[lc],press[lc],z[lc]); */
    if ( 40.0 <= press[lc] && press[lc] < 500.0 ) { 
      if ( temper[lc] < temper[lc-1]) {
        if ( temper[lc] < temper[lc-2] ) {
          /* 2 layer above are warmer */
          if ( i_CPT_2 == NOT_DEFINED_INTEGER)
            i_CPT_2 = lc;
          if ( temper[lc] < temper[lc-3] ) {
            /* three layers above are warmer */
            if ( i_CPT_3 == NOT_DEFINED_INTEGER)
              i_CPT_3 = lc;
            if ( temper[lc] < temper[lc-4] ) {
              /* four layers above are warmer, that is OK in each case */
              i_CPT = lc;
	      //20120816ak i_layers is not use, commented
	      //              i_layers=4;
              break;
            }
          }
        }
      }
    }
  }

  /* if we did not find a good inversion with 3 warmer layers above we check for smaller inversions */
  if ( i_CPT == NOT_DEFINED_INTEGER ) {
    if ( i_CPT_3 != NOT_DEFINED_INTEGER ) {
      /* we take second best guess */
      i_CPT = i_CPT_3;
      //20120816ak i_layers is not use, commented
      //      i_layers=3;
    }
    else if ( i_CPT_2 != NOT_DEFINED_INTEGER ) {
      /* we take second best guess */
      i_CPT = i_CPT_2;
      //20120816ak i_layers is not use, commented
      //      i_layers=2;
    }
    else {
      if ( !quiet )
        fprintf (stderr, "*** Warning: no cold point tropopause (CPT) found in the profile\n");
      return -1;
    }
  }

  /*   fprintf (stderr, " ... minimum temperature in the profile: T[%3d]=%f ( T[%3d]=%f, T[%3d]=%f, T[%3d]=%f )\n",  */
  /*                          i_CPT, temper[i_CPT], i_CPT-1, temper[i_CPT-1], i_CPT-2, temper[i_CPT-2], i_CPT-3, temper[i_CPT-3]); */

  if ( temper[i_CPT-1] < temper[i_CPT+1] )
    di_CPT = -1;
  else
    di_CPT = +1;
  /*   fprintf (stderr, " !!! T_min[%3d]=%f (z=%7.3f); T_min2[%3d]=%f (z=%7.3f)\n\n", */
  /*                         i_CPT,        temper[i_CPT],       z[i_CPT], */
  /*                         i_CPT+di_CPT, temper[i_CPT+di_CPT],z[i_CPT+di_CPT]); */

  if (di_CPT > 0) {
    i1=i_CPT+2;
    i2=i_CPT+1;
    i3=i_CPT;
    i4=i_CPT-1;
  }
  else {
    i1=i_CPT+1;
    i2=i_CPT;
    i3=i_CPT-1;
    i4=i_CPT-2;
  }


  /* intercept point of the linear regression line 1 through i1 and i2 */
  /* and regression line 2 through i3 and i4                           */
  *z_CPT_sea       = ( (z[i2]*temper[i1] - z[i1]*temper[i2])/ (z[i2] - z[i1])
		       -(z[i4]*temper[i3] - z[i3]*temper[i4])/ (z[i4] - z[i3])
		       ) /
    ( (temper[i4] - temper[i3]) / (z[i4] - z[i3])
      -(temper[i2] - temper[i1]) / (z[i2] - z[i1]) );

  /* if the double extrapolated value is not inbetween the innermost 2 layers */
  if (  !(z[i2] <= *z_CPT_sea && *z_CPT_sea <= z[i3]) ) {
    /* than the interpolation does not make sence                                              */
    /* we simply use the lowermost temperature level itself instead of any interpolated values */
    *z_CPT_sea = z[i_CPT];
    *T_CPT     = temper[i_CPT];
    *p_CPT     = press[i_CPT];
  }
  else {
    /* otherwise we use the interpolated z-value and calculate T and p */
    /* calculate temperature of the cold point tropopause via linear extrapolation */
    *T_CPT =   temper[i2]         + ( *z_CPT_sea - z[i2] ) / ( z[i2] - z[i1] ) * (temper[i2] - temper[i1]);
    /* calculate pressure of the cold point tropopause with logarithmic interpolation */
    *p_CPT = exp ( log(press[i2]) + ( *z_CPT_sea - z[i2] ) / ( z[i3] - z[i2] ) * ( log(press[i3]) - log(press[i2] ) ) );
  }

  /*   /\* check tropopause height once more *\/ */
  /*   if ( !(z[i2] <= *z_CPT_sea && *z_CPT_sea <= z[i3]) ) { */
  /*     fprintf (stderr, "*** Warning: strange tropopause height %8.4f %8.4f %8.4f\n", z[i2], *z_CPT_sea, z[i3]); */
  /*     fprintf (stderr, "*** i_layers = %3d\n", i_layers); */
  /*     fprintf (stderr, "*** z_above_sea = %8.4f %8.4f %8.4f %8.4f\n", */
  /*                      z[i1],z[i2],z[i3],z[i4]); */
  /*     fprintf (stderr, "*** T           = %8.4f %8.4f %8.4f %8.4f\n", */
  /*                      temper[i1],temper[i2],temper[i3],temper[i4]); */
  /*   } */

  if (verbose) {
    fprintf (stderr, " ... found cold point tropopause:\n");
    fprintf (stderr, "        z_CPT_above_sea = %8.4f, p_CPT = %f,  T_CPT = %10.5f\n",
	     *z_CPT_sea, *p_CPT, *T_CPT );
  }

  return status;
}


/****************************************************************/
/* compare_float_inv                                            */
/* small function, which is used by qsort (as function pointer) */
/*  in order to sort data e.g. altitudes in descending order    */
/****************************************************************/
static int compare_float_inv (void *ap, void *bp) { 
  float *a = (float *) ap;
  float *b = (float *) bp;
  
  if (*a == *b) 
    return 0;

  if (*a > *b) 
    return 1;
  else 
    return -1;
}


/**************************************************************/
/* Setup optical properties of trace gases.                   */
/**************************************************************/

int setup_gases (input_struct input, output_struct *output)
{
  char function_name[]="setup_gases";
  char file_name[]="atmosphere.c";
  static int first=1;
  int status=0;

  double bscar=0;


  double babso=0, babso_md=0, babso_amf=0, babso_amftot=0, vertical_column=0;
  double ***babso_mol;
  double bscartot=0;
  double babsotot=0;
  double ***babsotot_mol;
  double deltaz=0;
  int iq=0, is=0, isthis=0, iv=0, lc=0, ix=0, iy=0; /*ipa=0 not used */
  int nlev = output->atm.nlyr+1;
  int nlyr = output->atm.nlyr;
    
  int iq1=0, iq2=0, iq3=0, iq4=0;
  int nq=0, nq1=0, nq2=0, nq3=0, nq4=0;

  double factor=0;
  double weight[KATO_COMPONENTS];

  float *fu_rayleigh = calloc (nlyr, sizeof(float));
  float density=-999.0;
  int gas_nr;


  if (first!=0 || output->molecular3d) {
    first = 0;
   
    /* allocate memory for profiles */
    status =
      ASCII_calloc_float_3D(&output->atm.cdento, output->atm.Nxatm, output->atm.Nyatm, nlyr);
    if (status!=0){ 
      fprintf (stderr, "Error %d returned by ASCII_calloc_float3D in %s (%s) allocating output->atm.cdento \n", status, function_name, file_name); 
      return status;
    }

    //    output->atm.cdento   = (float *) calloc (nlyr, sizeof(float));
  
    output->atm.wght_r   = calloc (output->wl.nlambda_r, sizeof (double *));
    output->atm.nq_r     = calloc (output->wl.nlambda_r, sizeof (int));

    output->mc_photons_r = calloc (output->wl.nlambda_r, sizeof (double));


    /* copy ck weights */
    for (iv=output->wl.nlambda_rte_lower; iv<=output->wl.nlambda_rte_upper; iv++) {
      
      switch (input.ck_scheme) {
      case CK_FU:
        output->atm.nq_r[iv] = output->ck.wght[FU_WGHT][iv+1];
        break;
        
      case CK_KATO:
      case CK_KATO2:
      case CK_KATO2_96:
      case CK_KATO2ANDWANDJI:
        output->atm.nq_r[iv] = output->ck.wght[KATO_CO2][iv+1] *
	  output->ck.wght[KATO_O3] [iv+1] *
	  output->ck.wght[KATO_O2] [iv+1] *
	  output->ck.wght[KATO_H2O][iv+1];
	break;

      case CK_AVHRR_KRATZ:
        output->atm.nq_r[iv] = output->ck.wght[AVHRR_KRATZ_WGHT][iv+1];
        break;

      case CK_FILE:
        output->atm.nq_r[iv] = output->ck.wght[GENERIC_WGHT][iv+1];
        break;

      case CK_LOWTRAN:
        output->atm.nq_r[iv] = LOWTRAN_MAXINT;
        /* this is the maximum number, but can be less */
        /* this is optimised inside the solve_rte() ipa-loop, UH 2006-07-19 */
        break;

      case CK_CRS:
      case CK_REPTRAN:
      case CK_REPTRAN_CHANNEL:
	output->atm.nq_r[iv] = 1;
	break;

      case CK_RAMAN:
	output->atm.nq_r[iv] = 1;
	break;
	
      default:
	fprintf (stderr, "Error, unimplemented mol_abs_param scheme %d\n", input.ck_scheme);
	return -1;
	break;
      }

      output->atm.wght_r [iv] = calloc (output->atm.nq_r[iv], sizeof (double));

      /* allocate memory for photon fractions */
      output->mc_photons_r [iv] = calloc (output->atm.nq_r[iv], sizeof (double));


      /* copy photon fractions */
      switch(input.ck_scheme) {
      case CK_CRS:
      case CK_REPTRAN:
      case CK_REPTRAN_CHANNEL:
      case CK_RAMAN:
        
        iq=0; /* no loop over quadrature points required */
        output->mc_photons_r [iv][iq] = output->wl.pfraction[iv];
	
	break;
	
      case CK_KATO:
      case CK_KATO2:
      case CK_KATO2_96:
      case CK_KATO2ANDWANDJI:
        nq1 = output->crs_ck.profile[KATO_CO2-1][iv].ngauss;
        nq2 = output->crs_ck.profile[KATO_O3-1 ][iv].ngauss;
        nq3 = output->crs_ck.profile[KATO_O2-1 ][iv].ngauss;
        nq4 = output->crs_ck.profile[KATO_H2O-1][iv].ngauss;

        iq=0;  /* reset quadrature counter */
        /* loop over quadrature points */
        for (iq1=0; iq1<nq1; iq1++)  /* CO2 */
          for (iq2=0; iq2<nq2; iq2++)  /* O3 */
            for (iq3=0; iq3<nq3; iq3++)  /* O2 */
              for (iq4=0; iq4<nq4; iq4++)  /* H2O */
                switch (input.source) {
                case SRC_SOLAR: /* solar source */
                case SRC_LIDAR: /* BCA */
                case SRC_BLITZ: /* BCA */
                  output->mc_photons_r [iv][iq++] = ((double *****) output->ck.x_solar)  [iv+1][iq1+1][iq2+1][iq3+1][iq4+1];
                  break;
                case SRC_THERMAL: /* thermal source */
                  output->mc_photons_r [iv][iq++] = ((double *****) output->ck.x_thermal)[iv+1][iq1+1][iq2+1][iq3+1][iq4+1];
                  break;
                default:
                  fprintf (stderr, "Error, unknown source %d\n", input.source);
                  return -1;
                }
        break;
      
      case CK_FU:
      case CK_AVHRR_KRATZ:
      case CK_FILE:
      case CK_LOWTRAN:
        
	nq = output->crs_ck.profile[0][iv].ngauss;
        
	/* Need all three for LOWTRAN!     */
	/* We decide later (in solve_rte)  */
	/* if we need one or three; in the */
	/* further case, the three weights */
	/* are added.                      */
	if (input.ck_scheme == CK_LOWTRAN)
	  nq = LOWTRAN_MAXINT;
	
        for (iq=0; iq<nq; iq++)
          switch (input.source) {
          case SRC_SOLAR: /* solar source */
          case SRC_LIDAR: /* BCA */
          case SRC_BLITZ: /* BCA */
            output->mc_photons_r [iv][iq] = ((double **) output->ck.x_solar)  [iv+1][iq+1];
            break;
          case SRC_THERMAL: /* thermal source */
            output->mc_photons_r [iv][iq] = ((double **) output->ck.x_thermal)[iv+1][iq+1];
            break;
          default:
            fprintf (stderr, "Error, unknown source %d\n", input.source);
            return -1;
          }
        break;
        
      default:
        fprintf (stderr, "Error, unimplemented mol_abs_param scheme %d\n", input.ck_scheme);
        return -1;
        break;
      } /* end switch */
      
    } /* end wavelength loop */
    
    //fprintf(stderr, "3DAbs setup_gases() alloc nlyr %d \n", nlyr); 
    /* Rayleigh scattering optical thickness */
    if ((status = ASCII_calloc_double_5D_arylen (&output->atm.optprop.tau_rayleigh_r, 
						 output->atm.Nxatm, output->atm.Nyatm, 
						 nlyr, output->wl.nlambda_r, 
						 output->atm.nq_r)) != 0){
      fprintf (stderr, "Error %d returned by ASCII_calloc_double_5D_arylen in %s (%s) allocating output->atm.optprop.tau_rayleigh_r\n",
	       status, function_name, file_name); 
      return status;
    }
    /* Raman scattering optical thickness */
    /*    if ((status = ASCII_calloc_double_3D_arylen (&output->atm.optprop.tau_raman_t, 
	  nlyr, output->wl.nlambda_r, 
	  output->atm.nq_t)) != 0)  {
	  if (status!=0) {
	  fprintf (stderr, "Error %d returned by ASCII_calloc_double_3D_arylen in %s (%s)\n", status, function_name, file_name); 
	  return status;
	  }
	  }*/
    
    /* Molecular absorption optical thickness */
    //fprintf(stderr, "3DAbs setup_gases nlyr %d \n", nlyr);  
    if ((status = ASCII_calloc_double_5D_arylen (&output->atm.optprop.tau_molabs_r, 
						 output->atm.Nxatm, output->atm.Nyatm, 
						 nlyr, output->wl.nlambda_r, 
						 output->atm.nq_r)) != 0)   {
      fprintf (stderr, "Error %d returned by ASCII_calloc_double_5D_arylen in %s (%s) allocating output->atm.optprop.tau_molabs_r\n",
	       status, function_name, file_name); 
      return status;
    }
    
      
    /* Molecular absorption optical thickness for AMF calculations */
    if ((status = ASCII_calloc_double_5D_arylen (&output->atm.optprop.tau_molabs_md_r, 
						 output->atm.Nxatm, output->atm.Nyatm, 
						 nlyr, output->wl.nlambda_r, 
						 output->atm.nq_r)) != 0) {
      fprintf (stderr, "Error %d returned by ASCII_calloc_double_5D_arylen in %s (%s) allocating output->atm.optprop.tau_molabs_md_r\n",
	       status, function_name, file_name); 
      return status;
    } 
    
    /* Calculate column density of air (total density of the atmosphere). */
    /* Required to calculate Rayleigh optical depth.                      */

    if (output->atm.rayleigh==RAYLEIGH_CALC) {
      for (ix=0; ix<output->atm.Nxatm; ix++){
	for (iy=0; iy<output->atm.Nyatm; iy++){
	  for (lc=0; lc<nlyr; lc++) {
	    
	    /* The factor 1e5 converts z from km to cm.*/
	    deltaz = (output->atm.zd[lc] - output->atm.zd[lc+1]) * 1e5;
	    
	    /* average air density over a layer assuming  dens = exp (-c * z) */
	    output->atm.cdento[ix][iy][lc] = log_average
	      (output->atm.microphys.dens[MOL_AIR][ix][iy][lc],
	       output->atm.microphys.dens[MOL_AIR][ix][iy][lc+1]) * deltaz;
	  }
	}
      }
    }

  } /* End if first block */


  if (output->atm.microphys.denstab_id) {
    /* Find matrix profile corresponding to the sza this rte calculation is done for */
    /* This is airmass factor stuff */
    iv=0;  /* No solar zenith variations for a single amf calculation */
    for (is=0; is<output->atm.microphys.nsza_denstab; is++)
      if (fabs (output->atm.microphys.sza_denstab[is]-output->atm.sza_r[iv]) < 1.e-04)
	isthis = is;

    for (lc=0; lc<nlyr; lc++) { 
      output->atm.microphys.dens[output->atm.microphys.denstab_id][0][0][lc] = 
	output->atm.microphys.denstab[output->atm.microphys.denstab_id][isthis][lc];
    }
  }
  
  
  /***********************/
  /* Rayleigh scattering */
  /***********************/

  /* loop over wavelength */
  for (iv=output->wl.nlambda_rte_lower; iv<=output->wl.nlambda_rte_upper; iv++) {

    for (iq=0; iq<output->atm.nq_r[iv]; iq++) {

      /* loop over grid cells */
      for (ix=0; ix<output->atm.Nxatm; ix++){
	
	for (iy=0; iy<output->atm.Nyatm; iy++){
	  
	  
	  /* FU special: The Rayleigh optical depth depends on the subband */
	  if (input.ck_scheme == CK_FU) {
#if HAVE_FULIOU
	    status = ckdfuray (iv+1, iq+1, cos(output->atm.sza_r[iv]*PI/180.0), nlev,
			       fu_rayleigh);
	    
	    if (status!=0) {
	      fprintf (stderr, "Error, status %d returned by ckdfuray()\n", status);
	      return status;
	    }    
	  }
#else
     
	  fprintf (stderr, "Error, Fu and Liou not supported!\n");
	  return -1;
#endif  
	  
	  
	  
	  /* loop over layers */
	  for (lc=0; lc<nlyr; lc++) {
	    
	    switch (output->atm.rayleigh) {  /* if Rayleigh scattering is not switched off */
	    case RAYLEIGH_CALC:
	      if (input.ck_scheme != CK_FU)
		/* Rayleigh cross section depends only on wavelength, therefore ix=0, iy=0, iz=0 */
		bscar = output->atm.cdento[ix][iy][lc] * output->crs.crs[0][0][0][MOL_AIR][iv];
	      else 
		bscar = fu_rayleigh[lc];
	      break;
	  
	    case RAYLEIGH_FILE_MONO:
	      bscar = output->atm.optprop.tau_rayleigh_user[lc];
	      break;
	  
	    case RAYLEIGH_FILE_SPEC:
	      bscar = output->atm.optprop.tau_rayleigh_user_r[iv][lc];
	      break;
	  
	    case RAYLEIGH_NONE:
	      bscar = 0;
	      break;
	  
	    default:
	      fprintf (stderr, "Error, unknown Rayleigh option %d\n", output->atm.rayleigh);
	      return -1;
	    }
	    output->atm.optprop.tau_rayleigh_r[ix][iy][lc][iv][iq] = bscar;
	  }
	}
      }
    }
  } /* End wavelength loop */
  

  /***********************/
  /* Absorption by gases */
  /***********************/

  /* babso_mol = calloc (MOL_NN, sizeof(double)); */
  /* babsotot_mol = calloc (MOL_NN, sizeof(double)); */

  if ((status = ASCII_calloc_double_3D(&babso_mol, MOL_NN, 
				       output->atm.Nxatm, output->atm.Nyatm)) != 0) {
    
    fprintf (stderr, "Error %d returned by ASCII_calloc_double_3D in %s (%s) allocating babso_mol\n",
	     status, function_name, file_name); 
    return status;
  }
  
  if ((status = ASCII_calloc_double_3D(&babsotot_mol, MOL_NN, 
				       output->atm.Nxatm, output->atm.Nyatm)) != 0) {
    
    fprintf (stderr, "Error %d returned by ASCII_calloc_double_3D in %s (%s) allocating babsotot_mol\n",
	     status, function_name, file_name); 
    return status;
  }
  
 
  
  /* loop for indipendent columns */
  // CE ??? Nothing is changed in ipa loop, index ipa not used here !!
  //for (ipa=0; ipa<output->nipa; ipa++) {
    
  /* loop over wavelength */
  for (iv=output->wl.nlambda_rte_lower; iv<=output->wl.nlambda_rte_upper; iv++) {
    
    switch(input.ck_scheme) {
    case CK_CRS:
    case CK_RAMAN:
    case CK_REPTRAN:
    case CK_REPTRAN_CHANNEL:
     
      for (ix=0; ix< output->atm.Nxatm; ix++){ 
	for (iy=0; iy< output->atm.Nyatm; iy++){ 
	  iq=0; /* no loop over quadrature points required */
        
	  bscartot      = 0;
	  
	  babsotot      = 0;

	  for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
	    babsotot_mol[gas_nr][ix][iy]=0;
	  
	  babso_amftot  = 0;
	  
	  vertical_column=0;
     
	  /* loop over layers */
	  for (lc=0; lc<nlyr; lc++) {
	    
	    /* The factor 1.e+5 converts z from km to cm.*/
	    deltaz = (output->atm.zd[lc] - output->atm.zd[lc+1])*1e5;
	    
	    babso      = 0;
	    for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
	      babso_mol[gas_nr][ix][iy]=0;

	    babso_amf  = 0;
	    babso_md   = 0;
           
	    switch (output->atm.molabs) {
	    case MOLABS_CALC: 
	    case MOLABS_LOOKUP: 
	      //3DAbs not yet included for moltau file.
	      
	      for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) {
		
		density = output->atm.microphys.dens_avg[gas_nr][ix][iy][lc];
              
		babso_mol[gas_nr][ix][iy] += deltaz * 0.5 * (output->crs.crs[ix][iy][lc][gas_nr][iv] + output->crs.crs[ix][iy][lc+1][gas_nr][iv]) * density;
		// babso_mol[gas_nr] += deltaz * 0.5 * (output->crs.crs[lc][gas_nr][iv]*output->atm.microphys.dens[gas_nr][lc] + output->crs.crs[lc+1][gas_nr][iv]*output->atm.microphys.dens[gas_nr][lc+1]); // This line would give the same as WriteMolTau() in the ARTS radiative transfer package.
		babso += babso_mol[gas_nr][ix][iy];

		if (output->atm.microphys.denstab_id != gas_nr)  
		  babso_md += babso_mol[gas_nr][ix][iy];
           
		if (output->atm.microphys.denstab_id == gas_nr) {
		  /* This is airmass factor stuff */
		  /* Use column for matrix profile corresponding to the SZA this RTE calculation is done for */
		  switch (input.atm.interpol_method_gas[gas_nr]) {
		  case INTERP_METHOD_LINEAR:
		    density =         0.5 * (output->atm.microphys.denstab[gas_nr][isthis][lc] + 
					     output->atm.microphys.denstab[gas_nr][isthis][lc+1]);
		    break;
		  case INTERP_METHOD_LOG:
		    density = log_average   (output->atm.microphys.denstab[gas_nr][isthis][lc],
					     output->atm.microphys.denstab[gas_nr][isthis][lc+1]);
		    break;
		  case INTERP_METHOD_LINMIX:
		  case INTERP_METHOD_SPLINE:
		  case INTERP_METHOD_LOG_SPLINE: 
		    density = linmix_average(output->atm.microphys.denstab[gas_nr][isthis][lc],
					     output->atm.microphys.denstab[gas_nr][isthis][lc+1],
					     output->atm.microphys.dens[MOL_AIR][ix][iy][lc],
					     output->atm.microphys.dens[MOL_AIR][ix][iy][lc+1]);
		    break;
		  default:
		    fprintf (stderr, "Error, unkown z_interpol_gas option for denstab calculation\n");
		    return -1;
		    break;
		  }
		  
		  babso_amf += deltaz * 0.5 * (output->crs.crs[ix][iy][lc][gas_nr][iv] +
					       output->crs.crs[ix][iy][lc+1][gas_nr][iv]) *
		    density;
		  vertical_column += deltaz * density;
		}
	      }
            
	      output->atm.optprop.tau_molabs_r [ix][iy][lc][iv][iq] = babso;
	      if ( input.rte.solver == SOLVER_SDISORT) 
		output->atm.optprop.tau_molabs_md_r[ix][iy][lc][iv][iq] = babso_md;
	      break;
            
	    case MOLABS_FILE_MONO:
	      //3DAbs, user not yet included 
	      babso = output->atm.optprop.tau_molabs_user[lc];
	      output->atm.optprop.tau_molabs_r [ix][iy][lc][iv][iq] = babso;
	      break;
	  
	    case MOLABS_FILE_SPEC:
	      babso = output->atm.optprop.tau_molabs_user_r[iv][lc];
	      output->atm.optprop.tau_molabs_r [ix][iy][lc][iv][iq] = babso;
	      break;
            
	    case MOLABS_NONE:
	      babso = 0;
	      output->atm.optprop.tau_molabs_r [ix][iy][lc][iv][iq] = babso;
	      break;
             
	    default:
	      fprintf (stderr, "Error, unknown molecular absorption option %d\n", output->atm.molabs);
	      return -1;
	    }
	    
	    output->atm.wght_r[iv][iq] = 1;  /* default weight */
	    
           
	    if (input.verbose && !output->molecular3d) {
	      bscar=output->atm.optprop.tau_rayleigh_r[ix][iy][lc][iv][iq];
	      bscartot+=bscar;
          
	      babsotot      += babso;
	      
	      for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
		babsotot_mol[gas_nr][ix][iy]   += babso_mol[gas_nr][ix][iy];

	      // AMF not included for 3D trace gases
	      babso_amftot  += babso_amf;
          
	      if (lc == 0) {
		
		fprintf (stderr, "\n*** wavelength: iv = %d, %f nm\n", iv, output->wl.lambda_r[iv]);
		if (output->wl.lambda_r[iv] > 500.0 && output->atm.molabs==MOLABS_CALC)
		  fprintf (stderr, "\n*** WARNING: water vapour and ozone absorption are inexact for wavelength above 500nm!\n");
		fprintf(stderr,"*** Rayleigh scattering cross section: %e\n", output->crs.crs[ix][iy][0][MOL_AIR][iv]);
		fprintf(stderr,"*** setup_gases(), layer properties, layer range = from z-level written to the z-level above\n");
              
		fprintf(stderr, " ----------------------------------------------");
		for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
		  fprintf(stderr, "--------------");
		fprintf(stderr, "\n");
		
		fprintf(stderr, "   lc |      z[km] |  Rayleigh   |  Molecular\n");
		
		fprintf(stderr, "      |            |    dtau     |  absorption ");
		for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
		  fprintf(stderr, "|     %s    ", gas_number2string(gas_nr));
		fprintf(stderr, "\n");
		
		fprintf(stderr, " ----------------------------------------------");
		for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
		  fprintf(stderr, "--------------");
		fprintf(stderr, "\n");
		
	      }
	      
	      if (output->atm.zd[lc+1] >= output->alt.altitude) {
		fprintf (stderr, "%5d | %10.4f | %10.5e | %10.5e", lc, output->atm.zd[lc+1], bscar, babso);
		for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
		  fprintf (stderr, " | %10.5e", babso_mol[gas_nr][ix][iy]);
		fprintf(stderr, "\n");
	      }
	      
	    }
	  
	  } /* End layer loop */
	  
	  if (input.verbose) {
	    
	    fprintf(stderr, " ----------------------------------------------");
	    for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
	      fprintf(stderr, "--------------");
	    fprintf(stderr, "\n");
	    
	    fprintf (stderr, "%5s | %10.4f | %10.5e | %10.5e", "sum", 0.0/0.0, bscartot, babsotot);
	    for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
	      fprintf (stderr, " | %10.5e", babsotot_mol[gas_nr][ix][iy]);
	    fprintf(stderr, "\n");
	    
	    fprintf(stderr, " ----------------------------------------------");
	    for (gas_nr = 1; gas_nr < MOL_NN; gas_nr++) 
	      fprintf(stderr, "--------------");
	    fprintf(stderr, "\n");
	    
	    fprintf(stderr,"\n");
	  }
	} /* end iy loop*/
      } /* end ix loop*/

      break;
        
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
      nq1 = output->crs_ck.profile[KATO_CO2-1][iv].ngauss;
      nq2 = output->crs_ck.profile[KATO_O3-1 ][iv].ngauss;
      nq3 = output->crs_ck.profile[KATO_O2-1 ][iv].ngauss;
      nq4 = output->crs_ck.profile[KATO_H2O-1][iv].ngauss;
      
      if (input.verbose) {
	fprintf (stderr, "*** subbands:\n");
	fprintf (stderr, "***  CO2: %2d\n", nq1);
	fprintf (stderr, "***  O3:  %2d\n", nq2);
	fprintf (stderr, "***  O2:  %2d\n", nq3);
	fprintf (stderr, "***  H2O: %2d\n", nq4);
	fprintf (stderr, "\n");
      }
      
      for (ix=0; ix<output->atm.Nxatm; ix++){
	for (iy=0; iy<output->atm.Nyatm; iy++){
	  for (lc=0; lc<output->atm.nlyr; lc++) {
	    
	    /* The factor 1e5 converts z from km to cm.*/
	    deltaz = (output->atm.zd[lc] - output->atm.zd[lc+1]) * 1e5;
        
	    iq=0;  /* reset quadrature counter */
	    /* loop over components (4D correlated-k) */
	    for (iq1=0; iq1<nq1; iq1++) {  /* CO2 */
            
	      /* set CO2 cross section and weights */
	      output->crs.crs[ix][iy][lc][MOL_CO2][iv] = 
		output->crs_ck.profile[KATO_CO2-1][iv].crs[ix][iy][lc][iq1]; 
          
	      weight[KATO_CO2-1] = output->crs_ck.profile[KATO_CO2-1][iv].weight[iq1];
          
	      for (iq2=0; iq2<nq2; iq2++) { /* O3 */
            
		/* set O3 cross section and weights */
		output->crs.crs[ix][iy][lc][MOL_O3][iv]  = 
		  output->crs_ck.profile[KATO_O3-1][iv].crs[ix][iy][lc][iq2];
            
		weight[KATO_O3-1] = output->crs_ck.profile[KATO_O3-1][iv].weight[iq2];
            
		for (iq3=0; iq3<nq3; iq3++) { /* O2 */
              
		  /* set O2 cross section and weights */
		  output->crs.crs[ix][iy][lc][MOL_O2][iv]  =
		    output->crs_ck.profile[KATO_O2-1][iv].crs[ix][iy][lc][iq3];
              
		  weight[KATO_O2-1] = output->crs_ck.profile[MOL_O2][iv].weight[iq3];
              
		  for (iq4=0; iq4<nq4; iq4++) { /* H2O */
        	
		    /* set H2O cross section and weights */
		    output->crs.crs[ix][iy][lc][MOL_H2O][iv]  = 
		      output->crs_ck.profile[KATO_H2O-1][iv].crs[ix][iy][lc][iq4];
		    
		    weight[KATO_H2O-1] = output->crs_ck.profile[MOL_H2O][iv].weight[iq4];
        	
		    /* quadrature weight */
		    factor = weight[KATO_CO2-1]*weight[KATO_O3-1]*weight[KATO_O2-1]*weight[KATO_H2O-1];
                
		    /* sum absorption */
		    babso     = 0;

		    /* O3 absorption */
		    babso += deltaz * output->crs.crs[ix][iy][lc][MOL_O3][iv]  * output->atm.microphys.dens_avg[MOL_O3][ix][iy][lc];
		    /* babso += deltaz * output->crs.crs[lc][MOL_O3][iv]  *                                                 */
		    /*         0.5 * (output->atm.microphys.dens[MOL_O3][lc]  + output->atm.microphys.dens[MOL_O3][lc+1]);  */   /* before 2006-03-06, UH */

		    /* O2 absorption */
		    babso += deltaz * output->crs.crs[ix][iy][lc][MOL_O2][iv]  * output->atm.microphys.dens_avg[MOL_O2][ix][iy][lc]; 
		    /* babso += deltaz * output->crs.crs[lc][MOL_O2][iv]  *                                                 */                              
		    /*         0.5 * (output->atm.microphys.dens[MOL_O2][lc]  + output->atm.microphys.dens[MOL_O2][lc+1]);  */   /* before 2006-03-06, UH */

		    /* H2O absorption */
		    babso += deltaz * output->crs.crs[ix][iy][lc][MOL_H2O][iv] * output->atm.microphys.dens_avg[MOL_H2O][ix][iy][lc]; 
		    /* babso += deltaz * output->crs.crs[lc][MOL_H2O][iv] *                                                 */                              
		    /*         0.5 * (output->atm.microphys.dens[MOL_H2O][lc] + output->atm.microphys.dens[MOL_H2O][lc+1]); */   /* before 2006-03-06, UH */

		    /* CO2 absorption */
		    babso += deltaz * output->crs.crs[ix][iy][lc][MOL_CO2][iv] * output->atm.microphys.dens_avg[MOL_CO2][ix][iy][lc]; 
		    /* babso += deltaz * output->crs.crs[lc][MOL_CO2][iv] *                                                 */
		    /*         0.5 * (output->atm.microphys.dens[MOL_CO2][lc] + output->atm.microphys.dens[MOL_CO2][lc+1]); */   /* before 2006-03-06, UH */

		    /* copy absorption cross section and weight */
		    /* to final destination                     */
		    output->atm.wght_r[iv][iq] = factor;
                
		    output->atm.optprop.tau_molabs_r[ix][iy][lc][iv][iq]   = babso;

		    iq++;
		    
		  }  /* end of correlated-k (iq4) loop */
		}  /* end of correlated-k (iq3) loop */
	      }  /* end of correlated-k (iq2) loop */
	    }  /* end of correlated-k (iq1) loop */
	  } /* End layer loop */

	  if (input.verbose && !output->molecular3d) {
	    iq=0;
	    for (iq1=0; iq1<nq1; iq1++) { /* CO2 */
	      for (iq2=0; iq2<nq2; iq2++) { /* O3 */
		for (iq3=0; iq3<nq3; iq3++) { /* O2 */
		  for (iq4=0; iq4<nq4; iq4++) { /* H2O */
		    bscartot = 0.0;
		    babsotot = 0.0;

		    fprintf (stderr, "\n*** wavelength: iv = %d, lambda=%f nm, iq1(CO2)=%4d, iq2(O3)=%4d, iq3(O2)=%4d, iq4(H2O)=%4d\n",
			     iv, output->wl.lambda_r[iv], iq1, iq2, iq3, iq4);
		    fprintf(stderr,"*** setup_gases(), layer properties, layer range = from z-level written to the z-level above\n");
		    fprintf(stderr, " --------------------------------------------\n");
		    fprintf(stderr, "   lc |    z[km] |  Rayleigh   |   Molecular \n");
		    fprintf(stderr, "      |          |    dtau     |   absorption\n");
		    fprintf(stderr, " --------------------------------------------\n");

		    for (lc=0; lc<output->atm.nlyr; lc++) {
		      bscartot += output->atm.optprop.tau_rayleigh_r[ix][iy][lc][iv][iq];
		      babsotot += output->atm.optprop.tau_molabs_r  [ix][iy][lc][iv][iq];
		      if (output->atm.zd[lc+1] >= output->alt.altitude) {
			fprintf (stderr, "%5d | %8.2f | %11.6f | %11.6f |\n", lc, output->atm.zd[lc+1], 
				 output->atm.optprop.tau_rayleigh_r[ix][iy][lc][iv][iq], output->atm.optprop.tau_molabs_r[ix][iy][lc][iv][iq]);
		      }
		    }
		    fprintf(stderr, " --------------------------------------------\n");
		    fprintf (stderr, "%5s | %8.2f | %11.6f | %11.6f |\n", 
			     "sum", 0.0/0.0, bscartot, babsotot);
		    fprintf(stderr, " --------------------------------------------\n");
		    iq++;
		  }
		}
	      }
	    }
	  }
      	} /* end ix loop */ 
      } /* end iy loop */
      break;
    case CK_FU:
    case CK_AVHRR_KRATZ:
    case CK_LOWTRAN:
    case CK_FILE:
      
      nq = output->crs_ck.profile[0][iv].ngauss;
      if( input.rte.mc.spectral_is )
	nq = LOWTRAN_MAXINT;
        
      /* loop over quadrature points */
      
      for (iq=0; iq<nq; iq++) {
        
	for (ix=0; ix<output->atm.Nxatm; ix++){
	  for (iy=0; iy<output->atm.Nyatm; iy++){
	    /* set absorption coefficient */
	    for (lc=0; lc<output->atm.nlyr; lc++) {
	      
	      /* The factor 1e5 converts z from km to cm.*/
	      deltaz = (output->atm.zd[lc] - output->atm.zd[lc+1]) * 1e5;
      
	      /* crs_ck.profile are calculated in function setup_crs (in molecular.c) */    
	      output->atm.optprop.tau_molabs_r[ix][iy][lc][iv][iq] =
		deltaz/1e5 * output->crs_ck.profile[0][iv].crs[ix][iy][lc][iq];
	    }
          
	    /* quadrature weight */
	    factor = output->crs_ck.profile[0][iv].weight[iq];
	    output->atm.wght_r[iv][iq] = factor;
	  }
	}
      }
      
      break;
      
    default:
      fprintf (stderr, "Error, unimplemented mol_abs_param scheme %d\n", input.ck_scheme);
      return -1;
      break;
    }
  } /* End wavelength loop */
  
    //} /* End of indipendent pixel (ipa) loop */
  
  free(babso_mol);
  free(babsotot_mol);
  free(fu_rayleigh);
  
  return 0;  /* ok */
}



/**************************************************************/
/* Setup temperature at output levels.                        */
/**************************************************************/

int setup_temperature (input_struct input, output_struct *output)
{
  int status=0, lev=0, iz=0;
  float *zout=NULL;
  
  /* interpolate output temperature to zout */
  
  output->atm.temper_out  = (float *) calloc (output->atm.nlev, sizeof(float));
  
  /* copy temperature profile */
  for (lev=0; lev<output->atm.nlev; lev++) {
    output->atm.temper_out[lev] = output->atm.microphys.temper [0][0][lev];
  }
 
  /* Now we need to define a temporary zout grid;        */
  /* for rte_solver montecarlo it is possible to         */
  /* define an output altitude of -999 which stands      */
  /* for the surface altitude which may be a function    */
  /* of (x,y). This value is replaced by 0 which means   */
  /* that the surface temperature at the surface output  */
  /* level is set to temperature(z=0) which is not quite */
  /* correct if a 2D altitude grid has been defined.     */
  /* However, the only place where this information is   */
  /* currently used is the flexstor header which is      */
  /* required for photolysis calculations.               */

  zout = calloc (output->atm.nzout, sizeof(float));
  for (iz=0; iz<output->atm.nzout; iz++) {
    if (output->atm.zout_sea[iz] != ZOUT_MYSTIC && input.rte.solver!=SOLVER_SSLIDAR)
      zout[iz] = output->atm.zout_sea[iz];
    else {
      /* for MYSTIC surface marker set surface to 0 */
      /* would here not lowest level in the atmosphere file output->atm.zd[output->atm.nlev-1] be better than 0.0? UH 2008-02 */
      zout[iz] = 0.0; 
    }
  }
  
  /* interpolate zout_temperature to zout-grid with variable interpolation method */
  status = interpolate_profile (output->atm.zd , &(output->atm.temper_out), output->atm.nlev,
				zout, output->atm.nzout, input.atm.interpol_method_temper, 
				input.quiet);
  
  if (status != 0) {
    fprintf (stderr, "Error interpolating temperature\n");
    return status;
  }
  
  free(zout);


  /* Finally, copy the temperature of the lowest level to */
  /* output->temperature if the latter is not already set */
  /* ATTENTION: moved this code to redistribute.c because */
  /* we want the surface temperature to equal the         */
  /* z=altitude rather than z=0                           */

  /*
    if (output->surface_temperature < 0)
    output->surface_temperature = output->atm.microphys.temper[output->atm.nlev-1];
  */
  
  return 0;
}

 
 
 
/**************************************************************/
/* Read atmosphere file and write to atm_out_struct           */
/**************************************************************/
 
static int read_atmosphere (char *filename, atm_out_struct *out, 
			    int cdf, int quiet, int verbose) 
  
{
  int rows=0, min_columns=0, max_columns=0;
  int lev=0, status = 0;
  float **atm_data=NULL;
  
  int nlev=0;

  float tmpdens=0, idealdiff=0;
  int idealgas=0, firstnonideal=0;
  float *z_atm=NULL;
  
#if HAVE_LIBNETCDF
  int i=0;
  int ncid=0;
  int idd_nlev=0;
  int id_z=0, id_press=0, id_temper=0;
  int id_air=0, id_o3=0, id_o2=0, id_h2o=0, id_co2=0;
  
  size_t start[1] = {0};
  size_t count[1] = {0};
  size_t dimlen=0;
#endif
  
  if (verbose)
    fprintf (stderr, " ... reading atmosphere data file %s\n", filename);
  
  if (cdf==0) { /* classic ASCII atmospheric model */
    status = ASCII_file2float (filename, 
			       &rows, &max_columns, &min_columns, 
			       &atm_data);

    if (status!=0) {
      fprintf (stderr, "Error %d reading '%s'\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }
    
    if (max_columns!=min_columns) {
      fprintf (stderr, "Error, inconsistent number of columns in '%s'\n", filename);
      fprintf (stderr, "       (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      fprintf (stderr, "       min = %d, max =%d\n", min_columns, max_columns);
      return -1;
    }
    
    if (min_columns<3) {
      fprintf (stderr, "Error, too few columns in %s! Need at least %d.\n", filename, 3);
      fprintf (stderr, "       (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return -1;
    }
    
    nlev = rows;
  
    /* read data from atmosphere file into microphys.struct */
    if ( out->microphys.dens == NULL){
      
      /* allocate memory for density array */
      out->microphys.dens = calloc(MOL_NN, sizeof(float ***));
      if (out->microphys.dens == NULL) {
	fprintf (stderr, "Error allocating memory for 'dens'\n");
	fprintf (stderr, "       (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
	return -1;
      }
     
      if ( ASCII_calloc_float_3D (&out->microphys.press, 1, 1, nlev ) <0 ){
	fprintf( stderr, "Error allocating memory for pressure. \n");
	return -1; 
      }
      
      if ( ASCII_calloc_float_3D (&out->microphys.temper, 1, 1, nlev ) <0 ){
	fprintf( stderr, "Error allocating memory for pressure. \n");
	return -1; 
      }
    
      for (lev=0; lev<nlev; lev++){
	out->microphys.press[0][0][lev]  = atm_data[lev][ATMFILE_PRESS];
	out->microphys.temper[0][0][lev]  = atm_data[lev][ATMFILE_TEMPER];
      }
      
    }
    
    z_atm = ASCII_column_float (atm_data, nlev, ATMFILE_ZD);

    //fprintf(stderr, "3DAbs read_atm nlev %d \n", nlev);
    
    /**************************************/
    /* the following columns are optional */
    /**************************************/
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_AIR], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_O3], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_O2], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_H2O], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_CO2], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_NO2], 1, 1, nlev);
    
    for (lev=0; lev<nlev; lev++){
      if (min_columns>3) { /* air */
	out->microphys.dens[MOL_AIR][0][0][lev] = atm_data[lev][ATMFILE_AIR]; 
      }
      else{
	out->microphys.dens[MOL_AIR][0][0][lev] = out->microphys.press[0][0][lev]*1E-04/(BOLTZMANN * out->microphys.temper[0][0][lev]);
      }
      if (min_columns>4) /* o3   */
	out->microphys.dens[MOL_O3][0][0][lev]=atm_data[lev][ATMFILE_O3]; 
      if (min_columns>5) /* oxygen   */
	out->microphys.dens[MOL_O2][0][0][lev]=atm_data[lev][ATMFILE_O2];
      if (min_columns>6)  /* water vapour     */
	out->microphys.dens[MOL_H2O][0][0][lev]=atm_data[lev][ATMFILE_H2O];
      if (min_columns>7)  /* carbon dioxide   */
	out->microphys.dens[MOL_CO2][0][0][lev]=atm_data[lev][ATMFILE_CO2];
      if (min_columns>8)  /* nitrogen dioxide */
	out->microphys.dens[MOL_NO2][0][0][lev]=atm_data[lev][ATMFILE_NO2];
    }

    /* Allocate memory for species not in atmosphere file */
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_BRO], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_OCLO], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_HCHO], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_O4], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_SO2], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_CH4], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_N2O], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_CO], 1, 1, nlev);
    ASCII_calloc_float_3D(&out->microphys.dens[MOL_N2], 1, 1, nlev);

    ASCII_free_float(atm_data, rows);
        
  }
  
  
  else {

#if HAVE_LIBNETCDF

    /* open netcdf file */
    status = nc_open (filename, NC_NOWRITE, &ncid);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d opening netCDF file %s\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }

    /* get dimension id for "nlev" */
    status = nc_inq_dimid (ncid, "nlev", &idd_nlev);
  
    /* get dimension length for "nlev" */
    status = nc_inq_dimlen (ncid, idd_nlev, &dimlen);
    nlev = dimlen;

    start[0] = 0;
    count[0] = nlev;

    /* allocate memory for profiles */
    out->microphys.press  = calloc (nlev, sizeof (float));
    out->microphys.temper = calloc (nlev, sizeof (float));
    z_atm                 = calloc (nlev, sizeof (float));

    /* get variable id for "z" */
    status = nc_inq_varid (ncid, "z", &id_z);
    
    /* read altitude "z" */
    status = nc_get_vara_float (ncid, id_z, start, count, z_atm);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading z from '%s'\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }



    /* get variable id for "pressure" */
    status = nc_inq_varid (ncid, "pressure", &id_press);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error '%s' while getting id pressure from %s\n", nc_strerror(status), filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }
    
    /* read "pressure" */
    status = nc_get_vara_float (ncid, id_press, start, count, out->microphys.press[0][0]);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error '%s' while reading pressure from %s\n", nc_strerror(status), filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }


    /* get variable id for "temperature" */
    status = nc_inq_varid (ncid, "temperature", &id_temper);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error '%s' while getting id temperature from %s\n", nc_strerror(status), filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }
    
    /* read "temperature" */
    status = nc_get_vara_float (ncid, id_temper, start, count, out->microphys.temper[0][0]);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error '%s' while reading temperature from %s\n", nc_strerror(status), filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }


    /* allocate memory for density array */
    out->microphys.dens = calloc(MOL_NN, sizeof(float *));
    if (out->microphys.dens == NULL) {
      fprintf (stderr, "Error allocating memory for 'dens' \n");
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return -1;
    }

    for (i=0; i<MOL_NN; i++) {
      out->microphys.dens[i] = calloc(nlev, sizeof(float));
      if (out->microphys.dens[i] == NULL) {
  	fprintf (stderr, "Error allocating memory for dens[%d]\n", i);
	fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
  	return -1;
      }
    }


    /* get variable id for "density" */
    status = nc_inq_varid (ncid, "density", &id_air);
    
    /* read "density" */
    status = nc_get_vara_float (ncid, id_air, start, count, out->microphys.dens[MOL_AIR][0][0]);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading density from %s\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }


    /* get variable id for "o3" */
    status = nc_inq_varid (ncid, "o3", &id_o3);
    
    /* read "o3" */
    status = nc_get_vara_float (ncid, id_o3, start, count, out->microphys.dens[MOL_O3][0][0]);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading o3 from %s\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }


    /* get variable id for "o2" */
    status = nc_inq_varid (ncid, "o2", &id_o2);
    
    /* read "o2" */
    status = nc_get_vara_float (ncid, id_o2, start, count, out->microphys.dens[MOL_O2][0][0]);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading o2 from %s\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }


    /* get variable id for "h2o" */
    status = nc_inq_varid (ncid, "h2o", &id_h2o);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading h2o from %s\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }
    
    /* read "h2o" */
    status = nc_get_vara_float (ncid, id_h2o, start, count, out->microphys.dens[MOL_H2O][0][0]);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading h2o from %s\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }


    /* get variable id for "co2" */
    status = nc_inq_varid (ncid, "co2", &id_co2);
    
    /* read "co2" */
    status = nc_get_vara_float (ncid, id_co2, start, count, out->microphys.dens[MOL_CO2][0][0]);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading co2 from %s\n", status, filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }


    /* we don't have no2 in the netcdf file, hence we keep it 0;        */
    /* also we don't read n2o, co, and ch4 because we never used        */
    /* those before (should change that by introducing all relevant     */
    /* gases into dens; should also get rid of ndens  */
    /* and replace that by MOL_NN)                                      */

    nc_close (ncid);
    
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use the generic mol_abs_param option. Please get netcdf and rebuild. *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
    
  }

  /* it is not allowed, that any level is added to out->zd before this point */
  if (out->nlev != 0) {
    fprintf (stderr, "Error, this is a program bug!!!\n");
    fprintf (stderr, "       The common zd grid must be uninitialized at this point.\n");
    fprintf (stderr, "       This error can only occur after larger restructuring of uvspec.\n");
    fprintf (stderr, "       Please contact the programmers and add your input file. Thanx.\n");
    fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
    return -1;
  }
  set_common_z_grid (z_atm, nlev, &(out->zd), &(out->nlev));


  /* check hydrostatic equation of the input file */
  status = check_hydrostatic_equation ( z_atm, out->microphys.press[0][0], out->microphys.temper[0][0], nlev,
                                        filename, quiet, verbose );
  if (status!=0) {
    fprintf (stderr, "Error %d returned by 'check_hydrostatic_equation'\n", status);
    fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
    return status;
  }

  free(z_atm);

  /* Finally, check if profiles of pressure, temperature and air density are consistent */
  /* with ideal gas law over 6 orders of magnitude in pressure; check only first six    */
  /* orders of magnitude because for very small pressures roundoff errors due to too    */
  /* few digits cause problems                                                          */
  firstnonideal=1;

  //fprintf(stderr, "3DAbs read atm nlev %d\n", out->nlev); 

  for (lev=0; lev<out->nlev; lev++) {

    idealgas=1;


    if (out->microphys.press[0][0][lev]/out->microphys.press[0][0][nlev-1]>1e-6) {
      tmpdens = out->microphys.press[0][0][lev]*1E-04 / (BOLTZMANN*out->microphys.temper[0][0][lev]);
      
      if ((tmpdens==0 && out->microphys.dens[MOL_AIR][0][0][lev]!=0) || 
	  (tmpdens!=0 && out->microphys.dens[MOL_AIR][0][0][lev]==0))
	idealgas=0;
      
      if (out->microphys.dens[MOL_AIR][0][0][lev]>0)
	if ((idealdiff = fabs ((tmpdens - out->microphys.dens[MOL_AIR][0][0][lev]) / out->microphys.dens[MOL_AIR][0][0][lev])) > 0.01)
	  idealgas=0;
      
      if (!idealgas && firstnonideal) {
	firstnonideal=0;
	fprintf (stderr, "\n");
	fprintf (stderr, "*** ERROR: pressure, temperature, and air density not\n");
	fprintf (stderr, "*** consistent with ideal gas law at the following levels (read_atmosphere)\n");
	fprintf (stderr, "*** (assuming a Boltzmann constant %.6e):\n", BOLTZMANN);
	fprintf (stderr, "***    z [km]  density [cm-3]  p/kT [cm-3]  difference  press[hPa]  T[K]\n");
      }
      
      if (!idealgas) {
	fprintf (stderr, "***   %7.2f  %.4e      %.4e    %6.1f%%  %10.5f  %8.3f\n", 
		 out->zd[lev], out->microphys.dens[MOL_AIR][0][0][lev], tmpdens, idealdiff*100.0, 
                 out->microphys.press[0][0][lev], out->microphys.temper[0][0][lev]);
      }
    }
  }
  
  if (!firstnonideal) {
    fprintf (stderr, "*** Please check pressure, temperature, and air density in\n");
    fprintf (stderr, "*** %s for consistency!\n", filename);
    fprintf (stderr, "\n");
    return -1;
  }
  
  return 0;
}

/***********************************************************************************/
/* Function: check_hydrostatic_equation                                            */
/* Description:                                                                    */
/*  This function checks a profile, if the height grid z is                        */
/*  consistent with the profiles of the pressure and temperature                   */
/*  regarding the hydrostatic equation                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*  z      z-grid to be check                                                      */
/*  p,T    pressure and temperature                                                */
/*  nlev   number of levels                                                        */
/*                                                                                 */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    atmosphere.c                                                          */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Feb 2008   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/
int check_hydrostatic_equation ( float *z, float *p, float *T, int nlev, char *profile_name, int quiet, int verbose )
{
  int status=0;
  int lev = NOT_DEFINED_INTEGER;

  int first=0;
  float *z_test=NULL;
  float dz=NOT_DEFINED_FLOAT, dz_test=NOT_DEFINED_FLOAT;

  /* calculate z from T and p using hydrostatic equation (in ancillary.c) */
  status = calculate_z_from_p_and_T (p, T, &(z_test), nlev, 
                                     z[nlev-1], verbose);
  if (status!=0) {
    fprintf (stderr, "Error %d returned by function 'calculate_z_from_p_and_T'\n", status);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }
  /* check */
  for (lev=0; lev<nlev-1; lev++) {
    dz      = z[lev] - z[lev+1];
    dz_test = z_test[lev]  - z_test[lev+1];
    if ( fabs((dz_test - dz)/dz*100.0) > 3.35 && z[lev+1] < 100.0 ) { /* deviation larger than 3.35% and z below 100km */
      if (!quiet) {
        if (first==0) {
	  fprintf (stderr, "\n");  
          fprintf (stderr, "*** WARNING: The atmosphere defined in %s\n" , profile_name);
          fprintf (stderr, "*** is not in hydrostatic equilibrium (assuming the graviational acceleration of Earth).\n");
          fprintf (stderr, "*** This might not be what you want and could affect your interpretation of the results\n");
          fprintf (stderr, "*** because several relationships which are based on the implicit assumption of\n");
	  fprintf (stderr, "*** hydrostatic equilibrium are not true in this case. E.g. heating rates in uvspec are\n");
	  fprintf (stderr, "*** calculated from dE_net / dz and not from dE_net / dp which would give different\n");
          fprintf (stderr, "*** results. Also, layer optical thicknesses are calculated from pressure and layer\n");
	  fprintf (stderr, "*** thickness. The radiative transfer result is still correct - just make sure that\n");
	  fprintf (stderr, "*** you interpret it correctly! In particular the deviations are:\n");
          fprintf (stderr, " -----------------------------------------------\n");
          fprintf (stderr, "  z_atm1   z_atm2   dz_atm  dz_hydro   deviation\n");
          fprintf (stderr, "   [km]    [km]      [km]     [km]     per cent \n");
          fprintf (stderr, " -----------------------------------------------\n");
          first = 1;
        }
        fprintf (stderr, " %7.3f  %7.3f  %7.3f  %7.3f  %8.3f\n", 
		 z[lev+1], z[lev], dz, dz_test, (dz_test - dz)/dz*100.0);
      }
    }
  }
  if (first==1 && !quiet)
    fprintf (stderr,"\n");
  
  free(z_test);
  return status;
}

/**************************************************************/
/* Extrapolate atmospheric profiles for altitudes,            */
/* which are below the levels given by the atmosphere_file    */
/**************************************************************/

static int extrapolate_atmosphere (float altitude, atm_out_struct *out, 
                                   int quiet, int verbose)
{
  int status = 0;
  int gas    = NOT_DEFINED_INTEGER;
  float   g         = NOT_DEFINED_FLOAT;
  float   dz        = NOT_DEFINED_FLOAT;
  float   mix_ratio = NOT_DEFINED_FLOAT, mix_ratio1 = NOT_DEFINED_FLOAT, mix_ratio2 = NOT_DEFINED_FLOAT;
  /* float   rel_hum   = NOT_DEFINED_FLOAT; */

  /* realloc z-grid */
  if ((out->zd = realloc (out->zd, (out->nlev+1) * sizeof(float))) == NULL) {
    fprintf (stderr,"Error reallocating memory for output->atm.zd\n");
    return -1;
  }
  out->zd[out->nlev] = altitude;


  /* realloc temperature */
  if ((out->microphys.temper[0][0] = realloc (out->microphys.temper[0][0], (out->nlev+1) * sizeof(float))) == NULL) {
    fprintf (stderr,"Error reallocating memory for output->atm.zd\n");
    return -1;
  }
  out->microphys.temper[0][0][out->nlev] = (out->microphys.temper[0][0][out->nlev-1]-out->microphys.temper[0][0][out->nlev-2])/
    (out->zd[out->nlev-1] - out->zd[out->nlev-2] ) * 
    (out->zd[out->nlev] - out->zd[out->nlev-1]) 
    + out->microphys.temper[0][0][out->nlev-1];

  /* realloc pressure */
  if ((out->microphys.press[0][0] = realloc (out->microphys.press[0][0], (out->nlev+1) * sizeof(float) )) == NULL) {
    fprintf (stderr,"Error reallocating memory for output->atm.zd\n");
    return -1;
  }

  g = G_SURFACE * (R_EARTH * R_EARTH) / ((R_EARTH+out->zd[out->nlev-1]*1000.)*(R_EARTH+out->zd[out->nlev-1]*1000.));
  /* hydrostatic equation assuming linear temperature gradient and */ 
  /* NO variation of g within the last layer */
  /* dz = (out->zd[out->nlev] - out->zd[out->nlev-1]) * 1000.0; */
  /* or WITH variation of g within the last layer */
  dz = ( R_EARTH*R_EARTH/(R_EARTH + 1000.0*out->zd[out->nlev-1]) - R_EARTH*R_EARTH/(R_EARTH + 1000.0*out->zd[out->nlev]) );

  out->microphys.press[0][0][out->nlev] = out->microphys.press[0][0][out->nlev-1] * 
    exp( -(g*dz) / (R_AIR * 0.5 * (out->microphys.temper[0][0][out->nlev] + out->microphys.temper[0][0][out->nlev-1])));  

  /* realloc gas densities */
  for (gas=0;gas<MOL_NN;gas++) {
    if ((out->microphys.dens[gas][0][0] = realloc (out->microphys.dens[gas][0][0], (out->nlev+1) * sizeof(float) )) == NULL) {
      fprintf (stderr,"Error allocating memory for dens\n");
      return -1;
    }
  }
  /* calculate air number density from pressure and temperature */

  out->microphys.dens[MOL_AIR][0][0][out->nlev] = out->microphys.press[0][0][out->nlev] * 1E-04 / (BOLTZMANN*out->microphys.temper[0][0][out->nlev]);

  /* for all other gases use constant mixing ratio */
  for (gas=0;gas<MOL_NN;gas++) {
    switch (gas) {
    case MOL_AIR:
    case MOL_O4:  /* Calculated from O2 profile, ak 20110404 */
      /* do nothing */
      break;
      /* this is for constant relative humidity */
      /* case MOL_H2O: */  
      /*  rel_hum = out->microphys.dens[gas][out->nlev-1] / vapor_pressure(out->microphys.temper[out->nlev-1]) * 100.0; */
      /*  out->microphys.dens[gas][out->nlev] = rel_hum / 100.0 * vapor_pressure(out->microphys.temper[out->nlev]); */
      /*  fprintf (stderr,"rel_hum = %5.2f, %e, %e\n",rel_hum, */
      /*                     vapor_pressure(out->microphys.temper[out->nlev-1]), */
      /*                     vapor_pressure(out->microphys.temper[out->nlev]) ); */
      /*  break; */
      /*  here we assume a constant gradient of the mixing ratio */
    case MOL_H2O:
    case MOL_O3:
    case MOL_O2:
    case MOL_CO2:
    case MOL_NO2:
    case MOL_BRO:
    case MOL_OCLO: 
    case MOL_HCHO:
    case MOL_SO2:
    case MOL_CH4:
    case MOL_N2O:
    case MOL_CO:
    case MOL_N2: 
      mix_ratio2 = out->microphys.dens[gas][0][0][out->nlev-2] / out->microphys.dens[MOL_AIR][0][0][out->nlev-2];
      mix_ratio1 = out->microphys.dens[gas][0][0][out->nlev-1] / out->microphys.dens[MOL_AIR][0][0][out->nlev-1];
      mix_ratio  = mix_ratio1 + (mix_ratio1-mix_ratio2) / (out->zd[out->nlev-1] - out->zd[out->nlev-2] ) * 
	(out->zd[out->nlev] - out->zd[out->nlev-1]);
      out->microphys.dens[gas][0][0][out->nlev] = mix_ratio * out->microphys.dens[MOL_AIR][0][0][out->nlev]; 
      break;
      
    default:
      fprintf (stderr, "Error, unkown gas %d in setup_gases\n", gas);
      return -1;
      break;
    }
  }

  out->nlev = out->nlev + 1;
  
  return status;
}

 



/******************************************************************/
/* Read radiosonde file                                           */
/* and convert p and T data to a z-grid                           */
/* additional option to manipulate atmosphere profiles UH 2006-06 */
/******************************************************************/

static int read_radiosonde_file (input_struct input, output_struct *output)

{
  int rows=0, max_columns=0, min_columns=0;
  float **atm_data = NULL;
  float  *p_rs     = NULL;
  float  *T_rs     = NULL;  
  float **dens_rs  = NULL;
  float  *dens_air_atm = NULL;
  int     gas_nr   = NOT_DEFINED_INTEGER;
  int     lc       = NOT_DEFINED_INTEGER;
  int     status   = NOT_DEFINED_INTEGER;

  float   g        = NOT_DEFINED_FLOAT;
  float   dz       = NOT_DEFINED_FLOAT;
  int     i        = NOT_DEFINED_INTEGER;
  int     p_zero   = FALSE;
  char   *gas      = NULL;
  int     n_comb   = NOT_DEFINED_INTEGER;
  int     n_above  = NOT_DEFINED_INTEGER;
  float  *z_comb   = NULL;
  float  *p_comb   = NULL;
  float  *T_comb   = NULL;  
  float **dens_comb= NULL;
  int     n_outside= NOT_DEFINED_INTEGER;   /* number of points, that are in the radiosonde data, but not in the atmosphere_file */
  float scale_factor = NOT_DEFINED_FLOAT;
 
  if (input.atm.rs_source == RS_FROM_FILE) {
    /* read radiosonde ASCII file */

    if (input.verbose)
      fprintf (stderr, " ... reading radiosonde data file %s\n", input.filename[FN_RADIOSONDE]);

    status = ASCII_file2float (input.filename[FN_RADIOSONDE], 
                               &rows, &max_columns, &min_columns, 
                               &atm_data);

    if (status!=0) {
      fprintf (stderr, "Error %d reading radiosonde file: %s\n", status, input.filename[FN_RADIOSONDE]);
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }

  }
  else if (input.atm.rs_source == RS_FROM_ECMWF) {
    /* read ECMWF data file*/

    if (input.verbose)
      fprintf (stderr, " ... reading ECMWF atmosphere file %s\n", input.filename[FN_ECMWF]);
        
    status = read_ECMWF_atmosphere (input.latitude, input.longitude, input.UTC, input.atm.time_interpolate, 
                                    input.filename[FN_ECMWF], input.atm.ECMWF_new_format,
                                    &atm_data, &rows, output->alt.altitude, &(output->atm.z_model_layer), &(output->atm.nz_model_layer),
                                    &(output->atm.microphys.z_cpt_sea), &(output->atm.microphys.press_cpt), &(output->atm.microphys.temper_cpt),
                                    input.verbose, input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d reading ECMWF atmosphere file '%s'\n", status, input.filename[FN_ECMWF]);
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }
    
    /* dummy arguments for columns for checks later */
    max_columns=4;
    min_columns=4;
  }
  else if (input.atm.rs_source == RS_FROM_ECHAM) {
    /* read ECHAM data file */
    if (input.verbose)
      fprintf (stderr, " ... reading ECHAM data file %s\n", input.filename[FN_ECHAM]);
    
    status = read_ECHAM_atmosphere (input.latitude, input.longitude, input.UTC,
                                    input.atm.time_interpolate, 
                                    input.filename[FN_ECHAM],
                                    &atm_data, &rows, output->alt.altitude,
                                    &(output->atm.z_model_layer),
                                    &(output->atm.nz_model_layer),
                                    input.verbose, input.quiet);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading ECHAM file: %s\n", status, input.filename[FN_ECHAM]);
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }
      
    /* dummy arguments for columns for checks later */
    max_columns=3;
    min_columns=3;
  }
  
  /* check consistent number of columns */
  if (max_columns!=(2+input.atm.n_rs_gas) || min_columns!=(2+input.atm.n_rs_gas)) {
    fprintf (stderr, "Error, inconsistent number of columns in %s\n", input.filename[FN_RADIOSONDE]);
    fprintf (stderr, "       min = %3d, max =%3d\n", min_columns, max_columns);
    fprintf (stderr, "       Radiosonde data is expected to have 2 + %d columns\n", input.atm.n_rs_gas);
    fprintf (stderr, "       First pressure [hPa], second temperature [K], and specified gas densities\n");
    fprintf (stderr, "       (line %d, function '%s' in '%s') \n", __LINE__, __func__, __FILE__);
    return -1;
  }

  output->atm.nz_model_level = rows;

  /* copy radiosonde p-data to p_rs */
  p_rs = ASCII_column_float (atm_data, output->atm.nz_model_level, 0);

  /* copy radiosonde T-data to T_rs */
  T_rs = ASCII_column_float (atm_data, output->atm.nz_model_level, 1);

  /* if uppermost pressure level is 0, than replace it by uppermost background atmosphere z-height */
  if (p_rs[0] == 0.0) {
    p_zero=TRUE;
    p_rs[0] = output->atm.microphys.press[0][0][0];
  }

  /* calculate z from T and p using hydrostatic equation (in ancillary.c) */
  status = calculate_z_from_p_and_T (p_rs, T_rs, &(output->atm.z_model_level), output->atm.nz_model_level, output->alt.altitude, input.verbose);
  if (status!=0) {
    fprintf (stderr, "Error %d returned by calculate_z_from_p_and_T (line %d, function '%s' in '%s') \n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /*   for (lc=0; lc<output->atm.nz_model_level; lc++) */
  /*     fprintf (stderr," ### %3d atm p=%e T=%e, z=%e\n", lc, p_rs[lc], T_rs[lc], output->atm.z_model_level[lc] ); */

  /* correct uppermost level to something inside the range of the background atmosphere, if above */
  if ( p_zero == TRUE ) {
    if (input.verbose) 
      fprintf (stderr, " ... replace uppermost ECMWF level p=0 with uppermost level of the atmosphere\n");

    /* replace uppermost height level */
    output->atm.z_model_level[0] = output->atm.zd[0];

    /* find appropriate pressure for new z_top */
    /* hydrostatic equation assuming linear temperature gradient and variation of g within the last layer */
    dz = ( R_EARTH*R_EARTH/(R_EARTH + 1000.0* output->atm.z_model_level[1]) - R_EARTH*R_EARTH/(R_EARTH + 1000.0*output->atm.z_model_level[0]) );
    /* replace uppermost pressure level */
    p_rs[0] = p_rs[1] * exp( -(G_SURFACE*dz) / (R_AIR * 0.5 * (T_rs[0] + T_rs[1])));  

    if (input.verbose) 
      fprintf (stderr, "     z_TOA = %10.6f, p_TOA = %10.8f \n", output->atm.z_model_level[0], p_rs[0]);
  }

  /* count radiosonde points which are above the highest level of the atmosphere_file */
  n_outside = 0;
  for (lc=output->atm.nz_model_level-1; lc>=0;lc--) {
    if (output->atm.z_model_level[lc] > output->atm.zd[0] ) {
      n_outside = lc+1;
      break;
    }
  }

  /* throw away radiosonde data above of the highest atmosphere_file level and give warning */
  if (n_outside != 0 ) {
    if (!input.quiet) {
      fprintf (stderr, "*** Warning, %3d levels of the radiosonde data are outside\n", n_outside);
      fprintf (stderr, "*** the range of the atmosphere file:\n");
      for (lc=0; lc<n_outside; lc++)
        fprintf (stderr, "     (output->atm.z_model_level[%3d] = %8.4f) > (max(z_atm) = %8.4f)\n", lc, output->atm.z_model_level[lc], output->atm.zd[0]);
    }
    for (lc=0; lc<output->atm.nz_model_level-n_outside;lc++) {
      p_rs[lc] = p_rs [lc+n_outside];
      T_rs[lc] = T_rs [lc+n_outside];
      output->atm.z_model_level[lc] = output->atm.z_model_level [lc+n_outside];
      for ( i=0; i < input.atm.n_rs_gas+2; i++ )
        atm_data[lc][i] = atm_data[lc+n_outside][i];
    }
    output->atm.nz_model_level = output->atm.nz_model_level - n_outside;
  }

  /* Allocate dens_rs */
  if ((dens_rs = (float **) calloc (MOL_NN, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for dens_rs\n");
    return -5;
  }

  /* Allocate on the radiosonde_file grid */
  if ((dens_rs[MOL_AIR] = (float *) calloc (output->atm.nz_model_level, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for dens_rs[air]\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -6;
  }

  /* calculate air density of the radiosonde_file with ideal gas law */
  for (lc=0; lc < output->atm.nz_model_level ; lc ++ ) {
    dens_rs[MOL_AIR][lc] = p_rs[lc] * 1.0E-04 / (BOLTZMANN * T_rs[lc]);
    /* fprintf (stderr,"lc = %3d, n_rs = %e, p_rs = %10.5f, T_rs = %8.4f, z_rs = %7.3f\n", */
    /*                  lc, dens_rs[MOL_AIR][lc], p_rs[lc], T_rs[lc], output->atm.z_model_level[lc] ); */
  }


  /* calculate air number density of the atmosphere_file profile on the radiosonde grid     */
  /* for vertical interpolation of the background atmosphere, we have to take this density  */
  /* as the mixing ratio is kept constant, but number density would change with air density */
  /* and therefor column would change, which is not wanted                                  */

  if ((dens_air_atm = (float *) calloc (output->atm.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for dens_air_atm\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -7;
  }

  for (lc=0; lc<output->atm.nlev;lc++)
    dens_air_atm[lc] = output->atm.microphys.dens[MOL_AIR][0][0][lc];

  /* AIR (this is not 100% consistent with interpolation of p and T)*/
  status = interpolate_profile (output->atm.zd, &(dens_air_atm), output->atm.nlev, 
				output->atm.z_model_level, output->atm.nz_model_level, input.atm.interpol_method_gas[MOL_AIR], input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating dens_rs[air]\n", status);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }


  /* now (p, T and z from the radiosonde) are combined with */
  /* gas profiles from atmosphere file on rs-grid           */
  for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++) {
    switch (gas_nr) {
    case MOL_AIR:
      /* has to be done before, as interpolation of the other gases does need dens_rs */
      break;
    case MOL_O3:
    case MOL_O2:
    case MOL_H2O:
    case MOL_CO2:
    case MOL_NO2:
    case MOL_BRO:
    case MOL_OCLO:
    case MOL_HCHO:
    case MOL_SO2:
    case MOL_O4:
    case MOL_CH4:
    case MOL_N2O:
    case MOL_CO:
    case MOL_N2:

      /* allocation on atmosphere_file grid */
      if ((dens_rs[gas_nr] = (float *) calloc (output->atm.nlev, sizeof (float))) == NULL) {
        fprintf (stderr,"Error allocating memory for dens_rs[%d]\n",gas_nr);
        fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
        return -6;
      }

      /* copy number density */
      for (lc=0; lc<output->atm.nlev;lc++)
        dens_rs[gas_nr][lc] = output->atm.microphys.dens[gas_nr][0][0][lc]; 

      if ( input.atm.well_mixed_gas[gas_nr] == YES ) {
 
        /* Interpolate dens to the vertical grid of the radiosonde. */
        /* here the dens_rs[MOL_AIR] is used to keep mixing ratios  */
        /* constant in one height                                   */

	status = interpolate_density (output->atm.zd, &(dens_rs[gas_nr]), output->atm.nlev, 
				      output->atm.z_model_level, output->atm.nz_model_level, input.atm.interpol_method_gas[gas_nr],
				      output->atm.microphys.dens[MOL_AIR][0][0], dens_rs[MOL_AIR], input.quiet);
        if (status!=0) {
          fprintf (stderr, "Error %d interpolating dens_rs[%d]\n", status, gas_nr);
          fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
          return status;
        }
      }
      else {

        /* Interpolate dens to the vertical grid of the radiosonde.     */
        /* In usually case it would be wrong to use dens_air_atm, but   */
        /* for not well mixed gases we like to keep the number density  */
        /* constant in one leveland not the mixing ratio                */

        status = interpolate_density (output->atm.zd, &(dens_rs[gas_nr]), output->atm.nlev, 
				      output->atm.z_model_level, output->atm.nz_model_level, input.atm.interpol_method_gas[gas_nr],
				      output->atm.microphys.dens[MOL_AIR][0][0], dens_air_atm, input.quiet);
        if (status!=0) {
          fprintf (stderr, "Error %d interpolating dens_rs[%d]\n", status, gas_nr);
          fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
          return status;
        }
      }

      break;
    default:
      fprintf (stderr, "Error, unknown indentifier MOL_** = %d !\n", gas_nr);
      fprintf (stderr, "       This is a program bug, please contact the programmers!\n");
      fprintf (stderr, "       Please include the used input-file!\n");
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return -8;
      break;
    }
  }

  free(dens_air_atm);

  /* if there are any gas profiles provided by the radiosonde file  */
  /* take them to replace the standard profile                      */
  /* unit (e.g. massmixratio, volmixratio, rel.hum.) is converted   */ 
  /* to number density [cm-3]                                       */
  for (i=0; i < input.atm.n_rs_gas; i++) {

    gas_nr = input.atm.rs_gas[i];
    gas = gas_number2string(gas_nr);
    if (strcmp(" -1 ",gas)==0) {
      fprintf (stderr, "Error wrong gas identifier in atmophere.c\n");
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return -1;
    }

    if (input.verbose)
      fprintf (stderr, "     Reading %s profile from radiosonde data, unit=%d (1=cm-3, 2=m-3, 3=MMR, 4=VMR, 5=rel.hum.)\n",
	       gas, input.atm.rs_unit[gas_nr]);

    free(gas);

    free(dens_rs[gas_nr]);
    dens_rs[gas_nr] = ASCII_column_float (atm_data, output->atm.nz_model_level, i + 2);

    /* for (lc=0; lc<output->atm.nz_model_level;lc++) */
    /*   fprintf (stderr, "z=%6.3f, p=%12.6f, T=%8.3f, gas=%12.6e\n",output->atm.z_model_level[lc],p_rs[lc],T_rs[lc],dens_rs[gas_nr][lc]); */

    /* convert units of the dens_file to cm-3 */
    status =  convert_profile_units (&(dens_rs[gas_nr]), output->atm.z_model_level, output->atm.nz_model_level, gas_nr,
                                     dens_rs[MOL_AIR], T_rs, output->atm.z_model_level, output->atm.nz_model_level,
                                     input.atm.interpol_method_gas[MOL_AIR], input.atm.interpol_method_temper,
                                     input.atm.rs_unit[gas_nr], input.atm.mol_mass, input.quiet);
    
    if (status!=0) {
      fprintf (stderr, "Error %d returned by function 'convert_profile_unit'\n", status);
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }
  }

  for (lc=0; lc<output->atm.nz_model_level;lc++)
    free(atm_data[lc]);
  free(atm_data);

  /**************************************************************/
  /* now we have a complete new set of variables on the rs grid */
  /**************************************************************/

  /* /\* additional verbose output *\/ */
  /* if (input.verbose) { */
  /*   fprintf(stderr,"(radiosonde p, T and z) with densities from atmosphere file\n"); */
  /*   fprintf (stderr, " lc    z           p       T         air           O3            O2            "); */
  /*   fprintf (stderr,"H2O            CO2           NO2           BRO          OCLO         HCHO           O4\n"); */
  /*   for (lc=0; lc<output->atm.nz_model_level;lc++) */
  /*     fprintf (stderr,"%3d %7.3f %11.5f %7.2f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n", */
  /*                       lc, output->atm.z_model_level[lc], p_rs[lc], T_rs[lc], */
  /*                       dens_rs[0][lc],dens_rs[1][lc],dens_rs[2][lc],dens_rs[3][lc],dens_rs[4][lc], */
  /*                       dens_rs[5][lc],dens_rs[6][lc],dens_rs[7][lc],dens_rs[8][lc],dens_rs[9][lc]); */
  /* } */

  /************************************************************/
  /* combine radiosonde profiles and background profile above */
  /************************************************************/

  if (input.atm.rs_add_upper_levels) {

    /* count levels of the atmosphere file above the upper most level of the radiosonde */
    for (lc=0;lc<output->atm.nlev;lc++)
      if ( output->atm.zd[lc] <= output->atm.z_model_level[0] )
        break;

    n_above = lc;
    n_comb  = n_above + output->atm.nz_model_level;

    /* Allocate memory for combined profiles */
    if ((z_comb = (float *) calloc (n_comb, sizeof(float))) == NULL) {
      fprintf (stderr,"Error allocating memory for z_comb\n");
      return -7;
    }
    if ((p_comb = (float *) calloc (n_comb, sizeof(float))) == NULL) {
      fprintf (stderr,"Error allocating memory for z_comb\n");
      return -8;
    }
    if ((T_comb = (float *) calloc (n_comb, sizeof(float))) == NULL) {
      fprintf (stderr,"Error allocating memory for z_comb\n");
      return -9;
    }

    if ((dens_comb = (float **) calloc (MOL_NN, sizeof (float *))) == NULL) {
      fprintf (stderr,"Error allocating memory for dens_comb\n");
      return -10;
    }

    for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++) {
      if ((dens_comb[gas_nr] = (float *) calloc (n_comb, sizeof (float))) == NULL) {
        fprintf (stderr,"Error allocating memory for dens_comb[%d]\n", gas_nr);
        return -11;
      }
    }

    /* fprintf (stderr,"combine profiles, n_comb = %3d\n", n_comb); */
    /* fprintf (stderr, " lc    z           p       T         air           O3            O2         "); */
    /* fprintf (stderr, "   H2O            CO2           NO2           BRO          OCLO         HCHO           O4\n"); */

    /* points inside the radiosonde profile, starting from the ground (lc = n_comb-1) */
    for (lc=n_comb-1;lc>(n_comb-1-output->atm.nz_model_level);lc--) {
      p_comb[lc] = p_rs[lc-(n_above)];
      T_comb[lc] = T_rs[lc-(n_above)];
      z_comb[lc] = output->atm.z_model_level[lc-(n_above)];
      for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++)
        dens_comb[gas_nr][lc] = dens_rs[gas_nr][lc-(n_above)];
      /* fprintf (stderr,"%3d %7.3f %11.5f %7.2f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n", */
      /*                  lc, z_comb[lc], p_comb[lc], T_comb[lc], */
      /*                  dens_comb[0][lc],dens_comb[1][lc],dens_comb[2][lc],dens_comb[3][lc],dens_comb[4][lc], */
      /*                  dens_comb[5][lc],dens_comb[6][lc],dens_comb[7][lc],dens_comb[8][lc],dens_comb[9][lc]); */
    }

    if (n_comb-1-output->atm.nz_model_level >= 0) {
  
      /* transition layer to upper part of the atmosphere, given by background atmosphere */
      /* fprintf (stderr,"tansition layer %3d\n", n_comb-1-output->atm.nz_model_level); */
      lc=n_comb-1-output->atm.nz_model_level;
      z_comb[lc] = output->atm.zd [n_above-1];
      T_comb[lc] = output->atm.microphys.temper[0][0][n_above-1];

      /* g is treated as constant inside the next layer */
      g = G_SURFACE * (R_EARTH * R_EARTH) / ((R_EARTH+output->atm.z_model_level[lc]*1000.)*(R_EARTH+output->atm.z_model_level[lc]*1000.));

      p_comb[lc] = p_comb[lc+1] * exp( - (g * 1000.0 * (z_comb[lc]-z_comb[lc+1])) / (R_AIR * 0.5 * (T_comb[lc]+T_comb[lc+1])) );

      /* fprintf (stderr,"between:lc=%3d, z1=%7.3f, z2=%7.3f, T_rs=%7.2f, T_atm=%7.2f, p_rs=%12.6f, p_atm=%12.6f, p_new=%14.8f\n", */
      /*          lc,z_comb[lc+1],z_comb[lc],T_rs[0],output->atm.microphys.temper[n_above-1], */
      /*          p_comb[lc+1], output->atm.microphys.press[n_above-1], p_comb[lc]); */
        

      /* adjustment of the new pressure-grid */
      scale_factor = p_comb[lc] / output->atm.microphys.press[0][0][n_above-1];
      scale_pressure (&(output->atm.microphys.press),
		      &(output->atm.microphys.dens),output->atm.nlev,
		      scale_factor, input.atm.well_mixed_gas, n_above-1, FALSE);

      for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++) 
        dens_comb[gas_nr][lc] = output->atm.microphys.dens[gas_nr][0][0][n_above-1];

      /* fprintf (stderr,"%3d %7.3f %11.5f %7.2f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n", */
      /*                  lc, z_comb[lc], p_comb[lc], T_comb[lc], */
      /*                  dens_comb[0][lc],dens_comb[1][lc],dens_comb[2][lc],dens_comb[3][lc],densé_comb[4][lc], */
      /*                  dens_comb[5][lc],dens_comb[6][lc],dens_comb[7][lc],dens_comb[8][lc]); */


      /* points above the radiosonde data are copied from the background atmosphere */
      for (lc=n_above-2;lc>=0;lc--) {
        p_comb[lc] = output->atm.microphys.press[0][0][lc];
        T_comb[lc] = output->atm.microphys.temper[0][0][lc];
        z_comb[lc] = output->atm.zd[lc];
        for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++)
          dens_comb[gas_nr][lc] = output->atm.microphys.dens[gas_nr][0][0][lc];
	/* fprintf (stderr,"%3d %7.3f %11.5f %7.2f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n", */
	/*                 lc, z_comb[lc], p_comb[lc], T_comb[lc], */
	/*                 dens_comb[0][lc],dens_comb[1][lc],dens_comb[2][lc],dens_comb[3][lc],dens_comb[4][lc], */
	/*                 dens_comb[5][lc],dens_comb[6][lc],dens_comb[7][lc],dens_comb[8][lc]); */
      }
    }

    /* radiosonde arrays have done their job */
    free(p_rs);
    free(T_rs);
    /* free(output->atm.z_model_level); */
    for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++)
      free(dens_rs[gas_nr]);
    free(dens_rs);

  } /* if (input.atm.rs_add_upper_levels) { ... */
  else {

    /* only radiosonde/model levels are from now on atmosphere levels */
    n_comb = output->atm.nz_model_level;

    if ((z_comb = (float *) calloc (n_comb, sizeof(float))) == NULL) {
      fprintf (stderr,"Error allocating memory for z_comb\n");
      return -7;
    }

    /* copy z-levels of the radiosonde/model */
    for (lc=0;lc<n_comb;lc++)
      z_comb[lc] = output->atm.z_model_level[lc];

    /* attention, this are pointers! */
    p_comb = p_rs;
    T_comb = T_rs;
    dens_comb = dens_rs;
    /* p_rs, T_rs and dens_rs are free'ed at the same time, when p_comb, T_comb, and dens_comb will be free'ed */

  }

  /*--------------------------------*/
  /* copy data to final destination */
  /*--------------------------------*/

  /* reallocate output->atm.zd */
  free(output->atm.zd);
  free(output->atm.microphys.press[0][0]);
  free(output->atm.microphys.temper[0][0]);
  for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++)
    free(output->atm.microphys.dens[gas_nr][0][0]);

  if ((output->atm.zd               = (float *) calloc (n_comb, sizeof(float))) == NULL) {
    fprintf (stderr, "Error allocating memory for output->atm.zd\n");
    return -12;
  }

  if ((output->atm.microphys.press[0][0]  = (float *) calloc (n_comb, sizeof(float))) == NULL) {
    fprintf (stderr, "Error allocating memory for output->atm.microphys.press\n");
    return -13;
  }

  if ((output->atm.microphys.temper[0][0] = (float *) calloc (n_comb, sizeof(float))) == NULL) {
    fprintf (stderr, "Error allocating memory for output->atm.microphys.temper\n");
    return -14;
  }

  for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++)
    if ((output->atm.microphys.dens[gas_nr][0][0] = (float *) calloc (n_comb, sizeof(float))) == NULL) {
      fprintf (stderr, "Error allocating memory for output->atm.zd\n");
      return -15;
    }

  /* copy all to final destination */

  output->atm.nlev = n_comb;

  for (lc=0;lc<n_comb;lc++) {
    output->atm.zd                      [lc] = z_comb[lc];
    output->atm.microphys.press[0][0]         [lc] = p_comb[lc];
    output->atm.microphys.temper[0][0]        [lc] = T_comb[lc];
    for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++)
      output->atm.microphys.dens[gas_nr][0][0][lc] = dens_comb[gas_nr][lc];
  }

  /* combined arrays have done their job */
  free(z_comb);
  free(p_comb);
  free(T_comb);
  for (gas_nr = 0; gas_nr < MOL_NN; gas_nr++)
    free(dens_comb[gas_nr]);
  free(dens_comb);

  return status;
}



/**************************************************************/
/* Read ECHAM atmosphere file                                 */
/*                                                            */
/* This is basically the read_ECHAM_atmosphere function       */
/* adapted for ECHAM data.                                    */
/*                                                            */ 
/* ECHAM file must be in netCDF data format and must contain  */
/* p (pressure profile)                                       */
/* p_full (full layer pressure)                               */
/* T (temperature profile)                                    */
/* Q (H2O mass mixing ratio)                                  */
/**************************************************************/

static int read_ECHAM_atmosphere (float lat, float lon, struct tm UTC, 
                                  int time_interpolate, char *filename, 
                                  float ***atm_data, int *rows, float altitude,
                                  float **z_model_layer, int *nz_model_layer,
                                  int verbose, int quiet)
{
  int status=0;

#if HAVE_LIBNETCDF

  int ncid=NOT_DEFINED_INTEGER;

  int lc=NOT_DEFINED_INTEGER; 
  int t=NOT_DEFINED_INTEGER;
  int nt=-1;
  long ilat=NOT_DEFINED_INTEGER,ilon=NOT_DEFINED_INTEGER,itime=NOT_DEFINED_INTEGER; 

  int itime1 = -1, itime2 = -1;
  float dt = NOT_DEFINED_FLOAT;
  
  size_t nlat=0;
  double *ECHAM_lat=NULL;
  size_t nlon=0;
  double *ECHAM_lon=NULL;
  size_t nlev=39; /* ECHAM4 */
 
  float *Tsurf=NULL;
  float Tsurf_av=NOT_DEFINED_FLOAT;
  float **p=NULL;
  float *p_av=NULL;
  float **p_full=NULL;
  float *p_full_av=NULL;
  float **T=NULL;
  float *T_av=NULL;
  float **Q=NULL;
  float *Q_av=NULL;
  size_t start3D[3] = {0,0,0};  
  int id_data=0;
  
  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d opening ECHAM netCDF file %s\n", status, filename);
    return status;
  }
    
  /* read latitude array */
  alloc_and_read_netCDF_1D_double(ncid,"lat", &nlat, "lat", &(ECHAM_lat));

  /* search correct latitude index */
  ilat = get_grid_index(lat, ECHAM_lat, nlat, FALSE);
  if (ilat < 0) {
    fprintf (stderr, "Error -1 finding index for lat=%5.3f in %s\n", lat, filename);
    return -1;
  }

  /* read longitude */
  alloc_and_read_netCDF_1D_double(ncid,"lon", &nlon, "lon", &(ECHAM_lon));

  /* search correct longitude index */
  ilon = get_grid_index(lon, ECHAM_lon, nlon, TRUE);
  if (ilon < 0) {
    fprintf (stderr, "Error -2 finding index for lon=%5.2f in %s\n", lon, filename);
    return -2;
  }

  /* get time index */
  status = get_time_index (ncid, UTC, time_interpolate, 
                           &(nt), &(itime1), &(itime2), &(dt),
                           verbose, quiet);
  if (status != 0) {
    fprintf (stderr, "Error %d returned by get_time_index()\n", status);
    return status;
  }  
  
  /* alloc one or two timesteps for surface temperature */
  if ((Tsurf = calloc (nt, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for surface temperature\n");
    return -10;
  }

  /* alloc one or two timesteps for pressure */
  if ((p = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for pressure\n");
    return -10;
  }

  /* alloc one or two timesteps for pressure */
  if ((p_full = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for full level pressure\n");
    return -10;
  }
 
  /* alloc one or two timesteps for temperature */
  if ((T = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for temperature\n");
    return -10;
  }

  /* alloc one or two timesteps for mass mixing ratio water vapour */
  if ((Q = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for water vapour\n");
    return -10;
  }

  /* alloc variables for time interpolated p, p_full, T, Q */
  p_av      = calloc (nlev, sizeof (float));
  p_full_av = calloc (nlev, sizeof (float));
  T_av      = calloc (nlev, sizeof (float));
  Q_av      = calloc (nlev, sizeof (float));
  
  /* read nt (= 1 or 2) time steps */
  for (t=0;t<=nt-1;t++) {
    
    if (t == 0)
      itime = itime1;
    if (t == 1)
      itime = itime2;
    
    /* Read 2D fields */
    
    /* index for 2D fields */
    start3D[0]=itime;
    start3D[1]=ilat;
    start3D[2]=ilon;
   
    /* get variable id for "Tsurf" (surface temperature) */
    status = nc_inq_varid (ncid, "Tsurf", &id_data);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d, while getting id for surface temperature\n", status);
      return status;
    }
    /* read surface temperature */
    status = nc_get_var1_float (ncid, id_data, start3D, &Tsurf[t]);
    
    /* alloc and read half level pressure (level)  */
    status = alloc_and_read_netCDF_column_float(ncid, "p", &p[t], nlev, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading pressure from netCDF file %s\n", status, filename);
      return status;
    }
    
    /* alloc and read full level pressure (layer average)*/
    status = alloc_and_read_netCDF_column_float(ncid, "pfull", &p_full[t], nlev, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading pressure (full level) from netCDF file %s\n", status, filename);
      return status;
    } 
    
    /* allocate and read temperature */
    status = alloc_and_read_netCDF_column_float(ncid, "T", &T[t], nlev, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading temperature from netCDF file %s\n", status, filename);
      return status;
    }
    
    /* allocate and read h2o mass mixing ratio */
    status = alloc_and_read_netCDF_column_float(ncid, "Q", &Q[t], nlev, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading humidity from netCDF file %s\n", status, filename);
      return status;
    }
    
  }
  
  if (verbose) {
    fprintf (stderr, "     found %zd x %zd x %zd (lat x lon x lev) data points\n", nlat, nlon, nlev);
    fprintf (stderr, "     reading pixel lat = %5.2f (%4ld), lon = %5.2f (%4ld)\n",lat/3600.0,ilat,lon/3600.0,ilon);
  }
  
  /* allocate data for "mimikry" radiosonde data */
  if ((status = ASCII_calloc_float (atm_data,nlev,3)) != 0) {
    fprintf (stderr, " Error allocating memory for atm_data\n");   
    return status;
  }
  if (nt == 1) {
    for (lc=0; lc<nlev-1; lc++) {
      p_av[lc]=  p[0][lc]/100.;
      p_full_av[lc] = p_full[0][lc]/100.;
      T_av[lc]=  T[0][lc];
      Q_av[lc]=  Q[0][lc];
    }
    Tsurf_av= Tsurf[0];
  }
  else {
    /* time interpolation between the two timesteps */
    Tsurf_av = (1.0-dt)* Tsurf[0] + dt* Tsurf[1];
    /*   (*albedo_av) = (1.0-dt)* albedo[0] + dt* albedo[1]; */
    /*       (*sza_av) = (1.0-dt)* sza[0] + dt* sza[1]; */
    for (lc=0; lc<nlev-1; lc++) {
      p_av[lc] = ((1.0-dt)* p[0][lc]  +  dt* p[1][lc])/100.;
      T_av[lc] = (1.0-dt)* T[0][lc]  +  dt* T[1][lc];
      Q_av[lc] = (1.0-dt)* Q[0][lc]  +  dt* Q[1][lc];
      p_full_av[lc]= ((1.0-dt)* p_full[0][lc]  + 
                      dt* p_full[1][lc])/100.;
      
    }
  }
  (*rows)=nlev;
 
  
  /* Top level pressure is always 5 hPa, TOA properties taken as layer properties from layer below*/
  (*atm_data)[0][0]= 5.0;
  (*atm_data)[0][1]= T_av[0];
  (*atm_data)[0][2]= Q_av[0];

  float  *pressure_layer, *temperature_layer;
  
  pressure_layer=calloc(nlev, sizeof (float)); 
  temperature_layer=calloc(nlev, sizeof (float)); 
  (*nz_model_layer)=nlev;

  pressure_layer[0]=5.0; 
  temperature_layer[0]= T_av[0];

  for (lc=1; lc<nlev-1; lc++) {
    pressure_layer[lc]=p_av[lc-1];
    (*atm_data)[lc][0]=pressure_layer[lc];
    
    temperature_layer[lc]=(T_av[lc-1]*p_full_av[lc-1]*(p_full_av[lc]-p_av[lc-1])+
                           T_av[lc]*p_full_av[lc]*(p_av[lc-1]-p_full_av[lc-1]))/
      ( p_av[lc-1]*(p_full_av[lc]-p_full_av[lc-1]) );
    
    (*atm_data)[lc][1]=temperature_layer[lc];
    
    (*atm_data)[lc][2]=
      (Q_av[lc-1]*p_full_av[lc-1]*(p_full_av[lc]-p_av[lc-1]) 
       + Q_av[lc]*p_full_av[lc]*(p_av[lc-1]-p_full_av[lc-1]))/
      ( p_av[lc-1]*(p_full_av[lc]-p_full_av[lc-1]) ) ;  
  }
  
  /* surface properties */
  (*atm_data)[nlev-1][0]=p_av[lc-1];
  (*atm_data)[nlev-1][1]= Tsurf_av;
  (*atm_data)[nlev-1][2]= Q_av[nlev-2];
  
  pressure_layer[nlev-1]=p_av[lc-1];
  temperature_layer[nlev-1]=Tsurf_av;
  
  

  status = calculate_z_from_p_and_T (pressure_layer, temperature_layer, z_model_layer,
                                     (*nz_model_layer), 
                                     altitude, verbose);
  if (status!=0) {
    fprintf (stderr, "Error %d during calculate_z_from_p_and_T\n", status);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }
  /* additional verbose output 
     for (lc=0; lc<nlev; lc++) {
     fprintf(stderr, "atm_data: %d, %g %g %g %g\n",lc,(*z_model_layer)[lc], 
     (*atm_data)[lc][0], 
     (*atm_data)[lc][1],  (*atm_data)[lc][2] );
     } 
  */
  
  for (t=0;t<=nt-1;t++) {
    free(p[t]);
    free(T[t]);
    free(Q[t]);
    free(p_full[t]);
  }
  free(p);
  free(T);
  free(Q);
  free(p_av);
  free(T_av);
  free(Q_av);
  free(p_full_av);
  
#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use the ECHAM data file option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ******************************************************************\n");
  return -1;
#endif

  return status;
}



/****************************************************************/
/* Read ECMWF atmosphere file                                   */
/*                                                              */
/* Purpose:                                                     */
/* read p, T, Q, O3 from an ECMWF file                          */
/* ECMWF file must be in netCDF data format and must contain    */
/* hyai, hybi (hybrid coefficients at layer interfaces)         */
/* SP (surface pressure) or LNSP (logarithm of SP)              */
/* Q  (H2O mass mixing ratio)                                   */
/* O3 (Ozone mass mixing ratio)                                 */
/*                                                              */
/*  input                                                       */
/*  ------                                                      */
/*  lat        latitude                                         */
/*  lon        longitude                                        */
/*  UTC        universal time correlated                        */
/*  filename   name of the netCDF file including the wind map   */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  atm_data  [0]=p, [1]=T, [2]=O3, [3]=Q                       */
/*  December 2006  by Ulrich Hamann                             */
/****************************************************************/


static int read_ECMWF_atmosphere (float lat, float lon, struct tm UTC, int time_interpolate, char *filename, int ECMWF_new_format,
                                  float ***atm_data, int *rows, float altitude, float **z_model_layer, int *nz_model_layer,
                                  float *z_cpt_sea, float *press_cpt, float *temper_cpt, int verbose, int quiet)
{
  int status=0;

#if HAVE_LIBNETCDF

  int ncid    = NOT_DEFINED_INTEGER;
  int id_data = NOT_DEFINED_INTEGER;

  int lc=NOT_DEFINED_INTEGER, lc2=NOT_DEFINED_INTEGER, lc3=NOT_DEFINED_INTEGER;
  int t =NOT_DEFINED_INTEGER;
  int nt=-1;
  long ilat=NOT_DEFINED_INTEGER, ilon=NOT_DEFINED_INTEGER, itime=NOT_DEFINED_INTEGER;  /* index for lat, lon, time in netCDF file */

  size_t nlat  = 0;
  size_t nlon  = 0;

  double *ECMWF_lat  = NULL;
  double *ECMWF_lon  = NULL;

  int itime1 = -1, itime2 = -1;
  float dt = NOT_DEFINED_FLOAT;

  size_t nlev=0;
  size_t nlay=0;

  float **p_level = NULL;
  float **p_layer = NULL;
  float **T       = NULL;
  float **Q       = NULL;
  float **O3      = NULL;

  char lat_name[FILENAME_MAX]="";
  char lon_name[FILENAME_MAX]="";

  char z_name  [2]="";
  char T_name  [2]="";
  char O3_name [3]="";
  char Q_name  [2]="";

  int read_on_level_or_layer = 0;
  int level = 0;
  int layer = 1;

  float dx = NOT_DEFINED_FLOAT;

  /* float *z_layer=NULL; */
  float *press_layer=NULL;
  float *temper_layer=NULL;

  /* float kappa = (C_P_DRY_STD-C_V_DRY_STD) / C_P_DRY_STD;  /\* == R / cp *\/ */

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, filename);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }

  if ( ECMWF_new_format == TRUE ) {
    strcpy (z_name,   "z");
    strcpy (T_name,   "t");
    strcpy (O3_name,  "o3");
    strcpy (Q_name,   "q");
  }
  else {
    strcpy (z_name,   "Z");
    strcpy (T_name,   "T");
    strcpy (O3_name,  "O3");
    strcpy (Q_name,   "Q");
  }

  /* check format */
  status = nc_inq_varid (ncid, "lat", &id_data);
  if (status==NC_NOERR) {
    strcpy (lat_name, "lat");
    strcpy (lon_name, "lon");
  }
  else {
    status = nc_inq_varid (ncid, "latitude", &id_data);
    if (status==NC_NOERR) {
      strcpy (lat_name, "latitude");
      strcpy (lon_name, "longitude");
    }
    else {
      fprintf (stderr, "Error, neither 'lat' nor 'latitude' in %s\n", filename);
      fprintf (stderr, "      (line %d, function %s in %s)\n",  __LINE__, __func__, __FILE__);
      return status;
    }
  }

  /* read latitude array */
  status = alloc_and_read_netCDF_1D_double(ncid, lat_name, &nlat, lat_name, &(ECMWF_lat));
  if (status != 0) {
    fprintf (stderr, "Error %d reading latitude\n", status);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }

  /* search correct latitude index */
  ilat = get_grid_index(lat, ECMWF_lat, nlat, FALSE);
  if (ilat < 0) {
    fprintf (stderr, "Error -1 finding index for lat=%5.2f in %s\n", lat, filename);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  /* read longitude */
  status = alloc_and_read_netCDF_1D_double(ncid, lon_name, &nlon, lon_name, &(ECMWF_lon));
  if (status != 0) {
    fprintf (stderr, "Error %d reading longitude\n", status);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }

  /* search correct longitude index */
  ilon = get_grid_index(lon, ECMWF_lon, nlon, TRUE);
  if (ilon < 0) {
    fprintf (stderr, "Error -2 finding index for lon=%5.2f in %s\n", lon, filename);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -2;
  }

  /* get time index */ 
  status = get_time_index (ncid, UTC, time_interpolate, 
                           &(nt), &(itime1), &(itime2), &(dt),
                           verbose, quiet);
  if (status != 0){
    fprintf (stderr, "Error %d returned by function 'get_time_index'\n", status);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }  

  /* alloc nt timesteps for pressure (level) */
  if ((p_level = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for pressure\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -10;
  }

  /* alloc nt timesteps for pressure (layer) */
  if ((p_layer = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for pressure\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -10;
  }

  /* alloc nt timesteps for temperature */
  if ((T = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for temperature\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -10;
  }

  /* alloc nt timesteps for mass mixing ratio water vapour */
  if ((Q = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for water vapour\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -10;
  }

  /* alloc nt timesteps for mass mixing ratio Ozone */
  if ((O3 = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error allocating memory for Ozone\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -10;
  }


  /* read nt (= 1 or 2) time steps */
  for (t=0;t<=nt-1;t++) {

    if (t == 0)
      itime = itime1;
    if (t == 1)
      itime = itime2;

    /* alloc and read pressure */
    /* requires that SP/LNSP, hyai and hybi are in the ECMWF netCDF file */ /* in ancillary.c */
    status = alloc_and_read_ECMWF_netCDF_pressure(ncid, &(p_level[t]), &(nlev),&(p_layer[t]), &(nlay), itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading pressure from netCDF file %s\n", status, filename);
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }

    /* allocate and read temperature, defined on layers */
    status = alloc_and_read_netCDF_column_float(ncid, T_name, &T[t], nlay, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading temperature '%s' from netCDF file %s\n", status, T_name, filename);
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }

    /* allocate and read h2o mass mixing ratio, defined on layers */
    status = alloc_and_read_netCDF_column_float(ncid, Q_name, &Q[t], nlay, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading water vapour mass mixing ratio '%s' from netCDF file %s\n", status, Q_name, filename);
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }

    /* allocate and read h2o mass mixing ratio, defined on layers */
    status = alloc_and_read_netCDF_column_float(ncid, O3_name, &O3[t], nlay, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading ozone '%s' from netCDF file %s\n", status, O3_name, filename);
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }
  }

  if (verbose) {
    fprintf (stderr, "     found %zd x %zd x %zd (lat x lon x lev) data points\n", nlat, nlon, nlev);
    fprintf (stderr, "     reading pixel lat = %5.2f (%4ld), lon = %5.2f (%4ld)\n", lat, ilat, lon, ilon);
  }

  /* correct automatically the errors in the ECMWF files !!! */
  for (t=0;t<=nt-1;t++) {
    /* correct water vapour mixing ratio */
    if (Q[t][0] < 0.0 )
      Q[t][0] = 2.0e-06;
    for (lc=0; lc<nlay; lc++) {
      if (Q[t][lc] < 0.0) {
        /* search for next resonable entry */
        for (lc2=lc; lc2<nlay; lc2++) {
          if (0.0 < Q[t][lc2] && Q[t][lc2] < 20.0*Q[t][lc-1] ) {
            break;
          }
        }
        fprintf(stderr,"\n          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        fprintf(stderr,  "          !!! WARNING AUTOMATIC CORRECTION OF NEGATIV ECMWF INPUT DATA: WATER VAPOUR MIXING RATIO !!!\n");
        fprintf (stderr, "          !!! at pixel lat = %5.2f (%4ld), lon = %5.2f (%4ld)                                   !!!\n", lat, ilat, lon, ilon);
        fprintf(stderr,  "          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        if (lc2 < nlev) {
          for (lc3=lc; lc3<lc2; lc3++) {
            fprintf (stderr, "        original:  Q[lc=%3d]=%e, Q[lc2=%3d]=%e, Q[lc3=%3d]=%e\n", 
		     lc-1, Q[t][lc-1], lc2, Q[t][lc2], lc3, Q[t][lc3]);
            Q[t][lc3] = Q[t][lc-1] + (Q[t][lc2] - Q[t][lc-1])/( lc2-(lc-1) ) * (lc3 - (lc-1) );
            fprintf (stderr, "        corrected: Q[lc=%3d]=%e, Q[lc2=%3d]=%e, Q[lc3=%3d]=%e\n\n", 
		     lc-1, Q[t][lc-1], lc2, Q[t][lc2], lc3, Q[t][lc3]);
          }
        }
        else {
          fprintf (stderr, "          original:  Q[lc=%3d]=%e, Q[lc3=%3d]=%e\n", lc-1, Q[t][lc-1], lc, Q[t][lc]);
          Q[t][lc] = Q[t][lc-1];
          fprintf (stderr, "          corrected: Q[lc=%3d]=%e, Q[lc3=%3d]=%e\n\n", lc-1, Q[t][lc-1], lc, Q[t][lc]);
        }
      }
    }
    /* correct ozone mixing ratio */
    if (O3[t][0] < 0.0 )
      O3[t][0] = 2.0e-06;
    for (lc=0; lc<nlay; lc++) {
      if (O3[t][lc] < 0.0) {
        /* search for next resonable entry */
        for (lc2=lc; lc2<nlay; lc2++) {
          if (0.0 < O3[t][lc2] && O3[t][lc2] < 20.0*O3[t][lc-1] ) {
            break;
          }
        }
        fprintf(stderr,"\n          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        fprintf(stderr,  "          !!! WARNING AUTOMATIC CORRECTION OF NEGATIV ECMWF INPUT DATA: OZONE MIXING RATIO !!!\n");
        fprintf (stderr, "          !!! at pixel lat = %5.2f (%4ld), lon = %5.2f (%4ld)                            !!!\n", lat, ilat, lon, ilon);
        fprintf(stderr,  "          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        if (lc2 < nlev) {
          for (lc3=lc; lc3<lc2; lc3++) {
            fprintf (stderr, "          original:  O3[lc=%3d]=%e, O3[lc2=%3d]=%e, O3[lc3=%3d]=%e\n", lc-1, O3[t][lc-1], lc2, O3[t][lc2], lc3, O3[t][lc3]);
            O3[t][lc3] = O3[t][lc-1] + (O3[t][lc2] - O3[t][lc-1])/( lc2-(lc-1) ) * (lc3 - (lc-1) );
            fprintf (stderr, "          corrected: O3[lc=%3d]=%e, O3[lc2=%3d]=%e, O3[lc3=%3d]=%e\n\n", lc-1, O3[t][lc-1], lc2, O3[t][lc2], lc3, O3[t][lc3]);
          }
        }
        else {
          fprintf (stderr, "          original:  O3[lc=%3d]=%e, O3[lc3=%3d]=%e\n", lc-1, O3[t][lc-1], lc, O3[t][lc]);
          O3[t][lc] = O3[t][lc-1];
          fprintf (stderr, "          corrected: O3[lc=%3d]=%e, O3[lc3=%3d]=%e\n\n", lc-1, O3[t][lc-1], lc, O3[t][lc]);
        }
      }
    }
  }

  if (nt == 1) {
    /* no time interpolation needed, do nothing */
  }
  else {
    /* write time interpolated data into the zero'th entry */
    for (lc=0; lc<nlev; lc++) {
      p_level[0][lc] = (1.0-dt)* p_level[0][lc]  +  dt* p_level[1][lc];
    }
    for (lc=0; lc<nlay; lc++) {
      p_layer[0][lc] = (1.0-dt)* p_layer[0][lc]  +  dt* p_layer[1][lc];
      T      [0][lc] = (1.0-dt)*       T[0][lc]  +  dt*       T[1][lc];
      O3     [0][lc] = (1.0-dt)*      O3[0][lc]  +  dt*      O3[1][lc];
      Q      [0][lc] = (1.0-dt)*       Q[0][lc]  +  dt*       Q[1][lc];
      /* (*atm_data)[lc][0] = p_level[0][lc]; */
    }
  }

  /*   /\* additional verbose output data on layers *\/ */
  /*   fprintf (stderr, " ADDITIONAL_VERBOSE: ECMWF data at layer midpoints\n"); */
  /*   for (lc=0; lc<nlay; lc++) { */
  /*     fprintf (stderr,"lc=%3d, p_layer=%10.4f, p_level=%10.4f, T=%10.5f, O3=%12.6e, Q=%12.6e\n", */
  /*                      lc, p_layer[0][lc], p_level[0][lc], T[0][lc], O3[0][lc], Q[0][lc]); */
  /*   } */
  /*   fprintf (stderr,"lc=%3d, p_layer=%10.4f, p_level=%10.4f, T=%10.5f, O3=%12.6e, Q=%12.6e\n", */
  /*            nlev-1, 0.0/0.0, p_level[0][nlev-1], 0.0/0.0, 0.0/0.0, 0.0/0.0); */

  if (read_on_level_or_layer == level) {

    (*rows)=nlev;
    /* allocate data for "mimikry" radiosonde data */
    if ((status = ASCII_calloc_float (atm_data, (*rows), 4)) != 0) {
      fprintf (stderr," Error allocating memory for 'atm_data' \n"); 
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }

    /* inter- and extra-polate from LAYER to LEVEL values for of T, H2O and O3 on the log(p/p0) grid */
    /* interpolation */
    for (lc=1; lc<nlev-1; lc++) {
      dx = ( log ( p_level[0][lc] / p_level[0][nlev-1]) - log ( p_layer[0][lc-1] / p_level[0][nlev-1])) /
	( log ( p_layer[0][lc] / p_level[0][nlev-1]) - log ( p_layer[0][lc-1] / p_level[0][nlev-1]) );
      (*atm_data)[lc][0] = p_level[0][lc];
      (*atm_data)[lc][1] = (1.0-dx) *  T[0][lc-1] + dx *  T[0][lc];
      (*atm_data)[lc][2] = (1.0-dx) * O3[0][lc-1] + dx * O3[0][lc];
      (*atm_data)[lc][3] = (1.0-dx) *  Q[0][lc-1] + dx *  Q[0][lc];
    }
    /* extrapolate last layer midpoint to surface */
    dx = ( log ( p_level[0][nlev-1] / p_level[0][nlev-1]) - log ( p_layer[0][nlay-1] / p_level[0][nlev-1]) ) / 
      ( log ( p_layer[0][nlay-1] / p_level[0][nlev-1]) - log ( p_layer[0][nlay-2] / p_level[0][nlev-1]) );
    (*atm_data)[nlev-1][0] =  p_level[0][nlev-1];
    (*atm_data)[nlev-1][1] = (1.0+dx) *  T[0][nlay-1] - dx *  T[0][nlay-2];
    (*atm_data)[nlev-1][2] = (1.0+dx) * O3[0][nlay-1] - dx * O3[0][nlay-2];
    (*atm_data)[nlev-1][3] = (1.0+dx) *  Q[0][nlay-1] - dx *  Q[0][nlay-2];
    /* assume constant mixing ratio between last layer midpoint and TOA */
    (*atm_data)[0][0] =  p_level[0][0];
    (*atm_data)[0][1] =  T [0][0];
    (*atm_data)[0][2] =  O3[0][0];
    (*atm_data)[0][3] =  Q [0][0];
  }
  else if (read_on_level_or_layer == layer) {

    (*rows)=nlay+1;

    /* allocate data for "mimikry" radiosonde data */
    if ((status = ASCII_calloc_float (atm_data, (*rows), 4)) != 0) {
      fprintf (stderr," Error allocating memory for 'atm_data'\n"); 
      fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
      return status;
    }
    /* copy all layer values */
    for (lc=0; lc<nlay; lc++) {
      (*atm_data)[lc][0] = p_layer[0][lc];
      (*atm_data)[lc][1] = T [0][lc];
      (*atm_data)[lc][2] = O3[0][lc];
      (*atm_data)[lc][3] = Q [0][lc];
    }
    /* extrapolate last layer midpoint to surface */
    dx = ( log ( p_level[0][nlev-1] / p_level[0][nlev-1]) - log ( p_layer[0][nlay-1] / p_level[0][nlev-1]) ) / 
      ( log ( p_layer[0][nlay-1] / p_level[0][nlev-1]) - log ( p_layer[0][nlay-2] / p_level[0][nlev-1]) );
    (*atm_data)[nlay][0] =  p_level[0][nlev-1];
    (*atm_data)[nlay][1] = (1.0+dx) *  T[0][nlay-1] - dx *  T[0][nlay-2];
    (*atm_data)[nlay][2] = (1.0+dx) * O3[0][nlay-1] - dx * O3[0][nlay-2];
    (*atm_data)[nlay][3] = (1.0+dx) *  Q[0][nlay-1] - dx *  Q[0][nlay-2];
  }
  else {
    fprintf (stderr, "Error, unknown reading mode for ECMWF data\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }

  /* extrapolation might produce negative mixing ratios / number densities  */
  /* That's not very nice, but we cannot do better than correct it to zero. */
  if ((*atm_data)[nlev-1][2] < 0.0)
    (*atm_data)[nlev-1][2] = 0.0;
  if ((*atm_data)[nlev-1][3] < 0.0)
    (*atm_data)[nlev-1][3] = 0.0;

  if (verbose) {
    fprintf (stderr, " *** ECMWF atmosphere data \n");
    fprintf (stderr, "# lc |  Pressure  | Temp.  |   Ozone    | Water vap.  |\n");
    fprintf (stderr, "#    |   [hPa]    |  [K]   |   [mmr]    |   [mmr]     |\n");
    fprintf (stderr, "# -----------------------------------------------------\n");
    for (lc=0;lc<(*rows);lc++)
      fprintf (stderr, "%5d  %10.5f   %6.2f   %.5e   %.5e\n", 
	       lc, (*atm_data)[lc][0], (*atm_data)[lc][1],(*atm_data)[lc][2],(*atm_data)[lc][3]);
  }

  /* we look for the cold point troposphere dealing with layer data        */
  /* we do that with LAYER data, this is more accurate than to do it later */
  /* with data, which is alread interpolated from layers to levels         */
  if (verbose)
    fprintf (stderr, " ... search ECMWF cold point tropopause\n");

  if ((temper_layer = calloc (nlay+1, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for 'temper_layer'\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -1;
  }
  if ((press_layer = calloc (nlay+1, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for 'press_layer'\n");
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  /* copy layer data everywhere */
  for (lc=0; lc<nlay; lc++) {
    temper_layer [lc] = T[0][lc];
    press_layer  [lc] = p_layer[0][lc];
  }
  /* last entry is the surface level */
  temper_layer[nlay] = (*atm_data)[nlev-1][1];
  press_layer [nlay] = p_level[0][nlev-1];

  /* first calculate the height of the layer midpoints */
  (*nz_model_layer) = nlay+1;
  status = calculate_z_from_p_and_T (press_layer, temper_layer, z_model_layer, (*nz_model_layer), 
                                     altitude, verbose);
  if (status!=0) {
    fprintf (stderr, "Error %d during calculate_z_from_p_and_T\n", status);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }
  /* we skip the last entry, that is the surface (was only there for height integration) */
  (*nz_model_layer)= nlay;

  status = find_cold_point_tropopause ( press_layer, temper_layer, (*z_model_layer), (*nz_model_layer), 
                                        z_cpt_sea, press_cpt, temper_cpt, quiet, verbose);
  if ( status != 0 ) {
    if (!quiet)
      fprintf (stderr, "*** Warning, no ECMWF cold point tropopause was found.\n");
    /* we continue even if we do not find an cold point tropopause */
    status=0;
  }

  free(press_layer);
  free(temper_layer);

  for (t=0;t<=nt-1;t++) {
    free(p_level[t]);
    free(p_layer[t]);
    free(T[t]);
    free(Q[t]);
    free(O3[t]);
  }
  free(p_level);
  free(p_layer);
  free(T);
  free(Q);
  free(O3);

#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use the ECMWF data file option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ******************************************************************\n");
  return -1;
#endif

  return status;
}



/**************************************************************/
/* Read density file                                          */
/* Combine with background atmosphere and                     */
/* Interpolate on the common grid: z_all                      */
/**************************************************************/

static int read_and_combine_density (char *filename, atm_out_struct *out,
				     int dens_id, int *interpol_method_gas, int interpol_method_temper, 
                                     int unit_profile, float *mol_mass, 
				     float *dens_air_all, float *temper_all, float *z_all, int nlev_all,
                                     float *dens_air_atm, float *z_atm, int nlev_atm,
                                     int solver, int quiet, int verbose)

{
  int rows=NOT_DEFINED_INTEGER, min_columns=NOT_DEFINED_INTEGER, max_columns=NOT_DEFINED_INTEGER, nlev=NOT_DEFINED_INTEGER;
  int status = 0;
  float **atm_data = NULL;
  int lc           = NOT_DEFINED_INTEGER;
  int i            = NOT_DEFINED_INTEGER; 
  float *zd        = NULL;
  float **dens     = NULL;
  static int first_dens = 1;
  float *z_comb         = NULL;
  float **dens_comb     = NULL;
  int nlev_comb         = NOT_DEFINED_INTEGER;
  float *dens_air_comb  = NULL;
  int iv=0, is=0, isthis=0;
  int rows_in_file, lc1, icolumn;

  if ((strlen(filename) == 0) && (out->nlev == nlev_all)) {
    /* nothing to do, everything is fine */
  }
  else { /* something to do */
    
    /* densfile specified in the input file */
    if (strlen(filename) > 0) {

      if (verbose)
        fprintf (stderr, " ... reading mol_file %s\n", filename);

      /* read densfile */
      status = ASCII_file2float (filename, 
                                 &rows, &max_columns, &min_columns, 
                                 &atm_data);
  
      if (status!=0) {
        fprintf (stderr, "Error %d reading %s\n", status, filename);
        return status;
      }

      /* remember the number of rows (it is required for deallocation of atm_data) */
      rows_in_file=rows;

      /* check if there is an overlap between mol_file and atmosphere_file */ 
      if((atm_data[rows-1][0]>=z_all[0])||(atm_data[0][0]<=z_all[nlev_all-1])){
        fprintf(stderr,"Error, no overlap between mol_file %s and atmosphere_file.\n", filename);
        return -1;
      }
      
      /* Go through the levels in mol_file because levels not covered by the z_all range need a 
       * special treatment: They are interpolated to boundaries of the z_all range or sorted out. */
      for (lc=0; lc<rows; lc++) 
        if(atm_data[lc][0]>z_all[0] || atm_data[lc][0]<z_all[nlev_all-1]){ /* test if current mol_file level is outside range from z_all[0] to z_all[nlev_all-1] */

          /* linear interpolation of mol_file level to upper boundary in z_all */
          if(atm_data[lc][0]>z_all[0] && atm_data[lc+1][0]<z_all[0]){ /* test if z_all[0] is between current mol_file level and next mol_file level */
            fprintf(stderr," ... Linear interpolation of data from mol_file %s to upper boundary in z_all (%ekm).\n", filename, z_all[0]);
            for (icolumn=1; icolumn<max_columns; icolumn++){
              if((atm_data[lc][icolumn]<0)||(atm_data[lc+1][icolumn]<0)){
                fprintf(stderr,"Error, density=-1 in mol_file %s not supported in combination with interpolation at boundary of z_all.\n", filename);
                return -1;
              }
              atm_data[lc][icolumn]=atm_data[lc][icolumn]*(atm_data[lc+1][0]-z_all[0])+atm_data[lc+1][icolumn]*(z_all[0]-atm_data[lc][0]);
              atm_data[lc][icolumn]=atm_data[lc][icolumn]/(atm_data[lc+1][0]-atm_data[lc][0]);
            }
            atm_data[lc][0]=z_all[0];
          }
          /* linear interpolation of mol_file level to lower boundary in z_all */
          else if(atm_data[lc][0]<z_all[nlev_all-1] && atm_data[lc-1][0]>z_all[nlev_all-1]){ /* test if z_all[nlev_all-1] is between current mol_file level and previous mol_file level */
            fprintf(stderr," ... Linear interpolation of data from mol_file %s to lower boundary in z_all (%ekm).\n", filename, z_all[nlev_all-1]);
            for (icolumn=1; icolumn<max_columns; icolumn++){
              if((atm_data[lc-1][icolumn]<0)||(atm_data[lc][icolumn]<0)){
                fprintf(stderr,"Error, density=-1 in mol_file %s not supported in combination with interpolation at boundary of z_all.\n", filename);
                return -1;
              }
              atm_data[lc][icolumn]=atm_data[lc-1][icolumn]*(atm_data[lc][0]-z_all[nlev_all-1])+atm_data[lc][icolumn]*(z_all[nlev_all-1]-atm_data[lc-1][0]);
              atm_data[lc][icolumn]=atm_data[lc][icolumn]/(atm_data[lc][0]-atm_data[lc-1][0]);
            }
            atm_data[lc][0]=z_all[nlev_all-1];
          }
          /* sort out mol_file level */
          else{
            if(!quiet) 
              fprintf(stderr," ... Ignoring level %e km from mol_file %s because it is outside altitude grid z_all.\n", atm_data[lc][0], filename);
            for (lc1=lc; lc1<rows-1; lc1++)
              for (icolumn=0; icolumn<max_columns; icolumn++)
                atm_data[lc1][icolumn]=atm_data[lc1+1][icolumn]; /* overwrite level from mol_file with the levels above */
            rows--;
            lc--;
          }
	}

      if (max_columns==2) {   /* 2-columns: 1 altitude, 2 density */
    
        /* reading data from file */
        if (max_columns!=min_columns) {
          fprintf (stderr, "Error, inconsistent number of columns in %s\n", filename);
          fprintf (stderr, "     min = %d, max =%d\n", min_columns, max_columns);
          fprintf (stderr, "     Should have been %d. First altitude, second density.\n", 2);
          return -1;
        }
 
        nlev  = rows;

        /* copy dens-file z-grid to zd */
        zd   = ASCII_column_float (atm_data, nlev, 0);

        /* allocate dens-array */
        status  = ASCII_calloc_float(&dens, 1, nlev);
        if (status!=0) {
	  fprintf (stderr, "Error %d allocating memory for dens\n", status); 
	  return status; 
        }

        /* lc is the 2nd index */
        for (lc=0; lc<nlev; lc++) dens[0][lc] = atm_data[lc][1];
        
        /* convert units of the dens_file to cm-3 */
        status =  convert_profile_units (&(dens[0]), zd, nlev, dens_id,
					 dens_air_all, temper_all, z_all, nlev_all,
                                         interpol_method_gas[MOL_AIR],
					 interpol_method_temper, 
                                         unit_profile, mol_mass, quiet);
	
        if (status!=0) {
          fprintf (stderr, "Error %d returned by convert_profile_unit()\n", status);
          return status;
        }

        /* combine file data and background */
        status = combine_profiles (dens, zd, nlev, 1,
                                   out->microphys.dens[dens_id][0][0], dens_air_all, out->zd, out->nlev,
                                   interpol_method_gas[dens_id], interpol_method_gas[MOL_AIR],
                                   &(dens_comb), &(z_comb), &(nlev_comb), quiet);

        if (status!=0) {
          fprintf (stderr, "Error %d returned by combine_profiles()\n", status);
          return status;
        }

        /* check, if dens = 0 in the uppermost level, for log-interpolation types */
        if (interpol_method_gas[dens_id] == INTERP_METHOD_LOG || 
	    interpol_method_gas[dens_id] == INTERP_METHOD_LOG_SPLINE) {
          for (i=0; i<1; i++)
            if (dens_comb[i][0] == 0.0) {
              fprintf (stderr, "Error, dens(%5.1f km) = 0.0 cannot be interpolated with log assumption\n",z_comb[0]);
              fprintf (stderr, "  ***  By removing this line from the file %s, values for\n", filename);
              fprintf (stderr, "  ***  heights above the last entry in densfile will be\n");
              fprintf (stderr, "  ***  filled with the values given by the background atmosphere.\n");
              fprintf (stderr, "       in read_and_combine_density (atmosphere.c)\n");
              return -1;
            }
        }
        
        /* first determine air_dens profile on z_comb for linmix interpolation */
        if (interpol_method_gas[dens_id] == INTERP_METHOD_LINMIX) {
          dens_air_comb = (float *) calloc (nlev_comb, sizeof(float));
          status = arb_wvn (nlev_all, z_all,  dens_air_all, 
			    nlev_comb,  z_comb, dens_air_comb, 
			    interpol_method_gas[MOL_AIR], 1);
	  
	  if (status!=0) {
	    fprintf (stderr, "Error %d returned by arb_wvn()\n", status);
	    return status;
	  }
	}
	
        /* interpolate combined data on the z_all grid */
        status = interpolate_density (z_comb, &(dens_comb[0]), nlev_comb, 
				      z_all, nlev_all, interpol_method_gas[dens_id],
				      dens_air_comb, dens_air_all, quiet);

        if (status!=0) {
          fprintf (stderr, "Error %d returned by interpolate_density()\n", status);
          return status;
        }

        free(dens_air_comb);

        /* put result on final array */
        free(out->microphys.dens[dens_id][0][0]);
        out->microphys.dens[dens_id][0][0] = (float *) calloc (nlev_all, sizeof (float));
    
        for (lc=0; lc<nlev_all; lc++)
          out->microphys.dens[dens_id][0][0][lc] = dens_comb[0][lc];

        status  = ASCII_free_float(dens, 1);

      }

      else if (max_columns > 2) {   /* density matrix for airmass calculations */

        if (max_columns!=min_columns) {
          fprintf (stderr, "Error, inconsistent number of columns in %s\n", filename);
          fprintf (stderr, "     min = %d, max =%d\n", min_columns, max_columns);
          fprintf (stderr, "     in read_and_combine_density (atmosphere.c)\n");
          return -1;
        }

        if (solver != SOLVER_SDISORT) {
          fprintf (stderr, "Error, dens_file %s has the format of a matrix.\n", filename);
          fprintf (stderr, "            This input type is only possible for solver SDISORT.\n");
          fprintf (stderr, "            in read_and_combine_density (atmosphere.c)\n");
          return -1;
        }

    
        /* first row/column contains sza/altitude */
        out->microphys.nsza_denstab = max_columns-1;
        nlev = rows-1;

        zd  = (float *) calloc (nlev, sizeof (float));

        if (out->microphys.denstab_id > 0) {
          fprintf (stderr, "Error, dens_file with matrix instead of profile\n");
          fprintf (stderr, "       already specified. The matrix may be used\n");
          fprintf (stderr, "       only for one specie. Profiles for as many as you like\n");
          fprintf (stderr, "       in read_and_combine_density (atmosphere.c)\n");
          return -1;
        }

        if (first_dens) {
          out->microphys.denstab_id  = dens_id;
          out->microphys.sza_denstab = (float *) calloc (out->microphys.nsza_denstab, sizeof (float));

          status  = ASCII_calloc_float(&out->microphys.denstab_amf,
                                       out->microphys.nsza_denstab, nlev);
          if (status!=0) {
            fprintf (stderr, "Error %d allocating memory for out->denstab_amf\n", status);
            return status;
          }
        }

        /* copy first row to solar zenith angle array  */
        for (i=1; i<=out->microphys.nsza_denstab; i++)
          out->microphys.sza_denstab[i-1] = atm_data[0][i];
    
        /* copy first column to altitude array zd */
        for (lc=1; lc<=nlev; lc++)
          zd[lc-1] = atm_data[lc][0];


        /* index of lc has to be at the last position */
        for (i=1; i<=out->microphys.nsza_denstab; i++)
          for (lc=1; lc<=nlev; lc++)
            out->microphys.denstab_amf [i-1][lc-1] = atm_data[lc][i];


        /* convert units of the dens_file to cm-3 */
        for (i=0; i<out->microphys.nsza_denstab; i++) {
          status =  convert_profile_units (&(out->microphys.denstab_amf[i]), zd, nlev, dens_id, 
                                           dens_air_all, temper_all, z_all, nlev_all,
                                           interpol_method_gas[MOL_AIR], interpol_method_temper, 
                                           unit_profile, mol_mass, quiet);
          if (status!=0) {
            fprintf (stderr, "Error %d returned by convert_profile_unit()\n", status);
            return status;
          }
        }

        /* combine file data and background */
        status = combine_profiles (out->microphys.denstab_amf, zd, nlev, out->microphys.nsza_denstab,
                                   out->microphys.dens[dens_id][0][0], dens_air_all, out->zd, out->nlev,
                                   interpol_method_gas[dens_id], interpol_method_gas[MOL_AIR],
                                   &(dens_comb), &(z_comb), &(nlev_comb), quiet);

        if (status!=0) {
          fprintf (stderr, "Error %d returned by combine_profiles()\n", status);
          return status;
        }

	/* 	for (lc=1; lc<=nlev_all; lc++) */
	/*           fprintf (stderr, "lc = %d, z_comb = %5.2f, dens_comb %11.7e\n", lc, z_comb[lc], dens_comb[0][lc]); */

        /* check, if dens = 0 in the uppermost level, for log-interpolation types */
        if (interpol_method_gas[dens_id] == INTERP_METHOD_LOG || 
	    interpol_method_gas[dens_id] == INTERP_METHOD_LOG_SPLINE) {
          for (i=0; i<out->microphys.nsza_denstab; i++)
            if (dens_comb[i][0] == 0.0) {
              fprintf (stderr, "Error, dens(%5.1f km) = 0.0 cannot be interpolated with log assumption\n",z_comb[0]);
              fprintf (stderr, "  ***  By removing this line from the file %s, values for\n", filename);
              fprintf (stderr, "  ***  hights above the last entry in densfile will be\n");
              fprintf (stderr, "  ***  filled with the values given by the background atmosphere.\n");
              fprintf (stderr, "       in read_and_combine_density (atmosphere.c)\n");
              return -1;
            }
        }

        /* first determine air_dens profile on z_comb for linmix interpolation */
        if (interpol_method_gas[dens_id] == INTERP_METHOD_LINMIX) {
          dens_air_comb = (float *) calloc (nlev_comb, sizeof(float));
          status = arb_wvn (nlev_all, z_all,  dens_air_all, 
			    nlev_comb,  z_comb, dens_air_comb, 
			    interpol_method_gas[MOL_AIR], 1);
	  if (status!=0) {
	    fprintf (stderr, "Error %d returned by arb_wvn()\n", status);
	    return status;
	  }
	}

        /* interpolate combined data on the z_all grid */
        for (i=0; i<out->microphys.nsza_denstab; i++) {
          status = interpolate_density (z_comb, &(dens_comb[i]), nlev_comb, 
					z_all, nlev_all, interpol_method_gas[dens_id],
					dens_air_comb, dens_air_all, quiet);

          if (status!=0) {
            fprintf (stderr, "Error %d interpolating denstab (dens_comb) to z_all grid\n", status); 
            return status;
          }
        }

        free(dens_air_comb);

        /* reallocate denstab_amf*/
        ASCII_free_float(out->microphys.denstab_amf, out->microphys.nsza_denstab);
        status  = ASCII_calloc_float(&out->microphys.denstab_amf,
				     out->microphys.nsza_denstab, nlev_all);
        if (status!=0) {
          fprintf (stderr, "Error %d allocating memory for out->denstab_amf\n", status);
          return status;
        }

        /* copy combined data on final place  */
        for (i=0; i<out->microphys.nsza_denstab; i++)
          for (lc=0;lc<nlev_all;lc++)
            out->microphys.denstab_amf[i][lc]=dens_comb[i][lc];
           
        /* callocate denstab array */
        if (first_dens) {
          status  = ASCII_calloc_float_3D (&out->microphys.denstab, MOL_NN, out->microphys.nsza_denstab, nlev_all);
          if (status!=0) {
            fprintf (stderr, "Error %d allocating memory for out->denstab\n", status);
            return status;
          }
          status  = ASCII_calloc_float_3D (&out->microphys.denstab_avg, MOL_NN, out->microphys.nsza_denstab, nlev_all);
          if (status!=0) {
            fprintf (stderr, "Error %d allocating memory for out->denstab_avg\n", status);
            return status;
          }
        }

        /* copy the combined 2-dim double array to denstab */
        for (i=0; i<out->microphys.nsza_denstab; i++)
          for (lc=0; lc<nlev_all; lc++)
            out->microphys.denstab [dens_id][i][lc] = out->microphys.denstab_amf[i][lc];

        isthis = NOT_DEFINED_INTEGER;
        /* replace the profile in the atmosphere-data according to sza */
        iv=0;  /* No solar zenith variations for a single amf calculation */
        for (is=0; is<out->microphys.nsza_denstab; is++)
          if (fabs (out->microphys.sza_denstab[is]-out->sza_r[iv]) < 1.e-04)
	    isthis = is;

        if (isthis < 0) {
          /* sza not found in the denstab-file */
          if (!quiet) {
            fprintf (stderr, "*** Warning: requested sza is not specified in the denstab-file.\n");
            fprintf (stderr, "*** Using the trace gas profile from the atmosphere_file.\n");
          }

          /* interpolate data of the background-atmosphere */
          /* from atm-file grid onto the z_all-grid  */
          status = interpolate_density (z_atm, &(out->microphys.dens[dens_id][0][0]),
					nlev_atm, 
					z_all, nlev_all, interpol_method_gas[dens_id],
					dens_air_atm, dens_air_all, quiet);
          if (status!=0) {
            fprintf (stderr, "Error %d interpolating background denstab to z_all grid\n", status);
            return status;
          }
        }
        else {
          /* normal case, sza found in the denstab file */
          if (verbose)
            fprintf (stderr, " ... found denstab profile for sza = %6.2f\n", out->sza_r[0]);

          free(out->microphys.dens[dens_id]);
          out->microphys.dens[dens_id][0][0] = (float *) calloc (nlev_all, sizeof (float));

          for (lc=0; lc<nlev_all; lc++) {
            out->microphys.dens[dens_id][0][0][lc] = out->microphys.denstab[dens_id][isthis][lc];
          }
        }

        first_dens = 0;

      } /* if (max_columns==2) */

      free(zd);
      free(z_comb);
      free(dens_comb);
  
      ASCII_free_float(atm_data, rows_in_file);

    }       /* if (strlen(filename) > 0) */
    else {  /* (there is no densfile, but z-grids are different)    */

      /* interpolation from atm-file grid onto the z_all-grid  */
      status = interpolate_density (z_atm, &(out->microphys.dens[dens_id][0][0]), nlev_atm, 
				    z_all, nlev_all, interpol_method_gas[dens_id],
				    dens_air_atm, dens_air_all, quiet);
      if (status!=0) {
	fprintf (stderr, "Error %d returned by interpolate_density (gas =%3d)\n", status, dens_id); 
	return status;
      }     
      
    } /* if (strlen(filename) > 0) */
    
  } /* else { ( something to do ) */ 

  return status;
}


/*******************************************/
/* Scale pressure and all well mixed gases */
/*******************************************/

static int scale_pressure (float ****pressure, float *****dens, int nlev, float scale_factor, int *well_mixed_gas, int ialt, int verbose)
{
  int lc     = NOT_DEFINED_INTEGER;
  int gas_nr = NOT_DEFINED_INTEGER;
  char *gas  = NULL;

  for (lc=0;lc<nlev;lc++)
    (*pressure)[0][0][lc] *= scale_factor;
  if (verbose)
    fprintf (stderr, " ... scaling pressure    from  %11.6f hPa   to %11.6f hPa\n", 
	     (*pressure)[0][0][ialt]/scale_factor, (*pressure)[0][0][ialt]);

  for (gas_nr=0;gas_nr<MOL_NN;gas_nr++)
    if ( well_mixed_gas[gas_nr] == YES ) { 
      gas = gas_number2string(gas_nr);
      if (strcmp(" -1 ",gas)==0) {
        fprintf (stderr, "Error wrong gas identifier in atmophere.c\n");
        return -1;
      }

      for (lc=0;lc<nlev;lc++) 
        (*dens)[gas_nr][0][0][lc] *= scale_factor;
      if (verbose)
        fprintf (stderr, " ... scaling %s density from %11.5g cm-3  to %11.5g cm-3\n", 
		 gas, (*dens)[gas_nr][0][0][ialt]/scale_factor, (*dens)[gas_nr][0][0][ialt]);
    }
  
  return 0;
}



/**************************************************************/
/* Scale a density profile from it own column to input.column */
/**************************************************************/

static int scale_density_profile (float **dens, float **dens_avg, int alloc_avg, int i,
				  float input_column, float *output_column, 
				  float *scale_factor, int interpol_method_gas, 
				  input_struct input, output_struct *output)
{
  char *gas;
  int  lc = NOT_DEFINED_INTEGER;
  int status=0;

  gas = gas_number2string(i);
  if (strcmp(" -1 ",gas)==0) {
    fprintf (stderr, "Error wrong gas identifier in atmophere.c\n");
    return -1;
  }

  if (i == MOL_AIR) {
    fprintf (stderr, "Error you can't use the function scale_density_profile() for\n");
    fprintf (stderr, "      AIR, please have a look into setup_atmosphere() (atmosphere.c)\n");
    fprintf (stderr, "     -----\n");
    return -2;
  }

  if (i == MOL_H2O) {
    fprintf (stderr, "Error when using function scale_density_profile for\n");
    fprintf (stderr, "      WATER VAPOUR, saturate water vapour inside cloud is not considered.\n");
    fprintf (stderr, "     -------------- That's why, it is forbitten to use it for water.\n");
    fprintf (stderr, "                    Please have a look into setup_atmosphere() (atmosphere.c)\n");
    return -3;
  }

  status = column ((*dens), output->atm.zd, output->atm.nlev, output->alt.altitude, 
                   interpol_method_gas, output->atm.microphys.dens[MOL_AIR][0][0],
                   &(*dens_avg), alloc_avg, &(*output_column));

  if (status!=0) {
    fprintf (stderr, "Error %d calculating %s column\n", status, gas);
    return status;
  }

  if (input_column >= 0.0) {
    switch (input.atm.unit_column[i]) {
    case MOL_UNIT_DU:
      /* Everything ok, do nothing */
      if (!input.quiet)
        fprintf (stderr, " ... scaling %s column from %11.5g DU    to ", 
	         gas, (*output_column));
      break;
    case MOL_UNIT_CM_2:
      if (!input.quiet)
        fprintf (stderr, " ... scaling %s column from %11.5g cm^-2 to ", 
	         gas, (*output_column)*DU2particle);
      input_column /= DU2particle; /* not changed outside this function as input is not a pointer, UH */
      break;
    default:
      fprintf (stderr, "Error, unknown unit for %s column\n", gas);
      return -4;
    }
    
    /* determine scale factor */
    if ((*output_column) > 0.0)
      (*scale_factor) = input_column / (*output_column);
    else
      (*scale_factor) = 0;

    /* scale profile */
    for (lc = 0; lc<output->atm.nlev; ++lc)
      (*dens)[lc] *= (*scale_factor);

    status = column ((*dens), output->atm.zd, output->atm.nlev, output->alt.altitude, 
                     interpol_method_gas, output->atm.microphys.dens[MOL_AIR][0][0],
                     &(*dens_avg), NO, output_column);

    if (status!=0) {
      fprintf (stderr, "Error %d calculating %s column!\n", status, gas);
      return status;
    }

    switch (input.atm.unit_column[i]) {
    case MOL_UNIT_DU:
      if (!input.quiet) fprintf (stderr, "%11.5g DU\n", (*output_column));
      break;
    case MOL_UNIT_CM_2:
      if (!input.quiet) fprintf (stderr, "%11.5g cm^-2\n", (*output_column)*DU2particle);
      break;
    default:
      fprintf (stderr, "Error, unknown unit for %s column\n", gas);
      return -4;
    }
  }

  free(gas);

  return status;

}


/******************************************/
/* Converts the unit of a profile to cm-3 */
/******************************************/

static int convert_profile_units    (float **ptr_dens, float *z_dens, int nlev, int dens_id,   
                                     float *dens_air_atm, float *temper_atm, float *z_atm, int nlev_atm,
                                     int interpol_method_air, int interpol_method_temper,
                                     int unit_profile, float *mol_mass, int quiet)

{
  int    lc=0;
  float  mol_mass_ratio=0.0; 
  float *dens_air_filegrid = NULL; 
  float *temper_filegrid   = NULL;
  int    status=0;

  switch (unit_profile) {
  case CM_3: 
    /* do nothing, particles unit is already fine, particles per cubic centimeter */
    break;
  case M_3: /* particles per cubic meter */
    for (lc=0; lc<nlev; lc++)
      if ((*ptr_dens)[lc] != -1.0)
	(*ptr_dens)[lc] *= 1.E-6; /* m-3 -> cm-3 */
    break;
  case MMR: /* mass mixing ratio   */
  case VMR: /* volume mixing ratio */

    /* multiply mixing ratio with dens_air (and if nessesary molecular mass ratio) */

    /* interpolate n_air to z-grid of the dens file */
    dens_air_filegrid = (float *) calloc (nlev_atm, sizeof (float));
    if (dens_air_filegrid == NULL) {
      fprintf (stderr, "Error allocating memory for dens_air_filegrid\n");
      return -1;
    }
    for (lc=0; lc<nlev_atm; lc++) dens_air_filegrid[lc] = dens_air_atm[lc];
    
    status = interpolate_profile (z_atm, &(dens_air_filegrid), nlev_atm, 
				  z_dens, nlev, interpol_method_air, quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d interpolating density of air to z-grid of the dens-file\n", status);
      return status;
    }

    /* determine molecular mass ratio in order to convert mass mixing ratios */
    if (unit_profile == MMR)
      mol_mass_ratio = mol_mass[MOL_AIR]/mol_mass[dens_id];
    else if (unit_profile == VMR) 
      mol_mass_ratio = 1.0;
    else { 
      fprintf (stderr, "Error. Unknown unit for dens_file = %d!\n", unit_profile);
      return -1;
    }

    for (lc=0; lc<nlev; lc++)
      if ((*ptr_dens)[lc] != -1.0)
	(*ptr_dens)[lc] *= (dens_air_filegrid[lc] * mol_mass_ratio);

    free(dens_air_filegrid);

    break;
  case RH: /* relative humidity, only for water vapour */
    if (dens_id != MOL_H2O) {
      fprintf (stderr, "Error, unit RH is only allowed for water vapour!\n");
      return -2;
    }

    /* interpolate Temperature to z-grid of the dens file */
    temper_filegrid = (float *) calloc (nlev_atm, sizeof (float));
    if (temper_filegrid == NULL) {
      fprintf (stderr, "Error allocating memory for temper_filegrid\n");
      return -1;
    }
    for (lc=0; lc<nlev_atm; lc++) temper_filegrid[lc] = temper_atm[lc];
      
    status = interpolate_profile (z_atm, &(temper_filegrid), nlev_atm, 
				  z_dens, nlev, interpol_method_temper, quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d interpolating temperature to z-grid of the dens-file\n", status);
      return status;
    }

    for (lc=0; lc<nlev; lc++)
      if ((*ptr_dens)[lc] != -1.0) {
	if ( 0.0 <= (*ptr_dens)[lc] && (*ptr_dens)[lc] < 250.0)
	  (*ptr_dens)[lc] *=  vapor_pressure(temper_filegrid[lc]) / 100.0;
	else {
	  fprintf (stderr, "Error, relative humidity is outside of accepted range [0.0 ; 250]\n");
	  fprintf (stderr, "       relative humidity = %e\n", (*ptr_dens)[lc]);
	  return -3;
	}
      }

    free(temper_filegrid);

    break;
  default:
    fprintf (stderr, "Error. Unknown unit for dens_file = %d!\n", unit_profile);
    fprintf (stderr, "This is a program bug, please contact the programmers!\n");
    fprintf (stderr, "Please include the used input-file!\n");
    return -3;
    break;
  }
  return 0;
}


/**************************************/
/* Combine two profiles               */
/**************************************/

static int combine_profiles (float **dens, float *z, int nlev, int nsza, 
                             float *dens_atm, float *dens_air_atm, float *z_atm, int nlev_atm, 
                             int interpol_method_dens, int interpol_method_air,
                             float ***dens_comb, float **z_comb, int *nlev_comb, int quiet)

{  /* _atm means z-grid of the atmosphere files */
   /* comb means combined                */
   /* output are dens_comb and z_comb    */

  int lc=0,i=0;
  int n_up=-666, n_dn=-666;
  float *dens_atm_filegrid=NULL;
  float *dens_air_filegrid=NULL;
  int    status = 0;

  /* z-range of dens-file is larger than those of the background atmosphere file */
  /* reduce the data to the range of the background atmosphere file */
  if (z[0] > z_atm[0] || z[nlev-1] < z_atm[nlev_atm-1]) {
    fprintf (stderr, "Error in function combine_profiles (atmosphere.c)\n");  
    fprintf (stderr, "     Range of the densfile z = [%f,%f] is larger than\n",z[nlev-1],z[0]);
    fprintf (stderr, "     those of the atmosphere file z_atm = [%f,%f]\n",z_atm[nlev_atm-1],z_atm[0]);
    fprintf (stderr, "     Can't do interpolation, as information of dens_air on filegrid is missing\n");
    return -1;
  }


  /* determine AIR density on background atmosphere grid and filegrid */
  /* copy air density of the background atmosphere */
  dens_air_filegrid = calloc (nlev_atm, sizeof (float));
  for (lc=0; lc<nlev_atm; lc++)
    dens_air_filegrid[lc] = dens_air_atm[lc];

  /* interpolate air density to the z-grid of the dens_file */  
  status = interpolate_profile (z_atm, &(dens_air_filegrid), nlev_atm, 
				z, nlev, interpol_method_air, quiet);

  if (status!=0) {
    fprintf (stderr, "Error %d interpolating air density of air to z-grid of the dens-file\n", status);
    return status;
  }

  /* copy density of the background atmosphere */
  dens_atm_filegrid = calloc (nlev_atm, sizeof (float));
  for (lc=0; lc<nlev_atm; lc++)
    dens_atm_filegrid[lc] = dens_atm[lc];

  /* interpolate density to the z-grid of the dens-file*/
  status = interpolate_density (z_atm, &(dens_atm_filegrid), nlev_atm, 
				z, nlev, interpol_method_dens,
				dens_air_atm, dens_air_filegrid, quiet);
  
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating dens of backgroud atmosphere to z-grid of the dens-file\n", status);
    return status;
  }

  if (z[0] == z_atm[0] &&  z[nlev-1] == z_atm[nlev_atm-1]) {
    /* upper and lower boundary are identical, no combination nessesary */
    *nlev_comb = nlev;

    *z_comb    = calloc (*nlev_comb, sizeof (float));
    *dens_comb = calloc (nsza,      sizeof (float *));
    for (i=0; i<nsza; i++)
      (*dens_comb)[i] = calloc (*nlev_comb, sizeof (float));

    for (lc=0; lc<*nlev_comb; lc++)
      (*z_comb)[lc]    = z[lc];

    for (i=0; i<nsza; i++) 
      for (lc=0; lc<*nlev_comb; lc++)
        if (dens[i][lc] != -1.0)
          (*dens_comb)[i][lc] = dens[i][lc];
	else
          (*dens_comb)[i][lc] = dens_atm_filegrid[lc];
  }
  else {

    /* get z_atm index above and below the dens-data */
    n_up = 0;
    for  (lc=0; lc<nlev_atm; lc++) {
      if(z_atm[lc] <= z[0]) {        
        n_up = lc-1;     /* last index of z_atm above the z-grid-range */
        break;
      }
    }
    if (n_up < -1) {
      fprintf (stderr, "Error in function combine_profiles!\n");  
      fprintf (stderr, "upper index not found, n_up = %d\n", n_up);   
      return -2;
    }
    
    n_dn = nlev_atm-1;   
    for  (lc=(nlev_atm-1); lc>=0; lc--)
      if(z_atm[lc] >= z[nlev-1]) {
        n_dn = lc;       /* last index of z_atm within the z-grid-range */
	break;
      }

    if (n_dn < 0) {
      fprintf (stderr, "Error in function combine_profiles!\n");  
      fprintf (stderr, "lower index not found, n_down = %d\n", n_dn);   
      return -3;
    }

    (*nlev_comb) = (n_up + 1) + nlev + ((nlev_atm - 1) - n_dn);
    
    /* combine the two profiles, dens (on z) has priority */
    /* Allocation */
    *z_comb    = calloc (*nlev_comb, sizeof (float  ));
    *dens_comb = calloc (nsza,       sizeof (float *));
    for (i=0; i<nsza; i++)
      (*dens_comb)[i] = calloc (*nlev_comb, sizeof (float));

    /* combine z-grid */
    for (lc=0; lc<=n_up; lc++)
      (*z_comb)[lc]          =    z_atm[lc];
    for (lc=0; lc<nlev; lc++)
      (*z_comb)[lc+(n_up+1)] =    z[lc];
    for (lc=1; lc<(nlev_atm-n_dn);lc++)
      (*z_comb)[lc+(n_up+nlev)] =    z_atm[lc+n_dn];

    /* combine data */
    for (i=0; i<nsza; i++) {
      for (lc=0; lc<=n_up; lc++)
        (*dens_comb)[i][lc] = dens_atm[lc];
      for (lc=0; lc<nlev; lc++)
        if (dens[i][lc] != -1.0)
          (*dens_comb)[i][lc+(n_up+1)] = dens[i][lc];
	else
          (*dens_comb)[i][lc+(n_up+1)] = dens_atm_filegrid[lc];
      for (lc=1; lc<(nlev_atm-n_dn);lc++)
        (*dens_comb)[i][lc+(n_up+nlev)] = dens_atm[lc+n_dn];
    }
  }
  
  /*   for (lc=0; lc<*nlev_comb;lc++) */
  /*     fprintf (stderr,"z_comb[%3d] = %7.3f, dens_comb[%3d] = %11.7e\n", lc, (*z_comb)[lc], lc, (*dens_comb)[0][lc]); */

  free(dens_atm_filegrid);
  free(dens_air_filegrid);

  return 0;
}



/**************************************************************/
/* Function: calculate_or_read_refind                         */
/*                                                            */
/* Read refraction index file or calculate refractive index.  */
/* Combine with background atmosphere and                     */
/* interpolate on the common grid: z_all                      */
/*                                                            */
/* The parameterization of the refractive index is from       */
/* R. Penndorf, J. Opt. Soc. Am., Vol. 46, No.2, 1957         */
/*                                                            */
/* Parameters: refind   (refractive index)-1                  */
/* Retun value:                                               */
/* Example:                                                   */
/* Files:                                                     */
/* Known bugs:                                                */
/* Author: Bernhard Mayer                                     */
/*         Claudia Emde (included spectral dependence)        */
/**************************************************************/

static int calculate_or_read_refind  (char *filename,
				      float ***refind, float *z_all, int nlev_all,
				      int interpol_method_refind, 
                                      float *t, float *p, int nlambda, float *lambda,
                                      int quiet)

{
  int status=0, nlev=0, lc=0, iv=0;

  float *z=NULL;
  float *refind_all=NULL;
  float **refind_comb=NULL;
  float *refind_file=NULL; 
  float *z_comb=NULL;
  int nlev_comb=-666;
  float t_0 = 15.0;     /* normal temperature in celsius      */
  float p_0 = 1013.25;  /* normal pressure in mbar            */
  float ar  = 0.00366;  /* parameter for calculating refind   */
  float *dens_air_dummy = NULL;
  float nu = 0.0;       /* frequency                          */ 
  float n_0;            /* refractive index for standard pressure,
                           temperature */
  
  /* read refractive index from file, assume that this refractive
     index is for all wavelengths */
  if (strlen(filename)>0) {
    
    fprintf (stderr, " ... reading refractive index from %s\n", filename);

    status = read_2c_file_float (filename, &z, &refind_file, &nlev);   
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }

    if (z[0] == z_all[0] && z[nlev-1] == z_all[nlev_all-1]) {
      /* interpolate refind profile to common atmospheric grid */
      status = interpolate_profile (z, &refind_file, nlev,
                                    z_all, nlev_all,
                                    interpol_method_refind, quiet);
      if (status!=0) {
        fprintf (stderr, "Error %d interpolating refractive index to common atmospheric grid\n", status);
        return status;
      }
      /* use the same refractive index for all wavelengths*/
      for (iv=0; iv<nlambda; iv++)
        (*refind)[iv]=refind_file; 
      
    }
    else {
      /* combine refind_file and atmosphere-file data*/
      
      /* calculate refind of the atmosphere-file data */
      refind_all = (float *) calloc (nlev_all, sizeof(float));      
      if (refind_all == NULL) {
        fprintf (stderr, "Error allocating memory for refind_all\n");
        return -1;
      }
      /* Calculate n0 according to formula by Edlen, 1953*/
      
      for (iv=0; iv<nlambda; iv++){
        
        nu=1.0/(lambda[iv]*1e-3); /* conversion from nm to um */ 
        n_0=(6432.8+2949810/(146.0-nu*nu)+25540.0/(41.0-nu*nu))*1.0e-8;
        
        for (lc = 0; lc < nlev_all; lc++)
          refind_all[lc] = n_0 * (1.0 + ar*t_0) /
            (1.0 + ar * (t[lc]-273.15)) * p[lc] / p_0; 
        
        /* air dens dummy */
        dens_air_dummy = (float *) calloc (nlev_all, sizeof(float)); 
        
        /* combine refind-file data and atmosphere-file refind-data  */
        status = combine_profiles (refind[iv], z, nlev, 1,
                                   refind_all, dens_air_dummy, z_all, nlev_all, /* dens_air_dummy is dummy argument*/
                                   interpol_method_refind, 1,                   /* 1              is dummy argument*/
                                   &(refind_comb), &(z_comb), &(nlev_comb), quiet);
        
        if (status!=0) {
          fprintf (stderr, "Error %d returned by combine_profiles(refind)\n", status);
          return status;
        }
        
        free(dens_air_dummy);
        free(*refind[iv]);
        
        (*refind)[iv] = (float *) calloc (nlev_all, sizeof(float));
        
        status = arb_wvn (nlev_comb, z_comb, refind_comb[0],
                          nlev_all,  z_all, (*refind)[iv],
                          interpol_method_refind, 1);
        
        if (status!=0) {
          fprintf (stderr, "Error %d returned by arb_wvn()\n", status);
          return status;
        }
        
        free(z_comb);
        free(refind_comb[0]);
        free(refind_comb);
        free(refind_all);
      }
      free(z);
    }
  }
  else {  /* if (strlen(filename)>0) => no refind-input-file */
    
    /* allocate memory */
    ASCII_calloc_float(&(*refind), nlambda, nlev_all);
    if (*refind == NULL) {
      fprintf (stderr, "Error allocating memory for refind\n");
      return -1;
    }
    /* Calculate n0 according to formula by Edlen, 1953*/
    
    for (iv=0; iv<nlambda; iv++){
      
      nu=1.0/(lambda[iv]*1e-3); /* conversion from nm to um */ ; 
      
      /* The formula by Edlen (1953) was fitted to measurements starting at 185nm. */
      /* This formula has a singularity around 156nm. */
      /* To overcome this limitation the refractive index at 185nm is used also for lower wavelengths. */
      if (nu > 1.0/0.185) 
        nu=1.0/0.185;         

      n_0=(6432.8+2949810/(146.0-nu*nu)+25540.0/(41.0-nu*nu))*1.0e-8;

      for (lc = 0; lc < nlev_all; lc++)
        (*refind)[iv][lc] = 
          n_0 * (1.0 + ar*t_0) / (1.0 + ar * (t[lc]-273.15)) * p[lc] / p_0;
    }
  }
  return 0;   /* o.k. */
}



/***********************************************************************************/
/*  Read the molecular absorption file which can either be a two-column file       */
/*  containing altitude [km] and optical thickness of each layer or a matrix       */
/*  tau[wavelength, z]. Memory for dtau[] is allocated automatically.              */ 
/***********************************************************************************/

static int read_molecular_absorption (char* filename, 
				      float **zd, int *nlyr,
				      int quiet,
				      double **dtau_mono, 
				      float ***dtau_spec, float **lambda, int *nlambda,
				      int *monochromatic, int wl_start_index, int wl_end_index) 
{
  int rows=0, min_columns=0, max_columns=0, status=0;
  int lu=0, iv=0;

  int fulfilled=0;

  float **data=NULL;


#if HAVE_LIBNETCDF
  int file_with_nstk=0;
  /* try to open the CDF file */

  char cdf_filename[FILENAME_MAX]="";
  int ncid=0, iw=0;

  int idd_nlyr=0, idd_nwvl=0, idd_nstk=0;
  int id_wvl=0, id_z=0, id_tau=0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  size_t taustart[2] = {0,0};
  size_t taucount[2] = {0,0};

  size_t taustart_stk[4] = {0,0,0,0};
  size_t taucount_stk[4] = {0,0,0,0};

  double *tmp_zd=NULL, *tmp_lambda=NULL;

  size_t dimlen=0;

  strcpy (cdf_filename, filename);
  strcat (cdf_filename, ".cdf");

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) {
    if (!quiet) {
      fprintf (stderr, "*** %s is not a netcdf file, trying to open as ASCII instead.\n", filename);
      fprintf (stderr, "*** It is recommended to use the netCDF format because access is much faster.\n");
    }
  }
  else {
    if (!quiet)
      fprintf (stderr, " ... reading profiles from netcdf file %s\n", filename);

    /* get dimension id for "nlyr" */
    status = nc_inq_dimid (ncid, "nlyr", &idd_nlyr);
    
    /* get dimension length for "nlyr" */
    status = nc_inq_dimlen (ncid, idd_nlyr, &dimlen);
    *nlyr = dimlen;

    /* get dimension id for "nwvl" */
    status = nc_inq_dimid (ncid, "nwvl", &idd_nwvl);
    
    /* get dimension length for "nwvl" */
    status = nc_inq_dimlen (ncid, idd_nwvl, &dimlen);
    *nlambda = dimlen;

    /* read levels */
    tmp_zd = calloc (*nlyr+1, sizeof(double));
    *zd    = calloc (*nlyr+1, sizeof(float));

    /* get variable id for "z" */
    status = nc_inq_varid (ncid, "z", &id_z);

    count[0] = *nlyr+1;

    /* check if there is a dimension 'nstk' (it was introduced in ARTS-2.1) */
    status = nc_inq_dimid (ncid, "nstk", &idd_nstk);
    if (status==NC_NOERR){
      file_with_nstk=1;
      status = nc_inq_dimlen (ncid, idd_nstk, &dimlen);
      if (dimlen!=1) {
        fprintf (stderr, "Error: Only nstk=1 is implemented in uvspec, but %s uses different a different nstk.\n", cdf_filename);
        return -1;
      }
    }
    else
      file_with_nstk=0;

    /* read "z" */
    status = nc_get_vara_double (ncid, id_z, start, count, tmp_zd);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, cdf_filename);
      return status;
    }
    
    /* read wavelengths */
    tmp_lambda = calloc (*nlambda, sizeof(double));
    *lambda    = calloc (*nlambda, sizeof(float));
    /* aky20091018 *lambda    = calloc ( wl_end_index-wl_start_index+1, sizeof(float)); */

    /* get variable id for "wvl" */
    status = nc_inq_varid (ncid, "wvl", &id_wvl);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, cdf_filename);
      return status;
    }
    
    /* read "wvl" */
    count[0] = *nlambda;

    status = nc_get_vara_double (ncid, id_wvl, start, count, tmp_lambda);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, cdf_filename);
      return status;
    }


    /* get variable id for "tau" */
    status = nc_inq_varid (ncid, "tau", &id_tau);
 
    /* read optical thickness */
    if (file_with_nstk){
      taustart_stk[1] = 0;  /* start with first optical thickness */
      taucount_stk[0] = 1;  /* read one layer at a time           */
      taucount_stk[2] = 1; 
      taucount_stk[3] = 1; 
    }
    else{
      taustart[1] = 0;  /* start with first optical thickness */
      taucount[0] = 1;  /* read one layer at a time           */
    }

    *dtau_spec = calloc (*nlyr, sizeof(float *));
    for (lu=0; lu<*nlyr; lu++) {
      (*dtau_spec)[lu] = calloc (*nlambda, sizeof(float));

      if (file_with_nstk){
        taustart_stk[0] = lu;
        taucount_stk[1] = *nlambda;
        status = nc_get_vara_float (ncid, id_tau, taustart_stk, taucount_stk, (*dtau_spec)[lu]);
      }
      else{
        taustart[0] = lu;
        taucount[1] = *nlambda;
        status = nc_get_vara_float (ncid, id_tau, taustart, taucount, (*dtau_spec)[lu]);
      }

      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading %s\n", status, cdf_filename);
	return status;
      }
    }

    /* copy data to final destination and free temporary memory */
    for (lu=0; lu<=*nlyr; lu++)
      (*zd)[lu] = tmp_zd[lu];

    iw=0;
    /* aky20091018    for (iv=wl_start_index; iv<=wl_end_index; iv++) { */
    for (iv=0; iv<*nlambda; iv++) 
      (*lambda)[iw++] = tmp_lambda[iv];

    free (tmp_zd); free (tmp_lambda);

    /* this is obviously a wavelength-dependent cross section */
    *monochromatic=0;

    return 0;
  }

#endif

  if (!quiet)
    fprintf (stderr, " ... reading profiles from ASCII file %s\n", filename);

  /* read altitude and optical depth from file */
  status = ASCII_file2float (filename,   
			     &rows,        
			     &max_columns, 
			     &min_columns, 
			     &data);

  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }


  /* a rough format check */
  if (min_columns<2) {
    fprintf (stderr, "Error, found less than 2 columns in %s\n", filename);
    return -1;
  }

  if (min_columns!=max_columns) {
    fprintf (stderr, "Error, inconsistent number of columns in %s, %d - %d\n",
	     filename, min_columns, max_columns);
    return -1;
  }


  /* check if first column sorted in descending order;    */
  /* if yes, we assume old format: two columns z and dtau */
  fulfilled=1;
  for (lu=0; lu<rows-1; lu++)
    if (data[lu][0] < data[lu+1][0])
      fulfilled = 0;
  
  if (fulfilled) {

    /* here we assume that the file comes in the old format,  */
    /* two columns z and dtau; in this case, the first column */
    /* must equal the atmospheric profile                     */

    *monochromatic=1;

    *nlyr = rows-1;
    *zd        = calloc (*nlyr+1, sizeof (float));
    *dtau_mono = calloc (*nlyr,   sizeof (double));

    for (lu=0; lu<*nlyr+1; lu++)
      (*zd)[lu] = data[lu][0];

    for (lu=0; lu<*nlyr; lu++)
      (*dtau_mono)[lu] = data[lu+1][1];
  }
  else {
    /* now we need to assume that we deal with a new format */
    /* where the first row must be the altitude levels      */
    
    /* check if first column sorted in ascending order; */
    /* if yes, we assume wavelength-dependent format    */
    fulfilled=1;
    for (iv=1; iv<rows-1; iv++)
      if (data[iv][0] > data[iv+1][0])
	fulfilled = 0;

    if (!fulfilled) {
      fprintf (stderr, "Error reading %s\n", filename);
      return -1;
    }
    
    /* now we know that the format is correct, and the only thing   */
    /* we need to do is to copy the data to their final destination */
    
    *monochromatic=0;

    *nlyr    = min_columns-1;
    *nlambda = rows-1;

    /* allocate memory */
    *zd = calloc (*nlyr+1, sizeof (float));
    *lambda    = calloc (*nlambda, sizeof (float));
    *dtau_spec = calloc (*nlyr,    sizeof (float *));
    for (lu=0; lu<*nlyr; lu++)
      (*dtau_spec)[lu] = calloc (*nlambda, sizeof(float));


    for (lu=0; lu<=*nlyr; lu++)
      (*zd)[lu] = data[0][lu];

    for (iv=0; iv<*nlambda; iv++) {
      (*lambda)[iv] = data[iv+1][0];
      
      for (lu=0; lu<*nlyr; lu++)
	(*dtau_spec)[lu][iv] = data[iv+1][lu+1];
      
    }
  }
  
  /* free memory */
  ASCII_free_float (data, rows);

  return 0;
}


/***********************************************************************************/
/*  Read the molecular absorption file which can either be a two-column file       */
/*  containing altitude [km] and optical thickness of each layer or a matrix       */
/*  tau[wavelength, z].                                                            */
/***********************************************************************************/

static int read_molecular_absorption_z (char* filename, int quiet,
					float **zd, int *nlyr) 
{
  int rows=0, min_columns=0, max_columns=0, status=0;
  int lu=0, iv=0;

  int fulfilled=0;

  
  float **data=NULL;
  

#if HAVE_LIBNETCDF
  /* try to open the CDF file */

  char cdf_filename[FILENAME_MAX]="";
  int ncid=0;

  int idd_nlyr=0;
  int id_z=0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  size_t dimlen=0;

  double *z=NULL;

  strcpy (cdf_filename, filename);
  strcat (cdf_filename, ".cdf");

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {
    /* get dimension id for "nlyr" */
    status = nc_inq_dimid (ncid, "nlyr", &idd_nlyr);
    
    /* get dimension length for "nlyr" */
    status = nc_inq_dimlen (ncid, idd_nlyr, &dimlen);
    *nlyr = dimlen;

    /* read levels */
    z   = calloc (*nlyr+1, sizeof(double));
    *zd = calloc (*nlyr+1, sizeof(float));

    /* get variable id for "z" */
    status = nc_inq_varid (ncid, "z", &id_z);

    count[0] = *nlyr+1;

    /* read "z" */
    status = nc_get_vara_double (ncid, id_z, start, count, z);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, cdf_filename);
      return status;
    }

    /* copy profile to final destination and free temporary memory */
    for (lu=0; lu<=*nlyr; lu++)
      (*zd)[lu] = z[lu];
    free (z);

    nc_close (ncid);
    
    return 0;
  }

#endif

  if (!quiet)
    fprintf (stderr, " ... reading altitudes from ASCII file %s\n", filename);

  /* read altitude and optical depth from file */
  status = ASCII_file2float (filename,   
			     &rows,        
			     &max_columns, 
			     &min_columns, 
			     &data);

  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }


  /* a rough format check */
  if (min_columns<2) {
    fprintf (stderr, "Error, found less than 2 columns in %s\n", filename);
    return -1;
  }

  if (min_columns!=max_columns) {
    fprintf (stderr, "Error, inconsistent number of columns in %s, %d - %d\n",
	     filename, min_columns, max_columns);
    return -1;
  }


  /* check if first column sorted in descending order;    */
  /* if yes, we assume old format: two columns z and dtau */
  fulfilled=1;
  for (lu=0; lu<rows-1; lu++)
    if (data[lu][0] < data[lu+1][0])
      fulfilled = 0;
  
  if (fulfilled) {
    *nlyr = rows-1;
    *zd = calloc (rows, sizeof (float));

    for (lu=0; lu<rows; lu++)
      (*zd)[lu] = data[lu][0];
  }
  else {
    /* now we need to assume that we deal with a new format */
    /* where the first row must be the altitude levels      */
    
    /* check if first column sorted in ascending order; */
    /* if yes, we assume wavelength-dependent format    */
    fulfilled=1;
    for (iv=1; iv<rows-1; iv++)
      if (data[iv][0] > data[iv+1][0])
	fulfilled = 0;

    if (!fulfilled) {
      fprintf (stderr, "Error reading %s\n", filename);
      return -1;
    }
    
    *nlyr = min_columns-1;
    *zd = calloc (min_columns, sizeof (float));
    for (lu=0; lu<min_columns; lu++)
      (*zd)[lu] = data[0][lu];
  }

  /* free memory */
  ASCII_free_float (data, rows);

  return 0;
}




/***********************************************************************************/
/*  Read the molecular absorption file which can either be a two-column file       */
/*  containing altitude [km] and optical thickness of each layer or a matrix       */
/*  tau[wavelength, z]. Check if the levels equal those defined in zd[]. Memory    */
/*  for dtau[] is allocated automatically.                                         */ 
/***********************************************************************************/

int read_molecular_absorption_lambda (char* filename, int quiet,
				      float **lambda, int *nlambda, int *monochromatic) 
{
  int rows=0, min_columns=0, max_columns=0, status=0;
  int iv=0;

  int fulfilled=0;

  float **data=NULL;
  
#if HAVE_LIBNETCDF
  /* try to open the CDF file */

  char cdf_filename[FILENAME_MAX]="";
  int ncid=0;

  int idd_nwvl=0;
  int id_wvl=0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  size_t dimlen=0;

  double *wvl=NULL;

  strcpy (cdf_filename, filename);
  strcat (cdf_filename, ".cdf");

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {
    /* get dimension id for "nwvl" */
    status = nc_inq_dimid (ncid, "nwvl", &idd_nwvl);
    
    /* get dimension length for "nwvl" */
    status = nc_inq_dimlen (ncid, idd_nwvl, &dimlen);
    *nlambda = dimlen;

    /* read levels */
    wvl     = calloc (*nlambda, sizeof(double));
    *lambda = calloc (*nlambda, sizeof(float));

    /* get variable id for "wvl" */
    status = nc_inq_varid (ncid, "wvl", &id_wvl);

    count[0] = *nlambda;

    /* read "wvl" */
    status = nc_get_vara_double (ncid, id_wvl, start, count, wvl);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, cdf_filename);
      return status;
    }

    /* copy profile to final destination and free temporary memory */
    for (iv=0; iv<*nlambda; iv++)
      (*lambda)[iv] = wvl[iv];
    free (wvl);

    *monochromatic=0;

    nc_close (ncid);

    return 0;
  }

#endif

  if (!quiet)
    fprintf (stderr, " ... reading wavelengths from ASCII file %s\n", filename);

  /* read altitude and optical depth from file */
  status = ASCII_file2float (filename,   
			     &rows,        
			     &max_columns, 
			     &min_columns, 
			     &data);

  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }


  /* a rough format check */
  if (min_columns<2) {
    fprintf (stderr, "Error, found less than 2 columns in %s\n", filename);
    return -1;
  }

  if (min_columns!=max_columns) {
    fprintf (stderr, "Error, inconsistent number of columns in %s, %d - %d\n",
	     filename, min_columns, max_columns);
    return -1;
  }


  /* check if first column sorted in descending order;    */
  /* if yes, we assume old format: two columns z and dtau */
  fulfilled=1;
  for (iv=0; iv<rows-1; iv++)
    if (data[iv][0] < data[iv+1][0])
      fulfilled = 0;
  
  if (fulfilled){
    *monochromatic=1;
    *nlambda=1;
    *lambda = calloc (1, sizeof (float));
  }
  else {
    /* now we need to assume that we deal with a new format */
    /* where the first column must be the wavelengths       */
    
    *monochromatic=0;
    
    if(!quiet){
      fprintf (stderr, " ... now testing as wavelength-dependent absorption file\n");
      fflush (stderr);
    }
    
    /* check if first column sorted in ascending order; */
    /* if yes, we assume wavelength-dependent format    */
    fulfilled=1;
    for (iv=1; iv<rows-1; iv++)
      if (data[iv][0] > data[iv+1][0])
	fulfilled = 0;

    if (!fulfilled) {
      fprintf (stderr, "Error reading %s\n", filename);
      return -1;
    }
    
    *nlambda = rows-1;
    *lambda = calloc (rows-1, sizeof (float));
    for (iv=0; iv<*nlambda; iv++)
      (*lambda)[iv] = data[iv+1][0];
  }

  /* free memory */
  ASCII_free_float (data, rows);

  return 0;
}



/**************************************************************/
/* Calculate atmospheric column in DU.                        */
/* OUTPUT: dens_avg, column */
/**************************************************************/

int column (float *dens, float *zd, int nlev, float altitude, int interpol_method, float *dens_air, 
            float **dens_avg, int allocate, float *column)
{

  float deltaz=0;
  int lc=0, ialt=0;
  int status=0;
  
  
  for (ialt=0; ialt<nlev; ialt++)
    /*    if (altitude == zd[ialt])  which one is correct ??? */
    if (altitude >= zd[ialt])
      break;
    
  if (ialt == nlev) {
    fprintf (stderr, "Error, user-defined altitude %f not found in altitude grid;\n", altitude);
    fprintf (stderr, "exiting column (in atmosphere.c)!\n");
    return -1;
  }

  status = average_dens(dens, dens_air, zd, nlev, interpol_method, dens_avg, allocate); /* in ancillary.c */
  if (status!=0) {
    fprintf (stderr, "Error %d calculating average concentration (column, in atmosphere.c)\n", status);
    return status;
  }

  *column=0;
  for (lc=0; lc<ialt; ++lc) {
    deltaz = (zd[lc] - zd[lc+1]) * 1.0E5;  /* km to cm */
    *column += (*dens_avg)[lc] * deltaz;   /* cm-3 * cm = cm-2 */
  }

  *column /= DU2particle;   /* column [particles per cm^2] to DU */
  
  return 0;
}



/**************************************************************/
/* Convert relative humidity [%] to concentration [1/cm3]     */
/**************************************************************/

double vapor_pressure (double t)
{
  /*
    double c1=-6096.9385;
    double c2=16.635794;
    double c3=-2.711193e-2;
    double c4=1.673952e-5;
    double c5=2.433502;
  */
  
  double tzero = 273.15;
  double a = tzero/t;
  double mh2o = M_H2O/AVOGADRO;


  /* Sonntag, Meteorol. Zeitschr, n.f. 3, 51-66, 1994                                 */
  /* return 100.0 / BOLTZMANN / t * exp (c1/t +c2 +c3*t +c4*t*t + c5*log(t)) / 1.0e6; */
  
  /* SBDART, -> CRC handbook */
  return a * exp (18.916758 - a * (14.845878 + a*2.4918766)) / mh2o / 1.0e6;
}

/***********************************************************************************/
/* Function: vapor_pressure_over_ice                                               */
/* Description:                                                                    */
/*  This routine reuturns the water vapour pressure over ice                       */
/*  according to the Goff Gratch equation,                                         */
/*  (Smithsonian Met. Tables,  5th ed., pp. 350, 1984):                            */
/*                                                                                 */
/* Parameters:                                                                     */
/*     double t    Temperature in Kelvin                                           */
/* Return value:                                                                   */
/*     n           number density of water vapour in cm**-3, when                  */
/*                 air is saturated with respect to (water) ice                    */
/*                 NaN, if T > 0 degree C                                          */
/*                                                                                 */
/* Example:  p = vapor_pressure_over_ice (255)                                     */
/* Files:    atmosphere.c                                                          */
/* Known bugs: -                                                                   */
/* Authors:                                                                        */
/*    taken from cdf2c (http://cires.colorado.edu/~voemel/vp.html                  */
/*    Mar 2006   U. Hamann     Implemented into libRadtran                         */
/*                                                                                 */
/***********************************************************************************/

double vapor_pressure_over_ice (double t)
{

  double c1=-9.09718;
  double c2=-3.56654;
  double c3= 0.876793;
  double c4= 6.1071;
  
  double tzero = 273.16;
  double a = tzero/t;

  if (t > 273.16)
    return 0.0/0.0;
  else 
    /* n = p/ kT, p[hPa] = pow(10.0, c1*(a-1) + c2*log10(a) + c3*(1-1/a) + log10(c4)) */
    return 100.0 * pow(10.0, c1*(a-1) + c2*log10(a) + c3*(1-1/a) + log10(c4)) / (BOLTZMANN * t)  / 1.0e6;
  /*    (100.0 == hPa -> Pa)                                                                     1.0e6   ==   1m^3 -> 1/cm^3  */
   
}



/* converts number of interpolation method to human readable string*/
static char* interpol_number2string(int method_number) {

  char *method = (char *) calloc (strlen("spline             ")+1, sizeof (char));

  switch (method_number) {
  case INTERP_METHOD_SPLINE:
    strcpy (method, "spline             ");
    break;
  case INTERP_METHOD_LINEAR:
    strcpy (method, "linear             ");
    break;
  case INTERP_METHOD_LOG:
    strcpy (method, "logarithmic        ");
    break;
  case INTERP_METHOD_LINMIX:
    strcpy (method, "linear mixing ratio");
    break;
  case INTERP_METHOD_LOG_SPLINE:
    strcpy (method, "logarithmic spline ");
    break;
  default:
    fprintf (stderr, "Error. Unknown interpolation method z_interpolate = %d!\n", method_number);
    fprintf (stderr, "This is a program bug, please contact the programmers!\n");
    fprintf (stderr, "Please include the used input-file!\n");
    return "       -1          ";
    break;
  }

  return method;
}

/* converts number of MOL_** to human readable string*/
char* gas_number2string(int gas_number) {

  char *gas=NULL;
  /* all gas names should have 4 Characters */
  gas = (char *) calloc (strlen("Air ")+1, sizeof (char));

  switch (gas_number) {
  case MOL_AIR:
    strcpy (gas, "Air ");
    break;
  case MOL_O3:
    strcpy (gas, "O3  ");
    break;
  case MOL_O2:
    strcpy (gas, "O2  ");
    break;
  case MOL_H2O:
    strcpy (gas, "H2O ");
    break;
  case MOL_CO2:
    strcpy (gas, "CO2 ");
    break;
  case MOL_NO2:
    strcpy (gas, "NO2 ");
    break;
  case MOL_BRO:
    strcpy (gas, "BRO ");
    break;
  case MOL_OCLO:
    strcpy (gas, "OCLO");
    break;
  case MOL_HCHO:
    strcpy (gas, "HCHO");
    break;
  case MOL_O4:
    strcpy (gas, "O4  ");
    break;
  case MOL_SO2:
    strcpy (gas, "SO2 ");
    break;
  case MOL_CH4:
    strcpy (gas, "CH4 ");
    break;
  case MOL_N2O:
    strcpy (gas, "N2O ");
    break;
  case MOL_CO:
    strcpy (gas, "CO  ");
    break;
  case MOL_N2:
    strcpy (gas, "N2  ");
    break;
  default:
    fprintf (stderr, "Error. Unknown indentifier for gas species MOL_** = %d !\n", gas_number);
    fprintf (stderr, "Bug detected in gas_number2string (atmosphere.c).\n");
    fprintf (stderr, "Please contact us (the programmers) and include the used input-file!\n");
    return " -1 ";
    break;
  }

  return gas;
}



/****************************************************************/
/* read_dtheta_dx_from_ECMWF_file                               */
/*                                                              */
/* Purpose:                                                     */
/* read the horizontal gradient of the potential temperature    */
/* (theta) from an ECMWF netCDF file and interpolate this to    */
/* the zout levels                                              */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  input structure                                             */
/*  output structure                                            */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  output->atm.microphys.dtheta_dx                             */
/*            zonal gradient of theta                           */
/*                                                              */
/*  October 2007  by Ulrich Hamann                              */
/****************************************************************/

static int read_dtheta_dx_from_ECMWF_file (input_struct input, output_struct *output)
{

  int status=0;

#if HAVE_LIBNETCDF

  float   kappa = (C_P_DRY_STD-C_V_DRY_STD) / C_P_DRY_STD;  /* == R / cp */

  int ncid=NOT_DEFINED_INTEGER;

  int nt=-1;
  long ilat=NOT_DEFINED_INTEGER, ilon=NOT_DEFINED_INTEGER;  /* index for lat and lon in netCDF file */

  size_t nlat  = 0;
  size_t nlon  = 0;
  size_t nlev  = 0;

  size_t ilon_minus_one = 0;
  size_t ilon_plus_one  = 0;

  double *ECMWF_lat  = NULL;
  double *ECMWF_lon  = NULL;

  int itime1 = -1, itime2 = -1;
  float dt = NOT_DEFINED_FLOAT;

  size_t nlay=0;

  float *p_layer = NULL;
  float *T_layer = NULL;

  float *theta_minus_one = NULL;
  float *theta_plus_one  = NULL;

  float *z_layer=NULL;

  int lc=NOT_DEFINED_INTEGER;
  float dx=NOT_DEFINED_FLOAT;

  float *dtheta_dx=NULL;

  if (input.verbose)
    fprintf (stderr, " ... read 'dT/dx' from netCDF file %s\n", input.filename[FN_ECMWF]);

  status = get_all_netCDF_indices ( input.filename[FN_ECMWF], input.latitude, input.longitude,
                                    &(ncid), &(ilat), &(nlat), &(ilon), &(nlon), 
                                    &(ECMWF_lat), &(ECMWF_lon),
                                    input.UTC, input.atm.time_interpolate,
                                    &(nt), &(itime1), &(itime2), &(dt),
                                    input.verbose, input.quiet);

  if (status!=0) {
    fprintf (stderr, "Error %d, while get_all_netCDF_indices from %s\n", status, input.filename[FN_ECMWF]);
    return status;
  }

  if (ilon != nlon -1 ) 
    ilon_plus_one = ilon + 1; /* centered difference */
  else
    ilon_plus_one = nlon -1; /* backward difference at eastern boundary */

  if (ilon != 0 ) 
    ilon_minus_one = ilon - 1; /* centered difference */
  else
    ilon_minus_one = 0;     /* forward  difference at western  boundary */

  /* get temperature profile at ilon_minus_one */
  status = read_p_T_z_from_ECMWF_file( ncid, ilat, nlat, ilon_minus_one, nlon, nt, itime1, itime2, dt,
                                       &(nlay), &(nlev), output->alt.altitude,
                                       &(p_layer), &(theta_minus_one), &(z_layer),
                                       input.verbose, input.quiet) ;
  
  if (status!=0) {
    fprintf (stderr, "Error %d, while read_p_T_z_from_ECMWF_file (1) from %s\n", status, input.filename[FN_ECMWF]);
    return status;
  }

  for (lc=0; lc<nlay; lc++)
    theta_minus_one[lc] = theta_minus_one[lc] * pow((1000./p_layer[lc]),kappa);

  free(p_layer);
  free(z_layer);

  /* get temperature profile at ilon_plus_one */
  status = read_p_T_z_from_ECMWF_file( ncid, ilat, nlat, ilon_plus_one, nlon, nt, itime1, itime2, dt,
                                       &(nlay), &(nlev), output->alt.altitude,
                                       &(p_layer), &(theta_plus_one), &(z_layer),
                                       FALSE, input.quiet) ;

  if (status!=0) {
    fprintf (stderr, "Error %d, while read_p_T_z_from_ECMWF_file (2) from %s\n", status, input.filename[FN_ECMWF]);
    return status;
  }

  for (lc=0; lc<nlay; lc++)
    theta_plus_one[lc] = theta_plus_one[lc] * pow((1000./p_layer[lc]),kappa);

  free(p_layer);
  free(z_layer);

  /* get z profile at ilon */
  status = read_p_T_z_from_ECMWF_file( ncid, ilat, nlat, ilon, nlon, nt, itime1, itime2, dt,
                                       &(nlay), &(nlev), output->alt.altitude,
                                       &(p_layer), &(T_layer), &(z_layer),
                                       FALSE, input.quiet) ;

  if (status!=0) {
    fprintf (stderr, "Error %d, while read_p_T_z_from_ECMWF_file (3) from %s\n", status, input.filename[FN_ECMWF]);
    return status;
  }

  free(p_layer);
  free(T_layer);

  dx = (ECMWF_lon[ilon_plus_one]- ECMWF_lon[ilon_minus_one])/360 * (2*PI*R_EARTH) * cos (2.0*PI/360.0*ECMWF_lat[ilat]);

  if (input.verbose)
    fprintf (stderr, "     ilon-1=%zd (%6.2f), ilon+1=%zd (%6.2f), dx = %6.3f m\n", 
	     ilon_minus_one, ECMWF_lon[ilon_minus_one], ilon_plus_one, ECMWF_lon[ilon_plus_one], dx );

  if ((dtheta_dx = calloc (nlay, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for dtheta_dx\n");
    return -1;
  }

  for (lc=0; lc<nlay; lc++) {
    dtheta_dx[lc] = (theta_plus_one[lc] - theta_minus_one[lc]) / dx;
    /* fprintf (stderr, " %6.3f %8.3f %8.3f %13.6e\n", z_layer[lc], theta_minus_one[lc], theta_plus_one[lc], dtheta_dx[lc] ); */
  }

  if ((output->atm.microphys.dtheta_dx = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for dtheta_dx\n");
    return -2;
  }

  /* interpolation of  dtheta_dx to zout levels */
  status = arb_wvn (nlay, z_layer, dtheta_dx, 
                    output->atm.nzout, output->atm.zout_sea, output->atm.microphys.dtheta_dx, INTERP_METHOD_LINEAR, 1);

  if (status!=0) {
    fprintf (stderr, "Error %d interpolating d_theta_dx\n", status);
    return status;
  }


#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use the ECMWF data file option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ******************************************************************\n");
  return -1;
#endif

  return status;

}




/****************************************************************/
/* read_dtheta_dy_from_ECMWF_file                               */
/*                                                              */
/* Purpose:                                                     */
/* read the horizontal gradient of the potential temperature    */
/* (theta) from an ECMWF netCDF file and interpolate this to    */
/* the zout levels                                              */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  input structure                                             */
/*  output structure                                            */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  output->atm.microphys.dtheta_dy                             */
/*            meridional gradient of theta                      */
/*                                                              */
/*  October 2007  by Ulrich Hamann                              */
/****************************************************************/

static int read_dtheta_dy_from_ECMWF_file (input_struct input, output_struct *output)
{

  int status=0;

#if HAVE_LIBNETCDF

  float   kappa = (C_P_DRY_STD-C_V_DRY_STD) / C_P_DRY_STD;  /* == R / cp */

  int ncid=NOT_DEFINED_INTEGER;

  int nt=-1;
  long ilat=NOT_DEFINED_INTEGER, ilon=NOT_DEFINED_INTEGER;  /* index for lat and lon in netCDF file */

  size_t nlat  = 0;
  size_t nlon  = 0;
  size_t nlev  = 0;

  size_t ilat_minus_one = 0;
  size_t ilat_plus_one  = 0;

  double *ECMWF_lat  = NULL;
  double *ECMWF_lon  = NULL;

  int itime1 = -1, itime2 = -1;
  float dt = NOT_DEFINED_FLOAT;

  size_t nlay=0;

  float *p_layer = NULL;
  float *T_layer = NULL;

  float *theta_minus_one = NULL;
  float *theta_plus_one  = NULL;

  float *z_layer=NULL;

  int lc=NOT_DEFINED_INTEGER;
  float dy=NOT_DEFINED_FLOAT;

  float *dtheta_dy=NULL;

  if (input.verbose)
    fprintf (stderr, " ... read 'dT/dy' from netCDF file %s\n", input.filename[FN_ECMWF]);

  status = get_all_netCDF_indices ( input.filename[FN_ECMWF], input.latitude, input.longitude,
                                    &(ncid), &(ilat), &(nlat), &(ilon), &(nlon), 
                                    &(ECMWF_lat), &(ECMWF_lon),
                                    input.UTC, input.atm.time_interpolate,
                                    &(nt), &(itime1), &(itime2), &(dt),
                                    input.verbose, input.quiet);

  if (status!=0) {
    fprintf (stderr, "Error %d, while get_all_netCDF_indices from %s\n", status, input.filename[FN_ECMWF]);
    return status;
  }

  if (ilat != nlat -1 ) 
    ilat_plus_one = ilat + 1; /* centered difference */
  else
    ilat_plus_one = nlat -1; /* backward difference at northern boundary */

  if (ilat != 0 ) 
    ilat_minus_one = ilat - 1; /* centered difference */
  else
    ilat_minus_one = 0;     /* forward  difference at southern  boundary */

  /* get temperature profile at ilat_minus_one */
  status = read_p_T_z_from_ECMWF_file( ncid, ilat_minus_one, nlat, ilon, nlon, nt, itime1, itime2, dt,
				       &(nlay), &(nlev), output->alt.altitude,
				       &(p_layer), &(theta_minus_one), &(z_layer),
				       input.verbose, input.quiet) ;

  if (status!=0) {
    fprintf (stderr, "Error %d, while read_p_T_z_from_ECMWF_file (1) from %s\n", status, input.filename[FN_ECMWF]);
    return status;
  }

  for (lc=0; lc<nlay; lc++)
    theta_minus_one[lc] = theta_minus_one[lc] * pow((1000./p_layer[lc]),kappa);

  free(p_layer);
  free(z_layer);

  /* get temperature profile at ilat_plus_one */
  status = read_p_T_z_from_ECMWF_file( ncid, ilat_plus_one, nlat, ilon, nlon, nt, itime1, itime2, dt,
				       &(nlay), &(nlev), output->alt.altitude,
				       &(p_layer), &(theta_plus_one), &(z_layer),
				       FALSE, input.quiet) ;

  if (status!=0) {
    fprintf (stderr, "Error %d, while read_p_T_z_from_ECMWF_file (2) from %s\n", status, input.filename[FN_ECMWF]);
    return status;
  }

  for (lc=0; lc<nlay; lc++)
    theta_plus_one[lc] = theta_plus_one[lc] * pow((1000./p_layer[lc]),kappa);

  free(p_layer);
  free(z_layer);

  /* get z profile at ilat */
  status = read_p_T_z_from_ECMWF_file( ncid, ilat, nlat, ilon, nlon, nt, itime1, itime2, dt,
				       &(nlay), &(nlev), output->alt.altitude,
				       &(p_layer), &(T_layer), &(z_layer),
				       FALSE, input.quiet) ;

  if (status!=0) {
    fprintf (stderr, "Error %d, while read_p_T_z_from_ECMWF_file (3) from %s\n", status, input.filename[FN_ECMWF]);
    return status;
  }

  free(p_layer);
  free(T_layer);

  dy = (ECMWF_lat[ilat_plus_one]- ECMWF_lat[ilat_minus_one])/360 * (2*PI*R_EARTH);

  if (input.verbose)
    fprintf (stderr, "     ilat-1=%zd (%6.2f), ilat+1=%zd (%6.2f), dy = %6.3f m\n", 
	     ilat_minus_one, ECMWF_lat[ilat_minus_one], ilat_plus_one, ECMWF_lat[ilat_plus_one], dy );

  if ((dtheta_dy = calloc (nlay, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for dtheta_dy\n");
    return -1;
  }

  for (lc=0; lc<nlay; lc++) {
    dtheta_dy[lc] = (theta_plus_one[lc] - theta_minus_one[lc]) / dy;
    /* fprintf (stderr, " %6.3f %8.3f %8.3f %13.6e\n", z_layer[lc], theta_minus_one[lc], theta_plus_one[lc], dtheta_dy[lc] ); */
  }

  if ((output->atm.microphys.dtheta_dy = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for dtheta_dy\n");
    return -2;
  }

  /* interpolation of  dtheta_dy to zout levels */
  status = arb_wvn (nlay, z_layer, dtheta_dy, 
                    output->atm.nzout, output->atm.zout_sea, output->atm.microphys.dtheta_dy, INTERP_METHOD_LINEAR, 1);

  if (status!=0) {
    fprintf (stderr, "Error %d interpolating d_theta_dy\n", status);
    return status;
  }


#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use the ECMWF data file option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ******************************************************************\n");
  return -1;
#endif

  return status;

}

/****************************************************************/
/* calculate_dtheta_dz                                          */
/*                                                              */
/* Purpose:                                                     */
/* calculate the vertical gradient of the potential temperature */
/* from T and p of the common atmosphere profiles and           */
/* interpolate it to the zout levels                            */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  input structure                                             */
/*  output structure                                            */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  output->atm.microphys.dtheta_dz                             */
/*            vertical gradient of theta                        */
/*                                                              */
/*  October 2007  by Ulrich Hamann                              */
/****************************************************************/

static int calculate_dtheta_dz (input_struct input, output_struct *output) 
{
  int status=0;

  float *theta     = NULL;
  float *dtheta_dz = NULL;

  int iz=0;
  /* int lc=0; */

  float   kappa = (C_P_DRY_STD-C_V_DRY_STD) / C_P_DRY_STD;  /* == R / cp */

  /* allocate potential temperature */
  if ((theta = calloc (output->atm.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for theta\n");
    return -1;
  }

  /* calculate potential temperature */
  for (iz=0; iz<output->atm.nlev; iz++)
    theta[iz] = output->atm.microphys.temper[0][0][iz] * pow((1000./output->atm.microphys.press[0][0][iz]),kappa);

  /* allocate vertical gradient of the potential temperature */
  if ((dtheta_dz = calloc (output->atm.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for dtheta_dz\n");
    return -2;
  }

  /* calculate dtheta/dz  0.001 == (dz in km) -> (dz in m) */
  /* forward difference in first layer */
  dtheta_dz[       0          ] = 0.001 * (theta         [        0         ] - theta         [        1         ])/
    (output->atm.zd[        0         ] - output->atm.zd[        1         ]);
  /* backward difference in first layer */
  dtheta_dz[output->atm.nlev-1] = 0.001 * (theta         [output->atm.nlev-2] - theta         [output->atm.nlev-1])/
    (output->atm.zd[output->atm.nlev-2] - output->atm.zd[output->atm.nlev-1]);  
  /* centered difference in all other layers */
  for (iz=1; iz<output->atm.nlev-1; iz++)
    dtheta_dz[iz]               = 0.001 * (theta         [      iz-1        ] - theta         [      iz+1        ])/
      (output->atm.zd[      iz-1        ] - output->atm.zd[      iz+1        ]);

  /* for (lc=0; lc<output->atm.nlev; lc++) */
  /*   fprintf (stderr, " lc=%3d, z=%7.3f, theta=%10.4f, d_theta_dz=%8.4f\n", */
  /*                      lc, output->atm.zd[lc], theta[lc], dtheta_dz[lc] );  */

  free(theta);

  /* allocate final result: vertical gradient of the potential temperature at zout levels */
  if ((output->atm.microphys.dtheta_dz = calloc (output->atm.nzout, sizeof (float))) == NULL) {
    fprintf (stderr,"Error allocating memory for dtheta_dz\n");
    return -2;
  }

  /* interpolation of  dtheta_dz to zout levels */
  status = arb_wvn (output->atm.nlev, output->atm.zd, dtheta_dz, 
                    output->atm.nzout, output->atm.zout_sea, output->atm.microphys.dtheta_dz, INTERP_METHOD_LINEAR, 1);

  free(dtheta_dz);

  /* for (iz=0; iz<output->atm.nzout; iz++) */
  /*   fprintf (stderr, " iz=%3d, z=%7.3f, d_theta_dz=%8.4f\n", */
  /*                      iz, output->atm.zout_sea[iz], output->atm.microphys.dtheta_dz[iz] ); */

  return status;
  
}
 
/****************************************************************/
/* read_ECMWF_ozone_climatology                                 */
/*                                                              */
/* Purpose:                                                     */
/* replace the ozone profile in the atmosphere data by the      */
/* ozone climatology, which is used by the ECMWF for            */
/* radiative transfer                                           */
/*                                                              */
/*  input:                                                      */
/*  ------                                                      */
/*  input structure                                             */
/*  output structure                                            */
/*                                                              */
/*  output:                                                     */
/*  -------                                                     */
/*  output->atm.microphys.dtheta_dz                             */
/*            vertical gradient of theta                        */
/*                                                              */
/*  October 2007  by Ulrich Hamann                              */
/****************************************************************/
static int  read_ECMWF_ozone_climatology (char *data_directory, float lat, struct tm UTC,
                                          float *dens_air_atm, int nlev_atm, float *press_atm, float **ozone_atm, 
                                          int verbose, int quiet )
{

#if HAVE_LIBNETCDF

  int status = 0;

  double *lat_grid = NULL;
  int i=0;
  int ilat[2]={0,0};
  int il=0;
  float dlat=0.0;
  size_t nlat = 0;
  int it=0;
  int itime[2]={0,0};
  float dt=0.0;  
  int nt=0;
  size_t nlev=0;

  int id_o3 = 0;

  double scale_factor = 1.0;
  double add_offset   = 0.0;

  float  ***tmp_ozone = NULL;
  double *tmp_press = NULL;

  float *log_p_atm=NULL;
  float *log_tmp_p=NULL;

  int ncid=0;
  char filename[FILENAME_MAX]="";

  size_t index3D[3] = {0,0,0};  /* time, lev, lat */
  int lc=0;

  int time_interpolate = TIME_INTERPOLATION;

  if (verbose)
    fprintf (stderr," ... replace ozone profile with climatology \n");

  /* create filename of ozone climatology according to data directory */
  strcpy (filename, data_directory );
  strcat (filename, "atmmod/ECMWF_ozone.nc");

  if ( lat < -90.0 || 90.0 < lat ) {
    fprintf (stderr, "Error, latitude %f outside range, maybe not specified (if == -999) in %s (%s)\n", lat, __func__, __FILE__);
    return -1;
  }

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &(ncid));
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s (line %d, function '%s' in '%s')\n", status, filename, __LINE__, __func__, __FILE__ );
    return status;
  }

  /* read latitude array */
  status = alloc_and_read_netCDF_1D_double(ncid,"lat", &(nlat), "lat", &(lat_grid));
  if (status != 0) {
    fprintf (stderr, "Error %d reading latitude (line %d, function '%s' in '%s')\n", status, __LINE__, __func__, __FILE__ );
    return status;
  }

  /* search correct latitude index */
  if (lat == lat_grid[0] ) {
    ilat[0]=0;
    ilat[1]=1;
  }
  else{ 
    for ( i=1; i < nlat; i++ ) {
      if( lat_grid[i-1] < lat  && lat <= lat_grid[i] ) {
        ilat[0]=i-1;
        ilat[1]=i;
        break;
      }
    }
  }
  dlat = (lat-lat_grid[ilat[0]])/(lat_grid[ilat[1]]-lat_grid[ilat[0]]);

  /* read pressure array */
  status = alloc_and_read_netCDF_1D_double(ncid, "press", &(nlev), "press", &(tmp_press));
  if (status != 0) {
    fprintf (stderr, "Error %d reading latitude (line %d, function '%s' in '%s')\n", status, __LINE__, __func__, __FILE__ );
    return status;
  }

  /* get time index */
  status = get_time_index (ncid, UTC, time_interpolate,
                           &(nt), &(itime[0]), &(itime[1]), &(dt),
                           verbose, quiet);
  if (status != 0){
    fprintf (stderr, "Error %d, during get_time_index (line %d, function '%s' in '%s')\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  if (verbose) {
    if ( UTC.tm_mday > 0 && nt == 2 )
      fprintf (stderr, "     read data at: lat =%7.2f/%7.2f (%7.2f), time = %s", 
	       lat_grid[ilat[0]], lat_grid[ilat[1]], lat, asctime(&UTC) );
    else
      fprintf (stderr, "     read data at: lat =%7.2f/%7.2f (%7.2f), \n", 
	       lat_grid[ilat[0]], lat_grid[ilat[1]], lat);

    if ( nt == 1 )
      fprintf (stderr, "     element: j=(%6d/%6d) (lat1/lat2), t=%6d (time)\n", ilat[0], ilat[1], itime[0]); 
    if ( nt == 2 ) 
      fprintf (stderr, "     element: j=(%6d/%6d) (lat1/lat2/%f), (%4d/%4d/%f) (time1/time2/dt)\n", ilat[0], ilat[1], dlat, itime[0], itime[1], dt);
  }

  /* get variable id for data */
  status = nc_inq_varid (ncid, "ozone", &id_o3);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error '%s', while getting id for %s (line %d, function '%s' in '%s')\n", 
	     nc_strerror(status), "ozone", __LINE__, __func__, __FILE__ );
    return status;
  }

  /* read attribute scale_factor */
  status = nc_get_att_double (ncid, id_o3, "scale_factor", &scale_factor);
  if (status != NC_NOERR) {
    /* no given scale_factor -> default scale_factor = 1.0 */
    /* that's OK, put status to 0 */
    status = 0;
  }
  else {
    /* if (verbose) */
    /*   fprintf (stderr,"     scale_factor = %f \n", scale_factor); */
  }

  /* read attribute add_offset */
  status = nc_get_att_double (ncid, id_o3, "add_offset", &add_offset);
  if (status != NC_NOERR) {
    /* no given offset -> default offset = 0.0 */
    /* that's OK, put status to 0 */
    status = 0;
  }
  else {
    /* if (verbose) */
    /*   fprintf (stderr,"     add_offset = %f \n", add_offset); */
  }

  /* allocate temporary data array, 2 times, 2 lats, nlev levels */
  if ( (status = ASCII_calloc_float_3D(&tmp_ozone,2,2,nlev)) != 0 ) {
    fprintf (stderr, "Error allocating memory for 'tmp_ozone' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  /* read ozone values */
  for (lc=0; lc<nlev; lc++) {
    index3D[1]=lc;
    for (it=0;it<2;it++) {
      index3D[0]=itime[it];
      for (il=0;il<2;il++) {
        index3D[2]=ilat[il];
        status = nc_get_var1_float  (ncid, id_o3, index3D, &(tmp_ozone[it][il][lc]));
        /* fprintf (stderr," it = %3d, il = %3d, lc = %3d, %s = %e \n", it, il, lc, "tmp_ozone", tmp_ozone[it][il][lc]); */
      }
      /* interpolate for different latitudes */
      tmp_ozone[it][0][lc] = (tmp_ozone[it][1][lc]-tmp_ozone[it][0][lc])*dlat + tmp_ozone[it][0][lc];
      /* fprintf (stderr," it = %3d, lc = %3d, %s = %e \n", it, lc, "tmp_ozone", tmp_ozone[it][0][lc]); */
    }
    /* time interpolation */
    tmp_ozone[0][0][lc] = (tmp_ozone[1][0][lc]-tmp_ozone[0][0][lc])*dt + tmp_ozone[0][0][lc];
    /* fprintf (stderr,"it = %3d, lc = %3d, %s = %e \n", it, lc, "tmp_ozone", tmp_ozone[0][0][lc]); */
    /* apply scale factor and offset */
    tmp_ozone[0][0][lc] = (tmp_ozone[0][0][lc] * scale_factor + add_offset);
  }

  /* adjust pressure units Pa -> hPa */
  for (lc=0; lc<nlev; lc++)
    tmp_press[lc]/=100.0;  
  
  /* adjust ozone units ppmv -> VMR */
  for (lc=0; lc<nlev; lc++)
    tmp_ozone[0][0][lc]/=1000000.0;
  
  /* if uppermost climatology pressure level is 0.0, than replace it with uppermost atm pressure */
  if ( tmp_press[0] == 0.0 ) {
    if ( press_atm[0] < 0.003 )     /* second level of the ozone climatology */
      tmp_press[0] = press_atm[0];
    else 
      tmp_press[0] = 0.00000001;    /* hm, any suggestions what to put in here? */
  }

  /* if lowermost climatology level has lower pressure than atmosphere, than replace it with atm pressure */
  if ( tmp_press[nlev-1] < press_atm[nlev_atm-1] )
    tmp_press[nlev-1] = press_atm[nlev_atm-1];

  /* allocate temporary data array */
  if ((log_p_atm = calloc(nlev_atm, sizeof(float)))==NULL) {
    fprintf (stderr, "Error: Allocation of 'log_p_atm' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -1;
  }
  
  for (lc=0; lc<nlev_atm; lc++) {
    log_p_atm[lc] = log(press_atm[lc]);
  }

  /* allocate temporary data array */
  if ((log_tmp_p = calloc(nlev, sizeof(float)))==NULL) {
    fprintf (stderr, "Error: Allocation of 'log_p' (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  for (lc=0; lc<nlev; lc++) {
    log_tmp_p[lc] = log(tmp_press[lc]);
  }

  if (verbose) {
    fprintf (stderr,"  lc     press      ozone[VMR]\n");
    fprintf (stderr,"-------------------------------\n");
    for (lc=0; lc<nlev; lc++)
      fprintf (stderr," %3d  %9.4f  %14.6e \n", lc, tmp_press[lc], tmp_ozone[0][0][lc]);
  }

  /* interpolate volumn mixing ratio to atm grid */
  status = arb_wvn (nlev, log_tmp_p, tmp_ozone[0][0], nlev_atm, log_p_atm, (*ozone_atm), INTERP_METHOD_LINEAR, ASCENDING);
  if (status!=0) {
    fprintf (stderr, "Error %d returned by arb_wvn during interpolation of 'ozone climatology'\n", status);
    fprintf (stderr, "      (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status;
  }

  /* write ozone to final place and convert Volumn Mixing Ratio to number density */
  for (lc=0; lc<nlev_atm; lc++) {
    (*ozone_atm)[lc] *= dens_air_atm[lc];
  }

  return status;

#else
  fprintf (stderr, " *******************************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot  use any           *\n");
  fprintf (stderr, " * netCDF option (here ECMWF_ozone_climatology). Please get netcdf and rebuild.*\n");
  fprintf (stderr, " *******************************************************************************\n");
  return -1;
#endif
    
}
 

  
