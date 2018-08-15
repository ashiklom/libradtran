/*--------------------------------------------------------------------
 * $Id: ancillary.h 3276 2017-07-04 14:16:36Z Claudia.Emde $
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

#ifndef __ancillary_h
#define __ancillary_h

#include <stdio.h>
#include "uvspec.h"

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#define CAOTH_MOL 0
#define CAOTH_AER 1
#define CAOTH_FIR 2

int setup_wlgrid (input_struct input, output_struct *output);

int setup_rte_wlgrid (input_struct input, output_struct *output);

int reptran_filename(input_struct input,
                    int i_mol,
                    char *filename);

int arb_wvn( int n_old_x, float *old_x, float *old_y,  
	     int n_new_x, float *new_x, float *new_y, 
             int linear, int descend);

int arb_wvn_double( int n_old_x, double *old_x, double *old_y,  
	     int n_new_x, double *new_x, double *new_y, 
             int linear, int descend);

int arb_wvn_zero (int n_old_x, float *old_x, float *old_y,  
		  int n_new_x, float *new_x, float *new_y, 
		  int linear, int descend);

int arb_wvn2 (int n_old_x, float *old_x, float *old_y,  
	      int n_new_x_lower, int n_new_x_upper, float *new_x, float *new_y, int linear);
int arb_wvn2_double (int n_old_x, double *old_x, double *old_y,  
		     int n_new_x_lower, int n_new_x_upper, double *new_x, double *new_y, int linear);

int interpolate_profile (float *x, float **y, int n,  
			 float *xnew, int nnew, int interpol_method, 
                         int quiet );

int interpolate_density (float *x, float **y, int n, 
                         float *xnew, int nnew, int interp_method, 
                         float *dens_air_old, float *dens_air, int quiet );
     
int interpolate_atmosphere (float *z, float ****p_p, float ****p_T, 
			    float *****p_dens, float ****p_Tavg, float *****p_densavg, 
			    int n, float *znew, int nnew, 
                            int interpol_method_press, int interpol_method_temper,
			    int *interpol_method_gas, 
                            int quiet);

int convolve                  (input_struct input, output_struct *output);

int spline_interpolate        (input_struct input, output_struct *output);

int internal_to_transmittance_grid (input_struct input, output_struct *output);

int interpolate_transmittance (input_struct input, output_struct *output);

int multiply_extraterrestrial (input_struct input, output_struct *output);

int optical_properties(input_struct input, output_struct *output, 
		       double  wvl, int ir, int iv1, int iv2, int iq, 
		       int verbose, int skip_optical_properties);

int processing1D      (input_struct input, output_struct *output);
int integrate1D       (input_struct input, output_struct *output);
int sum1D             (input_struct input, output_struct *output);

int processing3D      (input_struct input, output_struct *output);
int integrate3D       (input_struct input, output_struct *output);
int sum3D             (input_struct input, output_struct *output);

int double2float_integrated_values(input_struct input, output_struct *output);

int calloc_int_values (input_struct input, output_struct *output);

void error_calloc     (char *variable, char *function, int *status);

int scaling_integrated_values (input_struct input, output_struct *output);

int average_dens (float *dens, float *dens_air, float *zd, int nlev, int interpol_method, float **dens_avg, int allocate);

float log_average          (float d1,float d2);
double dlog_average        (double d1,double d2);
float linmix_average       (float m1,float m2,float n1,float n2);
float mass_weighted_average(float x1,float x2,float n1,float n2);
int   spline_average       (float *x, float *y, int n, float **y_average);

int scale_output (input_struct input,
		  float ***p_rfldir, float ***p_rfldn,  float ***p_flup, float ***p_albmed, float ***p_transmed,
		  float ***p_uavgso, float ***p_uavgdn, float ***p_uavgup,
		  float ***p_uavg, float ****p_u0u, float *****p_uu, 
                  float ***p_heat, float ***p_emis, float ***p_w_zout,
		  float ****p_down_flux, float ****p_up_flux, float ******p_down_rad, float ******p_up_rad,
		  float *****p_rfldir3d, float *****p_rfldn3d, float *****p_flup3d, float *****p_uavgso3d,
		  float *****p_uavgdn3d, float *****p_uavgup3d, float *******p_radiance3d, float *****p_absback3d, 
		  float *****p_rfldir3d_var, float *****p_rfldn3d_var, float *****p_flup3d_var, float *****p_uavgso3d_var,
		  float *****p_uavgdn3d_var, float *****p_uavgup3d_var, float ******p_radiance3d_var, float *****p_abs3d_var, float *****p_absback3d_var, /* **CK added float *****p_abs3d_var for forward mc_std */
		  int nzout, int Nx, int Ny, int Nz, int Nc, int *threed, int passback3D,
		  int islower, int isupper, int jslower, int jsupper, 
		  float *****p_abs3d,       
		  double ffactor, double rfactor, double hfactor, int iv);

int alloc_and_read_netCDF_1D_double (int ncid, char *dimension, 
                                     size_t *dim, char *variable, double **data);

int get_grid_index (float number, double *array, int n_array, int periodic);

int get_time_index (int ncid, struct tm UTC, int time_interpolate, 
                    int *nt, int *t1, int *t2, float *dt,
                    int verbose, int quiet);

int get_all_netCDF_indices ( char* filename, float lat, float lon,
                                int *ncid, long *ilat, size_t *nlat, long *ilon, size_t *nlon, 
                                double **lat_grid, double **lon_grid,
                                struct tm UTC, int time_interpolate,
                                int *nt, int *itime1, int *itime2, float *dt,
                                int verbose, int quiet);

int read_p_T_z_from_ECMWF_file ( int ncid, long ilat, size_t nlat, long ilon, size_t nlon,
                                 int nt, int itime1, int itime2, float dt,
                                 size_t  *nlay, size_t  *nlev, float altitude,
                                 float **p_layer, float **T_layer, float **z_layer,
                                 int verbose, int quiet );

int get_local_apparent_time_index (int ncid, struct tm LAT, int time_interpolate, 
                                   int *nt, int *t1, int *t2, float *dt,
                                   int verbose, int quiet);

int get_number_from_netCDF_map (float lat, float lon, struct tm UTC, int time_interpolate, char *filename, 
                                void *data, int external_type, char *variable_name, 
                                int verbose, int quiet);

int read_sat_geometry ( int pixel_x, int pixel_y, struct tm UTC, char *filename,
                        float *latitude, float *longitude, int *numu, int *maxumu, float **umu, int *nphi, int *maxphi, float **phi,
                        int solver, int verbose, int quiet );

int read_u10_from_map (float lat, float lon, struct tm UTC, int time_interpolate, char *filename,
                       float *u10, int verbose, int quiet);

int read_wind_from_ECMWF_file (float lat, float lon, struct tm UTC, int time_interpolate, char *filename, float altitude,
                               float **u, float **v, float **w, float **z_wind, size_t *nlev, int verbose, int quiet);

int alloc_and_read_ECMWF_netCDF_pressure (int ncid, float **p_level, size_t *nlev, float **p_layer, size_t *nlay,
                                          size_t itime, size_t ilat, size_t ilon, int verbose);

int read_ECMWF_surface_pressure ( int ncid, size_t itime, size_t ilat, size_t ilon, float *SP, int verbose );

int alloc_and_read_netCDF_column_float (int ncid, char *variable_name, float **data, 
                                       size_t nlev, size_t itime, size_t ilat, size_t ilon, int verbose);

int calculate_z_from_p_and_T (float *p,float *T, float **z, int n, float altitude, int verbose);

int weekday(int year, int month, int day);

char *strtrim(char *str, const char *trim);
char *strrtrim(char *str, const char *trim);
char *strltrim(char *str, const char *trim);

int specific_heat_capacity_moist_air (float T_in_K, float dens_air, float dens_wv, float *c_p, int quiet);
float specific_heat_capacity_dry_air      (float T_in_K, int quiet);
float specific_heat_capacity_water_vapour (float T_in_K, int quiet);

time_t my_timegm (struct tm *tm);

int crs_raman_N2 (double lambda, int n_transitions, float temper, 
		  double ***crs, int *iv, int verbose);
int crs_raman_O2 (double lambda, int n_transitions, float temper, 
		  double ***crs, int *iv, int verbose);

float dewpoint (float press_h2o);
float Tcon(float T, float TD);
float EPT(float T, float TD, float p, float N_AIR, float N_H2O);

#endif
