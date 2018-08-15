/*--------------------------------------------------------------------
 * $Id: miecalc.h 3107 2015-05-12 12:17:19Z Claudia.Emde $
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

#ifndef __miecalc_h
#define __miecalc_h

#include <stdio.h>

typedef struct {
  float re;
  float im;
} mie_complex;

#define USER     1
#define ICE      2
#define WATER    3
#define AEROSOL  4

#define MIEV0    1
#define BHMIE    2

#define READ_SPF_PHASES      1
#define READ_SPF_MOMENTS     2
#define READ_SPF_BOTH        3
#define READ_SPF_ALL_MOMENTS 4

/* Wavelength input grid STRUCTURE */

typedef struct {
  double   wl_start;
  double   wl_end;
  double   wl_step;
} mie_wl_inp_struct;

/* Wavelength output grid STRUCTURE */

typedef struct {
  int      nlambda;
  double*  lambda;
} mie_wl_out_struct;

typedef struct {
  mie_complex *s1, *s2;  /* Mie scattering amplitudes */
  int         anyang;    /* If true xmu any angles may be given in xmu */
  int         ipolzn;    /* Type of pmom */
  int         momdim;    /* First dimension of pmom */
  int         nmom;      /* Highest Legendre moment */
  int         numang;    /* No. of angles in xmu */
  int         perfct;    /* If true assume refractive index is infinite */
  int         prnt[2];   /* Print flag */
  float       mimcut;    /* Imaginary component zeroed below mimcut */
  float       *xmu;      /* Cosine of angles for s1 and s2 */
  int quiet;             /* shut up! */
} mie_inp_struct;

typedef struct {
  mie_complex crefin;    /* Complex refractive index */
  mie_complex sback;     /* Back scattering amplitude */
  mie_complex sforw;     /* Forward scattering amplitude */
  mie_complex tback[2];  /* See miev.doc */
  mie_complex tforw[2];  /* See miev.doc */
  float       qext;      /* Extinction efficiency factor */
  float       qsca;      /* Scattering efficiency factor */
  float       qback;     /* Backscattering efficiency factor */
  float       gsca;      /* Asymmetry factor */
  float       gqsc;      /* Asymmetry times scattering efficiency factor */
  float       spike;     /* Magnitude of smallest denominator of Mie coefficients */
  double      xxd;        /* Mie size parameter */
  float       **pmom;    /* Phase function moments */
} mie_out_struct;


double mom2phase (double x, double *f, int L);

int mie_calc (mie_inp_struct input, mie_out_struct *output, 
	      int program, int medium, mie_complex crefin, 
	      double wavelength, double radius, 
	      float temperature, int nstokes, mie_complex *ref);

int mie_calc_sizedist (mie_inp_struct input, mie_out_struct *output, 
		       int program, int medium, mie_complex crefin, 
		       double wavelength, float temperature, int nstokes,
		       double *x_size, double *y_size, int n_size,
		       double *beta, double *omega, double *g,
		       mie_complex *ref);

int phase_function (float *moment, int L);

int read_mie_table (char *filename, 
		    double *r0, double *dr, size_t *nreff, double **reff, 
		    double *wavelen, double *nre, double *nim,
		    double **extinc, double **albedo, double **f, int **nleg,
		    float ****legen, size_t *nphamat,
 		    int ***ntheta, float ****theta, double ****mu, float ****phase,
		    int *alloc_moments, int *alloc_explicit,
		    int aerosol, int nstokes, int nstrmax,
		    int read_scattering_phase_function,
		    int iws, int niw, int quiet);

int read_mie_table_lambda (char *filename, size_t *nlam,
			   double **wavelen, int quiet);

int nc_inq_dimid_err (int ncid, char *varstr, int *id_var, char *cdf_filename, int quiet );
int nc_inq_varid_err (int ncid, char *varstr, int *id_var, char *cdf_filename, int quiet );
int nc_inq_dimlen_err (int ncid, char *varstr, int id_var, char *cdf_filename, size_t *var );
int nc_get_vara_int_err (int ncid, char *varstr, int id_var, size_t *start, size_t *count, char *cdf_filename, int *var );
int nc_get_vara_float_err (int ncid, char *varstr, int id_var, size_t *start, size_t *count, char *cdf_filename, float *var );
int nc_get_vara_double_err (int ncid, char *varstr, int id_var, size_t *start, size_t *count, char *cdf_filename, double *var );

int optimize_theta_grid(double **grid_opt, double **data_opt, int *Nopt,
                        double *grid_fine, double *data_fine, int Npoints,
                        double acc, int nthetamax,
                        int *warnings_ntheta);

int optimize_theta_grid_phamat(double **grid_opt,
                               double ***data_opt,
                               int *Nopt, 
                               double *grid_fine,
                               double **data_fine,
                               int Npoints, 
                               double acc,
                               int nthetamax,
                               int nstokes,
                               int *warnings_ntheta);

int element(int a, int*b, int L);

int element_double(float a, float *b, int L);

int find_opt_index(int index_acc, int *indices, int Nind,
                   int decimal, int Npoints);

int int_cmp(const void *a, const void *b);

#endif

