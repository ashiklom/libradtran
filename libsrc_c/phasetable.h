/*--------------------------------------------------------------------
 * $Id: phasetable.h 2698 2012-04-20 14:20:31Z robert.buras $
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

#ifndef __phasetable_h
#define __phasetable_h

#include <stdio.h>


/* phase function defined by moments or explicitely */
#define PHASE_NONE     0
#define PHASE_MOMENTS  1
#define PHASE_EXPLICIT 2
#define PHASE_HYBRID   3


/* table of phase function and cumulative phase function;  */
/* between the specified data points the phase function    */
/* is assumed to vary linearely in mu; in consequence, the */
/* cumulative phase function varies quadratic in mu; the   */
/* coefficients A, B, and C are used for the interpolation */
/* of mu between two adjacent cumulative phase function    */
/* values (needed for sampling of a random direction)      */
typedef struct {
  double **mu;   /* cosine of scattering angle;                     */ 
                 /* arbitrary values but sorted in increasing order */

  double ***p;   /* scattering phase matrix, normalized to 2      */

  double **F;    /* cumulative phase matrix, 0-2!                 */

    /* CE??? what to do here for polarization ???? */
  double **A;    /* coefficients for the interpolation of mu        */
  double **B;    /* between two adjancent values of the cumulative  */
  double **C;    /* phase function F                                */

  int *n;        /* number of data points                           */
  int nscales;
  int nphamat;    /* number of phase matrix elements                 */

  int truncate;   /* truncation yes/no                                         */ 
  double mutrnc;  /* truncation angle                                          */
  double ptrnc;   /* phase function at truncation angle                        */
  double Ftrnc;   /* cumulative phase function at truncation angle             */
  double dtrnc;   /* increment for truncation: p = ptrnc + dtrnc * (mu-mutrnc) */

  double dscale;  /*TZ ds*/

  int is_linear; /* switch that tells whether table is linear in theta */
  int is_one_theta_grid; /* switch that tells whether theta table is identical for all ip */
  double theta_0; /* value of theta at index 0 */
  double dthetainv; /* grid width in terms of theta, inverted */
  
} pft;

typedef struct {
  pft **iphase;    /* phase function table for each effective radius reff[i] */
  double *ext;     /* extinction coefficient for each effective radius       */
  double *dext;    /* extinction coefficient increment                       */

  double *sca;     /* scattering coefficient for each effective radius       */
  double *dsca;    /* scattering coefficient increment                       */
  double *f;       /* external delta-scaling factor CE (FU or Baum)         */
  double *dscale;  /* internal delta-scaling factor                    TZ ds*/
  double r0;       /* reff[i] = r0 + i * dr; i = 0 ... n-1                   */
  double dr;
  int    n;
} phase_function_table;

typedef struct {
  /* the following quantities are defined for each wavelength, */
  /* on an equidistant grid of effective radii,                */
  /* r0, r0+dr, r0+2*dr, ..., r0 + (nreff-1)*dr                */

  /* or on a nonequidistant grid stored in hum (relative humidity), if this structure
     is used for aerosols. */        

  int nlambda;
  
  double **extinc;      /* extinction coefficient           [iv][ireff]                */
  double **albedo;      /* single scattering albedo         [iv][ireff]                */
  double **f;           /* delta scaling factor                [iv][ireff]             */

  /* either define Legendre moments of the phase function                              */
  float  ****legen;     /* Legendre moments                 [iv][ireff][iphamat][ileg] */
  int    **nleg;        /* number of Legendre moments       [iv][ireff]                */

  /* read phase table as function of angle                                             */
  int  ***ntheta;       /* Number of phase function angles  [iv][ireff][iphamat]       */
  float ****theta;      /* Phase function angles            [iv][ireff][iphamat][imu]  */
  float ****phase;      /* Phase function                   [iv][ireff][iphamat][imu]  */
  double ****mu;        /* Phase function angles            [iv][ireff][iphamat][imu]  */

  double  *r0;          /* start radius for each wavelength [iv]                       */
  double  *dr;          /* wavelength step                  [iv]                       */
  size_t  *nreff;       /* number of effective radii        [iv]                       */ 
  double  **reff;       /* effective radius/relative humidity grid  [iv][ireff]        */  
  size_t   nphamat;     /* number of phase matrix elements                             */

  int   alloc;          /* flag, if basic memory is allocated                          */
  int   alloc_explicit; /* flag, if memory for explicit phase function is allocated    */
  int   alloc_moments;  /* flag, if memory for phase function moments is allocated     */

  int type;             /* phase function or Legendre moments                          */
} ssprop_struct;


void free_pft (pft *p);
int setup_Legendre_table (phase_function_table *phase,
			  ssprop_struct *ssprop,
			  int iv,
			  double rmin, double rmax, 
			  double delta_scaling_mucut, double truncate,
			  int quiet); /*TZ ds*/

int create_iphase_from_HG (pft *iphase, double g1, int quiet);

int calc_cumulative_table (double **mu, float **p, int *n, int nphamat, double truncate, 
			   pft *iphase, double delta_scaling_mucut, int quiet); /*TZ ds*/

int calc_Legendre_phase (double **pmom, int nmom, int nphamat, double ***mu, float ***p, int **n);

int read_phase_function (char *filename, 
			 double *r0, double *dr, size_t *nreff, size_t *nphamat,
			 double **extinc, double **albedo, int **ntheta, 
			 float ***theta, float ****phase, int *type, int quiet);

void calloc_phase_function_table (phase_function_table * phase, int n, double r0, double dr);
void free_phase_function_table (phase_function_table *phase);

int sort_theta_and_mu (int n_in, int *nmu_in,
		       float **theta_in, float *theta_out,
		       double **mu_in, double *mu_out);

int interpolate_phase_weighted (int n_in, int *nmu_in, double **mu_in, float **phase_in,
				double *weight,
				int nmu_out, double *mu_out, float *phase_out,
				int returnmax);

void calc_F_and_normalize (double **mu, float **p, int *n, int nphamat);

void normalize_phase (double **mu, float **p, double *F, int *n, int nphamat, int verbose);

int sort_and_add_weighted_phase (int n_in, double *epsilon,
				 int **ntheta, float ***theta, double ***mu, float ***phase,
				 int **ntheta_new, float ***theta_new, double ***mu_new,
				 float ***phase_new, int nphamat, int returnmax,
				 int quiet);

int calloc_iphase (pft *iphase, int nscales, int nphamat, int *n);
int free_iphase (pft *iphase);

void calc_iphase_coeffs (pft *iphase, int SC_mode);

int create_phase_from_Rayleigh (float depol, int *ntheta, float **theta, double **mu, float **phase, int only_backdir, int ip);

int create_phase_from_HG (double g, int *ntheta, float **theta, double **mu, float **phase, int only_backdir);
int choose_HG_grid ( double g, int *ntheta, float **theta);

void calc_pmom (double *mu, double *phase, int ntheta, int nmom,
		double *pmom);



void get_phase_matrix_rayleigh ( double  mu,
				 double  depol,
				 int     np,
				 /* output */
				 double *phase_matrix );

int get_phase_matrix_pft ( pft    *iphase,
			   double  mu,
			   int     scaled,
			   int     np,
			   /* output */
			   double *p );

int get_phase_matrix_pft_interpol_reff ( double                mu,
					 double                reff,
					 phase_function_table *phase_cld,
					 int                   scaled,
					 int                   np,
					 /* Output */
					 double               *phase_matrix_cld );

#endif



