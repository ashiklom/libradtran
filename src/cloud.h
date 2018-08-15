/*--------------------------------------------------------------------
 * $Id: cloud.h 3242 2016-08-18 22:10:06Z bernhard.mayer $
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

#ifndef __cloud_h
#define __cloud_h

#include <stdio.h>
#include "uvspec.h"


#define MAX_DEFF_FU96 130.24 /* maximum effective diameter from FU96 in um */
#define MIN_DEFF_FU96  18.63 /* minimum effective diameter from FU96 in um */

#define MAX_DEFF_FU98 129.6  /* maximum effective diameter from FU98 in um */
#define MIN_DEFF_FU98  11.0  /* minimum effective diameter from FU98 in um */


int setup_all_caoth (input_struct input, output_struct *output);

int setup_caoth (input_struct       input,
		 output_struct     *output,
		 caoth_inp_struct   input_caoth,
		 caoth_out_struct  *output_caoth,
		 caoth_out_struct **output_caoth_ipa );

int calloc_caoth_out ( caoth_out_struct *out,
		       char             *name,
		       char             *fullname,
		       int               nlambda,
		       int               nlev,
		       int               nphamat,
		       int               calloc_ssprop );
int free_caoth_out (caoth_out_struct *out, int free_ssprop);
int calloc_optprop (optprop_struct *out, int nlambda, int nlev, int nphamat);
int free_optprop   (optprop_struct *optprop);
int calloc_ssprop_struct (ssprop_struct *out, int nlambda); 
int free_ssprop_struct (ssprop_struct *out); 

int cp_optprop ( optprop_struct *target,
		 optprop_struct *source, 
		 int             nlev,
		 int             alloc_moment,
		 int             quiet );

int cp_caoth_out ( caoth_out_struct *target,
		   caoth_out_struct  source,
		   int               copy_ssproc,
		   int               alloc_moment,
		   int               quiet );

int cp_caoth3d_out ( caoth_out_struct   *target,
		     caoth3d_out_struct  source,
		     int                 copy_ssprop,
		     int                 alloc_moment,
		     int                 quiet,
		     int                 iipa,
		     int                 jipa );

int interpolate_ssprop_in_lambda ( ssprop_struct *src,
				   ssprop_struct *tgt,
				   double        *src_lambda,
				   int            src_nlambda,
				   float         *tgt_lambda,
				   int            tgt_nlambda,
				   int            nlambda_lower,
				   int            nlambda_upper );

int read_and_convert_caoth_file ( input_struct      input,
				  char             *filename,
				  wl_out_struct     wl_out,
				  float             altitude,
				  atm_out_struct    atm_out,
				  caoth_inp_struct  input_caoth,
				  /* output */
				  caoth_out_struct *output_caoth,
				  cf_out_struct    *cf_out,
				  int              *nmom );

int ssprop2optprop (caoth_out_struct *caoth, int iv, int quiet);

int aer_ssprop2optprop (aer_out_struct *aer, int iv, int i_aer, float *rh,
			       int nlambda, int nlambda_lower, int nphamat, int verbose);

int read_caoth_prop (char *filename, float *lambda_r, int nlambda_r, 
		     int nlambda_lower, int nlambda_upper,
		     int interpolate, int aerosol,
		     caoth_out_struct *caoth, int caoth_properties,
		     aer_out_struct *aer, float *rh, int i_aer,
		     int *nmom, int nstokes, int nstrmax, int solver,
                     int disort_icm, int verbose, int quiet);

int read_caoth_file ( int               source,
		      char             *caoth_name,
		      char             *caoth_fullname,
		      char             *filename,
		      /* input specific for ECMWF and ECHAM */
		      char             *filename_ic_reff,
		      float             latitude,
		      float             longitude,
		      struct tm         UTC,
		      int               time_interpolate,
		      int               reff_prop,
		      int               properties,
		      float             reff_fixed,
		      float             altitude,
		      float            *press_atm,
		      int               cloud_overlap,
		      float            *z_model_layer, 
		      /* other input */
		      int               optical,
		      int               nlambda,
		      int               layer, 
		      int               verbose,
		      int               quiet,
		      /* output */
		      caoth_out_struct *caoth,
		      cf_out_struct    *cf );

int apply_user_defined_properties_to_caoth ( caoth_inp_struct  cldin, 
					     int               nlambda_r,
					     float            *lambda_r,
					     float             altitude,
					     /* output */
					     caoth_out_struct *cldout );

int interpolate_ssprop_in_reff (/* Input */
				double reff,
				int nreff, double r0, double dr,
				double *extinc, double *albedo, double *fdelta,
				int interpolate_moments, 
				float ***legen, int *nleg,
				int interpolate_explicit,
				int **ntheta, float ***theta, double ***mu, float ***phase,
				int nphamat,
				/* Output */
				float *ext, float *alb, float *f,
				float ***coeffc, int *ncoeffc,
				float ***phase_new, float ***theta_new, double ***mu_new,
				int **ntheta_new );

int ic_fu96 (float wavelength, int nlyr,
	     char *path,
	     float *iwc, float *reff, 
	     float *dtau, float *g, float *ssa, float *f,
	     float *zd, int layer, int unscaled);

int ic_fu98 (float wavelength, int nlyr,
	     char *path,
	     float *iwc, float *reff, 
	     float *dtau, float *g, float *ssa, 
	     float *zd, int layer);

int read_isccp_reflectivity (int type, float sza, float tau, char *path, int quiet,
			     float *ref);

int wc_echam4 (float wavelength, int nlyr,
	       float *lwc, float *reff, 
	       float *dtau, float *g, float *ssa,
	       float *zd, int layer);

int ic_echam4 (float wavelength, int nlyr,
	       float *iwc, float *reff, 
	       float *dtau, float *g, float *ssa,
	       float *zd, int layer);

int copy_cloud_fraction(cf_out_struct *target, cf_out_struct source, int alloc);

int alloc_cloud_fraction (cf_out_struct *cf, int nlev);

int ic_yang (float wavelength, int nlyr, int habit,
	     char *path,
	     float *iwc, float *reff, 
	     float *dtau, 
	     float *f, float *g1, float *g2, float *ssa, 
	     float *zd, int layer, int newkey);

int caoth_prop_switch ( input_struct      input,
			caoth_inp_struct  input_caoth,
			wl_out_struct     wl_out,
			int               iv,
			/* Output */
			caoth_out_struct *output_caoth );

int read_raytracing_file (char *filename, crystal_prop_struct **raytracing_prop, int *n_raytracing_prop, int quiet);

#endif



