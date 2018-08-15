/*--------------------------------------------------------------------
 * $Id: aerosol.c 3276 2017-07-04 14:16:36Z Claudia.Emde $
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
#include <math.h>
#include <glob.h>

#include "uvspec.h"
#include "ascii.h"
#include "numeric.h"
#include "miecalc.h"
#include "fortran_and_c.h"
#include "f77-uscore.h"
#include "cloud.h"

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#define EPSILON 1E-6

/*  this variable is only required to avoid warnings        */
/*  during compilation, like                                */
/*  warning: 'extrapolate_atmosphere' defined but not used  */
/*  Once the function is used, the variable may be deleted  */
/*  from this file.                                         */
/* #define EXTRAPOLATE_ATMOSPHERE */




/************************************/
/* prototypes of internal functions */
/************************************/

#ifdef EXTRAPOLATE_ATMOSPHERE
static int extrapolate_atmosphere (float altitude, atm_out_struct *out, 
                                   int quiet, int verbose);
#endif

static int calculate_relative_humidity_on_aerosol_grid (input_struct input, output_struct *output, float **tmp_rh);

static int read_aerosol_profile (char *filename, int nlyr, int nlambda_r, 
				 float *zd_common, float **dtau, int scale);

static int read_aerosol_profile_from_map (char *filename, float latitude, float longitude, struct tm UTC, int time_interpolate, 
                                          float *z_atm, float *press, int nlev_atm,
                                          int nlyr, int nlambda_r, float *zd_common, 
                                          float **dtau, int scale, int verbose, int quiet);

static int read_aerosol_moments (char *filename,
				 float **moment, int *nmom);

static int aerosol_mie_calculation ( input_struct input, output_struct *output );

static int read_aerosol_species    ( input_struct input, output_struct *output );

static int read_aerosol_species_from_netCDF_map ( input_struct input, output_struct *output, char *filename );

static int read_aerosol_species_from_ASCII_file ( char *filename, aer_out_struct *aer, 
                                                  int quiet, int verbose );

static int aerosol_species_to_optical_properties ( input_struct input, output_struct *output );

static int average_aerosol_prop_opac ( aer_out_struct *aer, int nlambda, int nlambda_lower, int nlambda_upper );

static int scale_aer_tau (float **dtau, float *lambda_r, int nlambda_r, int nlyr,
			  float tau, float *zd, float altitude);

static int scale_aer_angstrom (float **dtau, float *lambda_r, int nlambda_r, int nlyr,
			       float alpha, float beta, float *zd, float altitude);
			       
static int scale_aer_king_byrne (float **dtau, float *lambda_r, int nlambda_r, int nlyr,
			       float alpha_0,float alpha_1, float alpha_2, float *zd, float altitude);			       
			       
static int aerosol_set_visibility (aer_inp_struct inp, char *path, float *visib, int verbose);





/**************************************************************/
/* Extrapolate atmospheric profiles for altitudes,            */
/* which are below the levels given by the atmosphere_file    */
/**************************************************************/

#ifdef EXTRAPOLATE_ATMOSPHERE
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
    fprintf (stderr,"Error, reallocating memory for output->atm.zd in extrapolate_atmosphere() (in atmosphere.c)\n");
    return -1;
  }
  out->zd[out->nlev] = altitude;


  /* realloc temperature */
  if ((out->microphys.temper = realloc (out->microphys.temper, (out->nlev+1) * sizeof(float))) == NULL) {
    fprintf (stderr,"Error, reallocating memory for output->atm.zd in extrapolate_atmosphere() (in atmosphere.c)\n");
    return -1;
  }
  out->microphys.temper[out->nlev] = (out->microphys.temper[out->nlev-1]-out->microphys.temper[out->nlev-2])/
                                     (out->zd[out->nlev-1] - out->zd[out->nlev-2] ) * 
                                     (out->zd[out->nlev] - out->zd[out->nlev-1]) 
                                     + out->microphys.temper[out->nlev-1];

  /* realloc pressure */
  if ((out->microphys.press = realloc (out->microphys.press, (out->nlev+1) * sizeof(float) )) == NULL) {
    fprintf (stderr,"Error, reallocating memory for output->atm.zd in extrapolate_atmosphere() (in atmosphere.c)\n");
    return -1;
  }

  g = g_surface * (r_earth * r_earth) / ((r_earth+out->zd[out->nlev-1]*1000.)*(r_earth+out->zd[out->nlev-1]*1000.));
  /* hydrostatic equation assuming linear temperature gradient and */ 
  /* NO variation of g within the last layer */
  /* dz = (out->zd[out->nlev] - out->zd[out->nlev-1]) * 1000.0; */
  /* or WITH variation of g within the last layer */
  dz = ( r_earth*r_earth/(r_earth + 1000.0*out->zd[out->nlev-1]) - r_earth*r_earth/(r_earth + 1000.0*out->zd[out->nlev]) );
  out->microphys.press[out->nlev] = out->microphys.press[out->nlev-1] * 
    exp( -(g*dz) / (R_air * 0.5 * (out->microphys.temper[out->nlev] + out->microphys.temper[out->nlev-1])));  

  /* realloc gas densities */
  for (gas=0;gas<MOL_NN;gas++) {
    if ((out->microphys.dens[gas] = realloc (out->microphys.dens[gas], (out->nlev+1) * sizeof(float) )) == NULL) {
      fprintf (stderr,"Error, allocating memory for dens in extrapolate_atmosphere() (in atmosphere.c)\n");
      return -1;
    }
  }
  /* calculate air number density from pressure and temperature */
  out->microphys.dens[MOL_AIR][out->nlev] = out->microphys.press[out->nlev] * 1E-04 / (boltzmann*out->microphys.temper[out->nlev]);
  /* for all other gases use constant mixing ratio */
  for (gas=0;gas<MOL_NN;gas++) {
    switch (gas) {
    case MOL_AIR:
    case MOL_O4:
      /* do nothing */
      break;
    /* this is for constant relative humidity */
    /* case MOL_H2O: */  
    /*  rel_hum = out->microphys.dens[gas][out->nlev-1] / vapor_pressure(out->microphys.temper[out->nlev-1]) * 100.0; */
    /*  out->microphys.dens[gas][out->nlev] = rel_hum / 100.0 * vapor_pressure(out->microphys.temper[out->nlev]); */
    /*  fprintf (stderr,"rel_hum = %5.2f, %e, %e \n",rel_hum, */
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
      mix_ratio2 = out->microphys.dens[gas][out->nlev-2] / out->microphys.dens[MOL_AIR][out->nlev-2];
      mix_ratio1 = out->microphys.dens[gas][out->nlev-1] / out->microphys.dens[MOL_AIR][out->nlev-1];
      mix_ratio  = mix_ratio1 + (mix_ratio1-mix_ratio2) / (out->zd[out->nlev-1] - out->zd[out->nlev-2] ) * 
	                              (out->zd[out->nlev] - out->zd[out->nlev-1]);
      out->microphys.dens[gas][out->nlev] = mix_ratio * out->microphys.dens[MOL_AIR][out->nlev]; 
      break;

    default:
      fprintf (stderr, "Error: unkown gas %d in setup_gases \n", gas);
      return -1;
      break;
    }
  }

  out->nlev = out->nlev + 1;

  return status;
}
#endif

/**************************************************************/
/* Read file holding  filenames for each layer.               */
/* Setup atmospheric aerosol.                                 */
/**************************************************************/

int setup_aerosol (input_struct input, output_struct *output)
{
  int status=0, ipa=0, iv=0, lc=0, k=0, nlyr=0, ialt=0, ic=0;
  static int maxawvn = MAXAWVN;

  int linear=0;
  int refrac=0, sizedist=0;

  int tau_from_netCDF=FALSE;
  
  int rows=0, min_columns=0, max_columns=0, max_length=0;

  float *tmp_aerext=NULL, *tmp_aerabs=NULL, *tmp_aersym=NULL, *tmp_aerwvn=NULL;  
  float *tmp_aer_dtau=NULL, *tmp_aer_gg=NULL, *tmp_aer_ssa=NULL;
  float *tmp_rh     = NULL;

  float *lref=NULL, *re=NULL, *im=NULL;
  int nref=0;
  
  double *tautot = NULL, *lambda = NULL;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double tauwvl=0, scale_factor=0;
  double settau=0, settauwvl=0;

  float *moment=NULL;

  float *tmp_zd=NULL;
  int tmp_nlev=0;

  float *aer_pro=NULL;
  
  int lc_6=0, lc_alt=0;
  float *new_z=NULL, *new_prof=NULL; 

#if HAVE_LIBNETCDF
  int ncid=0;
#endif

  aer_out_struct *tmpaer = calloc (1, sizeof(aer_out_struct));
  
  float *zalt=NULL;

  void F77_FUNC (aerprof, AERPROF) (int *nlyr, float *zd, float *aer_pro, 
			  int *seasn, int *vulcan, float *visib);
  void F77_FUNC (aeropt, AEROPT)  (int *vulcan, float *aerext, float *aerabs, 
			  float *aersym, float *aerwvn, float *rh, 
			  int *nlyr, float *zd, int *haze, int *maxawvn); 
  void F77_FUNC (aerwvn, AERWVN)  (float *aerwvn, float *aerabs, float *aersym, float *aerext,
			  float *aer_pro, float *aer_dtau, float *aer_gg, 
			  float *aer_ssa, int *nlyr, float *lambda_r, float *zd, 
			  int *nlambda_r, int *maxawvn);

  char function_name[] = "setup_aerosol";
  char file_name[]     = "aerosol.c";

  
  if (input.aer.vulcan > 0) {

#if HAVE_LIBNETCDF
    
    /* test if input file is netCDF format */
    if (strlen(input.aer.filename[FN_AER_TAU]) > 0) {
      
      /* open netcdf file */
      status = nc_open (input.aer.filename[FN_AER_TAU], NC_NOWRITE, &ncid);
      if (status==NC_NOERR) {
        tau_from_netCDF = TRUE;
        nc_close(ncid); 
      }
    }

#endif

    /* first, determine the complete set of altitudes   */
    /* for the aerosol; to do that, combine the z-grids */
    /* of all profiles                                  */
    
    /* we want all atmospheric levels in the aerosol profile  */
    /* because the default aerosol (Shettle) is defined       */
    /* on the atmospheric grid; user-defined levels are       */
    /* required as well, hence we use the common grid as      */
    /* input which contains all atmospheric levels as well as */
    /* the user-defined ones                                  */

    output->aer.nlev = 0; 
    set_common_z_grid (output->atm.zd_common, output->atm.nlev_common, 
                       &tmp_zd, &tmp_nlev);

    status=0;
    if (strlen(input.aer.filename[FN_AER_EXPLICIT]) > 0)
      status += add_file2grid (&tmp_zd, &tmp_nlev, input.aer.filename[FN_AER_EXPLICIT]);
    
    if (strlen(input.aer.filename[FN_AER_TAU]) > 0) { 
      if (!tau_from_netCDF)
        status += add_file2grid          (&tmp_zd, &tmp_nlev, input.aer.filename[FN_AER_TAU]);
      else 
        status += add_sigma_nc_file2grid (&tmp_zd, &tmp_nlev, 
					  input.aer.filename[FN_AER_TAU],
					  input.latitude,
					  input.longitude, input.UTC, 
					  input.atm.time_interpolate,
					  output->atm.zd, 
					  output->atm.microphys.press[0][0], 
					  output->atm.nlev,
					  input.verbose, input.quiet);
    }
    
    if (strlen(input.aer.filename[FN_AER_GG]) > 0)
      status += add_file2grid (&tmp_zd, &tmp_nlev, input.aer.filename[FN_AER_GG]);
    
    if (strlen(input.aer.filename[FN_AER_SSA]) > 0) 
      status += add_file2grid (&tmp_zd, &tmp_nlev, input.aer.filename[FN_AER_SSA]);
    
    if (status!=0) {
      fprintf (stderr, "Error creating combined aerosol grid\n");
      return status;
    }

    /* allocate memory for profiles */
    status = calloc_aer_out (&(output->aer), output->wl.nlambda_r, tmp_nlev, 1);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by calloc_aer_out()\n", status);
      return status;
    }
    
    /* copy levels */
    for (lc=0; lc<output->aer.nlev; lc++)
      output->aer.zd[lc] = tmp_zd[lc];
    
    free(tmp_zd);

    for (ipa=0; ipa<output->nipa; ipa++) {
      
      zalt = calloc (output->aer.nlev, sizeof(float));

   
      for (lc=0; lc<output->aer.nlev; lc++){
        /* squeeze aerosol profile when altitude is set to non-zero, enabled using option aerosol_profile_modtran */
        if (input.aer.profile_modtran){
          zalt[lc] =  output->aer.zd[lc];
        }
        /* the Shettle aerosol profile is defined starting       */
        /* at the surface not at z=0; need to subtract altitude  */
        /* this is the default setting                           */
        else
          zalt[lc] = output->aer.zd[lc]- output->alt.altitude;
      }
      
      /* copy visibility to destination and change to a reasonable value if needed */
      status = aerosol_set_visibility (input.aer, input.filename[FN_PATH], &(output->aer.visibility), input.verbose);
      if (status!=0) {
	fprintf (stderr, "Error %d, while setting aerosol visibility in %s (%s)\n", 
                          status, function_name, file_name);
	return status;
      }

      /**** Get appropriate aerosol profile and interpolate profile to zd grid ****/

      nlyr = output->aer.nlev-1;

      /* allocate memory for aerosol property profile */
      aer_pro = (float *) calloc (output->aer.nlev, sizeof (float));

      F77_FUNC (aerprof, ARPROF) (&nlyr, zalt, aer_pro, 
    		         &input.aer.seasn, &input.aer.vulcan, &(output->aer.visibility));
     
      /* squeeze aerosol profile as done in MODTRAN */
      if (input.aer.profile_modtran){
        /* find 6 km altitude level index, below this altitude aerosol profile is sqeezed*/
        lc=0;
        while ( lc < output->aer.nlev){
          if(zalt[lc]<= 6){
            lc_6=lc;
            break;
          }
          lc++;
        }
        
        /* find level index at specified altitude */
        lc=0;
        while ( lc < output->aer.nlev){
          if(zalt[lc]<= output->alt.altitude){
            lc_alt=lc;
            break;
          }
          lc++;
        }
        
        /* Calculate squeezed altitude grid */
        for (lc=lc_6+1; lc<output->aer.nlev; lc++)
          zalt[lc]=zalt[lc-1]-(zalt[lc_6]-output->alt.altitude)/zalt[lc_6]*(output->aer.zd[lc-1]-output->aer.zd[lc]);
        
        /* Interpolate profile on squeezed grid */
        new_z=(float *) calloc(lc_alt+1, sizeof(float)); 
        for (lc=0; lc<=lc_alt; lc++)
          new_z[lc]=output->aer.zd[lc];
            
        new_prof=(float *) calloc(lc_alt+1, sizeof(float)); 
        arb_wvn(output->aer.nlev, zalt, aer_pro, lc_alt, new_z, 
                new_prof, 1, 1); 
        new_prof[lc_alt]=aer_pro[output->aer.nlev-1]; 
        
        for (lc=0; lc<=lc_alt; lc++){
          zalt[lc]= new_z[lc]-output->alt.altitude;
          aer_pro[lc]=new_prof[lc]; 
        }
        for (lc=lc_alt+1; lc<output->aer.nlev; lc++){
          zalt[lc]=output->aer.zd[lc]-output->alt.altitude;
          aer_pro[lc]=0.0; 
        }
      }
      
      
      /**** Get appropriate aerosol optical properties ****/
    
      tmp_aerext = (float *) calloc (maxawvn * output->aer.nlev, sizeof(float));
      tmp_aerabs = (float *) calloc (maxawvn * output->aer.nlev, sizeof(float));
      tmp_aersym = (float *) calloc (maxawvn * output->aer.nlev, sizeof(float));
      tmp_aerwvn = (float *) calloc (maxawvn * output->aer.nlev, sizeof(float));
    

      status = calculate_relative_humidity_on_aerosol_grid (input, output, &(tmp_rh) );
      if (status!=0) {
	fprintf (stderr, "Error %d, calculating relative humidity on aerosol grid in %s (%s)\n", 
                          status, function_name, file_name);
	return status;
      }

      F77_FUNC (aeropt, AEROPT) (&input.aer.vulcan, tmp_aerext, tmp_aerabs, 
	  	        tmp_aersym, tmp_aerwvn, tmp_rh, 
		        &nlyr, zalt, &input.aer.haze, &maxawvn); 
    
      free (tmp_rh);

      tmp_aer_dtau = (float *) calloc (output->wl.nlambda_r * nlyr, sizeof(float));
      tmp_aer_gg   = (float *) calloc (output->wl.nlambda_r * nlyr, sizeof(float));
      tmp_aer_ssa  = (float *) calloc (output->wl.nlambda_r * nlyr, sizeof(float));
    
      F77_FUNC (aerwvn, AERWVN) (tmp_aerwvn, tmp_aerabs, tmp_aersym, tmp_aerext, aer_pro, 
		        tmp_aer_dtau, tmp_aer_gg, tmp_aer_ssa, &nlyr,
		        output->wl.lambda_r, zalt, &output->wl.nlambda_r, &maxawvn);
    
      fortran2c_2D_float_ary_noalloc (output->wl.nlambda_r, nlyr, tmp_aer_gg,   output->aer.optprop.g1);
      fortran2c_2D_float_ary_noalloc (output->wl.nlambda_r, nlyr, tmp_aer_ssa,  output->aer.optprop.ssa);
      fortran2c_2D_float_ary_noalloc (output->wl.nlambda_r, nlyr, tmp_aer_dtau, output->aer.optprop.dtau);

      /* for (lc=0; lc<nlyr+1; lc++) */
      /*         fprintf(stderr, "%d %g g %g ssa %g dtau %g \n", lc, zalt[lc+1], output->aer.optprop.g1[iv][lc],  */
      /*                 output->aer.optprop.ssa[iv][lc], output->aer.optprop.dtau[iv][lc]);  */
      
      
      free(tmp_aer_gg);
      free(tmp_aer_ssa);
      free(tmp_aer_dtau);
    
      free(tmp_aerext);
      free(tmp_aerabs);
      free(tmp_aersym);
      free(tmp_aerwvn);

      free(aer_pro);
      free(zalt);

      
      /* Now we need to unscale the optical thickness, if the user said aerosol_scale_tau,    */
      /* that is, to multiply with 1 / scaling factor. This may sound weird but the reason    */
      /* is that we changed the visibility in aerosol_set_visibility() which already implies  */
      /* the scaling of the aerosol optical thickness. However, at the end of setup_aerosol   */
      /* we have to apply this scaling factor again because the user might e.g. have          */
      /* specified an aerosol_tau_file which should also be scaled. Therefore we unscale here */
      /* in order to be able to scale later without having to consider if the user specified  */
      /* one of many options ... confusing? Yes!                                              */
      
      /**** Unscale aerosol optical thickness ****/
      if (input.aer.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE] > 0.0) /* we don't care about the 0 here because the aerosol is scaled to 0 later anyway */
        for (lc=0; lc<output->aer.nlev-1; lc++)
	  for (iv=0;iv<output->wl.nlambda_r;iv++)
	    output->aer.optprop.dtau[iv][lc] = output->aer.optprop.dtau[iv][lc] / input.aer.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE];


      /**** Read file holding filenames of layer aerosol ****/
      /**** optical properties if specified              ****/

      if (strlen(input.aer.filename[FN_AER_EXPLICIT]) > 0) {
       
        /* determine number of levels; required for memory allocation */
        if ((status = ASCII_checkfile (input.aer.filename[FN_AER_EXPLICIT], 
                                       &rows, &min_columns, &max_columns, &max_length)) != 0) {
	  fprintf (stderr, "Error %d reading %s\n", status, input.filename[FN_AER_EXPLICIT]);
	  return status;
        }
       
        /* allocate memory for temporary array */
        status = calloc_aer_out (tmpaer, output->wl.nlambda_r, rows, 1);
        if (status!=0) {
	  fprintf (stderr, "Error %d allocating memory for aer_out struct\n", status);
	  return status;
        }
      
        /* then read optical properties */
        status = read_optprop_files (input.aer.filename[FN_AER_EXPLICIT],
				     output->wl.lambda_r, output->wl.nlambda_r, 
				     output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
				     tmpaer->zd, tmpaer->nlev, 
				     tmpaer->optprop.dtau, tmpaer->optprop.ssa, tmpaer->optprop.g1,
				     tmpaer->optprop.moment, tmpaer->optprop.nmom, input.verbose);

        if (status != 0) {
	  fprintf (stderr, "Error %d reading aerosol optical property files: %s\n", 
	  	   status, input.aer.filename[FN_AER_EXPLICIT]);
	  return status;
        }
        /* redistribute optical properties to common aerosol grid */
        status = redistribute_optprop (&(tmpaer->optprop), output->wl.nlambda_r, 
				       output->wl.nlambda_rte_lower, 
                                       output->wl.nlambda_rte_upper,
				       tmpaer->zd, tmpaer->nlev-1, 1, 
				       output->aer.zd, output->aer.nlev-1, 1);
  
        /* re-write z profile */
        free (tmpaer->zd);
        tmpaer->nlev = output->aer.nlev;
        tmpaer->zd = calloc (tmpaer->nlev, sizeof(float));
        for (lc=0; lc<output->aer.nlev; lc++)
	  tmpaer->zd[lc] = output->aer.zd[lc];

        if (status!=0) {
	  fprintf (stderr, "Error regridding aerosol parameters\n");
	  return status;
        }
      
        /* copy common aerosol grid to final destination */
        calloc_aer_out (&(output->aer), output->wl.nlambda_r, output->aer.nlev, 1);
        status = cp_aer_out (tmpaer, &(output->aer), output->wl.nlambda_r, 1, input.quiet);
        if (status!=0) {
	  fprintf (stderr, "Error %d copying temporary aerosol properties to\n", status);
	  fprintf (stderr, "output->aer\n");
	  return status;
        }
	status = free_aer_out (tmpaer); 
        output->atm.nmom=0;
        /* increase number of moments if necessary */
        for (iv=0; iv<output->wl.nlambda_r; iv++)
	  for (lc=0; lc<output->aer.nlev-1; lc++)
            if (output->aer.optprop.nmom[iv][lc]-1 > output->atm.nmom)
              output->atm.nmom = output->aer.optprop.nmom[iv][lc]-1;
        
      } /* if (strlen(input.aer.filename[FN_AER_FILES]) > 0) { */

    
      /****   Read aerosol optical depth file if specified ****/
      if (strlen(input.aer.filename[FN_AER_TAU]) > 0) {
        if (!tau_from_netCDF) {
          status = read_aerosol_profile (input.aer.filename[FN_AER_TAU], 
                                         output->aer.nlev-1, output->wl.nlambda_r, output->aer.zd,
				         output->aer.optprop.dtau, 1);
        }
        else {
          /* assuming Koch aerosol data set format */
          /* aerosol optical depth tau (time, z, lat, lon) */
          /* and vertical grid as sigma = p/p_surface coordiante */
	  
          /* add z-grid of the aerosol data set to output->atm.zd */
          status += add_sigma_nc_file2grid 
	    (&tmp_zd, &tmp_nlev, input.aer.filename[FN_AER_TAU],
	     input.latitude, input.longitude, input.UTC, input.atm.time_interpolate,
	     output->atm.zd, output->atm.microphys.press[0][0], output->atm.nlev,
	     input.verbose, input.quiet);
	  
          /* read netCDF tau file */
          status = read_aerosol_profile_from_map 
	    (input.aer.filename[FN_AER_TAU], input.latitude, input.longitude,  
	     input.UTC, input.atm.time_interpolate,
	     output->atm.zd, output->atm.microphys.press[0][0], output->atm.nlev,
	     output->aer.nlev-1, output->wl.nlambda_r, output->aer.zd, 
	     output->aer.optprop.dtau, 1, input.verbose, input.quiet);


          if (status != 0) {
            fprintf (stderr, "Error reading aerosol_tau from netCDF map in %s (%s)\n", function_name, file_name );
            return status;
          }

        }
        if (status != 0) {
	  fprintf (stderr, "Error %d reading aerosol tau file: %s\n", 
	  	   status, input.aer.filename[FN_AER_TAU]);
	  return status;
        }
      }


      /****   Read aerosol asymmetry parameter file if specified ****/
      if (strlen(input.aer.filename[FN_AER_GG]) > 0) {
        status = read_aerosol_profile (input.aer.filename[FN_AER_GG], output->aer.nlev-1,
			  	       output->wl.nlambda_r, output->aer.zd,
				       output->aer.optprop.g1, 0);
        if (status != 0) {
	  fprintf (stderr, "Error %d reading aerosol asymmetry parameter file: %s\n", 
		   status, input.aer.filename[FN_AER_GG]);
	  return status;
        }
      }


      /****   Read aerosol single scattering albedo file if specified ****/
      if (strlen(input.aer.filename[FN_AER_SSA]) > 0) {
        status = read_aerosol_profile (input.aer.filename[FN_AER_SSA], output->aer.nlev-1,
				       output->wl.nlambda_r, output->aer.zd,
				       output->aer.optprop.ssa, 0);
        if (status != 0) {
	  fprintf (stderr, "Error %d reading aerosol single scattering albedo file: %s\n", 
		   status, input.aer.filename[FN_AER_SSA]);
	  return status;
        }
      }
    
    
      /**** Read aerosol moments file if specified ****/
      /**** At present: independent of wavelength  ****/
      /**** but this can be easily changed!        ****/
      
      if (strlen(input.aer.filename[FN_AER_MOMENTS]) > 0) {
      
        status = read_aerosol_moments (input.aer.filename[FN_AER_MOMENTS], 
				       &moment, &rows);
        if (status != 0) {
	  fprintf (stderr, "Error %d reading aerosol moments file: %s\n", 
		   status, input.aer.filename[FN_AER_MOMENTS]);
	  return status;
        }
            
        /* copy to other wavelengths and layers */
        for (iv=0; iv<output->wl.nlambda_r; iv++){
	  for (lc=0; lc<output->aer.nlev-1; lc++){
	    output->aer.optprop.nmom[iv][lc] = rows;
            /* dummy for phase matrix elements */
            output->aer.optprop.moment[iv][lc] = calloc(1, sizeof(float *));
	    output->aer.optprop.moment[iv][lc][0] = 
              calloc (output->aer.optprop.nmom[iv][lc], sizeof (float));

	    for (k=0; k<output->aer.optprop.nmom[iv][lc]; k++){
	      output->aer.optprop.moment[iv][lc][0][k] = moment[k];
            }
          }
        }

        /* increase number of moments if more found */
        for (iv=0; iv<output->wl.nlambda_r; iv++)
	  for (lc=0; lc<output->aer.nlev-1; lc++)
	    if (output->aer.optprop.nmom[iv][lc]-1 > output->atm.nmom)
	      output->atm.nmom = output->aer.optprop.nmom[iv][lc]-1;
      
        free (moment);
      }


      /****   Read aerosol refractive index   ****/
      if (strlen(input.aer.filename[FN_AER_REF]) > 0) {

        /* read wavelength-dependent aerosol refractive index */
        status = read_3c_file_float (input.aer.filename[FN_AER_REF], 
			    	     &lref, &re, &im, &nref);
      
        if (status!=0) {
	  fprintf (stderr, "Error %d reading aerosol refractive index file\n", status);
	  return status;
        }

        /* linear interpolation of the refractive index to the model wavelengths */
        linear = 1;

        output->aer.re_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));
        output->aer.im_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));

        status = arb_wvn (nref, lref, re, 
			  output->wl.nlambda_r, output->wl.lambda_r, output->aer.re_r, 
			  linear, 0);

        if (status != 0) {
	  fprintf (stderr, "Error, status %d returned by arb_wvn()\n", status);
	  return status;
        }

        status = arb_wvn (nref, lref, im, 
			  output->wl.nlambda_r, output->wl.lambda_r, output->aer.im_r, linear, 0);

        if (status != 0) {
	  fprintf (stderr, "Error, status %d returned by arb_wvn()\n", status);
	  return status;
        }
      
        free(lref); free(re); free(im);

        refrac=1;
      }
      else {
        if (input.aer.re>0) {   /* constant refractive index */
	  output->aer.re_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));
	  output->aer.im_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));

	  for (iv=0; iv<output->wl.nlambda_r; iv++) { 
	    output->aer.re_r[iv] = input.aer.re;
	    output->aer.im_r[iv] = input.aer.im;
	  }
	  refrac=1;
        }
      }


      /****   Read aerosol size distribution  ****/
      if (strlen(input.aer.filename[FN_AER_SIZ]) > 0) {

        /* read size distribution */
        status = read_2c_file (input.aer.filename[FN_AER_SIZ], 
			     &(output->aer.xsiz), &(output->aer.ysiz), &(output->aer.nsiz));
      
        if (status!=0) {
	  fprintf (stderr, "Error %d reading aerosol size distribution file\n", status);
	  return status;
        }
      

        if (output->aer.xsiz[0]==0) {
	  output->aer.xsiz[0] = 0.1 * output->aer.xsiz[1];
	  if (!input.quiet)
	    fprintf (stderr, " ... aerosol size distribution: setting first data point to r=%g\n",
		     output->aer.xsiz[0]);
        }

        sizedist=1;
      }

      /* mie calculation of aerosol optical properties */
      if (refrac&&sizedist) {
        status = aerosol_mie_calculation ( input, output );
        if (status != 0) {
          fprintf (stderr, "Error %d during mie calculation in %s (%s)\n", status, function_name, file_name );
          return status;
        }
      }

      if ( strlen(input.aer.filename[FN_AER_SPECIES]) != 0 ) {
       
        /* read aerosol species (mass density profiles from several species) */
        status = read_aerosol_species ( input, output );
        if (status != 0) {
          fprintf (stderr, "Error %d in %s (%s) during read_aerosol_species from %s \n", 
                   status, function_name, file_name, input.aer.filename[FN_AER_SPECIES] );
          return status;
        }
        
        /* convert aerosol mass density profiles in profiles of optical properties */
        status = aerosol_species_to_optical_properties ( input, output );
        if (status != 0) {
          fprintf (stderr, "Error %d in %s (%s) during aerosol_species_to_optical_properties \n", 
                   status, function_name, file_name );
          return status;
        }

      }
      
      /**** Set aerosol asymmetry factor ****/
      if (input.aer.modify[MODIFY_VAR_GG][MODIFY_TYPE_SET] >= -1.0)
        for (lc=0; lc<output->aer.nlev-1; lc++)
	  for (iv=0;iv<output->wl.nlambda_r;iv++)
	    output->aer.optprop.g1[iv][lc] = input.aer.modify[MODIFY_VAR_GG][MODIFY_TYPE_SET];
    
      /**** Set aerosol single scattering albedo ****/
      if (input.aer.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SET] >= 0.0)
        for (lc=0; lc<output->aer.nlev-1; lc++)
	  for (iv=0; iv<output->wl.nlambda_r; iv++) {
	    output->aer.optprop.ssa[iv][lc] = input.aer.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SET];
            if (output->aer.optprop.ssa[iv][lc]>1.0)	output->aer.optprop.ssa[iv][lc] = 1.0;
	  }
    
      /**** Scale aerosol single scattering albedo ****/
      if (input.aer.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SCALE] >= 0.0)
        for (lc=0; lc<output->aer.nlev-1; lc++)
	  for (iv=0;iv<output->wl.nlambda_r;iv++) {
	    output->aer.optprop.ssa[iv][lc] = output->aer.optprop.ssa[iv][lc] * input.aer.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SCALE];
            if (output->aer.optprop.ssa[iv][lc]>1.0)	output->aer.optprop.ssa[iv][lc] = 1.0;
	  }

      /**** Set aerosol optical thickness ****/
      if (input.aer.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SET] >= 0.0) {
        status = scale_aer_tau (output->aer.optprop.dtau, output->wl.lambda_r, output->wl.nlambda_r, 
			        output->aer.nlev-1, input.aer.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SET], 
			        output->aer.zd, output->alt.altitude);

        if (status!=0) {
	  fprintf (stderr, "Error %d returned by scale_aer_tau()\n", status);
	  return status;

        }
      }

      /**** Scale aerosol optical thickness ****/
      if (input.aer.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE] >= 0.0) 
        for (lc=0; lc<output->aer.nlev-1; lc++)
	  for (iv=0;iv<output->wl.nlambda_r;iv++)
	    output->aer.optprop.dtau[iv][lc] = output->aer.optprop.dtau[iv][lc] * input.aer.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE];
	  
	  
      /*for debugging*/
     /* fprintf(stderr,"ang alpha %f, ang beta %f, kb alpha_0 %f  kb alpha_1 %f  kb alpha_2 %f \n ", input.aer.alpha,input.aer.beta, input.aer.alpha_0, input.aer.alpha_1, input.aer.alpha_2); */

      /****test if angstrom and king byrne are not used simultaneaous****/
      if (input.aer.beta >= 0.0 && input.aer.alpha_2>=-800) {
        fprintf (stderr, "Warning: Angstrom and King Byrne parameters are both defined./n Only using King Byrne\n");
      }

      /**** Set aerosol optical thickness according to Angstrom parameters ****/
      if (input.aer.beta >= 0.0 && input.aer.alpha_2<-800) {
        status = scale_aer_angstrom (output->aer.optprop.dtau, output->wl.lambda_r, output->wl.nlambda_r, 
				     output->aer.nlev-1, input.aer.alpha, input.aer.beta,
				     output->aer.zd, output->alt.altitude);

        if (status!=0) {
	  fprintf (stderr, "Error %d returned by scale_aer_angstrom()\n", status);
	  return status;

        }
      }


      /**** Set aerosol optical thickness according to King Byrne parameters ****/
      if (input.aer.alpha_2 >= -800.0) {
        status = scale_aer_king_byrne (output->aer.optprop.dtau, output->wl.lambda_r, output->wl.nlambda_r, 
				     output->aer.nlev-1, input.aer.alpha_0, input.aer.alpha_1,
				     input.aer.alpha_2, output->aer.zd, output->alt.altitude);

        if (status!=0) {
	  fprintf (stderr, "Error %d returned by scale_aer_king_byrne()\n", status);
	  return status;

        }
      }




      /* Set the aerosol optical thickness to a user-defined value at 550nm */
      if (input.aer.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET] >= 0.0 || input.aer.tau_wvl_tau >=0.0) {
        lambda = calloc (output->wl.nlambda_r, sizeof(double));
        tautot = calloc (output->wl.nlambda_r, sizeof(double));

        /* determine level number of 'altitude' level */
        for (ialt=0; ialt<output->aer.nlev; ialt++)
	  if (output->alt.altitude == output->aer.zd[ialt])
	    break;
      
        if (ialt == output->aer.nlev) {
	  fprintf (stderr, "Error, user-defined altitude %f not found in altitude grid\n", output->alt.altitude);
	  return -1;
        }

        for (iv=0; iv<output->wl.nlambda_r; iv++) {
	  lambda[iv] = (double) output->wl.lambda_r[iv];

	  for (lc=0; lc<=ialt-1; lc++)
	    tautot[iv] += output->aer.optprop.dtau[iv][lc];
        }
      
        status = linear_coeffc (lambda, tautot, output->wl.nlambda_r, &a0, &a1, &a2, &a3);
        if (status!=0) {
	  fprintf (stderr, "Error %d calculating interpolation coefficients for aerosol optical thickness.\n", 
	  	   status);
	  return status;
        }
      
	if (input.aer.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET] >=0.0 ) {
	  settau = input.aer.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET];
	  settauwvl = 550.;
	}
	else if (input.aer.tau_wvl_tau>=0.0) {
	  settau = input.aer.tau_wvl_tau;
	  settauwvl = input.aer.tau_wvl_lambda;
	}
        
        status = calc_splined_value (settauwvl, &tauwvl, lambda, output->wl.nlambda_r, a0, a1, a2, a3);
        if (status!=0) {
          fprintf (stderr, "Error %d interpolating aerosol optical thickness to %f nm.\n", status, settauwvl);
          fprintf (stderr, "In the current implementation of aerosol_set_tau_at_wvl, it is required\n");
          fprintf (stderr, "that the internal wavelength range includes %f nm.\n", settauwvl);
	  fprintf (stderr, "The range of the internal wavelength grid is from %f to %f nm. \n", lambda[0], lambda[output->wl.nlambda_r-1]); 
          fprintf (stderr, "Hope to change that in future!\n");
          return status;
        }
        
        if (tauwvl>0.0) {
	  scale_factor = settau / tauwvl;
          
	  for (iv=0; iv<output->wl.nlambda_r; iv++)
	      for (lc=0; lc<output->aer.nlev-1; lc++)
		  output->aer.optprop.dtau[iv][lc] *= scale_factor;
        }
      

        free(a0); free(a1); free(a2); free(a3);
        free(lambda); free(tautot);
      }

    } /* end of the ipa-loop */
  }
  else {
    /* No aerosols included in this run */
    if (input.verbose)
      fprintf (stderr, "     No aerosols included in this run\n");

    status = calloc_aer_out (&(output->aer), output->wl.nlambda_r, 2, 1);
    if (status!=0) {
      fprintf (stderr, "Error %d allocating memory for aer_out-struct\n", status);
      return status;
    }

    output->aer.zd[0] = output->atm.zd[0];
    output->aer.zd[1] = output->atm.zd[output->atm.nlev-1];
    output->aer.nlev = 2;
  }


  /* /\* additional verbose output *\/ */
  /* fprintf(stderr,"   copied aerosol\n"); */
  /* fprintf(stderr,"   --------------\n"); */
  /* iv=0; */

  free(tmpaer);

  /* add the aerosol levels to the common profile */
  set_common_z_grid (output->aer.zd, output->aer.nlev, 
		     &output->atm.zd_common, &output->atm.nlev_common);

#if HAVE_MYSTIC
  /*  Read scaling factors for aerosol importance sampling if file is specified  */
  output->mc.alis.Nc=1;
  if(strlen(input.rte.mc.filename[FN_MC_AERIS])>0){
    if (input.verbose)
      fprintf(stderr, " ... aerosol concentration importance sampling, scaling factors from %s \n ... scaling factors: ",
            input.rte.mc.filename[FN_MC_AERIS]); 
    
    status=read_1c_file(input.rte.mc.filename[FN_MC_AERIS], &(output->mc.alis.aer_scaling_factors), 
                        &(output->mc.alis.Nc)); 
    if (status!=0)
      return err_out ("Error %d returned by read_2c_file() when reading aerosol importance sampling input file\n",
		      status);
    
    if (input.verbose){
      for (ic=0; ic<output->mc.alis.Nc; ic++)
        fprintf(stderr, "%g ", output->mc.alis.aer_scaling_factors[ic]);
      fprintf(stderr, "\n");
    }
  }
#endif

  return 0;
}


/**************************************************************/
/* Read file holding  filenames for each layer.               */
/**************************************************************/

int read_optprop_files (char *filename,
			float *lambda_r, int nlambda_r, 
			int nlambda_rte_lower, int nlambda_rte_upper,
			float *zd, int nlev, 
			float **dtau, float **ssa, float **gg,
			float ****mom, int **nmom, int verbose) 
{
  int i=0, iv=0, k=0, lc=0, ip=0;
  int status=0;
  int rows=0, min_columns=0, max_columns=0, max_length=0;
  int arows=0, amin_columns=0, amax_columns=0;
  int tmp_nlambda=0, tmp_nmom=0;
  int linear=1;
  float *tmp_lambda=NULL;
  float *tmp_dtau=NULL, *tmp_ssa=NULL, *tmp_gg=NULL;
  float **tmp_mom=NULL;
  float **data=NULL;
  float *temp=NULL, **temp_mom=NULL;
  char ***string=NULL;
  char *dummy=NULL;
  int nonzero=0;

  if ((status = ASCII_checkfile (filename, &rows, &min_columns, &max_columns, &max_length)) != 0) 
    return status;

  if (min_columns<2) {
    fprintf (stderr, "Error, need at least 2 columns in optical properties file %s\n", filename);
    return -1;
  }

  if ((status = ASCII_calloc_string(&string, rows, max_columns, FILENAME_MAX+1)) != 0) 
    return status;

  if ((status = ASCII_readfile(filename, string)) != 0) 
    return status;

  if (nlev != rows) {
    fprintf (stderr, "Fatal programming error in read_optprop_files(); nlev=%d, rows=%d\n",
	     nlev, rows);
    return -1;
  }
  

  /* copy altitude grid */
  for (i=0; i<=nlev-1; i++)
    zd[i] = (float) strtod (string[i][0], &dummy);
  
  
  temp = calloc (nlambda_r, sizeof(float));
					   
  /* check if the optical properties file of the              */
  /* uppermost layer has zero optical thickness;              */
  /* we need that one to define the upper layer boundary for  */
  /* the uppermost aerosol/cloud layer; as an alternative     */
  /* NULL as filename simply says that the uppermost layer    */
  /* is empty                                                 */                                               

  if (strcasecmp(string[0][1], "NULL")) {  /* if the uppermost string is not NULL */
    
    if ((status = ASCII_file2float (string[0][1], &arows, &amax_columns,
				    &amin_columns, &data)) != 0) {
      fprintf (stderr, "Error %d reading file %s for layer %d.\n", 
	       status, string[0][1], lc);
      return status;
    }
    
    if (amin_columns<4) {
      fprintf (stderr, "Error, need at least 4 columns in optical properties file %s\n", 
	       string[0][1]);
      return -1;
    }
    
    for (iv=0; iv<arows; iv++)
      if (data[iv][1] > 0)  /* optical thickness is first column */
	nonzero++;
    
    if (nonzero) {
      fprintf (stderr, "Error, need a layer with zero optical thickness as uppermost\n");
      fprintf (stderr, "layer in %s. Please add one and retry!\n", string[0][1]);
      fprintf (stderr, "See examples/AERO_FILES as an example!\n");
      return -1;
    }
    
    ASCII_free_float (data, arows);
  }
  else {
    fprintf (stderr, "  ... found NULL as uppermost layer, fine :-)\n");
  }

  /* For each layer read the optical properties file.                  */
  /* Note that the wavelength resolution may be different from the one */
  /* in the extraterrestrial solar flux file.                          */
  for (lc=0; lc<nlev-1; lc++) {

    if (verbose)
      fprintf (stderr, " ... reading optical properties %s for layer [%.4f, %.4f]\n",
	       string[lc+1][1], zd[lc+1], zd[lc]);

    if ((status = ASCII_file2float (string[lc+1][1], &arows, &amax_columns,
				    &amin_columns, &data)) != 0) {
      fprintf (stderr, "Error %d reading file %s for layer %d.\n", 
	       status, string[lc+1][1], lc);
      return status;
    }
    
    if (amin_columns<4) {
      fprintf (stderr, "Error, need at least 4 columns in optical properties file %s\n", 
	       string[lc+1][1]);
      return -1;
    }
    
    tmp_nmom = amax_columns-3;
    tmp_nlambda = arows;
      
    tmp_lambda = calloc (tmp_nlambda, sizeof(float));
    tmp_dtau   = calloc (tmp_nlambda, sizeof(float));
    tmp_gg     = calloc (tmp_nlambda, sizeof(float));
    tmp_ssa    = calloc (tmp_nlambda, sizeof(float));
      
    if ((ASCII_calloc_float (&tmp_mom,  tmp_nmom, tmp_nlambda)) != 0)
      return status;

    for (iv=0; iv<tmp_nlambda; iv++)
      tmp_lambda[iv] = data[iv][0];

    for (iv=0; iv<tmp_nlambda; iv++) {
      tmp_dtau [iv] = data[iv][1];	 /* Still in units of km-1      */
      tmp_ssa  [iv] = data[iv][2];	
      for (k=0; k<tmp_nmom; k++)
	if (data[iv][k+3]==data[iv][k+3]) /* Trick to replace NaN by 0  */
	  tmp_mom[k][iv] = data[iv][k+3];
	else
	  tmp_mom[k][iv] = 0;

      if (tmp_nmom>=2)
	tmp_gg [iv] = tmp_mom[1][iv];
      else 
	tmp_gg [iv] = 0;
    }

    ASCII_free_float (data, arows);

    for (iv=0; iv<tmp_nlambda; iv++)
      tmp_dtau[iv] *= (zd[lc]-zd[lc+1]);
    
    /* Interpolate from optical property wavelengths to */
    /* the internal wavelength grid.                    */
    
    /* interpolate optical thickness */
    status = arb_wvn2 (tmp_nlambda, tmp_lambda, tmp_dtau,  
		       nlambda_rte_lower, nlambda_rte_upper,
		       lambda_r, temp, linear);

    if (status!=0) {
      fprintf (stderr, "Error %d interpolating dtau in read_optprop_files()\n", status);
      return status;
    }

    for (iv=0; iv<nlambda_r; iv++)
      dtau[iv][lc] = temp[iv];

    /* interpolate asymmetry parameter */
    status = arb_wvn2 (tmp_nlambda, tmp_lambda, tmp_gg,  
		       nlambda_rte_lower, nlambda_rte_upper,
		       lambda_r, temp, linear);

    if (status!=0) {
      fprintf (stderr, "Error %d interpolating gg in read_optprop_files()\n", status);
      return status;
    }

    for (iv=0; iv<nlambda_r; iv++) 
      gg[iv][lc] = temp[iv];
    

    /* interpolate single scattering albedo */
    status = arb_wvn2 (tmp_nlambda, tmp_lambda, tmp_ssa,
		       nlambda_rte_lower, nlambda_rte_upper,
		       lambda_r, temp, linear);
    
    if (status!=0) {
      fprintf (stderr, "Error %d interpolating ssa in read_optprop_files()\n", status);
      return status;
    }
    
    for (iv=0; iv<nlambda_r; iv++) 
      ssa[iv][lc] = temp[iv];

    

    ASCII_calloc_float (&temp_mom, tmp_nmom, nlambda_r);

    for (k=0; k<tmp_nmom; k++) {
      status = arb_wvn2 (tmp_nlambda, tmp_lambda, tmp_mom[k],
			 nlambda_rte_lower, nlambda_rte_upper,
			 lambda_r, temp_mom[k], linear);
      if (status!=0) {
	fprintf (stderr, "Error %d interpolating mom in read_optprop_files()\n", status);
	return status;
      }
    }

    /* get index of maximum nonzero moment for each wavelength */
    for (iv=0; iv<nlambda_r; iv++) {
      for (k=tmp_nmom-1; k>=0; k--)  
	if (temp_mom[k][iv]!=0)
	  break;
      
      nmom[iv][lc] = k+1;
      
      mom[iv][lc] = calloc (1, sizeof(float *));
      /* Fix, if Legendre polynomials should be read not only for the 
         phase function */
      for (ip=0; ip<1; ip++)
        mom[iv][lc][ip] = calloc (nmom[iv][lc], sizeof(float));
      
      for (k=0; k<nmom[iv][lc]; k++){
	mom[iv][lc][0][k] = temp_mom[k][iv];
      }
    }
    
    /* free memory */
    ASCII_free_float (temp_mom, tmp_nmom);
    ASCII_free_float (tmp_mom,  tmp_nmom);

    free (tmp_dtau);
    free (tmp_gg);
    free (tmp_ssa);
    free (tmp_lambda);
  }
  
  ASCII_free_string (string, rows, max_columns);
  free (temp);
  
  return 0;
}



/***********************************************************************************/
/* Function: calculate_relative_humidity_on_aerosol_grid                           */
/* Description:                                                                    */
/*  calculate relative humidity on aerosol grid                                    */
/*  interpolation is done for interpol_method_gas[MOL_H2O],                        */
/*  which is per default: linear interpolation of the mixing ratio                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: status, status<0 means not OK, status == 0 means OK               */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: Ulrich Hamann, 12.2007                                                  */
/*                                                                                 */
/***********************************************************************************/

static int calculate_relative_humidity_on_aerosol_grid (input_struct input, output_struct *output, float **tmp_rh)
{

  int status=0;

  int lc=0;

  float *tmp_air    = NULL;
  float *tmp_h2o    = NULL;
  float *tmp_temper = NULL;

  (*tmp_rh)     = calloc (output->aer.nlev, sizeof (float));

  /* interpolate air number density on aerosol grid */
  tmp_air    = calloc (output->aer.nlev, sizeof (float));
  status = arb_wvn (output->atm.nlev, output->atm.zd, output->atm.microphys.dens[MOL_AIR][0][0],
		    output->aer.nlev, output->aer.zd, tmp_air, input.atm.interpol_method_gas[MOL_AIR], 1);
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating air density to aerosol grid\n", status);
    return status;
  }

  /* calculate water vapour number density on aerosol grid */
  tmp_h2o    = calloc (output->atm.nlev, sizeof (float));  
  for (lc=0; lc<output->atm.nlev; lc++)
    tmp_h2o[lc] = output->atm.microphys.dens[MOL_H2O][0][0][lc]; 

  status += interpolate_density (output->atm.zd, &(tmp_h2o), output->atm.nlev, 
                                 output->aer.zd, output->aer.nlev, input.atm.interpol_method_gas[MOL_H2O],
                                 output->atm.microphys.dens[MOL_AIR][0][0], tmp_air, input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating water vapour to aerosol grid (in aerosol.c)\n", status);
    return status;
  }

  tmp_temper = calloc (output->aer.nlev, sizeof (float));

  status = arb_wvn (output->atm.nlev, output->atm.zd, output->atm.microphys.temper[0][0],
                    output->aer.nlev, output->aer.zd, tmp_temper, input.atm.interpol_method_temper, 1);
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating temperature to aerosol grid\n", status);
    return status;
  }

  for (lc=0; lc<output->aer.nlev; lc++)
    (*tmp_rh)[lc] = tmp_h2o[lc] / vapor_pressure (tmp_temper[lc]) * 100.0;

  free(tmp_air);
  free(tmp_h2o);
  free(tmp_temper);

  return status;

}

/**************************************************************/
/* Read aerosol profile.                                      */
/**************************************************************/

static int read_aerosol_profile (char *filename, int nlyr, int nlambda_r, 
				 float *zd_common, float **dtau, int scale) 
{
  int rows=0;
  int i=0, j=0, lc=0;
  int status = 0;
  float *aero_zd=NULL, *aero_tau=NULL, *aero_tau_common=NULL;


  /* read altitude and optical depth from aerosol file */
  status = read_2c_file_float (filename, 
			       &aero_zd, &aero_tau, &rows);

  if (status!=0) {
    fprintf (stderr, "Error %d reading aerosol file %s\n", status, filename);
    return status;
  }

  aero_tau_common = calloc (nlyr+1, sizeof (float));
  
  status = regrid_float (aero_zd, aero_tau,  rows,
			 zd_common, aero_tau_common, nlyr+1, scale);
  
  if (status!=0) {
    fprintf (stderr, "Error regridding %s to common aerosol grid\n", filename);
    for (lc=0; lc<nlyr+1; lc++)
      fprintf (stderr, "%d %f\n", lc, aero_tau_common[lc]);
    

    return status;
  }

  /* copy optical depth to final destination dtau[][] */
  for (i=0; i<nlambda_r; i++)
    for (j=0; j<nlyr; j++)
      dtau[i][j] = aero_tau_common[j+1];

  /* free memory */
  free (aero_zd); free (aero_tau);
  free (aero_tau_common);

  return 0;
}

/**************************************************************/
/* Read aerosol profile.                                      */
/**************************************************************/

static int read_aerosol_profile_from_map (char *filename, float latitude, float longitude, struct tm UTC, int time_interpolate, 
                                          float *z_atm, float *press, int nlev_atm,
                                          int nlyr, int nlambda_r, float *zd_common, 
                                          float **dtau, int scale, int verbose, int quiet) 
{
  int ncid=NOT_DEFINED_INTEGER;

  long ilat=NOT_DEFINED_INTEGER,ilon=NOT_DEFINED_INTEGER,itime=NOT_DEFINED_INTEGER;  /* index for lat, lon, time in netCDF file */
  size_t nlat  = 0;
  size_t nlon  = 0;
  double *lat_grid  = NULL;
  double *lon_grid  = NULL;
  int nt=-1;
  int itime1 = -1, itime2 = -1;
  float dt = NOT_DEFINED_FLOAT;

  size_t nnew=0;              /* number of layers in the netCDF file */
  double *sigma_grid=NULL;
  float *log_p_aero=NULL;
  float *log_p_atm =NULL;

  float *aero_zd=NULL;
  float tmp=NOT_DEFINED_FLOAT;

  int t =NOT_DEFINED_INTEGER;

  float **tau=NULL;

  size_t mlev=0;              /* number of layers in the netCDF file */
  double *mlev_grid=NULL;

  int i=NOT_DEFINED_INTEGER, j=NOT_DEFINED_INTEGER, lc=NOT_DEFINED_INTEGER;
  int status = 0;
  float *aero_tau=NULL, *aero_tau_common=NULL;

  char function_name[] = "read_aerosol_profile_from_map";
  char file_name[]     = "aerosol.c";

  status = get_all_netCDF_indices ( filename, latitude, longitude,
                                    &(ncid), &(ilat), &(nlat), &(ilon), &(nlon), 
                                    &(lat_grid), &(lon_grid),
                                    UTC, time_interpolate,
                                    &(nt), &(itime1), &(itime2), &(dt),
                                    verbose, quiet);

  /* read number of layers */
  status = alloc_and_read_netCDF_1D_double(ncid,"ilev", &nnew, "sigma", &(sigma_grid));
  if (status!=0) {
    fprintf (stderr, "Error %d reading 'ilev' from %s in %s (%s)\n", status, filename, function_name, file_name );
    return status;
  }

  if ((log_p_aero = calloc (nnew, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'log_p_aero' in %s (%s)\n", function_name, file_name);
    return -10;
  }

  for (lc=0; lc<nnew; lc++) {
    /* convert from sigma(=p/p0) -> log (sigma) */
    log_p_aero[lc] = log(sigma_grid[lc]);
  }

  free(sigma_grid);

  /* alloc one or two timesteps for pressure */
  if ((log_p_atm = calloc (nlev_atm, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'log_p_atm' in %s (%s)\n", function_name, file_name);
    return -10;
  }

  for (lc=0; lc<nlev_atm; lc++) {
    /* convert from sigma(=p/p0) -> log (sigma) */
    log_p_atm[lc] = log(press[lc]/press[nlev_atm-1]);
  }

  if ((aero_zd = calloc (nnew, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'aero_zd' in %s (%s)\n", function_name, file_name);
    return -10;
  }

  /* linear interpolation of z on the log(p/p0)-grid */
  status += arb_wvn (nlev_atm, log_p_atm,  z_atm, 
                     (int) nnew, log_p_aero, aero_zd, 
                     INTERP_METHOD_LINEAR, 0);
  if (status!=0) {
    fprintf (stderr,"Error, during interpolation of 'z' in %s (%s)\n", function_name, file_name);
    return status;
  }

  free(log_p_aero);
  free(log_p_atm);

  /* reverse order of aero_zd */
  for (lc=0; lc<nnew/2; lc++) {
    tmp             = aero_zd[nnew-1-lc];
    aero_zd[nnew-1-lc] = aero_zd[lc];
    aero_zd[lc]        = tmp;
  }

  /* correct small rounding errors for lowermost level */
  if (fabs(aero_zd[nnew-1]-zd_common[nlyr-1]) < 0.00000001 )
    aero_zd[nnew-1]=zd_common[nlyr-1];

  /* read number of layers */
  status = alloc_and_read_netCDF_1D_double(ncid,"mlev", &mlev, "mlev",  &(mlev_grid));
  if (status!=0) {
    fprintf (stderr, "Error %d reading 'mlev' from %s in %s (%s)\n", status, filename, function_name, file_name );
    return status;
  }

  if ((tau = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'tau' in %s (%s)\n", function_name, file_name);
    return -10;
  }

  /* read nt (= 1 or 2) time steps */
  for (t=0;t<=nt-1;t++) {

    if (t == 0)
      itime = itime1;
    if (t == 1)
      itime = itime2;

    /* allocate and read temperature, defined on layers */
    status = alloc_and_read_netCDF_column_float ( ncid, "tau", &(tau[t]), mlev, itime, ilat, ilon, TRUE );
    if (status != 0) {
      fprintf (stderr, "Error %d reading temperature 'tau' from netCDF file\n", status);
      return status;
    }
  }


  /* alloc aero_tau */
  if ((aero_tau = calloc (mlev+1, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for pressure in read_ECHAM_atmosphere() (in atmosphere.c)\n");
    return -10;
  }

  if (nt == 1) {
    for (lc=0; lc<mlev; lc++) {
      aero_tau[mlev-lc] = tau[0][lc];
    }
  }
  else {
    for (lc=0; lc<mlev; lc++) {
      aero_tau[mlev-lc] = (1.0-dt)* tau[0][lc]  +  dt* tau[1][lc];
    }
  }

  if (verbose) {
    fprintf (stderr," lc    z_aero       tau_aero \n");
    fprintf (stderr,"------------------------------\n");
    for (lc=0; lc<mlev+1; lc++)
      fprintf (stderr,"%3d  %11.7f  %14.7e\n", lc, aero_zd[lc], aero_tau[lc]);
  }

  aero_tau_common = calloc (nlyr+1, sizeof (float));  /* number of layer of the common grid */
  
  status = regrid_float (aero_zd, aero_tau,  mlev,
			 zd_common, aero_tau_common, nlyr+1, scale);
  
  if (status!=0) {
    fprintf (stderr, "Error regridding %s to common aerosol grid\n", filename);
    for (lc=0; lc<nlyr+1; lc++)
      fprintf (stderr, "%3d %13.6e\n", lc, aero_tau_common[lc]);
    return status;
  }

  /* copy optical depth to final destination dtau[][] */
  for (i=0; i<nlambda_r; i++)
    for (j=0; j<nlyr; j++)
      dtau[i][j] = aero_tau_common[j+1];

  /* free memory */
  free (aero_zd); free (aero_tau);
  free (aero_tau_common);

  return 0;
}



/**************************************************************/
/* Read moments of the aerosol phase function; memory for     */
/* float *moment is allocated automatically.                  */
/**************************************************************/

static int read_aerosol_moments (char *filename, 
	                         float **moment, int *nmom) 
{
  int i=0;
  int status = 0;
  double *aero_mom=NULL;
  
  /* read moments of the aerosol phase function from file */
  status = read_1c_file (filename, &aero_mom, nmom);
  
  if (status!=0) {
    switch(status)  {
    case ASCIIFILE_NOT_FOUND: 
      fprintf (stderr, "Aerosol tau file %s not found\n", filename);
      return status;
      break;
    case ASCII_NO_MEMORY: 
      fprintf (stderr, "Not enough memory for reading aerosol moments file %s\n",
	       filename);
      return status;
      break;
    default:
      fprintf (stderr, "Error %d reading aerosol file %s\n",  
	       status, filename);
      return status;
      break;
    }
  }


  /* normalize and copy moments to final destination */

  if (aero_mom[0]<=0) {
    fprintf (stderr, "Error, first moment cannot be zero or negative\n");
    return -1;
  }

  *moment = calloc (*nmom, sizeof(float ));
  
  for (i=0;i<*nmom;i++){
    (*moment)[i] = aero_mom[i] / aero_mom[0];
  }

  /* free memory */
  free (aero_mom);

  return 0;
}




/**************************************************************/
/* Scale the total aerosol optical thickness with Angstrom    */
/* parameters alpha and beta.                                 */
/**************************************************************/

static int scale_aer_angstrom (float **dtau, float *lambda_r, int nlambda_r, int nlyr,
			       float alpha, float beta, float *zd, float altitude)
{
  int lc=0, iv=0, ialt=0;
  float fact=0, depth=0, ext=0;
  
  /* determine level number of 'altitude' level */
  for (ialt=0; ialt<=nlyr; ialt++)
    if (altitude == zd[ialt])
      break;
  
  if (ialt == nlyr+1) {
    fprintf (stderr, "Error, user-defined altitude %f not found in altitude grid\n", altitude);
    return -1;
  }

  /* Function Body */
  for (iv=0; iv<nlambda_r; iv++) {
    
    ext = 0.0;
    for (lc=0; lc<=ialt-1; lc++)
      ext += dtau[iv][lc];
    
    depth = beta * pow(lambda_r[iv]/1000.0, -alpha);
    fact = depth / ext;
   

    for (lc=0; lc<nlyr; lc++)
      dtau[iv][lc] *= fact;
  }
  
  return 0;
} 



/***********************************************************************************/
/*Function: scale_aer_king_byrne                                                   */
/* Description: Scale the total aerosol optical thickness with the King-Byrne      */ 
/*               parameters alpha_0,  alpha_1 and alpha_2.                         */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: status                                                            */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: Jethro Betcke                                                           */  
/* Date:   Sept 2014                                                               */
/* Based on: scale_aer_angstrom  by anonymous   only one line, and declaration     */
/*           changed                                                               */
/* Note:     It would be more elegant to combine this function with                */
/*           scale_aer_angstrom since the equation only differs by one factor      */
/*           However, this could easily lead to errors due to the sign difference  */
/*           between alpha in the Angstrom Equation and alpha_1 in the King-Byrne  */
/*           equation.                                                             */             
/***********************************************************************************/

static int scale_aer_king_byrne (float **dtau, float *lambda_r, int nlambda_r, int nlyr,
				 float alpha_0,float alpha_1, float alpha_2, float *zd, float altitude)
{
  int lc=0, iv=0, ialt=0;
  float fact=0, depth=0, ext=0;
  
  /* determine level number of 'altitude' level */
  for (ialt=0; ialt<=nlyr; ialt++)
    if (altitude == zd[ialt])
      break;
  
  if (ialt == nlyr+1) {
    fprintf (stderr, "Error, user-defined altitude %f not found in altitude grid\n", altitude);
    return -1;
  }


  /* Function Body */
  for (iv=0; iv<nlambda_r; iv++) {
    
    ext = 0.0;
    for (lc=0; lc<=ialt-1; lc++)
      ext += dtau[iv][lc];
    
    depth = exp(alpha_0) * pow(lambda_r[iv]/1000.0, alpha_1) * pow(lambda_r[iv]/1000.0, alpha_2*log(lambda_r[iv]/1000.0));
    fact = depth / ext;
    
    for (lc=0; lc<nlyr; lc++)
      dtau[iv][lc] *= fact;
  }
  
  return 0;
} 









/***********************************************************************************/
/* Function: aerosol_mie_calculation                                               */
/* Description:                                                                    */
/*  preparation of the mie calculation for sizedistribution of aerosols            */
/*  using complex refraction index and mie-theory                                  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: status                                                            */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: unknown                                                                 */
/*                                                                                 */
/***********************************************************************************/

static int aerosol_mie_calculation ( input_struct input, output_struct *output )
{

  int status=0;
  int iv=0, ip=0, lc=0, k=0;

  int nmom=1000;   /* ??? should tie this to nstr to save computational time ??? */

  mie_inp_struct *mie_inp=NULL;
  mie_out_struct *mie_out=NULL;

  mie_complex crefin = {0.0, 0.0}, ref = {0.0, 0.0};

  double temperature = 300.0;  /* The temperature is only used for wavelengths > 167 um */
  double beta=0, beta0=0, ext=0, omega=0, g=0;

  int kmax=0;

  if (!input.quiet)
    fprintf (stderr, " ... now doing aerosol Mie calculation\n");

  mie_inp = calloc (1, sizeof(mie_inp_struct));
  mie_out = calloc (1, sizeof(mie_out_struct));

  mie_out->pmom = calloc (4, sizeof(float *));
  for (ip=0; ip<4; ip++)
    mie_out->pmom[ip] = calloc (nmom+1, sizeof(float));
      

  /**********************************************************/
  /* The wavelength dependence of the optical thickness     */
  /* follows Mie theory; the absolute value of the optical  */
  /* thickness is defined by its value at this place of the */
  /* code; that is, the value defined by 'visibility' or    */
  /* (much better) by aerosol_set_tau. Later, it may still  */
  /* be scaled with the Angstrom parameters (if anybody     */
  /* wants that.			*/
  /* ??? I guess there are much better ways to do that ???  */
  /* yes, King Byrne!                                       */
  /**********************************************************/

  for (iv=0; iv<output->wl.nlambda_r; iv++) {

    mie_inp->anyang      = 0;
    mie_inp->ipolzn      = 0;
    mie_inp->nmom        = nmom;
    mie_inp->momdim      = nmom+1;
    mie_inp->numang      = 0;
    mie_inp->mimcut      = 1.0e-08;
    mie_inp->perfct      = 0;
    mie_inp->prnt[0]     = 0;
    mie_inp->prnt[1]     = 0;
    mie_inp->xmu         = (float *)   calloc (mie_inp->numang+1, sizeof(float));
    mie_inp->s1          = (mie_complex *) calloc (mie_inp->numang+1, sizeof(mie_complex));
    mie_inp->s2          = (mie_complex *) calloc (mie_inp->numang+1, sizeof(mie_complex));


    crefin.re = output->aer.re_r[iv];
    crefin.im = output->aer.im_r[iv];

    status = mie_calc_sizedist (*mie_inp, mie_out, 
                     	        MIEV0, USER, crefin, 
        			(output->wl.lambda_r[iv])/1000.0, temperature, 1,
        			output->aer.xsiz, output->aer.ysiz, output->aer.nsiz, 
        			&beta, &omega, &g, &ref);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by mie_calc_sizedist()\n", status);
      return status;
    }

    /* determine maximum non-zero moment */
    kmax=0;
    for (k=1; k<=nmom; k++)
      if (mie_out->pmom[0][k] != 0)
        kmax=k;
	  
    for (lc=0; lc<output->aer.nlev-1; lc++)
      output->aer.optprop.nmom[iv][lc] = kmax+1;

    /* copy moments to final destination */
    for (lc=0; lc<output->aer.nlev-1; lc++) {
      output->aer.optprop.moment[iv][lc]= calloc (1, sizeof(float *));
      for (ip=0; ip<1; ip++){
        output->aer.optprop.moment[iv][lc][ip] = calloc (output->aer.optprop.nmom[iv][lc], sizeof(float));
        for (k=0; k<kmax; k++)
          /* FIXME: Normalization should be done only for phase function */
          output->aer.optprop.moment[iv][lc][ip][k] = mie_out->pmom[ip][k]/mie_out->pmom[ip][0];
      }
    }
	

    /* extinction at first wavelength; required for scaling the optical thickness */
    if (iv==0)
      beta0 = beta;


    ext = 0.0;
	  
    /* overwrite optical thickness, single-scattering albedo, and asymmetry parameter */
    for (lc=0; lc<output->aer.nlev-1; lc++) {
      output->aer.optprop.dtau [iv][lc] *= beta/beta0;
      output->aer.optprop.ssa  [iv][lc]  = omega;
      output->aer.optprop.g1   [iv][lc]  = g;
	  
      ext += output->aer.optprop.dtau[iv][lc];
    }
	

    if (!input.quiet) {
      if (iv==0)
	fprintf (stderr, "mie  lambda  nmom       tau       ssa       asy\n");
	  
      fprintf (stderr, "mie %7.2f  %4d  %g  %g  %g\n",
	       output->wl.lambda_r[iv], kmax, ext, omega, g);
    }
  }

  /* free Mie memory */
  for (ip=0; ip<4; ip++)
    free (mie_out->pmom[ip]);
  free(mie_out->pmom);

  free (mie_inp);
  free (mie_out);

  return status;

}

/***********************************************************************************/
/* Function: read_aerosol_species                                                  */
/* Description:                                                                    */
/*  read aerosol mass densities profiles for several optical species from netCDF   */
/*  map. Then look into aerosol_library, where for each species a netCDF file      */
/*  with optical properties is expected. Adept optical properties to relative      */
/*  humidity (if nessesary) and calculate total optical properties.                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: Ulrich Hamann, 17.12.2007                                               */
/*                                                                                 */
/***********************************************************************************/

static int read_aerosol_species ( input_struct input, output_struct *output )
{

  int status=0;
  char aerosol_species_name_list[FILENAME_MAX]="";
  char ***aerosol_species=NULL;

  int i_aer=NOT_DEFINED_INTEGER;
  int max_columns=NOT_DEFINED_INTEGER, min_columns=NOT_DEFINED_INTEGER;

  char momfilename[FILENAME_MAX]="";
  glob_t momfiles;

  char aerosol_species_file[FILENAME_MAX]="";
  int species_from_netCDF = FALSE;

  int lev=0;

  char function_name[] = "read_aerosol_species";
  char file_name[]     = "aerosol.c";

#if HAVE_LIBNETCDF
  int ncid=0;
#endif

  if (( output->aer.aerosol_library = (char *) calloc (FILENAME_MAX+1, sizeof (char))) == NULL) {
    fprintf (stderr,"Error, allocating memory for aerosol_library in %s (%s)\n", function_name, file_name);
    return -2;
  }
  
  /* replace library keywords with data folders */
  strcpy(output->aer.aerosol_library, input.aer.filename[FN_AER_SPECIES_LIB]);
  if (strcmp("OPAC",output->aer.aerosol_library) == 0 || strcmp("opac",output->aer.aerosol_library) == 0) {

    strcpy(output->aer.aerosol_library, input.filename[FN_PATH]);
    strcat(output->aer.aerosol_library,"aerosol/OPAC/optprop/");

    /* directory 'aerosol/OPAC/mie' was renamed to 'aerosol/OPAC/optprop/' in June 2012 */
    /* check whether netcdf files are available at new location */
    /* if not, try to find netcdf files at old location */
    strcpy (momfilename, output->aer.aerosol_library);
    strcat (momfilename, "*.*.cdf");
    status = glob (momfilename, 0, NULL, &momfiles);
    if ( momfiles.gl_pathc == 0 ) {
      strcpy (momfilename, output->aer.aerosol_library);
      strcat (momfilename, "../mie/*.*.cdf");
      status = glob (momfilename, 0, NULL, &momfiles);
      if ( momfiles.gl_pathc > 0 ) {
	if ( !input.quiet ) {
	  fprintf (stderr, "*** Warning, the location of the OPAC optical properties was changed\n");
	  fprintf (stderr, "*** to aerosol/OPAC/optprop/ in June 2012. On this machine it is still in\n");
	  fprintf (stderr, "*** the old location aerosol/OPAC/mie/ - this might not work anymore in\n");
	  fprintf (stderr, "*** future.\n");
	}
	strcpy(output->aer.aerosol_library, input.filename[FN_PATH]);
	strcat(output->aer.aerosol_library,"aerosol/OPAC/mie/");
      }
    }

  }

  if (input.verbose) 
    fprintf (stderr, " ... calculating aerosol optical properties using library: %s\n", output->aer.aerosol_library);
  
  if ( input.aer.n_species != NOT_DEFINED_INTEGER ) {
    
    /* aerosol species defined in the input, copy input structure to output structure */
    output->aer.n_species = input.aer.n_species;
    
    /* allocating memory for aerosol species names */
    status = ASCII_calloc_char(&output->aer.species_names,output->aer.n_species,FILENAME_MAX+1);
    if(status!=0){
      fprintf (stderr,"Error, allocating memory for aer.species_names in %s (%s)\n", function_name, file_name);
      return status;
    }
    
    /* copy to final place */
    for (i_aer=0;i_aer<output->aer.n_species;i_aer++) {
      strcpy ( output->aer.species_names[i_aer], input.aer.species_names[i_aer] );
    }
  }
  else {
    
    /* read aerosol species names from */
    strcpy (aerosol_species_name_list, output->aer.aerosol_library);
    strcat (aerosol_species_name_list, "aerosol_species_names.dat");
    
    if (input.verbose) 
      fprintf (stderr, "     reading aerosol species name list: %s \n", aerosol_species_name_list);
    
    status = ASCII_file2string (aerosol_species_name_list, &(output->aer.n_species), &(max_columns), &(min_columns), &(aerosol_species));
    if (status!=0) {
      fprintf (stderr, "Error %d in %s (%s) reading %s\n", status, function_name, file_name, aerosol_species_name_list);
      return status;
    }
    
    /* allocating memory for aerosol species names */
    status = ASCII_calloc_char(&output->aer.species_names,output->aer.n_species,FILENAME_MAX+1);
    if(status!=0){
      fprintf (stderr,"Error, allocating memory for aer.species_names in %s (%s)\n", function_name, file_name);
      return status;
    }
    
    /* copy to final place */
    for (i_aer=0;i_aer<output->aer.n_species;i_aer++)
      strcpy (output->aer.species_names[i_aer], aerosol_species[i_aer][0]);
    
    ASCII_free_string(aerosol_species, output->aer.n_species, max_columns);
    
  }

  if (input.verbose) 
    fprintf (stderr, "     library contains %4d aerosol species \n", output->aer.n_species);
  
  strcpy(aerosol_species_file, input.aer.filename[FN_AER_SPECIES]);
  
#if HAVE_LIBNETCDF
  
  /* test if input file is netCDF format */
  if (strlen(aerosol_species_file) > 0) {
    
    /* open netcdf file */
    status = nc_open (aerosol_species_file, NC_NOWRITE, &ncid);
    if (status==NC_NOERR) {
      species_from_netCDF = TRUE;
      nc_close(ncid); 
    }
  }

#endif
  

  if (species_from_netCDF) {
    /* read aerosol species from netCDF file */
    status = read_aerosol_species_from_netCDF_map ( input, output, aerosol_species_file );
  }
  else {
    /* read aerosol species from ASCII file */
    status = read_aerosol_species_from_ASCII_file ( aerosol_species_file, &(output->aer), 
                                                    input.quiet, input.verbose); 
  }
  if (status != 0) {
    fprintf (stderr, "Error %d while reading aerosol species (line %d, function %s in %s)\n", 
             status, __LINE__, __func__, __FILE__ );
    return status;
  }
  
  /* need to shift by one layer as OPAC aerosols are automatically considered layer properties */
  for ( i_aer=0; i_aer<output->aer.n_species; i_aer++ ){
    for (lev=0; lev<output->aer.nlev-1; lev++) {
      output->aer.massdens[i_aer][lev] = output->aer.massdens[i_aer][lev+1];
    }
  }

  /* Arve Kylling 20120516: aerosol_species_file did not work with altitude option.  */
  /* Solve this by shifting aerosol profile upwards by output->alt.altitude. This is */
  /* consistent with the documentation for the altitude option.                      */
  if ( output->alt.altitude ) {
    for (lev=0; lev<output->aer.nlev; lev++) {
      output->aer.zd[lev] = output->aer.zd[lev] + output->alt.altitude;
    }
  }


  /* verbose output */
  if (input.verbose) {
    /* fist line: physical properties */
    fprintf (stderr, "\n#    z         ");
    for ( i_aer=0; i_aer<output->aer.n_species; i_aer++ )
      fprintf (stderr, " %-9s   ", output->aer.species_names[i_aer] );
    fprintf (stderr, "\n");
    /* second line: units */
    fprintf (stderr, "#  [km]     ");
    for ( i_aer=0; i_aer<output->aer.n_species; i_aer++ )
      fprintf (stderr, "   [g/m-3]   ");
    fprintf (stderr, "\n");
    /* third line: grid */
    fprintf (stderr, "#-----------");
    for ( i_aer=0; i_aer<output->aer.n_species; i_aer++ )
      fprintf (stderr, "-------------");
    fprintf (stderr, "\n");
    /* following lines: data values */
    for (lev=0; lev<output->aer.nlev; lev++) {
      fprintf (stderr, " %11.6f  ", output->aer.zd[lev]);
      for (i_aer=0; i_aer<output->aer.n_species; i_aer++)
        fprintf (stderr, "%.5e  ", output->aer.massdens[i_aer][lev]);
      fprintf (stderr, "\n");
    }
    /* last 2 lines: grid */
    fprintf (stderr, "#-----------");
    for ( i_aer=0; i_aer<output->aer.n_species; i_aer++ )
      fprintf (stderr, "-------------");
    fprintf (stderr, "\n\n");
  }

  return status;

}

/***********************************************************************************/
/* Function: read_aerosol_species_from_netCDF_map                                  */
/* Description:                                                                    */
/*  read aerosol mass densities profiles from an netCDF file                       */
/*  from position input.latitude, input.longitude                                  */
/*  at time input.UTC                                                              */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: Ulrich Hamann, 03.01.2008                                               */
/*                                                                                 */
/***********************************************************************************/

int read_aerosol_species_from_netCDF_map ( input_struct input, output_struct *output, char *filename )
{
  int status=0;

  int ncid=NOT_DEFINED_INTEGER;

  int nt=-1;
  long ilat=NOT_DEFINED_INTEGER, ilon=NOT_DEFINED_INTEGER;  /* index for lat and lon in netCDF file */

  double *aer_lat = NULL;
  double *aer_lon = NULL;

  size_t nlat  = 0;
  size_t nlon  = 0;
  size_t aer_nlev = 0;
  size_t aer_nlay = 0;

  int itime = -1, itime1 = -1, itime2 = -1;
  float dt = NOT_DEFINED_FLOAT;

  double *hyai = NULL;
  double *hybi = NULL;
  double *hyam = NULL;
  double *hybm = NULL;

  float *p_aer = NULL;
  float *p_aer_m = NULL;

  float *log_p_aer = NULL;
  float *log_p_aer_m = NULL;

  float *log_p_atm = NULL;

  float ***aerosol_mass_dens=NULL;

#if HAVE_LIBNETCDF
  int id_data = 0;
  char  data_unit[50] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
#endif

  float *rho_air_aero=NULL;
  int t=0;

  int lc=NOT_DEFINED_INTEGER;

  int aer_nlayer = NOT_DEFINED_INTEGER;
  int i_aer = NOT_DEFINED_INTEGER;

  if (input.verbose)
    fprintf (stderr, " ... read aerosol mass profiles from netCDF file %s \n", input.filename[FN_ECMWF]);

  status = get_all_netCDF_indices ( filename, input.latitude, input.longitude,
                                    &(ncid), &(ilat), &(nlat), &(ilon), &(nlon), 
                                    &(aer_lat), &(aer_lon),
                                    input.UTC, input.atm.time_interpolate,
                                    &(nt), &(itime1), &(itime2), &(dt),
                                    input.verbose, input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d, during get_all_netCDF_indices from '%s'\n", status, filename );
    fprintf (stderr, "          (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__ );
    return status;
  }

  /* alloc pressure for atmosphere levels */
  if ((log_p_atm = calloc (output->atm.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'log_p_atm' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  /* calculate log-p grid for atmosphere */
  for (lc=0; lc<output->atm.nlev; lc++) {
    log_p_atm[lc] = log(output->atm.microphys.press[lc][0][0]);
  }

  /* read hybrid A coefficient at layer interfaces */
  status = alloc_and_read_netCDF_1D_double(ncid,"ilev", &(aer_nlev), "hyai", &(hyai));
  if (status != 0) {
    fprintf (stderr, "Error %d reading 'hyai' from %s \n", status, filename );
    fprintf (stderr, "         (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__ );
    return status;
  }
  /* read hybrid B coefficient at layer interfaces */
  status = alloc_and_read_netCDF_1D_double(ncid,"ilev", &(aer_nlev), "hybi", &(hybi));
  if (status != 0) {
    fprintf (stderr, "Error %d reading 'hybi' from %s \n", status, filename );
    fprintf (stderr, "         (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__ );
    return status;
  }

  /* read hybrid A coefficient at layer midpoints */
  status = alloc_and_read_netCDF_1D_double(ncid,"mlev", &(aer_nlay), "hyam", &(hyam));
  if (status != 0) {
    fprintf (stderr, "Error %d reading 'hyam' from %s \n", status, filename );
    fprintf (stderr, "         (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__ );
    return status;
  }
  /* read hybrid B coefficient at layer midpoints */
  status = alloc_and_read_netCDF_1D_double(ncid,"mlev", &(aer_nlay), "hybm", &(hybm));
  if (status != 0) {
    fprintf (stderr, "Error %d reading 'hybm' from %s \n", status, filename );
    fprintf (stderr, "         (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__ );
    return status;
  }

  output->aer.nlev = aer_nlev;  /* convert: int <-size_t */
  aer_nlayer       = aer_nlay;  /* convert: int <-size_t */

  /* alloc pressure for aerosol levels */
  if ((p_aer = calloc (output->aer.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'p_aer' (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__ );
    return -1;
  }
  /* alloc pressure for aerosol levels */
  if ((log_p_aer = calloc (output->aer.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'log_p_aer' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  /* alloc pressure for aerosol layers */
  if ((p_aer_m = calloc (aer_nlayer, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'p_aer_m' (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__ );
    return -1;
  }
  /* alloc pressure for aerosol layers */
  if ((log_p_aer_m = calloc (aer_nlayer, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'log_p_aer_m' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  /* calculate pressure from hybrid coefficients and surface pressure */
  for ( lc=0; lc<output->aer.nlev; lc++ ) {
    p_aer[lc] = 0.01*hyai[lc] + hybi[lc]*output->atm.microphys.press[output->atm.nlev-1][0][0]; /* 0.01 == Pa -> hPa */

    if (p_aer[lc] == 0.0)
      log_p_aer[lc]=log_p_atm[0];        /* uppermost layer is equal to uppermost layer of the atmosphere */
    else 
      log_p_aer[lc]=log(p_aer[lc]);
  }

  free(hyai);
  free(hybi);

  /* calculate pressure from hybrid coefficients and surface pressure */
  for ( lc=0; lc<aer_nlayer; lc++ ) {
    p_aer_m[lc] = 0.01*hyam[lc] + hybm[lc]*output->atm.microphys.press[output->atm.nlev-1][0][0]; /* 0.01 == Pa -> hPa */
    log_p_aer_m[lc]=log(p_aer_m[lc]);
  }

  free(hyam);
  free(hybm);

  free(p_aer);

  /* calculate density of air on aer-z-grid */
  if (( rho_air_aero = (float  *) calloc (aer_nlayer, sizeof(float))) == NULL) {
    fprintf (stderr, "Error allocating memory for 'rho_air_aero' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
    return -1;
  }
  
  /* interpolation of temper on the zd-grid to the aerosol levels */
  status = arb_wvn (output->atm.nlev, log_p_atm, output->atm.microphys.dens[MOL_AIR][0][0],
                    aer_nlayer, log_p_aer_m, rho_air_aero,
                    input.atm.interpol_method_temper, 0);
  if (status!=0) {
    fprintf (stderr, "Error %d during interpolation of 'rho_air_aero' (line %d, function %s in %s)\n", 
             status, __LINE__, __func__, __FILE__ );
    return status;
  }

  for ( lc=0; lc<aer_nlayer; lc++ ) {
                    /* 1.e+6: convert from cm-3 to m-3; 1.e-3: convert g -> kg */
    rho_air_aero[lc] *= 1.e+6 * 1.e-3 * input.atm.mol_mass[MOL_AIR] / AVOGADRO;  
  }

  free(p_aer_m);

  /* aerosol z grid */
  if (( output->aer.zd = calloc (output->aer.nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'z_aer' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  /* linear interpolation of z on the log-p grid to the aerosol levels */
  status = arb_wvn (output->atm.nlev, log_p_atm, output->atm.zd, 
                    output->aer.nlev, log_p_aer, output->aer.zd, 
                    INTERP_METHOD_LINEAR, 0);
  if (status!=0) {
    fprintf (stderr, "Error %d during interpolation of 'z_aer' (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__ );
    return status;
  }

  /* correction of interpolation uncertainty */
  if ( fabs(output->aer.zd[output->aer.nlev-1] - output->atm.zd[output->atm.nlev-1]) < 0.0000001 )
    output->aer.zd[output->aer.nlev-1] = output->atm.zd[output->atm.nlev-1];

  free(log_p_atm);
  free(log_p_aer);

  /* alloc aerosol mass density with size (timesteps,n_species,nlev) */
  if ((status = ASCII_calloc_float_3D(&(aerosol_mass_dens),nt,output->aer.n_species, aer_nlayer)) != 0) {
    fprintf (stderr,"Error, allocating memory for 'aerosol_mass_dens' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  /* read nt (= 1 or 2) time steps */
  for (t=0;t<=nt-1;t++) {

    if (t == 0)
      itime = itime1;
    if (t == 1)
      itime = itime2;

    for ( i_aer=0; i_aer<output->aer.n_species; i_aer++ ) {
      /* alloc and read aerosol mass density profile */
      status = alloc_and_read_netCDF_column_float( ncid, output->aer.species_names[i_aer], &aerosol_mass_dens[t][i_aer], 
                                                   aer_nlayer, itime, ilat, ilon, input.verbose);
      if (status != 0) {
        fprintf (stderr, "Error %d reading '%s' from netCDF file %s\n", status, output->aer.species_names[i_aer], filename);
        return status;
      }

#if HAVE_LIBNETCDF
      nc_inq_varid   ( ncid, output->aer.species_names[i_aer], &id_data ); /* status = */ 
      status = nc_get_att_text( ncid, id_data, "units", data_unit );                
      if (status == 0) {
        if (input.verbose)
          fprintf (stderr, " ... aerosol '%s', units: '%s'\n", output->aer.species_names[i_aer], data_unit );

        /* some information is given about the unit, try to interprete it */
        if (      strncasecmp( "g m**-3",   data_unit,  7) == 0 || strncasecmp( "g/m**3", data_unit,  6) == 0 ) {
          /* nothing to do, aerosols already in the right units */
        }
        else if ( strncasecmp( "kg kg**-1", data_unit, 10) == 0 || strncasecmp( "kg(aerosol)/kg(air)", data_unit, 19) == 0 ) {
          /* convert kg/kg -> g/m3 */
          for (lc=0; lc<aer_nlayer; lc++)                                     /*kg/kg->g/kg*/ /* g/kg -> g/m3 */
            aerosol_mass_dens[t][i_aer][lc] = aerosol_mass_dens[t][i_aer][lc]  *   1000.0   * rho_air_aero[lc];
        }
        else {
          fprintf (stderr," ... Error, unknown aerosol format %s (line %d, function %s in %s)\n", data_unit, __LINE__, __func__, __FILE__ );
          return -1;
        }

      }
      else /* no units given, but thats OK, we assume g/m3 */
        status = 0;
#else
      fprintf (stderr, "Error, uvspec has been compiled without netcdf which is\n");
      fprintf (stderr, "is required for the aerosol species options.\n");
      return -1;
#endif
    }
  }

  if ( (status = ASCII_calloc_float (&(output->aer.massdens), output->aer.n_species, output->aer.nlev)) != 0 ) {
    fprintf (stderr,"Error, allocating memory for 'aer.massdens' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -1;
  }

  if (nt == 1) {
    for ( i_aer=0; i_aer<output->aer.n_species; i_aer++ )
      for (lc=0; lc<aer_nlayer; lc++)
        output->aer.massdens[i_aer][lc] = aerosol_mass_dens[0][i_aer][lc]; 
        /* (lc+1 <- lc) uppermost level is empty, marks upper layer boundary */
  }
  else {
    for ( i_aer=0; i_aer<output->aer.n_species; i_aer++ )
      for (lc=0; lc<aer_nlayer; lc++)   /* (lc+1 <- lc) == uppermost level is empty, marks upper layer boundary */
        output->aer.massdens[i_aer][lc] = (1.0-dt)* aerosol_mass_dens[0][i_aer][lc]  +  dt* aerosol_mass_dens[1][i_aer][lc];
  }

  ASCII_free_float_3D(aerosol_mass_dens, nt, output->aer.n_species);

  return status;

}


/***********************************************************************************/
/* Function: read_aerosol_species_from_ASCII_file                                  */
/* Description:                                                                    */
/*  read aerosol mass densities profiles from an ASCII file                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: Claudia Emde, 17.12.2007                                                */
/*                                                                                 */
/***********************************************************************************/

static int read_aerosol_species_from_ASCII_file (char *filename, aer_out_struct *aer,
                                                 int quiet, int verbose) 
  
{
  int rows=0, min_columns=0, max_columns=0;
  int status = 0;
  float **aer_data=NULL;
  
  int i_aer=0;

  char function_name[] = "read_aerosol_species_from_ASCII_file";
  char file_name[]     = "aerosol.c";

  if (verbose)
    fprintf (stderr," ... reading aerosol species mass density profiles from %s\n", filename);
  
  status = ASCII_file2float (filename, 
                             &rows, &max_columns, &min_columns, 
                             &aer_data);

  if (status!=0) {
    fprintf (stderr, "Error %d in %s (%s) reading %s\n", status, function_name, file_name, filename);
    return status;
  }

  aer->nlev = rows;
  
  if ( max_columns != min_columns ) {
    fprintf (stderr, "Error in %s (%s): Inconsistent number of columns in %s\n", function_name, file_name, filename);
    fprintf (stderr, "    min_columns = %d, max_columns =%d\n", min_columns, max_columns);
    return -1;
  }

  if ( (min_columns-1) != aer->n_species ) { /* -1 as first column is z profile */
    fprintf (stderr, "Error in %s (%s): Number of aerosol species in %s \n", function_name, file_name, filename);
    fprintf (stderr, "   is different to the number of aerosol species specified in the library \n");
    return -1;
  }
  

  /* alloc and read height profile */
  aer->zd = (float *) calloc (aer->nlev, sizeof(float));
  aer->zd = ASCII_column_float (aer_data, aer->nlev, 0);
  

  /* allocate memory for density profiles */
  aer->massdens = calloc(aer->n_species, sizeof(float *));
  if (aer->massdens == NULL) {
    fprintf (stderr, "Error allocating memory for dens\n");
    return -1;
  }

  /* read aerosol species */
  for (i_aer=0;i_aer<aer->n_species;i_aer++)
    aer->massdens[i_aer] = ASCII_column_float(aer_data, aer->nlev, i_aer+1);

  ASCII_free_float(aer_data, rows);
  
  return 0;
}


/***********************************************************************************/
/* Function: aerosol_species_to_optical_properties                                 */
/* Description:                                                                    */
/*  convert mass density profiles from the aerosol_species_file to                 */
/*  optical properties with precalculated tables in netCDF format                  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: Ulrich Hamann , 02.01.2008 (initial version)                            */ 
/*         Claudia Emde , 22.02.2008                                               */
/*                                                                                 */
/***********************************************************************************/

static int aerosol_species_to_optical_properties ( input_struct input, output_struct *output )
{

  int status=0;
  int i_aer=0;
  
  char aerosol_opt_prop_file[FILENAME_MAX]="";
  
  float *rh;
  
  rh = calloc (output->aer.nlev, sizeof (float));

  /*   char function_name[] = "aerosol_species_to_optical_properties"; */
  /*   char file_name[]     = "aerosol.c"; */

  calculate_relative_humidity_on_aerosol_grid (input, output, &rh);
  
  /* Allocate memory for optical properties of all aerosol species */
  output->aer.optprop_opac = (optprop_struct *) calloc (output->aer.n_species, sizeof( optprop_struct ) );  
  
  for (i_aer=0; i_aer<output->aer.n_species; i_aer++) {
    
    if (status!=0) {
      fprintf (stderr, "Error %d allocating memory for aerosol optical properties\n", status);
      return status;
    }

    strcpy( aerosol_opt_prop_file, output->aer.aerosol_library);
    strcat( aerosol_opt_prop_file, output->aer.species_names[i_aer]);
    
    if (input.verbose) 
      fprintf (stderr, " ... reading aerosol optical properties from: %s \n", aerosol_opt_prop_file );
    
    /* Read aerosol properties file */
    status = read_caoth_prop (aerosol_opt_prop_file, output->wl.lambda_r, output->wl.nlambda_r, 
			      output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
			      1, 1,
			      NULL, 0,
			      &(output->aer), rh, i_aer,
			      &(output->atm.nmom), input.rte.polradtran[POLRADTRAN_NSTOKES],
			      input.rte.nstr,
			      input.rte.solver, input.rte.disort_icm,
			      input.verbose, input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d reading aerosol properties from OPAC data.\n", status);
    return status;
    }
    
  }
  
  
  status = average_aerosol_prop_opac (  &(output->aer), output->wl.nlambda_r, output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper );
  if (status!=0) {
    fprintf (stderr, "Error %d averaging aerosol properties (OPAC data).\n", status);
    return status;
    
  }
  return 0; 
  
}

/***********************************************************************************/
/* Function: average_aerosol_prop_opac                                             */
/* Description:                                                                    */
/*  The optical properties of the varios aerosol species are averaged.             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: Claudia Emde, 22.02.2008                                                */       
/*                                                                                 */
/***********************************************************************************/

static int average_aerosol_prop_opac ( aer_out_struct *aer, int nlambda, int nlambda_lower, int nlambda_upper )
{
  
  int i_aer=0, ip=0, im=0, iv=0, lc=0;
  float tau_tot=0.0, ssa_tau_tot=0.0;
  float **mom_ssa_tau_tot;
  int nphamat=0, maxmom=0; 
  int status=0; 

  double *interp_optprop_weight=NULL;
  int **tmp_ntheta=NULL;
  float ***tmp_theta=NULL, ***tmp_phase=NULL;
  double ***tmp_mu=NULL;
  int n_tot=0, n_in=0;

  /* copy nphamat */
  nphamat = aer->optprop_opac[0].nphamat;
  aer->optprop.nphamat = nphamat;

  /* alloc */
  /* this is probably not needed: */
  status = calloc_optprop(&aer->optprop, nlambda, aer->nlev, nphamat);
  if (status!=0) {
    fprintf (stderr, "Error %d allocating aer->optprop structure.\n", status);
    return status;
  }

  for (iv=nlambda_lower; iv<=nlambda_upper ; iv++) {
    for (lc=0; lc<aer->nlev-1; lc++) {
      
      tau_tot=0.0;
      ssa_tau_tot=0.0;

      /* perform sum over species from tau, ssa*tau, pmom*ssa*tau */
      for (i_aer=0; i_aer < aer->n_species; i_aer++)
        if (aer->massdens[i_aer][lc] > 0) {
          tau_tot +=  aer->optprop_opac[i_aer].dtau[iv][lc];
          ssa_tau_tot +=  aer->optprop_opac[i_aer].ssa[iv][lc] * aer->optprop_opac[i_aer].dtau[iv][lc];
        }

      /* tau summed over species */      
      aer->optprop.dtau[iv][lc] = tau_tot;
     

      /* average ssa */
      if(tau_tot>0)
        aer->optprop.ssa[iv][lc] = ssa_tau_tot/tau_tot; 
      
    }
  }

  /* moments */
  if (aer->ssprop.type==PHASE_HYBRID || aer->ssprop.type==PHASE_MOMENTS) {

    mom_ssa_tau_tot = calloc( nphamat, sizeof (float *));

    for (iv=nlambda_lower; iv<=nlambda_upper ; iv++) {
    
      for (lc=0; lc<aer->nlev-1; lc++) {
      
	maxmom=aer->optprop_opac[0].nmom[iv][lc];

	/* find maxmom */
	for (i_aer=1; i_aer < aer->n_species; i_aer++) 
	  if (aer->optprop_opac[i_aer].nmom[iv][lc] > maxmom) 
	    maxmom = aer->optprop_opac[i_aer].nmom[iv][lc];

	aer->optprop.nmom[iv][lc] = maxmom;

	/* calloc ; this is kind of umstaendlich */
	for (ip=0; ip<nphamat; ip ++){
	  mom_ssa_tau_tot[ip] = calloc ( maxmom, sizeof(float) ); 
	  aer->optprop.moment[iv][lc][ip]= calloc( maxmom, sizeof(float));
	}

	/* perform sum over species from tau, ssa*tau, pmom*ssa*tau */
	for (i_aer=0; i_aer < aer->n_species; i_aer++) { 

	  if (aer->massdens[i_aer][lc] > 0) {

	    for (ip=0; ip<nphamat; ip++)
	      for (im=0; im<aer->optprop_opac[i_aer].nmom[iv][lc]; im++)
		mom_ssa_tau_tot[ip][im] +=  
		  aer->optprop_opac[i_aer].moment[iv][lc][ip][im] *
		  aer->optprop_opac[i_aer].ssa[iv][lc] * 
		  aer->optprop_opac[i_aer].dtau[iv][lc];

	  }
	}
      
	/* average pmom */
	ssa_tau_tot = aer->optprop.dtau[iv][lc]*aer->optprop.ssa[iv][lc];
	if (ssa_tau_tot>0)
	  for (ip=0; ip<nphamat; ip++){
	    for (im=0; im<maxmom; im++)
	      aer->optprop.moment[iv][lc][ip][im] =  mom_ssa_tau_tot[ip][im]/ssa_tau_tot;

	    free ( mom_ssa_tau_tot[ip] ); 
	  }
      }
    }

    free ( mom_ssa_tau_tot );

    for (iv=nlambda_lower; iv<=nlambda_upper ; iv++) {
      for (lc=0; lc<aer->nlev-1; lc++) {
	for (ip=0; ip<nphamat; ip++){
	  maxmom=aer->optprop_opac[0].nmom[iv][lc];

	  /* find maxmom */
	  for (i_aer=1; i_aer < aer->n_species; i_aer++) 
	    if (aer->optprop_opac[i_aer].nmom[iv][lc] > maxmom) 
	      maxmom = aer->optprop_opac[i_aer].nmom[iv][lc];
	}
      }
    }
  }

  /* phase functions */
  if (aer->ssprop.type==PHASE_HYBRID || aer->ssprop.type==PHASE_EXPLICIT) {

    interp_optprop_weight = calloc ( aer->n_species, sizeof(double));

    for (iv=nlambda_lower; iv<=nlambda_upper ; iv++) {

      for (lc=0; lc<aer->nlev-1; lc++) {

	/* the following has strong resemblance with part of the function interpolate_ssprop_in_reff in cloud.c */

	/* The moments and phases need to be weighted with alb*ext */
	n_tot=0;
	ssa_tau_tot=0.0;
	for (i_aer=0; i_aer < aer->n_species; i_aer++)
	  if (aer->massdens[i_aer][lc] > 0) {
	    interp_optprop_weight[n_tot] = aer->optprop_opac[i_aer].ssa[iv][lc] * 
	      aer->optprop_opac[i_aer].dtau[iv][lc];
	    ssa_tau_tot += interp_optprop_weight[n_tot];
	    n_tot++;
	  }

	/* normalize */
	for (i_aer=0; i_aer < n_tot; i_aer++)
	  interp_optprop_weight[i_aer] /= ssa_tau_tot;

	/* allocate temporary theta */
	tmp_ntheta = calloc(nphamat, sizeof(int *));
	tmp_theta  = calloc(nphamat, sizeof(float **));
	tmp_mu     = calloc(nphamat, sizeof(double **));
	tmp_phase  = calloc(nphamat, sizeof(float **));

	for (ip=0; ip<nphamat; ip++){
	  tmp_ntheta[ip] = calloc(n_tot, sizeof(int));
	  tmp_theta [ip] = calloc(n_tot, sizeof(float *));
	  tmp_mu    [ip] = calloc(n_tot, sizeof(double *));
	  tmp_phase [ip] = calloc(n_tot, sizeof(float *));

	  n_in=0;
	  for (i_aer=0; i_aer < aer->n_species; i_aer++)
	    if (aer->massdens[i_aer][lc] > 0) {
	      tmp_ntheta[ip][n_in] = aer->optprop_opac[i_aer].ntheta[iv][lc][ip];
	      tmp_theta [ip][n_in] = aer->optprop_opac[i_aer].theta[iv][lc][ip];
	      tmp_mu    [ip][n_in] = aer->optprop_opac[i_aer].mu[iv][lc][ip];
	      tmp_phase [ip][n_in] = aer->optprop_opac[i_aer].phase[iv][lc][ip];
	      n_in++;
	    }
	}

	/* perform interpolation */
	status = sort_and_add_weighted_phase (n_tot, interp_optprop_weight,
					      tmp_ntheta, tmp_theta, tmp_mu, tmp_phase,
					      &(aer->optprop.ntheta[iv][lc]),
					      &(aer->optprop.theta [iv][lc]),
					      &(aer->optprop.mu    [iv][lc]),
					      &(aer->optprop.phase [iv][lc]),
					      nphamat, 0, 1 );

	if (status) {
	  fprintf(stderr,"Error! something (%d) went wrong in sort_and_add_weighted_phase \n",status);
	  return status;
	}

	for (ip=0; ip<nphamat; ip++){
	  free(tmp_ntheta[ip]);
	  free(tmp_theta[ip]);
	  free(tmp_mu[ip]);
	  free(tmp_phase[ip]);
	}
	free(tmp_ntheta);
	free(tmp_theta);
	free(tmp_mu);
	free(tmp_phase);

      }
    }

    free(interp_optprop_weight);

  } /* end phase functions */

  return 0;
  
}


/**************************************************************/
/* Scale the total aerosol optical thickness to a given value */
/**************************************************************/

static int scale_aer_tau (float **dtau, float *lambda_r, int nlambda_r, int nlyr,
			  float tau, float *zd, float altitude)
{
  int lc=0, iv=0, ialt=0;
  float ext=0, fact=0;
 
  /* determine level number of 'altitude' level */
  for (ialt=0; ialt<=nlyr; ialt++)
    if (altitude == zd[ialt])
      break;
  
  if (ialt == nlyr+1) {
    fprintf (stderr, "Error, user-defined altitude %f not found in altitude grid\n", altitude);
    return -1;
  }

  for (iv=0; iv<nlambda_r; iv++) {
 
    ext = 0.0;
    for (lc=0; lc<=ialt-1; lc++)
      ext += dtau[iv][lc];
    
    fact = tau / ext;
    
    for (lc=0; lc<nlyr; lc++)
      dtau[iv][lc] *= fact;
  }
  
  return 0;
} 


int aerosol_set_visibility (aer_inp_struct inp, char *path, float *visib, int verbose)
{
  int status=0;

  char shettle_filename[FILENAME_MAX] = "";

  double *t0=NULL, *t1=NULL, *t2=NULL, *t3=NULL;
  double *v0=NULL, *v1=NULL, *v2=NULL, *v3=NULL;
  double *shettle_vis=NULL, *shettle_tau550=NULL;
  int shettle_n=0;
  double tau550=0, visibility=0;

  double first=0, last=0;


  /* set visibility */
  *visib = inp.visibility;
    
  /* we may need to correct the visibility if the user fixed the optical */
  /* thickness and did not adjust the visibility properly;               */
  if (inp.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SET]   >= 0.0 || 
      inp.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE] >= 0.0 || 
      inp.beta      >= 0.0 ||
      inp.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET]    >= 0.0)  {

    if (verbose)
      fprintf (stderr, " ... need to adjust visibility because user changed the optical thickness:\n");


      /* Read table of visibility vs. integrated optical thickness at 550nm; */
      /* the extinction profile is determined by aerosol_season and          */
      /* aerosol_vulcan whereas aerosol_haze and relative humidity only      */
      /* affect asymmetry parameter and single scattering albedo.            */

    sprintf (shettle_filename, "%s/aerosol/shettle/tau550.%1d.%1d", 
	     path, inp.vulcan, inp.seasn);

    status = read_2c_file (shettle_filename, 
			   &shettle_vis, &shettle_tau550, &shettle_n);
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, shettle_filename);
      return status;
    }

    /* now determine the user-defined optical thickness at 550nm */

    if (inp.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SET] >= 0.0) {
      tau550 = inp.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SET];

      if (verbose)
	fprintf (stderr, " ...   aerosol_set_tau sets the optical thickness at 550nm to %f\n", 
		 tau550);
    }
      
    if (inp.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE] >= 0.0) {

      /* determine lower and upper visibility limits */
      first = shettle_vis[0];
      last  = shettle_vis[shettle_n-1];
      sort_double (&first, &last);

      if (*visib >= first && *visib <= last) {
	/* calculate coefficients for tau550(vis) */
	status = slinear_coeffc (shettle_vis, shettle_tau550, shettle_n, &t0, &t1, &t2, &t3);
	if (status!=0) {
	  fprintf (stderr, "Error %d calculating interpolation coefficients for %s\n", 
		   status, shettle_filename);
	  return status;
	}
	
	status = calc_splined_value (*visib, &tau550, shettle_vis, shettle_n, t0, t1, t2, t3);
	if (status!=0) {
	  fprintf (stderr, "Error %d interpolating optical thickness at visibility %f\n", 
		   status, *visib);
	  return status;
	}

	free(t0); free(t1); free(t2); free(t3);
      }
      else {
	if (*visib <= first)
	  tau550 = (shettle_tau550[0]>shettle_tau550[shettle_n-1]?shettle_tau550[0]:shettle_tau550[shettle_n-1]);

	if (*visib >= last)
	  tau550 = (shettle_tau550[0]<shettle_tau550[shettle_n-1]?shettle_tau550[0]:shettle_tau550[shettle_n-1]);
      }

      tau550 *= inp.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE];

      if (verbose)
	fprintf (stderr, " ...   aerosol_scale_tau sets the optical thickness at 550nm to %f\n", 
		 tau550);
    }

    if (inp.beta >= 0.0) {
      tau550 = inp.beta * pow(0.55, -inp.alpha);
	
      if (verbose)
	fprintf (stderr, " ...   aerosol_angstrom sets the optical thickness at 550nm to %f\n", 
		 tau550);
    }

    if (inp.alpha_2 >= -800.0) {
      tau550 = exp(inp.alpha_0) * pow(0.55, inp.alpha_1) * pow(0.55, inp.alpha_2*log(0.55));
	
      if (verbose)
	fprintf (stderr, " ...   aerosol_king_byrne sets the optical thickness at 550nm to %f\n", 
		 tau550);
    }    

      
    if (inp.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET] >= 0.0) {
      tau550 = inp.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET];

      if (verbose)
	fprintf (stderr, " ...   aerosol_set_tau550 sets the optical thickness at 550nm to %f\n", 
		 tau550);
    }



    /* determine lower and upper visibility limits */
    first = shettle_tau550[0];
    last  = shettle_tau550[shettle_n-1];
    sort_double (&first, &last);


    /* now determine the visibility which corresponds to tau550 */
    if (tau550 >= first && tau550 <= last) {

      /* calculate linear interpolation coefficients for vis(tau550) */
      status = slinear_coeffc (shettle_tau550, shettle_vis, shettle_n, &v0, &v1, &v2, &v3);
      if (status!=0) {
	fprintf (stderr, "Error %d calculating interpolation coefficients for %s\n", 
		 status, shettle_filename);
	return status;
      }

      status = calc_splined_value (tau550, &visibility, shettle_tau550, shettle_n, v0, v1, v2, v3);

      if (status!=0) {
	fprintf (stderr, "Error %d interpolating visibility at optical thickness %f\n", 
		 status, tau550);
	return status;
      }
	
      free(v0); free(v1); free(v2); free(v3);
    }
    else {
      if (tau550 <= first)
	visibility = (shettle_vis[0]>shettle_vis[shettle_n-1]?shettle_vis[0]:shettle_vis[shettle_n-1]);
	
      if (tau550 >= last)
	visibility = (shettle_vis[0]<shettle_vis[shettle_n-1]?shettle_vis[0]:shettle_vis[shettle_n-1]);
    }
      
    *visib = visibility;

    free(shettle_vis);
    free(shettle_tau550);


    if (verbose)
      fprintf (stderr, " ...   changing visibility to %f km\n", *visib);
  }

  return 0;
}



int calloc_aer_out (aer_out_struct *out, int nlambda, int nlev, int nphamat) 
{
  int  status=0;
  
  out->nlev    = nlev;
  out->zd      = (float *) calloc (nlev, sizeof (float));

  status = calloc_optprop (&(out->optprop), nlambda, nlev, nphamat);
  if (status!=0) {
    fprintf (stderr, "Error %d allocating memory for aerosol optical properties\n", status);
    return status;
  }

  return 0;

}

int cp_aer_out (aer_out_struct *source, aer_out_struct *target, int nlambda, int alloc_moment, int quiet) 
{
  int lc=0, nlev=0, status=0;
  nlev = (target->nlev < source->nlev ? target->nlev : source->nlev);

  for (lc=0; lc<nlev; lc++)
    target->zd[lc] = source->zd[lc];
  
  status = cp_optprop (&(target->optprop), &(source->optprop), nlev, alloc_moment, quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d copying aerosol optical properties\n", status);
    return status;
  }
  return 0;
}


int free_aer_out (aer_out_struct *out) 
{
  free (out->zd);

  free_optprop (&(out->optprop));
		     
  return 0;
}



