/*--------------------------------------------------------------------
 * $Id: albedo.c 3224 2016-06-03 18:36:57Z bernhard.mayer $
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
#include "fortran_and_c.h"
#include "f77-uscore.h"
#include "ascii.h"

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#define nigbp 20
#define nband 15
#define PI 3.14159265358979323846264338327


/************************************/
/* prototypes of internal functions */
/************************************/

static int read_albedo (char *filename, alb_out_struct *alb_out);
static int read_rpv    (char *filename, rpv_out_struct *rpv,
			float sigma, float t1, float t2, float scale);
static int interpolate_rossli (rossli_out_struct *rossli, float *lambda, 
			       int nlambda_lower, int nlambda_upper);
static int interpolate_hapke (hapke_out_struct *hapke, float *lambda, 
			      int nlambda_lower, int nlambda_upper);
static int interpolate_rpv (rpv_out_struct *rpv, float *lambda, 
			    int nlambda_lower, int nlambda_upper);
static int read_spectral_BRDFs (char *filename, surfaces_out_struct *surf,
				float *lambda, int nlambda,
				int nlambda_lower, int nlambda_upper,
				float sigma, float t1, float t2, float scale,
				int quiet);
static int read_spectral_albedos (char *filename, surfaces_out_struct *alb, float *lambda, int nlambda,
				  int nlambda_lower, int nlambda_upper,
				  int quiet);

void F77_FUNC (landspec, LANDSPEC) (int *igbp, float *u0, float *wv, float *salb15, float *bbalb);
void F77_FUNC (getemiss, GETEMISS) (int *igbp, float *ems12);
void F77_FUNC (swdnwgts, SWDNWGTS) (float *CSZ, float *PW, float *WGT);

static void number2surfacestring (int surface_nr, char *surface_string);


/**************************************************************/
/* Setup surface albedo.                                      */
/**************************************************************/

int setup_albedo (input_struct input, output_struct *output)
{
  int status=0, iv=0, linear=0;
  float albedo = NOT_DEFINED_FLOAT;

#if HAVE_LIBNETCDF
  int  ncid    = NOT_DEFINED_INTEGER;
  int  id_data = NOT_DEFINED_INTEGER;
#endif

  char *albedo_name = calloc (FILENAME_MAX+1, sizeof(char));

  char *filename       = calloc (FILENAME_MAX+1, sizeof(char));
  char *surface_string = calloc (30,             sizeof(char));

  /* Input Parameters (solar and thermal) */
  /* igbp - IGBP scene identification index [1 - 20] */
  int band;
  float *bands;
  float *subbands;

  /* Input Parameters (solar): */
  /*  u0   - Cosine of the solar zenith angle [0.0 - 1.0] */
  /*  wv   - PRECIPITABLE WATER (cm) Range (0 to~7) */
  float u0, wv;

  /* Output Parameters: */
  /* salb15 - spectral albedo in 15 shortwave bands */
  /* wgt15  - estimated weights for salb15 */
  /* bbalb  - solar broad band albedo */
  float *salb15   = NULL;
  float *wgt15    = NULL;
  float bbalb;

  /* Output Parameters: */
  /* ems12 - spectral emissivity in 12 longwave bands */
  float *ems12 = NULL;

  float sum_wgt=0;
  float sum_alb=0;
  int iq=0;

  float *file_lower_lambda_albedo=NULL;
  float *file_upper_lambda_albedo=NULL;

  char tmp_string[2];

  char function_name[]="setup_albedo";
  char file_name[]="albedo.c";

  /* copy albedo input to output structure*/
  output->alb.source        = input.alb.source;
  output->alb.surface       = input.alb.surface;
  output->alb.albedo        = input.alb.albedo;
  output->alb.surf_type_map = input.alb.surf_type_map;

  /* allocate required arrays */
  output->alb.albedo_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));

  output->rpv.rho0_r   = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->rpv.k_r      = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->rpv.theta_r  = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->rpv.scale_r  = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->rpv.sigma_r  = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->rpv.t1_r     = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->rpv.t2_r     = (float *) calloc (output->wl.nlambda_r, sizeof(float));

  output->rossli.iso_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->rossli.vol_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->rossli.geo_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));

  output->hapke.b0_r   = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->hapke.h_r    = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->hapke.w_r    = (float *) calloc (output->wl.nlambda_r, sizeof(float));

  /* if there is a surface map given                            */
  /* search surface type in the map according to lat, lon, time */
  /* used by ALBEDO_USER_LIBRARY, ALBEDO_IGBP_LIBRARY, */
  if (input.alb.surf_type_map == TRUE) {

    /* read surface type from netCDF map */
    status = get_number_from_netCDF_map (input.latitude, input.longitude, input.UTC, input.atm.time_interpolate,
                                         input.filename[FN_SURFACE_TYPE_MAP], &(output->alb.surface), TYPE_INT, input.alb.netCDF_surf_name,
                                         input.verbose, input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d reading surface type map from \n %s\n (line %d, function %s in %s) \n", 
                       status, input.filename[FN_SURFACE_TYPE_MAP], __LINE__, __func__, __FILE__);
      return status;
    }
  }


  switch(input.alb.source) {
  case ALBEDO_CONSTANT:
  case ALBEDO_FROM_ALBEDO_MAP:
  case ALBEDO_FROM_EMISSIVITY_MAP:

    if (input.alb.source == ALBEDO_FROM_ALBEDO_MAP || input.alb.source == ALBEDO_FROM_EMISSIVITY_MAP) {

#if HAVE_LIBNETCDF
      /* open netcdf file */
      status = nc_open (input.filename[FN_ALBEDO_MAP], NC_NOWRITE, &ncid);
      if (status != NC_NOERR) { 
        fprintf (stderr, "Error '%s' opening netCDF file '%s' (line %d, function %s in %s)\n", 
                 nc_strerror(status), input.filename[FN_ALBEDO_MAP], __LINE__, __func__, __FILE__);
        return -abs(status);
      }

      /* try to find the correct name for the albedo */
      if ( strlen(input.alb.netCDF_alb_name) > 0 )
        strcpy (albedo_name, input.alb.netCDF_alb_name);
      else if (input.alb.source == ALBEDO_FROM_ALBEDO_MAP) {
        if      ( (status = nc_inq_varid (ncid, "FAL", &id_data)) == NC_NOERR )  /* forecast albedo */
          strcpy (albedo_name, "FAL");
        else if ( (status = nc_inq_varid (ncid, "fal", &id_data)) == NC_NOERR )  /* forecast albedo */
          strcpy (albedo_name, "fal");
        else if ( (status = nc_inq_varid (ncid, "AL",  &id_data)) == NC_NOERR )   /* albedo */
          strcpy (albedo_name, "AL");
        else if ( (status = nc_inq_varid (ncid, "al",  &id_data)) == NC_NOERR )   /* albedo */
          strcpy (albedo_name, "al"); 
        else {
          fprintf (stderr, "Error '%s' while getting id for surface albedo from '%s' \n", 
                           nc_strerror(status), input.filename[FN_ALBEDO_MAP]);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return -1;
        }
      }
      else if (input.alb.source == ALBEDO_FROM_EMISSIVITY_MAP) {
        if      ( (status = nc_inq_varid (ncid, "EMIS", &id_data)) == NC_NOERR )   /* emissivity */
          strcpy (albedo_name, "EMIS");
        else if ( (status = nc_inq_varid (ncid, "emis", &id_data)) == NC_NOERR )   /* emissivity */
          strcpy (albedo_name, "emis"); 
        else {
          fprintf (stderr, "Error '%s' while getting id for surface emissivity from '%s' \n", 
                           nc_strerror(status), input.filename[FN_ALBEDO_MAP]);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return -1;
        }
      }
      nc_close (ncid);

      /* read surface type from netCDF map */
      status = get_number_from_netCDF_map (input.latitude, input.longitude, input.UTC, input.atm.time_interpolate,
                                           input.filename[FN_ALBEDO_MAP], &(albedo), TYPE_FLOAT, albedo_name,
                                           input.verbose, input.quiet);
      if (status!=0) {
        fprintf (stderr, "Error %d reading albedo map %s\n", status, input.filename[FN_ALBEDO_MAP]);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }

      albedo /= input.alb.scale_factor;

      if (input.alb.source == ALBEDO_FROM_EMISSIVITY_MAP)
        albedo = 1.0 - albedo;

#else
      fprintf (stderr, "Error, uvspec has been compiled without netcdf which is\n");
      fprintf (stderr, "is required for the emissivity map function.\n");
      return -1;
#endif

    }
    else {
      albedo = input.alb.albedo;
    }

    /* security check of albedo range */
    if ( albedo < 0.0 || 1.0 < albedo ) {
      fprintf (stderr, "\nError, albedo value = %f is out of range [ 0, 1 ] !!! \n", albedo ); 
      fprintf (stderr,   "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    if (input.verbose)
      fprintf (stderr, "     wavelength indenpendent albedo of %f \n", albedo);

    for (iv=0;iv<output->wl.nlambda_r;iv++)
      output->alb.albedo_r[iv] = albedo; 
    break;

  case ALBEDO_FROM_ALBEDO_FILE:
  case ALBEDO_USER_LIBRARY:

    if (input.alb.source == ALBEDO_FROM_ALBEDO_FILE)

      /* copy albedo file name */
      strcpy (filename, input.filename[FN_ALBEDO]);

    else if (input.alb.source == ALBEDO_USER_LIBRARY) { 

      /* check if surface_type was initialized */
      if (output->alb.surface == NOT_DEFINED_INTEGER || output->alb.surface != output->alb.surface) {
        fprintf (stderr, "\nError, surface type = %d is not initialsed.\n", output->alb.surface ); 
        fprintf (stderr,   "       Please fix your input file.\n"); 
        return -1;
      }

      /* security check for surface type range */
      if (output->alb.surface < 1 || 99 < output->alb.surface) {
        fprintf (stderr, "Error, albedo_surface %5d out of range \n", output->alb.surface);
        fprintf (stderr, "       only range [ 1, 99] is valid in %s (%s)\n", function_name, file_name);
        return -1;
      }

      /* create albedo file name */
      strcpy (filename, input.filename[FN_ALBEDO_LIB_PATH]);
      strcat (filename, "albedo_");
      if (output->alb.surface < 10)
        sprintf(tmp_string, "0%1d",output->alb.surface);
      else
        sprintf(tmp_string, "%2d",output->alb.surface);
      strcat (filename, tmp_string);
      strcat (filename, ".dat");

    }
    else {
      fprintf (stderr, "Error, unknown albedo source %d in %s (%s) \n", input.alb.source, function_name, file_name);
      return -1;
    }

    /* black underground, for safety set albedo_r = 0 (should already be 0) */
    if (input.alb.source == ALBEDO_USER_LIBRARY && output->alb.surface == 0) {
      for (iv=0; iv<output->wl.nlambda_r; iv++)
        output->alb.albedo_r[iv]=0;
    }
    else {
      status = read_albedo (filename, &output->alb);
      if (status != 0) {
        fprintf (stderr, "Error reading albedo file '%s' \n", input.filename[FN_ALBEDO]);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }
    
      /**** Interpolate albedo to internal wavelength grid ****/
    
      linear = 1;  /**** Use linear interpolation for the albedo, otherwise 
                         negative albedos or albedos larger than 1 may occur, 
                         I should know :-) AKy, ****/ 

      status = arb_wvn2 (output->alb.file_nlambda_albedo, output->alb.file_lambda_albedo,
                         output->alb.file_albedo, 
                         output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper, 
                         output->wl.lambda_r, output->alb.albedo_r, 
                         linear);

      if (status != 0) {
        fprintf (stderr, "Error %d returned by arb_wvn2()\n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }
    }
    
    break;
  case ALBEDO_IGBP_LIBRARY:

    /* check if surface_type was initialised */
    if (output->alb.surface == NOT_DEFINED_INTEGER || output->alb.surface != output->alb.surface) {
      fprintf (stderr, "\nError, surface type is not initialsed.\n"); 
      fprintf (stderr,   "       Please fix your input file using one of the options\n"); 
      fprintf (stderr,   "       'surface_type' or 'surface_type_map'.\n\n"); 
      fprintf (stderr,   "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    /* check, if 'surface_type' is inside the known range of surface_types of IGBP library */
    if (output->alb.surface < 0 || 20 < output->alb.surface) {
      fprintf (stderr, "Error: surface type nr %d for albedo_surface not known!\n", output->alb.surface);
      return -1;
    }

    /* get input values for albedo functions */

    /* Cosine of the solar zenith angle [0.0 - 1.0] */
    /* if there are clouds, then albedo is calculated as if the SZA is 53 degrees */
    if ( input.n_caoth>0 && ( input.i_ic!=-1 || input.i_wc!=-1) )
      /* cloud has been detected, IGBP wants sza=53 */
      u0 = cos(53./360.*2.*PI);
    else {
      /* no cloud detected */
      u0 = cos(output->atm.sza_r [iv]/360*2*PI); 
      if ( input.n_caoth>0 ) {
	/* ... but caoth detected! Need to warn user */
	fprintf(stderr," *** Warning! You are using the option 'input' and igbp-albedo!\n");
	fprintf(stderr," *** It is not clear whether your 'input' are clouds,\n");
	fprintf(stderr," *** Therefore I will assume that they are not.\n");
	fprintf(stderr," *** As a consequence, sza will not be set to 53degr\n");
	fprintf(stderr," *** for the albedo calculation (which is normally\n");
	fprintf(stderr," *** done in the presence of clouds).\n");
	fprintf(stderr," *** If you do want sza=53, try defining a cloud, or\n");
	fprintf(stderr," *** contact the libRadtran developers for solution.\n");
      }
    }

    /* Precipitable water in cm (== 10 kg/m2), range = [ 0.0, 7.0]  */
    wv = 0.1 * output->precip_water;  /* 1/10. == kg/m2 -> cm */ 
    if (wv > 7.0) { 
      wv = 7.0;  /* maximum range given by parameterization */
                 /* should libRadtran give a warning here?  */
    }

    if (output->alb.surface == 0) /* black surface */
      return 0;
    else {
      if (input.verbose) {
        number2surfacestring(output->alb.surface, surface_string);
        fprintf (stderr, "    IGBP surface type nr %d, %s\n", output->alb.surface, surface_string);
      }

      strcpy (input.filename[FN_ALBEDO_LIB_WVL], input.filename[FN_PATH]);
      strcat (input.filename[FN_ALBEDO_LIB_WVL], "/albedo/IGBP_map/wvl_albedo.dat");

      /* read center wavelength and band limits [wavenumbers] from file */
      status = read_5c_file_float (input.filename[FN_ALBEDO_LIB_WVL],
				   &(bands),
				   &(subbands),
				   &(output->alb.file_lambda_albedo),
                                   &(file_lower_lambda_albedo),
                                   &(file_upper_lambda_albedo),
				   &(output->alb.file_nlambda_albedo));
      if (status != 0) {
        fprintf (stderr, "Error %d, while reading file %s in setup_albedo (albedo.c) !\n", status, input.filename[FN_ALBEDO_LIB_WVL]);
        return -1;
      }

      for (band=0; band<output->alb.file_nlambda_albedo; band++) {
        file_lower_lambda_albedo[band]       *= 1000.0; /* 1000 == micro m -> nm */
        file_upper_lambda_albedo[band]       *= 1000.0; /* 1000 == micro m -> nm */
        output->alb.file_lambda_albedo[band] *= 1000.0; /* 1000 == micro m -> nm */
      }

      if ((ems12 = calloc (12, sizeof(float)) ) == NULL) {
        fprintf (stderr,"Error, allocating memory for 'ems12' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -17;
      }

      if ((salb15 = calloc (15, sizeof(float)) ) == NULL) {
        fprintf (stderr,"Error, allocating memory for 'salb15' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -17;
      }

      if ((wgt15 = calloc (15, sizeof(float)) ) == NULL) {
        fprintf (stderr,"Error, allocating memory for 'wgt15' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -17;
      }

      /* solar albedo */
      F77_FUNC (landspec, LANDSPEC) (&output->alb.surface, &u0, &wv, salb15, &bbalb);
      /* estimated weights for solar albedo */
      F77_FUNC (swdnwgts, SWDNWGTS) (&u0, &wv, wgt15);
      /* thermal albedo */
      F77_FUNC (getemiss, GETEMISS) (&output->alb.surface, ems12);

      output->alb.file_albedo = calloc (15+12, sizeof(float));

      for (band=0;band<(15+12);band++) {
        if (band<15)
          output->alb.file_albedo[band] = salb15[band];
	else
	  output->alb.file_albedo[band] = 1 - ems12[band-15];
      }

      free(ems12);
      free(salb15);

      if (input.verbose) { 
        fprintf(stderr,"\n    band    lambda     lambda      lambda    albedo    weight\n");
        fprintf(stderr,  "            lower      center      upper                     \n");
        fprintf(stderr,  "   -----------------------------------------------------------");
        for (band=0;band<(15+12);band++) {
          fprintf(stderr,"\n    %3d %10.2f %10.2f %11.2f   %8.6f ",
                  band+1, 
                  file_lower_lambda_albedo[band], 
                  output->alb.file_lambda_albedo[band], 
                  file_upper_lambda_albedo[band], 
		  output->alb.file_albedo[band]);
          if (band<15)
            fprintf(stderr," %8.6f", wgt15[band]);
        }
        fprintf(stderr,"\n\n");
      }

    }

    if (input.ck_scheme == CK_FU) {

      /* set wavelength grid */
      output->alb.file_nlambda_albedo = output->wl.nlambda_r;

      for (band = output->wl.nlambda_rte_lower; band <= output->wl.nlambda_rte_upper; band++) {
	if (band == 0) {
          sum_wgt=0;
          sum_alb=0;
          for (iq=0;iq<10;iq++) {
            sum_wgt += wgt15[iq];
            sum_alb += wgt15[iq]*output->alb.file_albedo[iq];
          }
          output->alb.albedo_r[band] = sum_alb / sum_wgt; 
	}
	else
          output->alb.albedo_r[band] = output->alb.file_albedo[band+9];
      }
    }
    else {

      /* spectral albedo constant inside one band */

      for (iv=0; iv<output->wl.nlambda_r; iv++) {
        output->alb.albedo_r[iv] = -999;
        for (band = 0; band < output->alb.file_nlambda_albedo; band ++) {
          if ((file_lower_lambda_albedo[band] <= output->wl.lambda_r[iv]) && 
              (output->wl.lambda_r[iv] <= file_upper_lambda_albedo[band])) {
	    /* fprintf (stderr, "band[%2d]=%f,lambda_r[%2d]=%f, band[%d]=%f \n",             */
            /*   band,  output->alb.file_lambda_albedo[band]  , iv, output->wl.lambda_r[iv], */
            /*   band+1,output->alb.file_lambda_albedo[band+1]);                            */
            output->alb.albedo_r[iv] = output->alb.file_albedo[band];
	    break;
	  }
	}
        if (output->alb.albedo_r[iv] == -999) {
          if (4000.0 <= output->wl.lambda_r[iv] && output->wl.lambda_r[iv] <= 4545.45) { 
            /* not defined in this region: make linear interpolation (is this OK????) */
            output->alb.albedo_r[iv] = output->alb.file_albedo[14] + 
	          (output->wl.lambda_r[iv]-4000.0) / (4545.45 - 4000.0) * 
                  (output->alb.file_albedo[15]-output->alb.file_albedo[14]);
	  }
          else {
            fprintf (stderr, "Error %d at setup albedo_IGBP! \n", -999);
            return -999;
	  }
	}
      }
    }
    free(wgt15);
    break;
  default:
    fprintf (stderr, "Error: unsupported albedo source in %s (%s)\n", function_name, file_name);
    return -1;
  }
 
  free(surface_string);

  /**** Check interpolated albedo values ****/

  for (iv=0; iv<output->wl.nlambda_r; iv++) {
/*      if (input.verbose) { */
/*        fprintf (stderr, "lambda   albedo"); */
/*        fprintf (stderr, "%10.3f %9.5f\n",output->wl.lambda_r[iv], output->alb.albedo_r[iv]); */
/*      } */
    if (output->alb.albedo_r[iv] > 1.0 || output->alb.albedo_r[iv] < 0.0)  {
      fprintf (stderr, "Albedo interpolated from %s out of range:\n", input.filename[FN_ALBEDO]);
      fprintf (stderr, "wavelength: %f, albedo %f\n", 
                 output->wl.lambda_r[iv], output->alb.albedo_r[iv]);
      fprintf (stderr, "Modify albedo file %s\n", input.filename[FN_ALBEDO]);
      return -1;
    } 
  }


  /* read RossLi BRDF */
  switch(input.rossli.source) {
  case ROSSLI_CONSTANT:
 
    /* direct input from input file, wavelength independent */
    for (iv=0;iv<output->wl.nlambda_r;iv++) {
      output->rossli.iso_r [iv] = input.rossli.rossli[BRDF_ROSSLI_ISO];
      output->rossli.geo_r [iv] = input.rossli.rossli[BRDF_ROSSLI_GEO];
      output->rossli.vol_r [iv] = input.rossli.rossli[BRDF_ROSSLI_VOL];
    }

    break;
  case ROSSLI_AMBRALS_CONSTANT:
 
    /* direct input from input file, wavelength independent */
    for (iv=0;iv<output->wl.nlambda_r;iv++) {
      output->rossli.iso_r [iv] = input.rossli.rossli[BRDF_ROSSLI_ISO] / PI;
      output->rossli.geo_r [iv] = input.rossli.rossli[BRDF_ROSSLI_GEO] / PI;
      output->rossli.vol_r [iv] = input.rossli.rossli[BRDF_ROSSLI_VOL] * 3. / 4.;
    }

    break;
  case ROSSLI_FROM_ROSSLI_FILE:
  case ROSSLI_FROM_AMBRALS_FILE:
    /* copy albedo file name */
    strcpy (filename, input.filename[FN_ROSSLI]);

    /* check, if filename is given */
    if (strlen(filename) == 0) {
      fprintf (stderr, "Error (program bug), no filename is given for RossLi input (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    }

    if (input.verbose)
      fprintf (stderr, " ... read BRDF function (RossLi) from file %s (line %d, function %s in %s)\n", filename, __LINE__, __func__, __FILE__);

    /* read file */
    status = read_4c_file_float (filename,
				 &(output->rossli.lambda),
				 &(output->rossli.iso),
				 &(output->rossli.vol),
				 &(output->rossli.geo),
				 &(output->rossli.n));

    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    if (input.rossli.source==ROSSLI_FROM_AMBRALS_FILE){
      for (iv=0; iv<output->rossli.n; iv++){
	output->rossli.iso[iv] /= PI;
	output->rossli.geo[iv] /= PI;
	output->rossli.vol[iv] /= 4. / 3.;
      }
    }

    /* interpolate rpv values to internal wavelength grid */
    status = interpolate_rossli (&(output->rossli), output->wl.lambda_r, 
				 output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by interpolate_rossli() (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
      return status;
    }

    break;
  default:
    fprintf (stderr, "Error (program bug!): unknown rossli source (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__);
    return -1;
  }


  /* read Hapke BRDF */
  switch(input.hapke.source) {
  case HAPKE_CONSTANT:
 
    /* direct input from input file, wavelength independent */
    for (iv=0;iv<output->wl.nlambda_r;iv++) {
      output->hapke.w_r  [iv] = input.hapke.hapke[BRDF_HAPKE_W];
      output->hapke.b0_r [iv] = input.hapke.hapke[BRDF_HAPKE_B0];
      output->hapke.h_r  [iv] = input.hapke.hapke[BRDF_HAPKE_H];
    }

    break;
  case HAPKE_FROM_HAPKE_FILE:
    /* copy albedo file name */
    strcpy (filename, input.filename[FN_HAPKE]);

    /* check, if filename is given */
    if (strlen(filename) == 0) {
      fprintf (stderr, "Error (program bug), no filename is given for Hapke input (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    }

    if (input.verbose)
      fprintf (stderr, " ... read BRDF function (Hapke) from file %s (line %d, function %s in %s)\n", filename, __LINE__, __func__, __FILE__);

    /* read file */
    status = read_4c_file_float (filename,
				 &(output->hapke.lambda),
				 &(output->hapke.w),
				 &(output->hapke.b0),
				 &(output->hapke.h),
				 &(output->hapke.n));

    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    /* interpolate rpv values to internal wavelength grid */
    status = interpolate_hapke (&(output->hapke), output->wl.lambda_r, 
				output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by interpolate_hapke() (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
      return status;
    }

    break;
  default:
    fprintf (stderr, "Error (program bug!): unknown hapke source (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__);
    return -1;
  }

  /* read RPV BRDF */
  switch(input.rpv.source) {
  case RPV_CONSTANT:
 
    /* direct input from input file, wavelength independent */
    for (iv=0;iv<output->wl.nlambda_r;iv++) {
      output->rpv.rho0_r  [iv] = input.rpv.rpv[BRDF_RPV_RHO0];    
      output->rpv.k_r     [iv] = input.rpv.rpv[BRDF_RPV_K];    
      output->rpv.theta_r [iv] = input.rpv.rpv[BRDF_RPV_THETA];    
      output->rpv.scale_r [iv] = input.rpv.rpv[BRDF_RPV_SCALE];    
      output->rpv.sigma_r [iv] = input.rpv.rpv[BRDF_RPV_SIGMA];    
      output->rpv.t1_r [iv]    = input.rpv.rpv[BRDF_RPV_T1];    
      output->rpv.t2_r [iv]    = input.rpv.rpv[BRDF_RPV_T2];
    }

    break;
  case RPV_FROM_RPV_FILE:
  case RPV_IGBP_LIBRARY:
  case RPV_USER_LIBRARY:

    /* if library is used -> create rpv-filename */

    if (input.rpv.source == RPV_FROM_RPV_FILE)

      /* copy albedo file name */
      strcpy (filename, input.filename[FN_RPV]);

    else if (input.rpv.source == RPV_IGBP_LIBRARY || input.rpv.source == RPV_USER_LIBRARY ) { 

      /* check if surface_type was initialised */
      if (output->alb.surface == NOT_DEFINED_INTEGER || output->alb.surface != output->alb.surface) {
        fprintf (stderr, "\nError, surface type = %d is not initialsed.\n", output->alb.surface ); 
        fprintf (stderr,   "       Please fix your input file.\n"); 
        return -1;
      }

      /* security check for surface type range */
      if (output->alb.surface < 1 || 99 < output->alb.surface) {
        fprintf (stderr, "Error, surface_type %5d is out of range \n", output->alb.surface);
        fprintf (stderr, "       only range [ 1, 99] is valid\n");
        return -1;
      }

      /* create rpv file name */
      strcpy (filename, input.filename[FN_RPV_LIB_PATH]);
      strcat (filename, "IGBP.");
      if (output->alb.surface < 10)
        sprintf(tmp_string, "0%1d",output->alb.surface);
      else
        sprintf(tmp_string, "%2d",output->alb.surface);
      strcat (filename, tmp_string);

      /* IGBP library is also dependent on NDVI, but here most common NDVI class is chosen */
      /* might be replaced by reading NDVI from netCDF map in future, UH 07-2007           */
      if ( input.rpv.source == RPV_IGBP_LIBRARY ) {
#if HAVE_MYSTIC
        switch(output->alb.surface) {
        case 1:
          strcat (filename, ".07");
          break;
        case 2:
          strcat (filename, ".08");
          break;
        case 3:
          strcat (filename, ".02");
          break;
        case 4:
          strcat (filename, ".07");
          break;
        case 5:
          strcat (filename, ".07");
          break;
        case 6:
          strcat (filename, ".07");
          break;
        case 7:
          strcat (filename, ".03");
          break;
        case 8:
          strcat (filename, ".08");
          break;
        case 9:
          strcat (filename, ".08");
          break;
        case 10:
          strcat (filename, ".04");
          break;
        case 11:
          strcat (filename, ".03");
          break;
        case 12:
          strcat (filename, ".07");
          break;
        case 13:
          strcat (filename, ".05");
          break;
        case 14:
          strcat (filename, ".07");
          break;
        case 15:
          strcat (filename, ".02");
          break;
        case 16:
          strcat (filename, ".03");
          break;
        case 17:
          strcat (filename, ".02");
          break;
        default:
          fprintf (stderr, "Error: unknown surface type %d while reading IGBP_RPV (line %d, function %s in %s)\n", 
                           output->alb.surface, __LINE__, __func__, __FILE__);
          return -1;
        }    
#else
        fprintf (stderr, " Error, the BRDF library for IGBP surface types is not availeble in this version \n");
#endif
      }
      strcat (filename, ".rpv");
    }
    else {
      fprintf (stderr, "Program bug during creation of the rpv filename (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }


    /* check, if filename is given */
    if (strlen(filename) == 0) {
      fprintf (stderr, "Error (program bug), no filename is given for rpv input (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    }

    if (input.verbose)
      fprintf (stderr, " ... read BRDF function (RPV) from file %s (line %d, function %s in %s)\n", filename, __LINE__, __func__, __FILE__);

    status = read_rpv (filename, &(output->rpv), input.rpv.rpv[BRDF_RPV_SIGMA], input.rpv.rpv[BRDF_RPV_T1], input.rpv.rpv[BRDF_RPV_T2], input.rpv.rpv[BRDF_RPV_SCALE]);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by read_rpv from %s (line %d, function %s in %s)\n", 
                       status, input.filename[FN_RPV], __LINE__, __func__, __FILE__);
      return status;
    }
    
    /* interpolate rpv values to internal wavelength grid */
    status = interpolate_rpv (&(output->rpv), output->wl.lambda_r, 
			      output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by interpolate_rpv() (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
      return status;
    }

    break;
  default:
    fprintf (stderr, "Error (program bug!): unknown rpv source (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__);
    return -1;
  }

  /* read surface information for MYSTIC */

  /* 2D spectral BRDF */
  if (strlen(input.rte.mc.filename[FN_MC_RPV_SPECTRAL]) > 0) {

    status = read_spectral_BRDFs (input.rte.mc.filename[FN_MC_RPV_SPECTRAL], &(output->rpv_surf),
				  output->wl.lambda_r, output->wl.nlambda_r,
				  output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
				  input.rpv.rpv[BRDF_RPV_SIGMA], input.rpv.rpv[BRDF_RPV_T1], input.rpv.rpv[BRDF_RPV_T2], input.rpv.rpv[BRDF_RPV_SCALE], input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by read_spectral_BRDFs from %s (line %d, function %s in %s)\n", 
                       status, input.rte.mc.filename[FN_MC_RPV_SPECTRAL], __LINE__, __func__, __FILE__);
      return status;
    }
  }
  
  /* 2D spectral albedo */
  if (strlen(input.rte.mc.filename[FN_MC_ALBEDO_SPECTRAL]) > 0) { 
    
    status = read_spectral_albedos (input.rte.mc.filename[FN_MC_ALBEDO_SPECTRAL], &(output->rpv_surf),
				    output->wl.lambda_r, output->wl.nlambda_r,
				    output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
				    input.quiet);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by read_spectral_albedos from %s (line %d, function %s in %s)\n", 
	       status, input.rte.mc.filename[FN_MC_ALBEDO_SPECTRAL], __LINE__, __func__, __FILE__);
      return status;
    }
  }


  free(albedo_name);
  free(filename);

  return 0;
}





/**************************************************************/
/* Read 2-column albedo file.                                 */
/**************************************************************/

static int read_albedo (char *filename, alb_out_struct *alb_out) 
{
  int status=0;
  
  status = read_2c_file_float (filename, 
			       &(alb_out->file_lambda_albedo),
			       &(alb_out->file_albedo),
			       &(alb_out->file_nlambda_albedo));

  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  return 0;
}


/**************************************************************/
/* Interpolate RossLi to internal wavelength grid.            */
/**************************************************************/

static int interpolate_rossli (rossli_out_struct *rossli, float *lambda, 
			       int nlambda_lower, int nlambda_upper)
{
  int status=0;

  int linear=1;   /* linear interpolation */

  /* interpolate to internal RTE grid */
  status = arb_wvn2 (rossli->n, rossli->lambda, rossli->iso,
		     nlambda_lower, nlambda_upper, lambda, rossli->iso_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (iso)\n", status);
    return status;
  }
  
  status = arb_wvn2 (rossli->n, rossli->lambda, rossli->vol,
		     nlambda_lower, nlambda_upper, lambda, rossli->vol_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (vol)\n", status);
    return status;
  }
  
  status = arb_wvn2 (rossli->n, rossli->lambda, rossli->geo,
		     nlambda_lower, nlambda_upper, lambda, rossli->geo_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (geo)\n", status);
    return status;
  }
  
  return 0;  /* if o.k. */
}


/**************************************************************/
/* Interpolate Hapke to internal wavelength grid.            */
/**************************************************************/

static int interpolate_hapke (hapke_out_struct *hapke, float *lambda, 
			      int nlambda_lower, int nlambda_upper)
{
  int status=0;

  int linear=1;   /* linear interpolation */

  /* interpolate to internal RTE grid */
  status = arb_wvn2 (hapke->n, hapke->lambda, hapke->w,
		     nlambda_lower, nlambda_upper, lambda, hapke->w_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (iso)\n", status);
    return status;
  }
  
  status = arb_wvn2 (hapke->n, hapke->lambda, hapke->b0,
		     nlambda_lower, nlambda_upper, lambda, hapke->b0_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (vol)\n", status);
    return status;
  }
  
  status = arb_wvn2 (hapke->n, hapke->lambda, hapke->h,
		     nlambda_lower, nlambda_upper, lambda, hapke->h_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (geo)\n", status);
    return status;
  }
  
  return 0;  /* if o.k. */
}


/**************************************************************/
/* Read 4- or 5- or 7- or 8-column RPV file.                  */
/**************************************************************/

static int read_rpv (char *filename, rpv_out_struct *rpv,
		     float sigma, float t1, float t2, float scale)
{
  int iv=0, status=0;

  status = read_8c_file_float (filename,
			       &(rpv->lambda),
			       &(rpv->rho0),
			       &(rpv->k),
			       &(rpv->theta),
			       &(rpv->sigma),
			       &(rpv->t1),
			       &(rpv->t2),
			       &(rpv->scale),
			       &(rpv->n));

  if (status!=0) {

    /* try 7 columns instead */
    status = read_7c_file_float (filename,
				 &(rpv->lambda),
				 &(rpv->rho0),
				 &(rpv->k),
				 &(rpv->theta),
				 &(rpv->sigma),
				 &(rpv->t1),
				 &(rpv->t2),
				 &(rpv->n));

    if (status!=0) {

      /* try 5 columns instead */
      status = read_5c_file_float (filename,
				   &(rpv->lambda),
				   &(rpv->rho0),
				   &(rpv->k),
				   &(rpv->theta),
				   &(rpv->scale),
				   &(rpv->n));

      if (status!=0) {

	/* try 4 columns instead */
	status = read_4c_file_float (filename,
				     &(rpv->lambda),
				     &(rpv->rho0),
				     &(rpv->k),
				     &(rpv->theta),
				     &(rpv->n));
    
	if (status!=0) {
	  fprintf (stderr, "Error %d reading %s\n", status, filename);
	  return status;
	}
      }
    }
  }

  /* sigma was not read, set to constant */
  if (rpv->sigma==NULL) {
    rpv->sigma = calloc (rpv->n, sizeof(float));
    for (iv=0; iv<rpv->n; iv++)
      rpv->sigma[iv] = sigma;
  }

  /* t1 was not read, set to constant */
  if (rpv->t1==NULL) {
    rpv->t1 = calloc (rpv->n, sizeof(float));
    for (iv=0; iv<rpv->n; iv++)
      rpv->t1[iv] = t1;
  }

  /* t2 was not read, set to constant */
  if (rpv->t2==NULL) {
    rpv->t2 = calloc (rpv->n, sizeof(float));
    for (iv=0; iv<rpv->n; iv++)
      rpv->t2[iv] = t2;
  }

  /* scale was not read, set to constant */
  if (rpv->scale==NULL) {
    rpv->scale = calloc (rpv->n, sizeof(float));
    for (iv=0; iv<rpv->n; iv++)
      rpv->scale[iv] = scale;
  }

  return 0; /* if o.k. */
}  



/**************************************************************/
/* Interpolate RPV to internal wavelength grid.               */
/**************************************************************/

static int interpolate_rpv (rpv_out_struct *rpv, float *lambda, 
			    int nlambda_lower, int nlambda_upper)
{
  int status=0;

  int linear=1;   /* linear interpolation */

  /* interpolate to internal RTE grid */
  status = arb_wvn2 (rpv->n, rpv->lambda, rpv->rho0,
		     nlambda_lower, nlambda_upper, lambda, rpv->rho0_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (rho0)\n", status);
    return status;
  }
  
  status = arb_wvn2 (rpv->n, rpv->lambda, rpv->k,
		     nlambda_lower, nlambda_upper, lambda, rpv->k_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (k)\n", status);
    return status;
  }
  
  status = arb_wvn2 (rpv->n, rpv->lambda, rpv->theta,
		     nlambda_lower, nlambda_upper, lambda, rpv->theta_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (theta)\n", status);
    return status;
  }
  
  status = arb_wvn2 (rpv->n, rpv->lambda, rpv->sigma,
		     nlambda_lower, nlambda_upper, lambda, rpv->sigma_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (sigma)\n", status);
    return status;
  }

  status = arb_wvn2 (rpv->n, rpv->lambda, rpv->t1,
		     nlambda_lower, nlambda_upper, lambda, rpv->t1_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (t1)\n", status);
    return status;
  }

  status = arb_wvn2 (rpv->n, rpv->lambda, rpv->t2,
		     nlambda_lower, nlambda_upper, lambda, rpv->t2_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (t2)\n", status);
    return status;
    }
  
  status = arb_wvn2 (rpv->n, rpv->lambda, rpv->scale,
		     nlambda_lower, nlambda_upper, lambda, rpv->scale_r,
		     linear);
  
  if (status != 0) {
    fprintf (stderr, "Error %d returned by arb_wvn2 (scale)\n", status);
    return status;
  }
  
  return 0;  /* if o.k. */
}
  
  

/*****************************************************************/
/* Read surface types from file, that is, label plus a filename. */
/* containing a wavelength-dependent RPV. These are assigned to  */
/* the labels used in the rpv_file. Please note that "CaM" is    */
/* interpreted as "Cox and Munk" and isn't assigned a RPV.       */
/* "Cox and Munk" is assigned index 0.                           */ 
/*****************************************************************/

static int read_spectral_BRDFs (char *filename, surfaces_out_struct *surf, float *lambda, int nlambda,
				int nlambda_lower, int nlambda_upper,
				float sigma, float t1, float t2, float scale,
				int quiet) 
{
  int il=0, status=0;
  int rows=0, min_columns=0, max_columns=0, max_length=0;
  char ***string=NULL;

  status = ASCII_checkfile (filename, 
			    &rows,
			    &min_columns,
			    &max_columns,
			    &max_length);
  
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  if (min_columns<2)
    fprintf (stderr, "Error, found less than 2 columns in %s\n", filename);
  
  if (min_columns != max_columns)
    fprintf (stderr, "*** Warning, inconsistent number of columns in %s\n", filename);

  status = ASCII_calloc_string (&string,
				rows,
				max_columns,
				max_length);

  if (status!=0) {
    fprintf (stderr, "Error %d allocating memory\n", status);
    return status;
  }

  status = ASCII_readfile (filename, string);

  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  surf->n = rows + 1;  /* plus 1 because 0 reserved for Cox and Munk */
  surf->label    = calloc (surf->n, sizeof (char *));
  surf->filename = calloc (surf->n, sizeof (char *));
  surf->rpv      = calloc (surf->n, sizeof (rpv_out_struct));

  for (il=0; il<surf->n; il++) {
    surf->label[il]    = calloc (max_length+1, sizeof(char));
    surf->filename[il] = calloc (max_length+1, sizeof(char));
  }

  strcpy (surf->label[0], "CaM");  /* Cox and Munk */

  for (il=1; il<surf->n; il++) {
    strcpy (surf->label[il], string[il-1][0]);
    strcpy (surf->filename[il], string[il-1][1]);

    /* make sure that CaM is not used as label by the user */
    if (!strcmp(surf->label[il],"CaM")) {
      fprintf (stderr, "Error, label \"CaM\" reserved for Cox and Munk! \n");
      fprintf (stderr, "Please choose a different label in %s\n", filename);
      return -1;
    }
  }
  
  ASCII_free_string (string, rows, max_columns);

  /* now read filenames and interpolate data to internal resolution */
  for (il=0; il<surf->n; il++) {
    /* allocate memory */
    surf->rpv[il].rho0_r   = (float *) calloc (nlambda, sizeof(float));
    surf->rpv[il].k_r      = (float *) calloc (nlambda, sizeof(float));
    surf->rpv[il].theta_r  = (float *) calloc (nlambda, sizeof(float));
    surf->rpv[il].scale_r  = (float *) calloc (nlambda, sizeof(float));
    surf->rpv[il].sigma_r  = (float *) calloc (nlambda, sizeof(float));
    surf->rpv[il].t1_r     = (float *) calloc (nlambda, sizeof(float));
    surf->rpv[il].t2_r     = (float *) calloc (nlambda, sizeof(float));

    /* read wavelength-dependent RPV data */
    if (il>0) {
      if (!quiet)
	fprintf (stderr, " ... reading type %s from %s\n",
		 surf->label[il], surf->filename[il]);
      status = read_rpv (surf->filename[il], &(surf->rpv[il]),
			 sigma, t1, t2, scale);
      if (status!=0) {
	fprintf (stderr, "Error %d returned by read_rpv(%s)\n", status, surf->filename[il]);
	return status;
      }

      /* interpolate to internal wavelength grid */
      status = interpolate_rpv (&(surf->rpv[il]), lambda, nlambda_lower, nlambda_upper);
      if (status!=0) {
	fprintf (stderr, "Error %d returned by interpolate_rpv()\n", status);
	return status;
      }
    }
  }

  return 0;
}


/***********************************************************************/
/* Read (MYSTIC) spectral albedo file, that is, label plus a filename. */
/* Each file contains a wavelength-dependent albedo. These are         */
/* assigned to the labels and interpolated to the internal wavelength  */
/* grid.                                                               */
/***********************************************************************/

static int read_spectral_albedos (char *filename, surfaces_out_struct *alb, float *lambda, int nlambda,
				  int nlambda_lower, int nlambda_upper,
				  int quiet)
{
  int il=0, iv=0, status=0;
  int rows=0, min_columns=0, max_columns=0, max_length=0;
  int nlambda_tmp=0;
  float *lambda_tmp=NULL, *albedo_tmp=NULL;
  char ***string=NULL;

  status = ASCII_checkfile (filename, 
			    &rows,
			    &min_columns,
			    &max_columns,
			    &max_length);
  
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  if (min_columns<2)
    fprintf (stderr, "Error, found less than 2 columns in %s\n", filename);
  
  if (min_columns != max_columns)
    fprintf (stderr, "*** Warning, inconsistent number of columns in %s\n", filename);

  status = ASCII_calloc_string (&string,
				rows,
				max_columns,
				max_length);

  if (status!=0) {
    fprintf (stderr, "Error %d allocating memory\n", status);
    return status;
  }

  status = ASCII_readfile (filename, string);

  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  alb->n = rows;  /* number of spectral albedos */
  alb->label    = calloc (alb->n, sizeof (char *));
  alb->filename = calloc (alb->n, sizeof (char *));
  alb->albedo_r = calloc (alb->n, sizeof (float *));


  for (il=0; il<alb->n; il++) {
    alb->label[il]    = calloc (max_length+1, sizeof(char));
    alb->filename[il] = calloc (max_length+1, sizeof(char));

    strcpy (alb->label[il], string[il][0]);
    strcpy (alb->filename[il], string[il][1]);
  }
  
  ASCII_free_string (string, rows, max_columns);

  /* now read filenames and interpolate data to internal resolution */
  for (il=0; il<alb->n; il++) {
    /* allocate memory */
    alb->albedo_r[il] = (float *) calloc (nlambda, sizeof(float));

    /* read wavelength-dependent RPV data */
    if (!quiet)
      fprintf (stderr, " ... reading type %s from %s\n", alb->label[il], alb->filename[il]);
    
    status = read_2c_file_float (alb->filename[il], &lambda_tmp, &albedo_tmp, &nlambda_tmp);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by read_2c_file (%s)\n", status, alb->filename[il]);
      return status;
    }
    
    /* interpolate linearly to internal wavelength grid */
    status = arb_wvn2 (nlambda_tmp, lambda_tmp, albedo_tmp,
		       nlambda_lower, nlambda_upper, lambda, alb->albedo_r[il],
		       1); 

    for (iv=0; iv<nlambda; iv++)
      fprintf (stderr, "AAA %f %f\n", lambda[iv], alb->albedo_r[il][iv]);
    
    if (status!=0) {
      fprintf (stderr, "Error %d returned by arb_wvn2()\n", status);
      return status;
    }
  }
  
  return 0;
}




static void number2surfacestring (int surface_nr, char *surface_string) {

  if      (surface_nr ==  0) strcpy (surface_string, "BLACK"); 
  else if (surface_nr ==  1) strcpy (surface_string, "EVERGREEN_NEEDLE_FOREST");  
  else if (surface_nr ==  2) strcpy (surface_string, "EVERGREEN_BROAD_FOREST"); 
  else if (surface_nr ==  3) strcpy (surface_string, "DECIDUOUS_NEEDLE_FOREST");
  else if (surface_nr ==  4) strcpy (surface_string, "DECIDUOUS_BROAD_FOREST");
  else if (surface_nr ==  5) strcpy (surface_string, "MIXED_FOREST");
  else if (surface_nr ==  6) strcpy (surface_string, "CLOSED_SHRUBS");
  else if (surface_nr ==  7) strcpy (surface_string, "OPEN_SHRUBS");
  else if (surface_nr ==  8) strcpy (surface_string, "WOODY_SAVANNA"); 
  else if (surface_nr ==  9) strcpy (surface_string, "SAVANNA"); 
  else if (surface_nr == 10) strcpy (surface_string, "GRASSLAND");
  else if (surface_nr == 11) strcpy (surface_string, "WETLAND"); 
  else if (surface_nr == 12) strcpy (surface_string, "CROPLAND");
  else if (surface_nr == 13) strcpy (surface_string, "URBAN");
  else if (surface_nr == 14) strcpy (surface_string, "CROP_MOSAIC"); 
  else if (surface_nr == 15) strcpy (surface_string, "ANTARCTIC_SNOW");
  else if (surface_nr == 16) strcpy (surface_string, "DESERT");
  else if (surface_nr == 17) strcpy (surface_string, "OCEAN_WATER"); 
  else if (surface_nr == 18) strcpy (surface_string, "TUNDRA"); 
  else if (surface_nr == 19) strcpy (surface_string, "FRESH_SNOW"); 
  else if (surface_nr == 20) strcpy (surface_string, "SEA_ICE");
  else {
    fprintf(stderr,"Unknown surface number %d\n", surface_nr); 
  }

}
