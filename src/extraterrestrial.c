/*--------------------------------------------------------------------
 * $Id: extraterrestrial.c 3194 2015-11-25 18:19:46Z svn-kylling $
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
#include "ascii.h"
#include "sun.h"

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif


/************************************/
/* prototypes of internal functions */
/************************************/

static int read_extraterrestrial(char* extraterrestrial_filename, 
				 wl_inp_struct *wl_inp, wl_out_struct *wl_out,
				 float wl_start, float wl_end, int ck_scheme,
                                 float wlband_start, float wlband_end);

static int cp_lambdat_lambdah (output_struct *output);



/**************************************************************/
/* Setup extraterrestrial irradiance.                         */
/**************************************************************/

int setup_extraterrestrial (input_struct input, output_struct *output)
{
  int status=0, iv=0;
  char temp[FILENAME_MAX]="";
  char filename[FILENAME_MAX]="";

  float *x_flt=NULL, *y_flt=NULL;
  int n_flt=0;
  char temp_filename[FILENAME_MAX+1]="";
  char *loc=NULL;
  float wlband_start = 0.0;
  float wlband_end   = 0.0;
  char function_name[]="setup_extraterrestrial";
  char file_name[]="extraterrestrial.c";
  float filter_sum=0;

  float weight_sum=0;
  int i_h, i_band, i;

  /* For special applications, pre-defined extraterrestrial */
  /* spectra are used.                                      */

  /* construct location of the extraterrestrial flux */
  strcpy (temp, input.filename[FN_PATH]);
  switch (input.ck_scheme) {
  case CK_KATO:
  case CK_KATO2:
  case CK_KATO2_96:
  case CK_KATO2ANDWANDJI:
      strcat (temp, "solar_flux/kato");
      break;
      
  case CK_FU:
      strcat (temp, "solar_flux/fu");
      break;
    
  case CK_AVHRR_KRATZ:
      strcat (temp, "solar_flux/kratz");
      break;
      
  case CK_CRS:
  case CK_RAMAN:
  case CK_LOWTRAN:
  case CK_FILE:
  case CK_REPTRAN:
  case CK_REPTRAN_CHANNEL:
      break;
      
  default:
    fprintf (stderr, "Error: unsupported correlated-k scheme %d in %s (%s)\n", input.ck_scheme, function_name, file_name );
      return -1;
  }
    

  switch (input.ck_scheme) {
  case CK_KATO:
  case CK_KATO2:
  case CK_KATO2_96:
  case CK_KATO2ANDWANDJI:
  case CK_FU:
  case CK_AVHRR_KRATZ:
      
      strcpy (filename, temp);

      /* check if we need to replace the extraterrestrial flux */
      if (strcmp(temp, input.filename[FN_EXTRATERRESTRIAL]))
	  if (!input.quiet)
	      fprintf (stderr, " ... replaced extraterrestrial flux by %s\n", filename);
      
      break;
      
  case CK_FILE:
      if (!input.quiet)
	  if (strlen(input.filename[FN_EXTRATERRESTRIAL]) > 0)
	      fprintf (stderr, " ... ignoring %s, reading extraterrestrial flux from %s\n",
		       input.filename[FN_EXTRATERRESTRIAL], input.ck_scheme_filename);
      
      strcpy (filename, input.ck_scheme_filename);
     
      break;
      
  case CK_CRS:
  case CK_RAMAN:
  case CK_LOWTRAN:
  case CK_REPTRAN:
  case CK_REPTRAN_CHANNEL:
      strcpy (filename, input.filename[FN_EXTRATERRESTRIAL]);
      break;
      
  default:
      fprintf (stderr, "Error: unsupported correlated-k scheme\n");
      return -1;
  }
    

  /* also in case of RGB conversion we need a special extraterrestrial */
  /* spectrum with a step width of 5nm                                 */
  if (input.processing == PROCESS_RGB || input.processing == PROCESS_RGBNORM) {
      strcat (temp, "solar_flux/rgb");
      
      strcpy (filename, temp);

      /* check if we need to replace the extraterrestrial flux */
      if (strcmp(temp, input.filename[FN_EXTRATERRESTRIAL]))
	  if (!input.quiet)
	      fprintf (stderr, " ... replaced extraterrestrial flux by %s\n", filename);
  }

  /* additional comments for the user (probably wrong solar spectrum file) */
  /* check kato spectrum */
  strcpy(temp_filename,"solar_flux/kato");
  loc = strstr(filename,temp_filename); /* search for spectrum_filename in filename */

  if (loc == NULL && (input.ck_scheme == CK_KATO || 
		      input.ck_scheme == CK_KATO2 || 
		      input.ck_scheme == CK_KATO2ANDWANDJI || 
		      input.ck_scheme == CK_KATO2_96) && !input.quiet)   /* no match found */
      fprintf (stderr, 
	       "\n*** Warning: correlated-k KATO_xxx specified, but another solar file: %s than standard is used!\n\n", 
	       temp_filename);
  if (loc != NULL && (input.ck_scheme != CK_KATO && 
		      input.ck_scheme != CK_KATO2 && 
		      input.ck_scheme != CK_KATO2ANDWANDJI &&
		      input.ck_scheme != CK_KATO2_96) && !input.quiet) /* match found */
      fprintf (stderr, 
	       "\n*** Warning: using solar file = %s, but not correlated-k KATO_xxx scheme! \n\n", 
	       temp_filename);

  /* check fu spectrum */
  strcpy(temp_filename,"solar_flux/fu");
  loc = strstr(filename,temp_filename); /* search for spectrum_filename in filename */

  if (loc == NULL && input.ck_scheme == CK_FU && !input.quiet)   /* no match found */
      fprintf (stderr, 
	       "\n*** Warning: correlated-k Fu specified, but another solar file: %s than standard is used!\n\n", 
	       temp_filename);
  if (loc != NULL && input.ck_scheme != CK_FU && !input.quiet) /* match found */
      fprintf (stderr, 
	       "\n*** Warning: using solar file = %s, but not correlated-k Fu scheme! \n\n", 
	       temp_filename);

  /* check AVHRR-Kratz spectrum */
  strcpy(temp_filename,"solar_flux/kratz");
  loc = strstr(filename,temp_filename); /* search for spectrum_filename in filename */

  if (loc == NULL && input.ck_scheme == CK_AVHRR_KRATZ && !input.quiet) /* no match found */
      fprintf (stderr, 
	       "\n*** Warning: correlated-k AVHRR-Kratz specified, but another solar file: %s than standard is used!\n\n", 
	       temp_filename);
  if (loc != NULL && input.ck_scheme != CK_AVHRR_KRATZ && !input.quiet) /* match found */
      fprintf (stderr, 
	       "\n*** Warning: using solar file = %s, but not correlated-k AVHRR-Kratz scheme! \n\n", 
	       temp_filename); 

  /* According to predefined spectra, set default for spectrum units */
  /* (if not explecitly specified by user in the input file)         */
  output->spectrum_unit = input.spectrum_unit;
  if (output->spectrum_unit == UNIT_NOT_DEFINED) {
    switch(input.ck_scheme) {
    case CK_FU:
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
    case CK_AVHRR_KRATZ:
    case CK_FILE:
      output->spectrum_unit = UNIT_PER_BAND;
      break;
    case CK_LOWTRAN:
    case CK_RAMAN:
    case CK_CRS:
      switch(input.source) {
      case SRC_SOLAR:
      case SRC_BLITZ:
      case SRC_LIDAR:
	output->spectrum_unit = UNIT_NOT_DEFINED;
	break;
      case SRC_THERMAL:
        if ( input.bandwidth != 1.0 || strlen(input.filename[FN_WLBANDS]) > 0 )
          output->spectrum_unit = UNIT_PER_BAND;
        else 
          output->spectrum_unit = UNIT_PER_CM_1;
	break;
      default:
	fprintf (stderr, "Error, unsupported source %d\n",
		 input.source);
	return -1;
      }
      break;
    case CK_REPTRAN:
    case CK_REPTRAN_CHANNEL:
      switch(input.source) {
      case SRC_SOLAR:
	output->spectrum_unit = UNIT_PER_NM;
	break;
      case SRC_THERMAL:
        output->spectrum_unit = UNIT_PER_CM_1;
	break;
      default:
	fprintf (stderr, "Error, unsupported source %d\n",
		 input.source);
	return -1;
      }
      break;
    default:
      fprintf (stderr, "Error, unsupported correlated-k scheme %d\n", input.ck_scheme);
      return -1;
      break;
    }
  }


  /* some warning and checks */
  switch(input.ck_scheme) {
    case CK_FU:
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
    case CK_AVHRR_KRATZ:
    case CK_FILE:
      /* no spline interpolation for fixed corr-k schemes */
      if (input.spline == 1) {
        fprintf (stderr, "Error, option 'spline ...' is not adequate for the chosen correlated-k scheme.\n");
        return -1;
      }
      break;
    case CK_LOWTRAN:
    case CK_RAMAN:
    case CK_REPTRAN:
    case CK_REPTRAN_CHANNEL:
    case CK_CRS:
      /* tell the user that he do something  which is not the usual way */
      if (input.source == SRC_SOLAR && (output->spectrum_unit!=UNIT_PER_NM && output->spectrum_unit!=UNIT_NOT_DEFINED) && !input.quiet) {
        fprintf (stderr, "*** Warning, solar_spectrum is usually not given in W/(m2 band), but in W/(m2 nm)\n"); 
        fprintf (stderr, "*** Please make sure that W/(m2 band) is really the unit of your spectrum!\n\n");
      }

      if (input.processing == PROCESS_SUM && !input.quiet && input.ck_scheme!=CK_REPTRAN) {
        fprintf (stderr, "*** Warning, 'output_process sum' might not be appropriate for spectral calculations\n");
        fprintf (stderr, "*** the sum of the output over wavelength is only the integrated flux \n");
        fprintf (stderr, "*** if there are no gaps between the (wavelength or wavenumber) bands.\n\n");
        fprintf (stderr, "*** Please make sure that this is what you really want!\n\n");
      }
      break;
    default:
      fprintf (stderr, "Error, unsupported correlated-k scheme %d\n", input.ck_scheme);
      return -1;
      break;
  }


  /* check: for heating rate calculation it is nessesary to know the output units */
  if ( input.heating != HEAT_NONE && output->spectrum_unit == UNIT_NOT_DEFINED ) {
    fprintf (stderr, "Error, for heating rate calculation, uvspec needs to know\n");
    fprintf (stderr, "the units of the solar spectrum. Please use:\n");
    fprintf (stderr, " 'source solar filename per_nm' for solar spectrums in W/(m^2 nm) or\n");
    fprintf (stderr, " 'source solar filename per_band' for bands in W/(m^2 band).\n");
    return -1;
  }

  /* can't convert undefined spectrum to a defined output */
  if ( output->spectrum_unit == UNIT_NOT_DEFINED && input.output_unit != UNIT_NOT_DEFINED ) {
    fprintf (stderr, "Error, in order to use 'output_process per_nm', 'output_process per_cm-1', or 'output_process per_band' \n");
    fprintf (stderr, " you have to specify the unit of the solar spectrum:'\n");
    fprintf (stderr, " 'source solar filename per_nm' for solar spectrums in W/(m^2 nm) or\n");
    fprintf (stderr, " 'source solar filename per_band' for bands in W/(m^2 band).\n");
    return -1;
  }

  /* if representative wavelengths are used, no solar_file is read if none was given by the user */
  if (strlen(filename)==0 && output->wl.use_reptran==1) 
    output->wl.ignore_solar_file=1;

  /* if no extraterrestrial file is defined, check if we really need it */
  if (strlen(filename)==0 && output->wl.use_reptran==0) {
    if (input.source != SRC_LIDAR) {
      switch (input.calibration) {
      case OUTCAL_ABSOLUTE:
        /* in this case we really need input from a file */
        fprintf (stderr, "Error, no extraterrestrial spectrum defined!\n");
        fprintf (stderr, "Please use 'source solar filename' and retry\n");
        return -1;
        break;
  
      case OUTCAL_TRANSMITTANCE:
      case OUTCAL_REFLECTIVITY:
        if (input.processing == PROCESS_SUM || input.processing == PROCESS_INT) {
          /* in this case we really need input from a file */
          fprintf (stderr, "Error, no extraterrestrial spectrum defined!\n");
          fprintf (stderr, "Please use 'source solar filename' and retry\n");
          return -1;
        }
        else{
          /* if we don't have an extraterrestrial spectrum, we set it to 1 */
          if (!input.quiet)
  	  fprintf (stderr, " ... no extraterrestrial spectrum defined, assuming 1\n");
  
          output->wl.ignore_solar_file=1; /* no need for an extraterrestrial file */
        }
        break;
  
      case OUTCAL_THERMAL:
      case OUTCAL_BRIGHTNESS:
        /* if we don't have an extraterrestrial spectrum, we set it to 1 */
        if (!input.quiet)
  	fprintf (stderr, " ... no extraterrestrial spectrum defined, assuming 1\n");
  
        output->wl.ignore_solar_file=1; /* no need for an extraterrestrial file */
        break;
  
      default:
        fprintf (stderr, "Error, unknown output calibration %d\n", 
  	       input.calibration);
        return -1;
      }
    }
    else {
      if (!input.quiet)
      fprintf (stderr, " ... no need for extraterrestrial spectrum, assuming 1\n");
      output->wl.ignore_solar_file=1; /* no need for an extraterrestrial file */
    }
  }

  /* we also don't want a solar file for the "real" k-distributions */
  /* in the thermal spectral range; here everything has been set    */
  /* when the wavelength bands were read; reading the solar_file    */
  /* could spoil that                                               */
  
  if (input.source == SRC_THERMAL) {
    switch (input.ck_scheme) {
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
    case CK_FU:
    case CK_AVHRR_KRATZ:
      output->wl.ignore_solar_file = 1;
      break;
      
    case CK_CRS:
    case CK_RAMAN:
    case CK_REPTRAN:
    case CK_REPTRAN_CHANNEL:
    case CK_LOWTRAN:
    case CK_FILE:
      break;
      
    default:
      fprintf (stderr, "Error: unsupported correlated-k scheme\n");
      return -1;
    }
  }


  wlband_start = 1.e7/output->wl.wvnmhi_r[0];
  wlband_end   = 1.e7/output->wl.wvnmlo_r[output->wl.nlambda_r-1];

  if (strlen(filename)>0 && output->wl.ignore_solar_file && !input.quiet)
    fprintf (stderr, " ... ignoring extraterrestrial spectrum %s\n", filename);

  if (!output->wl.ignore_solar_file) {
    /**** Read extraterrestrial irradiance ****/
    status = read_extraterrestrial(filename, 
				   &input.wl, &output->wl, 
				   output->wl.start, output->wl.end, input.ck_scheme,
                                   wlband_start, wlband_end);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }

    /**** Allocate memory for wavelength bands ****/
    output->wl.wvnmlo_h = calloc (output->wl.nlambda_h, sizeof(float));
    output->wl.wvnmhi_h = calloc (output->wl.nlambda_h, sizeof(float));

    /**** Set to default bandwidth, not very useful ****/
    if (output->bandwidth_unit == UNIT_PER_CM_1) {
      for (iv=0; iv<output->wl.nlambda_h; iv++) {
        output->wl.wvnmlo_h[iv] = 1.0E7 / output->wl.lambda_h[iv] - input.bandwidth / 2.0;
        output->wl.wvnmhi_h[iv] = output->wl.wvnmlo_h[iv] + input.bandwidth;
      }
    }
    else if (output->bandwidth_unit == UNIT_PER_NM) {
      for (iv=0; iv<output->wl.nlambda_h; iv++) {
        output->wl.wvnmlo_h[iv] = 1.0E7 / (output->wl.lambda_h[iv] + input.bandwidth / 2.0);
        output->wl.wvnmhi_h[iv] = 1.0E7 / (output->wl.lambda_h[iv] - input.bandwidth / 2.0);
      }
    }
    else {
      fprintf (stderr, "Error, unsupported bandwidth_unit %d in %s (%s)\n", output->bandwidth_unit, function_name, file_name);
      return -1;
    }
  }
  else {
    status = cp_lambdat_lambdah (output);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by cp_lambdat_lambdah()\n", status);
      return status;
    }
    /* in case of representative wavelengths, the extraterrestrial flux is calculated from fluxes stored in the parameterization file (output->wl.extra_reptran_r) */
    if(output->wl.use_reptran){  
      for (i_h=0; i_h<output->wl.nlambda_h; i_h++) {
        i_band = output->wl.reptran_band_t[output->wl.map_e2h[i_h]];  
        output->wl.fbeam[i_h] = 0;
        weight_sum = 0;
        for(i=0; i<output->wl.nlambda_in_reptran_band[i_band]; i++){
           output->wl.fbeam[i_h] += output->wl.weight_reptran_band[i_band][i] * output->wl.extra_reptran_r[output->wl.reptran_band[i_band][i]];
           weight_sum += output->wl.weight_reptran_band[i_band][i];
        }
        output->wl.fbeam[i_h] = output->wl.fbeam[i_h] / weight_sum;
      }
    }
  }

  /* allocate memory for filter function in any case */
  output->wl.filter = calloc (output->wl.nlambda_h, sizeof(float));
  
  /* set filter function to 1 (default) */
  for (iv=0; iv<output->wl.nlambda_h; iv++)
    output->wl.filter[iv] = 1;


  /**** Read filter function ****/
  if (strlen(input.filename[FN_FILTERFUNCTION])>0) {

    if (!input.quiet)
      fprintf (stderr, " ... reading filter function from  %s\n",
	       input.filename[FN_FILTERFUNCTION]);
    
    status = read_2c_file_float (input.filename[FN_FILTERFUNCTION], &x_flt, &y_flt, &n_flt); 
    if (status!=0) {
      fprintf (stderr, "Error %d reading filter function from %s\n", 
	       status, input.filename[FN_FILTERFUNCTION]);
      return status;
    }

    status = arb_wvn (n_flt, x_flt, y_flt,  
		      output->wl.nlambda_h, output->wl.lambda_h, output->wl.filter, 1, 0);

    if (status!=0) {
      fprintf (stderr, "Error %d interpolating filter function\n", status);
      return status;
    }
  }
  
  /* Normalize filter function such that its integral over wavelength is one */
  if(input.filter_function_normalize){
    for(iv=1; iv<output->wl.nlambda_h; iv++)
      filter_sum += 0.5*(output->wl.filter[iv]+output->wl.filter[iv-1])*(output->wl.lambda_h[iv]-output->wl.lambda_h[iv-1]);
    for(iv=0; iv<output->wl.nlambda_h; iv++)
      output->wl.filter[iv] = output->wl.filter[iv]/filter_sum;
  }

  /**** Correction factor for Earth-Sun distance ****/
  output->eccentricity_factor = 1.0;
  
  if (input.UTC.tm_yday > 0.0) {
    output->eccentricity_factor = eccentricity (input.UTC.tm_yday);
    for (iv=0; iv<output->wl.nlambda_h; iv++)
      output->wl.fbeam[iv] *= output->eccentricity_factor;
  }

  /* stupidity check */
  if (strlen(filename) != 0) {
    /* none of the specified wavelength inside the range of the solar spectrum */
    if (output->wl.nlambda_h == 0 && !input.quiet) {
      fprintf (stderr, "\n*** Warning, the solar spectrum does not include a value inside the specified\n");
      fprintf (stderr,   "*** wavelength range. Are you sure, that you use the correct solar spectrum?\n\n");
      fprintf (stderr,   "*** uvspec will probably not generate any output!\n\n");
    }
  }

  /*  BCA: Is this the appropriate spot for sun radius calculation? BR */
  if (!input.rte.mc.sun_radius && strlen(input.rte.mc.filename[FN_MC_SUNSHAPE_FILE]) >0) {
    /* Angular sun radius (0.26648??) as derived from mean sun radius as =6.96000e8 meter (mean radius of the */
    /* sun (visibe disk) in  Liou, Demtr??der Expphysik 1, Klett Formelsammlung */
    input.rte.mc.sun_radius = 0.26648 * sqrt(eccentricity (input.UTC.tm_yday));
  }
  output->mc.sample.sun_radius = input.rte.mc.sun_radius;

  return 0;
}



/**************************************************************/
/* Read extraterrestrial irradiance.                          */
/**************************************************************/

static int read_extraterrestrial (char *filename, 
				  wl_inp_struct *wl_inp, wl_out_struct *wl_out,
				  float wl_start, float wl_end, int ck_scheme,
                                  float wlband_start, float wlband_end) 
{
  int iv=0, nlambda=0, status=0, have_wl_start=0, have_wl_end=0;
  float *tmp_lambda=NULL, *tmp_fbeam=NULL;

  int is_cdf=0;

#if HAVE_LIBNETCDF
  double *wvnlo=NULL, *wvnhi=NULL, *wvl=NULL;
  int id_wvl=0, id_wvnlo=0, id_wvnhi=0;
  int ncid=0;
  int idd_nwvn=0;
  int id_extra=0;
    
  size_t n=0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  /* test if the file is as netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) 
    is_cdf=1;
#endif

  /* read extraterrestrial irradiance */
  if (is_cdf==0) {
    status = read_2c_file_float (filename, 
				 &tmp_lambda, &tmp_fbeam, &(wl_out->nlambda_h));
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
  }
  else {

#if HAVE_LIBNETCDF
    
    /* get dimension id for "nwvn" */
    status = nc_inq_dimid (ncid, "nwvn", &idd_nwvn);
    if (status!=NC_NOERR) 
      status = nc_inq_dimid (ncid, "nwvl", &idd_nwvn);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
  
    /* get dimension length for "nwvn" */
    status = nc_inq_dimlen (ncid, idd_nwvn, &n);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }

    wl_out->nlambda_h=n;

    start[0] = 0;
    count[0] = wl_out->nlambda_h;
    
    /* allocate memory for fields */
    tmp_fbeam  = calloc (wl_out->nlambda_h, sizeof(float));
    tmp_lambda = calloc (wl_out->nlambda_h, sizeof(float));
    wvnlo      = calloc (wl_out->nlambda_h, sizeof(double));
    wvnhi      = calloc (wl_out->nlambda_h, sizeof(double));
    wvl        = calloc (wl_out->nlambda_h, sizeof(double));
    

    /* read wavelengths from wvl; if not available, use wvnmlo and wvnmhi instead */

    /* get variable id for "wvl" */
    status = nc_inq_varid (ncid, "wvl", &id_wvl);
    
    if (status==NC_NOERR) {  /* read "wvl" */
      status = nc_get_vara_double (ncid, id_wvl, start, count, wvl);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading %s\n", status, filename);
	return status;
      }

      for (iv=0; iv<wl_out->nlambda_h; iv++)
	tmp_lambda[iv] = wvl[iv];

      free (wvl);
    }
    else {
      
      /* get variable id for "wvnlo" */
      status = nc_inq_varid (ncid, "wvnlo", &id_wvnlo);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading %s\n", status, filename);
	return status;
      }
      
      /* read "wvnlo" */
      status = nc_get_vara_double (ncid, id_wvnlo, start, count, wvnlo);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading %s\n", status, filename);
	return status;
      }
      
      /* get variable id for "wvnhi" */
      status = nc_inq_varid (ncid, "wvnhi", &id_wvnhi);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading %s\n", status, filename);
	return status;
      }
      
      /* read "wvnhi" */
      status = nc_get_vara_double (ncid, id_wvnhi, start, count, wvnhi);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading %s\n", status, filename);
	return status;
      }
      
      for (iv=0; iv<wl_out->nlambda_h; iv++)
	tmp_lambda[iv] = 1.0E7 / (0.5 * (wvnlo[iv] + wvnhi[iv]));

      free (wvnlo);
      free (wvnhi);
    }

    /* get variable id for "extra" */
    status = nc_inq_varid (ncid, "extra", &id_extra);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading extraterrestrial irradiance from %s\n", 
	       status, filename);
      return status;
    }
    
    /* read "extra" */
    status = nc_get_vara_float (ncid, id_extra, start, count, tmp_fbeam);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }

    nc_close(ncid);

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use the generic correlated_k option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif

  }

  /* count the extraterrestrial wavelengths within the */
  /* user-defined wavelength range                     */
  nlambda = 0;
  for (iv=0; iv<wl_out->nlambda_h; iv++)
    if ((tmp_lambda[iv]>=wl_start) && (tmp_lambda[iv]<=wl_end))
      nlambda++;

  wl_out->lambda_h = (float *) calloc (nlambda, sizeof(float));
  wl_out->fbeam    = (float *) calloc (nlambda, sizeof(float));

  wl_out->map_e2h  = (int *)   calloc (nlambda, sizeof(int));


  /* Check if the extraterrestrial spectrum covers */
  /* the required wavelength range  */          
  
  if (ck_scheme == CK_CRS || ck_scheme == CK_RAMAN || ck_scheme == CK_REPTRAN || ck_scheme == CK_REPTRAN_CHANNEL) { 
    if (wl_start != -999.0 && wl_start < tmp_lambda[0]) {
      fprintf (stderr, "Error: required start wavelength: %.3f nm outside extraterrestrial\n", wl_start);
      fprintf (stderr, "       spectrum range [%.3f - %.3f]\n", 
	       tmp_lambda[0], tmp_lambda[wl_out->nlambda_h-1]);
      return -1;
    }
    if (wl_end != -999.0 && tmp_lambda[wl_out->nlambda_h-1] < wl_end) {
      fprintf (stderr, "Error: required end wavelength: %.3f nm outside extraterrestrial\n", wl_end);
      fprintf (stderr, "       spectrum range [%.3f - %.3f]\n", 
	       tmp_lambda[0], tmp_lambda[wl_out->nlambda_h-1]);
      return -1;
    }
  }
  else {
    if (wl_start != -999.0 && wl_start < wlband_start) {
      fprintf (stderr, "Error: required start wavelength: %.3f nm outside extraterrestrial\n", wl_start);
      fprintf (stderr, "       spectrum range [%.3f - %.3f]\n", 
	       wlband_start, wlband_end);
      return -1;
    }
    if (wl_end  != -999.0 && wlband_end < wl_end) {
      fprintf (stderr, "Error: required end wavelength: %.3f nm outside extraterrestrial\n", wl_end);
      fprintf (stderr, "       spectrum range [%.3f - %.3f]\n", 
	       wlband_start, wlband_end);
      return -1;
    }
  }


  if (ck_scheme == CK_CRS || ck_scheme == CK_RAMAN || ck_scheme == CK_REPTRAN || ck_scheme == CK_REPTRAN_CHANNEL) { 

    /* Both wl_end and wl_start should be in the extraterrestrial solar file         */
    /* Otherwise problems may occurr below if a single wavelength which is not       */
    /* part of the solar flux file is specified with the wvn option. The user should */
    /* use the spline option to get arbitrary wavelengths or include the wanted      */
    /* wavelength in the extraterrestrial solar flux file.                           */
    for (iv=0; iv<wl_out->nlambda_h; iv++) {
      if (wl_start == tmp_lambda[iv]) 
	have_wl_start=1;
      if (wl_end   == tmp_lambda[iv]) 
	have_wl_end=1;
    }
    if (have_wl_start==0) {
      fprintf (stderr, "Error: required start wavelength: %.6f nm is not in extraterrestrial solar flux file.\n", 
	       wl_start);
      fprintf (stderr, "       Either include wavelength in solar flux file, or\n");
      fprintf (stderr, "       or use spline option.\n");      
      return -1;
    }
    if (have_wl_end==0) {
      fprintf (stderr, "Error: required end wavelength: %.6f nm is not in extraterrestrial solar flux file.\n", 
	       wl_end);
      fprintf (stderr, "       Either include wavelength in solar flux file, or\n");
      fprintf (stderr, "       or use spline option.\n");      
      return -1;
    }
  }


  nlambda = 0;

  for (iv=0; iv<wl_out->nlambda_h; iv++)
    if ((tmp_lambda[iv]>=wl_start) && (tmp_lambda[iv]<=wl_end)) {
      wl_out->lambda_h[nlambda] = tmp_lambda[iv];
      wl_out->fbeam[nlambda]    = tmp_fbeam[iv];
      wl_out->map_e2h[nlambda]  = iv;
      nlambda++;
    }

  wl_out->nlambda_h = nlambda;

  /* free memory */
  free(tmp_lambda);
  free(tmp_fbeam);

  return 0;
}



/* copy the transmission wavelength grid */
/* to the output grid                    */
static int cp_lambdat_lambdah (output_struct *output) 
{
  int iv=0, nlambda=0;

  output->wl.nlambda_h = output->wl.nlambda_t;
  output->wl.lambda_h  = calloc (output->wl.nlambda_h, sizeof(float));
  output->wl.fbeam     = calloc (output->wl.nlambda_h, sizeof(float));
  output->wl.map_e2h   = calloc (output->wl.nlambda_h, sizeof(int));

  output->wl.wvnmlo_h = calloc (output->wl.nlambda_h, sizeof(float));
  output->wl.wvnmhi_h = calloc (output->wl.nlambda_h, sizeof(float));
  
  for (iv=0; iv<output->wl.nlambda_h; iv++)
    if ((output->wl.use_reptran || ((output->wl.lambda_t[iv]>=output->wl.start) && (output->wl.lambda_t[iv]<=output->wl.end)))) {
      output->wl.lambda_h[nlambda] = output->wl.lambda_t[iv];
      output->wl.wvnmlo_h[nlambda] = output->wl.wvnmlo_t[iv];
      output->wl.wvnmhi_h[nlambda] = output->wl.wvnmhi_t[iv];
      output->wl.fbeam[nlambda]    = 1.0;
      output->wl.map_e2h[nlambda]  = iv;
      nlambda++;
    }
  
  output->wl.nlambda_h = nlambda;
  
  return 0;
}
