/*--------------------------------------------------------------------
 * $Id: cloud.c 3276 2017-07-04 14:16:36Z Claudia.Emde $
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
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <glob.h>
#include <float.h>

#include "uvspec.h"
#include "ckdfu.h"
#include "cloud.h"
#include "ipa.h"
#include "ascii.h"
#include "numeric.h"
#include "miecalc.h"
#include "f77-uscore.h"
#include "wcloud3d.h"
#include "errors.h"

#if HAVE_KEY56
#include "yang56.h"
#endif
 
#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#define EPSILON 1E-6

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/************************************/
/* prototypes of internal functions */
/************************************/

int read_cloud_fraction_file ( char           *filename,
			       atm_inp_struct  atm_in,
			       atm_out_struct *atm_out,
			       cf_out_struct  *cf,
			       int             quiet );

int read_ECMWF_clouds ( char            *filename,
			char            *file_ic_reff,
			float            lat,
			float            lon,
			struct tm        UTC,
			int              time_interpolate, 
			char            *cloud, 
			int              reff_prop,
			int              ic_prop,
			float            reff_fixed,
			float            altitude,
			float           *press_atm,
			int              cloud_overlap, 
			int              verbose,
			int              quiet,
			int             *rows,
			int             *max_columns,
			int             *min_columns, 
			float         ***cld_data,
			cf_out_struct   *cf );

int read_ECHAM_clouds ( char            *filename,
			float            lat,
			float            lon,
			struct tm        UTC,
			int              time_interpolate, 
			char            *cloud,
			float           *zd,
			int              verbose,
			int              quiet,
			int             *rows,
			int             *max_columns,
			int             *min_columns,
			float         ***cld_data,
			cf_out_struct   *cf );

static int setup_cloudfraction ( int                cloud_overlap,
				 char              *caoth_name,
				 char              *caoth_fullname,
				 caoth_out_struct   cld,
				 int                nlambda,
				 int                cld_source,
				 int                verbose,
				 int                quiet,
				 /* in/output */
                                 int               *nipa,
				 float             *cloudcover,
				 cf_out_struct     *cf,
				 /* output */
				 cf_out_struct    **cfipa,
				 caoth_out_struct **cldipa,
				 float            **ipaweight );

static int get_total_column_caoth ( caoth_out_struct *cld,
				    int               nipa,
				    float            *ipaweight,
				    caoth_out_struct *cldipa, 
				    cf_out_struct     cf,
				    int               verbose );


int ic_geom (float wavelength, int nlyr,
	     float *iwc, float *reff, float *zd,
	     float *dtau, float *ssa);

int read_raytracing_file (char *filename, crystal_prop_struct **raytracing_prop, int *n_raytracing_prop, int quiet);


/***********************************************************************************/
/* Function: setup_all_caoth                                              @62_30i@ */
/* Description:                                                                    */
/*  Run loop over all caoth and set them up. Also get cloud fraction.              */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int setup_all_caoth (input_struct input, output_struct *output)
{
  int status=0, isp=0;
  int iscloudcover=0;

  if ( strlen(input.filename[FN_CLOUD_FRACTION])>0 ) {
    status = read_cloud_fraction_file (input.filename[FN_CLOUD_FRACTION],
				       input.atm,
				       &(output->atm),
				       &(output->cf),
				       input.quiet);
    if (status)
      return fct_err_out ( status, "read_cloud_fraction_file", ERROR_POSITION );
  }

  output->caoth = calloc ((size_t) input.n_caoth, sizeof(caoth_out_struct));
  if (output->caoth==NULL)
    return mem_err_out ( "output->caoth", ERROR_POSITION );

  output->caoth_ipa = calloc ((size_t) input.n_caoth, sizeof(caoth_out_struct *));
  if (output->caoth_ipa==NULL)
    return mem_err_out ( "output->caoth_ipa", ERROR_POSITION );

  for (isp=0; isp<input.n_caoth; isp++) {
    status = setup_caoth ( input,
			   output,
			   input.caoth[isp],
			   &(output->caoth[isp]),
			   &(output->caoth_ipa[isp]) );
    if (status)
      return fct_err_out ( status, "setup_caoth", ERROR_POSITION );
  }

  /* give default value to total cloud fraction, even if we do not */
  /*   have any overlap assumption                                 */
  /* Attention: this must be after setup_cloudfraction is called   */
  /*            as a not initialzied cf.cf_total is used there as  */
  /*		marker for "cf not calculated"                     */

  if (input.n_caoth==0)
    output->cf.cf_total = 0.0;

  if (input.cloud_overlap == CLOUD_OVERLAP_OFF && input.rte.solver!=SOLVER_TWOMAXRND) {
    output->cf.cf_total = 0.0;
    for (isp=0; isp<input.n_caoth; isp++)
      if ( input.caoth[isp].cloudcover != NOT_DEFINED_INTEGER ) {
	iscloudcover=1;
	if ( input.caoth[isp].cloudcover > output->cf.cf_total )
	  output->cf.cf_total = input.caoth[isp].cloudcover;
      }

    if (iscloudcover==0)
      output->cf.cf_total = 1.0;

    if (input.verbose)
      fprintf(stderr,"     total cloud cover = %f\n", output->cf.cf_total);
  }

  return 0;    
}


/***********************************************************************************/
/* Function: read_cloud_fraction_file                                     @62_30i@ */
/* Description:                                                                    */
/*  Read get cloud fraction file.                                                  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int read_cloud_fraction_file ( char           *filename,
			       atm_inp_struct  atm_in,
			       atm_out_struct *atm_out,
			       cf_out_struct  *cf,
			       int             quiet )
{
  int lc=0;
  float *z_atm=NULL;
  int nlev_atm = 0;
  int status=0;

  /* read cloud fraction file if specified */
  
  if (strlen(filename)>0) {
    status = read_2c_file_float (filename,
				 &(cf->zd),
				 &(cf->cf),
				 &(cf->nlev));
    if (status!=0) {
      fprintf (stderr, "Error %d reading cloud fraction file %s\n", 
	       status, filename);
      return status;
    }
  }

  /* need to shift by one layer as cloud fractions are automatically
     considered layer properties */
  for (lc=0;lc<cf->nlev-1;lc++)
    cf->cf[lc] = cf->cf[lc+1];

  if (!quiet) {
    fprintf (stderr, " ... read %d layers from cloud fraction file %s:\n",
	     cf->nlev, filename);
  }

  /* copy z_grid of the atmosphere file */
  if ((z_atm = calloc(atm_out->nlev, sizeof(float)))==NULL) {
    fprintf (stderr, "Error: Allocation of z_atm (line %d, function %s in %s) \n",
	     __LINE__, __func__, __FILE__);
    return -1;
  }

  for (lc=0; lc<atm_out->nlev; lc++)
    z_atm[lc] = atm_out->zd[lc];

  nlev_atm = atm_out->nlev;

  /* add level to common z grid */
  set_common_z_grid (cf->zd, cf->nlev, 
		     &(atm_out->zd), &(atm_out->nlev));


  /*???? overkill? should be in redistribute*/
  /* interpolate the atmospheric profiles to the new z-grid */
  status = interpolate_atmosphere (z_atm,
				   &(atm_out->microphys.press),
				   &(atm_out->microphys.temper),
				   &(atm_out->microphys.dens),
				   &(atm_out->microphys.temper_avg),
				   &(atm_out->microphys.dens_avg),
				   nlev_atm,
				   atm_out->zd,
				   atm_out->nlev,
				   atm_in.interpol_method_press,
				   atm_in.interpol_method_temper,
				   atm_in.interpol_method_gas,
				   quiet);

  if (status != 0) {
    fprintf (stderr, "Error %d interpolating atmosphere (line %d, function %s in %s)\n",
	     status, __LINE__, __func__, __FILE__);
    return status; 
  }

  free(z_atm);

  return 0;
}

/***********************************************************************************/
/* Function: setup_caoth                                                  @62_30i@ */
/* Description:                                                                    */
/*  Setup single caoth (Cloud, Aerosol, or Other THing).                           */
/*  This function read the cloud/aerosol 1D profile and converts them to           */
/*  optical properties. Also, the cloud fraction is calculated.                    */
/* Parameters: input, output, input_caoth                                          */
/* Return values: output_caoth, output_caoth_ipa                                   */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int setup_caoth ( /* input */
		  input_struct       input,
		  output_struct     *output,
		  caoth_inp_struct   input_caoth,
		  /* output */
		  caoth_out_struct  *output_caoth,
		  caoth_out_struct **output_caoth_ipa )
{
  int status=0, lc=0, iv=0, ipa=0;
  int rows=0, min_columns=0, max_columns=0, max_length=0;

  float ipasum=0.0;

  /* if caoth microphysical properties profile specified */
  switch (input_caoth.source) {
  /* read caoth from one file */
  case CAOTH_FROM_1D:
  case CAOTH_FROM_ECMWF:
  case CAOTH_FROM_ECHAM: 
    
    /* Read microphysical properties      */
    /* and convert to optical properties. */

    if (input.verbose)
      fprintf (stderr, "     reading %s properties from %s\n", 
	       input_caoth.name, input_caoth.filename);

    /* Read caoth microphysical properties and apply */
    /* default conversion to optical properties      */
    status = read_and_convert_caoth_file (input,
					  input_caoth.filename,
					  output->wl,
					  output->alt.altitude,
					  output->atm,
					  input_caoth,
					  output_caoth,
					  &(output->cf),
					  &(output->atm.nmom) );
    if (status)
      return fct_err_out ( status, "read_and_convert_caoth_file", ERROR_POSITION );

    /* Copy cloudcover to output variable ... */
    output_caoth->cloudcover = input_caoth.cloudcover;

    /* ... and modify it if cloud_overlap assumption is used. */
    status = setup_cloudfraction ( input.cloud_overlap,
				   input_caoth.name,
				   input_caoth.fullname,
				   *output_caoth,
				   output->wl.nlambda_r, 
				   input_caoth.source,
				   input.verbose,
				   input.quiet,
				   &(output->nipa),
				   &(output_caoth->cloudcover),
				   &(output->cf),
				   &(output->cfipa),
				   output_caoth_ipa,
				   &(output->ipaweight) );
    if (status)
      return fct_err_out ( status, "setup_cloudfraction", ERROR_POSITION );

    break;

  case CAOTH_FROM_IPA_FILES:
    /* Read caoth files using caoth_ipa_files option */
    status = read_caoth_ipa ( input,
			      input_caoth,
			      output->alt.altitude,
			      output->wl,
			      output->atm,
			      &(output->nipa),
			      &(output->ipaweight),
			      output_caoth,
			      output_caoth_ipa,
			      &(output->cf) );
    if (status)
      return fct_err_out ( status, "read_caoth_ipa", ERROR_POSITION );

    break;

  case CAOTH_FROM_MOMENTS:
  
    /* Read file holding filenames of layer caoth optical properties */    

    if (input_caoth.properties != PROP_EXPLICIT) {
      fprintf (stderr, "Error, inconsistent properties for %s ",
	       input_caoth.fullname);
      fprintf (stderr, "(%s_files%s, but not explicit properties).\n",
	       input_caoth.name1, input_caoth.name2);
      fprintf (stderr, " This is not supposed to happen and indicates a coding error. ");
      fprintf (stderr, "(line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    /* just to make sure, check if filename has been defined */
    if (strlen(input_caoth.filename) == 0) {
      fprintf (stderr, "Error, inconsistent properties %s ",
	       input_caoth.fullname);
      fprintf (stderr, "(layer files, but no filename).\n");
      fprintf (stderr, " This is not supposed to happen and indicates a coding error. ");
      fprintf (stderr, "(line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    if (output->alt.altitude != 0) {
      fprintf (stderr, "Error, option 'altitude' currently not implemented in combination\n");
      fprintf (stderr, " with '%s_files' (line %d, function %s in %s, %s) \n",
	       input_caoth.name1, __LINE__, __func__, __FILE__,
	       input_caoth.fullname);
      return -1;
    }

    /* determine number of levels; required for memory allocation */
    if ((status = ASCII_checkfile (input_caoth.filename, 
				   &rows, &min_columns, &max_columns, &max_length)) != 0) {
      fprintf (stderr, "Error %d reading %s(line %d, function %s in %s) \n", 
	       status, input_caoth.filename, __LINE__, __func__, __FILE__);
      return status;
    }

    /* allocate memory */
    status = calloc_caoth_out ( output_caoth,
				input_caoth.name,
				input_caoth.fullname,
				output->wl.nlambda_r,
				rows,
				1,
				1 );
    if (status)
      return mem_err_out ( "output_caoth struct ", ERROR_POSITION );

    /* read optical properties profile, */
    /* allocate memory for moments      */

    status = read_optprop_files ( input_caoth.filename,
				  output->wl.lambda_r,
				  output->wl.nlambda_r, 
				  output->wl.nlambda_rte_lower,
				  output->wl.nlambda_rte_upper,
				  output_caoth->zd,
				  output_caoth->nlev,
				  output_caoth->optprop.dtau, 
				  output_caoth->optprop.ssa,
				  output_caoth->optprop.g1,
				  output_caoth->optprop.moment,
				  output_caoth->optprop.nmom,
				  input.verbose );

    if (status != 0) {
      fprintf (stderr, "Error %d reading %s optical property files '%s' ",
	       status, input_caoth.name, input_caoth.filename);
      fprintf (stderr, "(line %d, function %s in %s)\n", 
	       __LINE__, __func__, __FILE__);
      return status;
    }

    /* increase number of moments if more found */
    for (iv=0; iv<output->wl.nlambda_r; iv++)
      for (lc=0; lc<output_caoth->nlev-1; lc++)
	if (output_caoth->optprop.nmom[iv][lc]-1 > output->atm.nmom)
	  output->atm.nmom = output_caoth->optprop.nmom[iv][lc]-1;

    break;

  case CAOTH_FROM_3D:

    /* allocate memory, even if no caoth is defined; */
    /* it may be required by other functions;          */
    
    output_caoth->nlev = 2;
    output_caoth->optprop.nphamat = 1;
    status = calloc_caoth_out ( output_caoth,
				input_caoth.name,
				input_caoth.fullname,
				output->wl.nlambda_r,
				output_caoth->nlev,
				output_caoth->optprop.nphamat,
				1 );
    if (status)
      return mem_err_out ( "output_caoth struct ", ERROR_POSITION );
    
    output_caoth->zd[0] = output->atm.zd[0];
    output_caoth->zd[1] = output->atm.zd[output->atm.nlev-1];

    /* if ipa structure is needed, declare dummy ipa caoth */
    if (input_caoth.ipa) {

      /* allocate output_caoth_ipa */
      if (( *output_caoth_ipa = calloc (output->nipa, sizeof(caoth_out_struct)) ) == NULL) {
        fprintf (stderr, "Error, allocating dummy output ipa structure ");
	fprintf (stderr, "(line %d, function %s in %s) \n", __LINE__, __func__, __FILE__);
        return -1;
      }

      for (ipa=0; ipa<output->nipa; ipa++) {
        /* allocate caoth ipa */
        status = calloc_caoth_out ( &((*output_caoth_ipa)[ipa]),
				    input_caoth.name,
				    input_caoth.fullname,
				    output->wl.nlambda_r,
				    output_caoth->nlev,
				    output_caoth->optprop.nphamat,
				    1 );
	if (status)
	  return mem_err_out ("output_caoth_ipa struct ", ERROR_POSITION );

        /* copy dummy wc to dummy wcipa */
        status = cp_caoth_out ( &((*output_caoth_ipa)[ipa]),
				*output_caoth,
				FALSE,
				TRUE,
				input.quiet );
        if (status!=0) {
          fprintf (stderr, "Error %d copying %s to dummy output ipa structure %3d ",
		   status, input_caoth.fullname, ipa );
	  fprintf (stderr, "(line %d, function %s in %s)\n", 
		   __LINE__, __func__, __FILE__ );
          return status;
        }
      }
    }
    break;

  default:
    fprintf (stderr, "Error, unknown %s source %d (line %d, function %s in %s)\n", 
	     input_caoth.fullname, input_caoth.source, __LINE__, __func__, __FILE__);
    return -1;
  }

  /* overwrite these properties with user-defined optical thickness, ssa, etc */
  status = apply_user_defined_properties_to_caoth ( input_caoth,
						    output->wl.nlambda_r,
						    output->wl.lambda_r,
						    output->alt.altitude,
						    output_caoth);
  if (status)
    return fct_err_out ( status, "apply_user_defined_properties_to_caoth", ERROR_POSITION );

  /* Check that ipa weights add up to one */
  if (input_caoth.ipa)  {
    ipasum=0.0;
    for (ipa=0; ipa<output->nipa; ipa++)
      ipasum += output->ipaweight[ipa];
    if (ipasum < 1.0-EPSILON || ipasum > 1.0+EPSILON) {
      fprintf (stderr, "Error, ipa weights add up to %f. But the sum should be 1.0000 ! ",
	       ipasum);
      fprintf (stderr, "(line %d, function %s in %s)\n",
	       __LINE__, __func__, __FILE__);
      return -1;
    }
  }

  /* calculate: total column caoth */
  status = get_total_column_caoth ( output_caoth,
				    output->nipa,
				    output->ipaweight,
				    *output_caoth_ipa,
				    output->cf,
				    input.verbose );
  if (status)
    return fct_err_out ( status, "get_total_column_caoth", ERROR_POSITION );
  if (input.verbose)
    fprintf (stderr, "     total column %s = %e kg /m2 \n",
	     input_caoth.fullname, output_caoth->microphys.tcw );

  return 0;
}


/***********************************************************************************/
/* Function: setup_cloudfraction                                                   */
/* Description:                                                                    */
/*  calculates the column (total) cloud fraction (cf)                              */
/*  scales the cloud fraction in the cloudy column, accordingly to the total cf    */
/*                                                                                 */
/* Parameters:                                                                     */
/*   input                                                                         */
/*     int cloud_overlap                                                           */
/*     verbose                                                                     */
/*   output                                                                        */
/*     cf_out_struct *cf                                                           */
/*     float *cloudcover                                                           */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    cloud.c                                                               */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Jun 2008   U. Hamann     Created                                             */
/*    Apr 2009   U. Hamann     Move scaling of cloud water content                 */
/*                             from ECMWF cloud into this function                 */
/*                                                                                 */
/***********************************************************************************/

static int setup_cloudfraction ( int                  cloud_overlap,
				 char                *caoth_name,
				 char                *caoth_fullname,
				 caoth_out_struct     cld,
				 int                  nlambda,
				 int                  cld_source,
				 int                  verbose,
				 int                  quiet,
				 /* in/output */
                                 int                 *nipa,
				 float               *cloudcover,
				 cf_out_struct       *cf,
				 /* output */
				 cf_out_struct      **cfipa,
				 caoth_out_struct   **cldipa,
				 float              **ipaweight )
{

  int status=0;
  int lc=0;
  float cloudfree = NOT_DEFINED_FLOAT;
  float cc_max = 0.0;
  int ipa=0;
  static int done = FALSE;   /* calculate TCC only once and especially scale cf only once !!! */
  int iclear=1, icloud=0;    /* marker for cloudy and clear column, for ipa calculations      */

  if ( cloud_overlap != CLOUD_OVERLAP_OFF ) { 

    if (verbose) 
      fprintf(stderr, " ... %s \n", __func__);

    /* check if there is a cloud fraction profile */
    if (cf->nlev==0) {
      fprintf(stderr, "Error, no cloudfraction levels given (line %d, function %s in %s)\n",
                       __LINE__, __func__, __FILE__ );
      return -1;
    }
  }
      
  /* calculate TOTAL CLOUD COVER */
  switch (cloud_overlap) {
  case CLOUD_OVERLAP_MAX:
    /* find maximum cloud fraction */
    cf->cf_total = 0.0;
    for (lc=0;lc <cf->nlev-1;lc++) {
      if ( cf->cf[lc] > cf->cf_total )
        cf->cf_total = cf->cf[lc];
    }
    break;
  case CLOUD_OVERLAP_MAXRAND:
    cloudfree = 1.0;
    cc_max = cf->cf[0]; /* initialisation of maximum cloud cover in one cloud */
    for ( lc=1; lc<cf->nlev-1; lc++ ) {
      if ( cf->cf[lc] > cc_max ) cc_max = cf->cf[lc]; /* find cc_max inside one cloud */
      if (( cf->cf[lc] == 0.0 && cf->cf[lc-1] > 0.0 ) || /* end of connected cloud */ 
          ( cf->cf[lc]  > 0.0 && lc == cf->nlev-1-1 )) { /* cloud reach until bottom */ 
        cloudfree = cloudfree * (1-cc_max);
        cc_max=0.0;
      }
    }
    cf->cf_total = 1.0 - cloudfree;
    break;
  case CLOUD_OVERLAP_RAND:
    cloudfree=1.0;
    for (lc=0;lc <cf->nlev-1;lc++) {
      cloudfree *= ( 1.0 - cf->cf[lc] );
    }
    cf->cf_total = 1.0 - cloudfree;
    break;
  case CLOUD_OVERLAP_OFF:
    if ((*cloudcover) == NOT_DEFINED_FLOAT)
      cf->cf_total = 1.0;
    else 
      cf->cf_total = (*cloudcover);
    break;
  default:
    fprintf (stderr, "Error, unknown cloud_overlap assumption %d. (line %d, function %s in %s)\n", 
             cloud_overlap, __LINE__, __func__, __FILE__  );
    return -1;
  }
  
  if (cf->cf_total == NOT_DEFINED_FLOAT) {
    fprintf (stderr, "Error, calculating total cloud cover (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
    return -1;
  }

  if ( verbose && cloud_overlap != CLOUD_OVERLAP_OFF ) 
    fprintf(stderr, "     total cloud fraction is %f \n", cf->cf_total );

  switch (cloud_overlap) {
  case CLOUD_OVERLAP_MAX:
  case CLOUD_OVERLAP_MAXRAND:
    
    /* the column is treated as one cloudy column with cf_total and one cloud free column */
    /* the cloud fraction in the cloudy column is now cf/cf_total                         */
      
    /* ipa weight of the cloud column is equal to the total cloud cover */
    /* this must be done for both water and ice clouds seperately       */
    (*cloudcover) = cf->cf_total;
    

    /* cloudfraction has only to be done once, that's why it might already be done (if we have both water and ice clouds) */
    if (done == FALSE ) {
      if ( (*nipa) != 2 && (*nipa) != 0 ) {
        fprintf(stderr, "Error, cloud overlap max or maxrand only possible with 2 columns. (nipa=%d)\n", (*nipa) );
        return -1;
      }

      (*nipa) = 2;   /* 2 independent pixels */

      /* allocate memory for IPA arrays */
      (*ipaweight) = calloc ((*nipa), sizeof(float));
      (*cfipa)     = calloc ((*nipa), sizeof(cf_out_struct));
      
      (*ipaweight)[iclear] = 1.0-(*cloudcover); /* cloudfree pixel */
      (*ipaweight)[icloud] =     (*cloudcover); /* overcast  pixel */

      if (!quiet)
        fprintf (stderr, " ... IPA with %2d columns, %6.2f%% cloudfree and %6.2f%% overcast\n", 
                 (*nipa), 100.0*(*ipaweight)[iclear], 100.0*(*ipaweight)[icloud]);
            
      /* alloc ipa cloud fraction */
      for (ipa=0; ipa<(*nipa); ipa++) {
        status = alloc_cloud_fraction(&((*cfipa)[ipa]),cf->nlev);
        if (status != 0) {
          fprintf (stderr, "Error, during alloc_cloud_fraction (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__ );
          return -1;
        }
      }

      /* first pixel is cloudfree; only levels need to be set */
      (*cfipa)[iclear].nlev = cf->nlev;
      for (lc=0; lc<(*cfipa)[0].nlev; lc++)
        (*cfipa)[iclear].zd[lc] = cf->zd[lc];
      /* we have clear sky inside the first column */      
      (*cfipa)[iclear].cf_total = 0.0;

      /* second entry [1] is the cloudy column */
      (*cfipa)[icloud].nlev = cf->nlev;
      for (lc=0;lc<(*cfipa)[1].nlev;lc++) {
        (*cfipa)[icloud].cf[lc] = cf->cf[lc];
        (*cfipa)[icloud].zd[lc] = cf->zd[lc];
      }

      /* scale cloud fraction inside cloudy column */
      if ( cf->cf_total != 0.0 )
        for ( lc=0; lc<(*cfipa)[icloud].nlev-1; lc++ ) {
          (*cfipa)[icloud].cf[lc] /= cf->cf_total;
        }
      /* we have full coverage inside the second column */      
      (*cfipa)[icloud].cf_total = 1.0;
      done = TRUE;  /* marker that tcc is calculated and that cf is already scaled to cloudy column */

      if (verbose) {
        fprintf(stderr,"\n#--- total column --  |");
        for (ipa=0; ipa<(*nipa); ipa++)
          fprintf(stderr," -- column %3d -- |", ipa);
        fprintf (stderr, "\n");
        fprintf(stderr,"# lc    zd      cf    |");
        for (ipa=0; ipa<(*nipa); ipa++)
          fprintf(stderr,"    zd       cf   |");
        fprintf (stderr, "\n");
        fprintf(stderr,"#---------------------|");
        for (ipa=0; ipa<(*nipa); ipa++)
          fprintf(stderr,"------------------|");
        fprintf (stderr, "\n");
        fprintf (stderr, " %3d %7.3f %8.6f |", lc, cf->zd[0], 0.0);
        for (ipa=0; ipa<(*nipa); ipa++)
          fprintf (stderr, " %7.3f %8.4f |", (*cfipa)[ipa].zd[0], 0.0 );
        fprintf (stderr, "\n");
        for (lc=1; lc<cf->nlev; lc++) {
          fprintf (stderr, " %3d %7.3f %8.6f |", lc, cf->zd[lc], cf->cf[lc-1]);                /* -1 == layer property */
          for (ipa=0; ipa<(*nipa); ipa++)
            fprintf (stderr, " %7.3f %8.4f |", (*cfipa)[ipa].zd[lc], (*cfipa)[ipa].cf[lc-1] ); /* -1 == layer property */
          fprintf (stderr, "\n");
        }
        fprintf (stderr, "\n");
      }
    }
    else {
      /* check if cloud fractions are identical */
      for ( lc=0; lc<(*cfipa)[icloud].nlev-1; lc++ ) {
        if ( fabs((*cfipa)[icloud].cf[lc] - (cf->cf[lc]/cf->cf_total)) > EPSILON ) {
          fprintf (stderr, "Error, different cloud_fraction for water and ice clouds (line %d, function %s in %s)\n",
                   __LINE__, __func__, __FILE__ );
          fprintf (stderr, "  wc_cf = %f, ic_cf = %f \n", (*cfipa)[icloud].cf[lc], (cf->cf[lc]/cf->cf_total) );
          return -1;
        }
      }
    }
        
    /* allocate */
    (*cldipa)    = calloc ((*nipa), sizeof(caoth_out_struct));

    for (ipa=0; ipa<(*nipa); ipa++) {
      status = calloc_caoth_out ( &((*cldipa)[ipa]),
				  caoth_name,
				  caoth_fullname,
				  nlambda,
				  cld.nlev,
				  cld.optprop.nphamat,
				  1 );
      if (status!=0) {
        fprintf (stderr, "Error %d allocating memory for caoth_out_struct %d (line %d, function %s in %s)\n", 
                 status, ipa, __LINE__, __func__, __FILE__);
        return status;
      }
    }

    /* first pixel is cloudless; only levels need to be set */
    for (lc=0; lc<cld.nlev; lc++)
      (*cldipa)[iclear].zd[lc] = cld.zd[lc]; 
    
    /* copy 1D caoth to 2nd pixel */
    status = cp_caoth_out (&((*cldipa)[icloud]), cld,     TRUE,        TRUE,     quiet);
    if (status!=0) {       /*  target,      source, copy_ssprop, alloc_moment  quiet  */
      fprintf (stderr, "Error %d copying cld to cloudy cldipa (line %d, function %s in %s)\n",
	       status, __LINE__, __func__, __FILE__ );
      return status;
    }
    
    /* scale water content to cloudy column */
    switch (cld_source) {
      /* read caoth from one file */
    case CAOTH_FROM_1D:
    case CAOTH_FROM_ECHAM:
      /* no need to change the cloud water content, the caoth file contains the cloud water of the cloudy column */
      break;
    case CAOTH_FROM_ECMWF:
      /* scale cloud water content to cloudy column */
      for (lc=0; lc<cld.nlev; lc++) {
        (*cldipa)[icloud].microphys.lwc      [lc] /= cf->cf_total;
        (*cldipa)[icloud].microphys.lwc_layer[lc] /= cf->cf_total;
      }
      break;
    case CAOTH_FROM_IPA_FILES:
    case CAOTH_FROM_3D:
    case CAOTH_FROM_MOMENTS:
      fprintf (stderr, "Error, cloud overlap only possible with clouds from file, ECMWF, or ECHAM (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__); 
      fprintf (stderr, "       clouds given from %d (ipa=%d, 3Dclouds=%d, cloudlayer=%d)\n", 
                       cld_source, CAOTH_FROM_IPA_FILES, CAOTH_FROM_3D, CAOTH_FROM_MOMENTS);
      return -1;
      break;
    default:
      fprintf (stderr, "Error, unknown cloud source %d (line %d, function %s in %s)\n", 
               cld_source, __LINE__, __func__, __FILE__);
      return -1;
    }

    if (verbose) {
      fprintf(stderr,"\n#---- total column -----  |");
      for (ipa=0; ipa<(*nipa); ipa++)
        fprintf(stderr," ---- column %3d ---- |", ipa);
      fprintf (stderr, "\n");
      fprintf(stderr,"# lc    zd      cwc       |");
      for (ipa=0; ipa<(*nipa); ipa++)
        fprintf(stderr,"    zd       cwc      |");
      fprintf (stderr, "\n");
      fprintf(stderr,"#-------------------------|");
      for (ipa=0; ipa<(*nipa); ipa++)
        fprintf(stderr,"----------------------|");
      fprintf (stderr, "\n");
      for (lc=0; lc<cld.nlev; lc++) {
        fprintf (stderr, " %3d %7.3f %e |",
                 lc, cld.zd[lc], cld.microphys.lwc[lc]); 
                 /* cld.microphys.lwc_layer[lc-1] == cld.microphys.lwc[lc] for wc_layer */
        for (ipa=0; ipa<(*nipa); ipa++)
          fprintf (stderr, " %7.3f %e |", (*cldipa)[ipa].zd[lc], (*cldipa)[ipa].microphys.lwc[lc]);
        fprintf (stderr, "\n");
      }
      fprintf (stderr, "\n");
    }


    break;
  case CLOUD_OVERLAP_RAND:
    /* For cloud_overlap == CLOUD_OVERLAP_RAND the overlap is simulated                   */
    /* by the effective optical thickness (see optical properties)                        */
  case CLOUD_OVERLAP_OFF:
    /* nothing to do here, only one column considered, */
    if ( (*cloudcover) >= 0) {
      (*nipa)=2;

      if (!quiet)
        fprintf (stderr, " ... IPA with %d columns, cloudless and overcast\n", (*nipa));
      
      /* allocate */
      (*cldipa)    = calloc ((*nipa), sizeof(caoth_out_struct));
      (*ipaweight) = calloc ((*nipa), sizeof(float));

      for (ipa=0; ipa<(*nipa); ipa++) {
        status = calloc_caoth_out ( &((*cldipa)[ipa]),
				    caoth_name,
				    caoth_fullname,
				    nlambda,
				    cld.nlev,
				    cld.optprop.nphamat,
				    1 );
        if (status!=0) {
          fprintf (stderr, "Error %d allocating memory for caoth_out_struct (line %d, function %s in %s)\n", 
                   status, __LINE__, __func__, __FILE__);
          return status;
        }
      }

      (*ipaweight)[iclear] = 1.0-(*cloudcover); /* cloudfree pixel */
      (*ipaweight)[icloud] =     (*cloudcover); /* overcast  pixel */
      
      if (!quiet)
        fprintf (stderr, " ... IPA with %2d columns, %6.2f%% cloudfree and %6.2f%% overcast\n", 
                 (*nipa), 100.0*(*ipaweight)[iclear], 100.0*(*ipaweight)[icloud]);

      /* first pixel is cloudless; only levels need to be set */
      for (lc=0; lc<cld.nlev; lc++)
        (*cldipa)[iclear].zd[lc] = cld.zd[lc]; 

      /* copy 1D cloud to 2nd pixel */
      status = cp_caoth_out (&((*cldipa)[icloud]), cld,      FALSE,        TRUE,      quiet);
      if (status!=0) {      /*  target,    source, copy_ssprop, alloc_moment  quiet  */
        fprintf (stderr, "Error %d copying cld to cloudy cldipa (line %d, function %s in %s)\n",
		 status, __LINE__, __func__, __FILE__ );
        return status;
      }

    }

    break;
  default:
    fprintf (stderr, "Error, unknown cloud_overlap assumption %d. (line %d, function %s in %s)\n", 
             cloud_overlap, __LINE__, __func__, __FILE__  );
    return -1;
  }
  
  return status;
}

static int calloc_caoth_microphys_struct (caoth_microphys_struct *out, int nlev) 
{
  out->lwc  = calloc (nlev, sizeof(float));
  out->effr = calloc (nlev, sizeof(float));

  /* for safety reasons allocate nlev even for layer quantities */
  out->lwc_layer  = calloc (nlev, sizeof(float));
  out->effr_layer = calloc (nlev, sizeof(float));

  out->alloc = 1;

  return 0;
}

/***********************************************************************************/
/* Function: get_total_column_caoth                                                */
/* Description:                                                                    */
/*  calculates the total column (cloud) water (tcw)                                */
/*                                                                                 */
/* Parameters:                                                                     */
/*   input                                                                         */
/*     int cloud_overlap      marker for cloud overlap assumption                  */
/*     cf_out_struct *cf      cloud fraction structure                             */
/*     cld_out_struct *cld    might be water or ice cloud strucutre                */
/*     verbose                                                                     */
/*   output                                                                        */
/*     float *cloudcover                                                           */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    cloud.c                                                               */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Apr 2009   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

static int get_total_column_caoth ( caoth_out_struct *cld,
				    int               nipa,
				    float            *ipaweight,
				    caoth_out_struct *cldipa, 
				    cf_out_struct     cf,
				    int               verbose )
{

  int status = 0;
  int ipa    = NOT_DEFINED_INTEGER;
  int lc     = NOT_DEFINED_INTEGER;

  if (verbose) fprintf(stderr," ... calculate total column ice/liquid water\n");

  /* intergrate total column (cloud) water */

  if (nipa == 1) {
    for (lc=0; lc<cld->nlev-1; lc++) {
      cld->microphys.tcw += cld->microphys.lwc_layer[lc]*(cld->zd[lc]-cld->zd[lc+1]); /* g/m3 * km -> kg/m2 */
      /* fprintf(stderr,"     z=%6.2f, dz=%6.4f, CWC=%e \n", cld->zd[lc], cld->zd[lc]-cld->zd[lc+1], cld->microphys.lwc_layer[lc] ); */
    }
  }
  else {
    for (ipa=0; ipa<nipa; ipa++) {
      /* fprintf(stderr,"     ipa column %3d of %3d, weight: %f, %d \n", ipa, nipa, ipaweight[ipa], cldipa[ipa].nlev); */
      for (lc=0; lc<cldipa[ipa].nlev-1; lc++)
        cldipa[ipa].microphys.tcw += cldipa[ipa].microphys.lwc_layer[lc] * ( cldipa[ipa].zd[lc] - cldipa[ipa].zd[lc+1] );
      cld->microphys.tcw += ipaweight[ipa] * cldipa[ipa].microphys.tcw;
      /* fprintf(stderr,"     CWC=%e, z=%f\n", cldipa[ipa].microphys.lwc_layer[lc], cldipa[ipa].zd[lc]); */
    }
  }

  return status;
}


static int free_caoth_microphys_struct (caoth_microphys_struct *out) 
{
  free (out->lwc_layer); 
  free (out->effr_layer);

  free (out->lwc); 
  free (out->effr);
  
  return 0;
}


int calloc_optprop (optprop_struct *out, int nlambda, int nlev, int nphamat) 
{
  int iv=0, lc=0, ip=0;

  /* BCA: nlev should be replaced by nlyr! Makes more sense */

  out->nlam = nlambda;
  out->nlev = nlev;
  out->nphamat = nphamat;

  out->dtau   = calloc (nlambda, sizeof (float *));
  out->ssa    = calloc (nlambda, sizeof (float *));
  out->f      = calloc (nlambda, sizeof (float *));

  out->dscale = calloc (nlambda, sizeof (float *));

  /* parameters for double HG */
  out->g1   = calloc (nlambda, sizeof (float *));
  out->g2   = calloc (nlambda, sizeof (float *));
  out->ff   = calloc (nlambda, sizeof (float *));
  
  for (iv=0; iv<nlambda; iv++) {
    out->dtau  [iv] = calloc (nlev, sizeof(float));
    out->ssa   [iv] = calloc (nlev, sizeof(float));
    out->f     [iv] = calloc (nlev, sizeof(float));
    out->dscale[iv] = calloc (nlev, sizeof(float));
    out->g1    [iv] = calloc (nlev, sizeof(float));
    out->g2    [iv] = calloc (nlev, sizeof(float));
    out->ff    [iv] = calloc (nlev, sizeof(float));

    for (lc=0; lc<nlev; lc++)
      out->ff[iv][lc] = 1.0;
  }  
  
  out->nmom   = calloc (nlambda, sizeof (int *));
  out->moment = calloc (nlambda, sizeof (float ***)); 
  
  for (iv=0; iv<nlambda; iv++) {
    out->nmom[iv]   = calloc (nlev-1, sizeof (int));
    out->moment[iv] = calloc (nlev-1, sizeof (float **));
    for (lc=0; lc<nlev-1; lc++) {
      out->nmom[iv][lc] = 0;
      out->moment[iv][lc] = calloc (nphamat, sizeof (float *));
    }
  }
  
  out->ntheta = calloc (nlambda, sizeof (int **));
  out->theta  = calloc (nlambda, sizeof (float ***)); 
  out->mu     = calloc (nlambda, sizeof (double ***)); 
  out->phase  = calloc (nlambda, sizeof (float ***)); 

  for (iv=0; iv<nlambda; iv++) {
    out->ntheta[iv]   = calloc (nlev-1, sizeof (int *));
    out->theta [iv]   = calloc (nlev-1, sizeof (float **));
    out->mu    [iv]   = calloc (nlev-1, sizeof (double **));
    out->phase [iv]   = calloc (nlev-1, sizeof (float **));
    for (lc=0; lc<nlev-1; lc++) {
      out->ntheta[iv][lc] = calloc (nphamat, sizeof (int));
      out->theta [iv][lc] = calloc (nphamat, sizeof (float *));
      out->mu    [iv][lc] = calloc (nphamat, sizeof (double *));
      out->phase [iv][lc] = calloc (nphamat, sizeof (float *));
      for (ip=0; ip<nphamat; ip++)
	out->ntheta[iv][lc][ip] = 0;
    }
  }

  out->alloc = 1;

  return 0;
}


int free_optprop (optprop_struct *optprop)
{
  int iv=0, lc=0, ip=0;
  
  for (iv=0; iv<optprop->nlam; iv++) {
    free (optprop->g1[iv]);
    free (optprop->g2[iv]);
    free (optprop->ff[iv]);
    free (optprop->f[iv]);
    free (optprop->dscale[iv]); /*TZ ds*/
    free (optprop->ssa[iv]);
    free (optprop->dtau[iv]);
    free (optprop->nmom[iv]);

    for (lc=0; lc<optprop->nlev-1; lc++) {
      for (ip=0; ip<optprop->nphamat; ip++) {
	if (optprop->moment[iv][lc]!=NULL)
	  if (optprop->moment[iv][lc][ip]!=NULL)
	    free(optprop->moment[iv][lc][ip]);
	if (optprop->theta[iv][lc]!=NULL)
	  if (optprop->theta[iv][lc][ip]!=NULL)
	    free(optprop->theta[iv][lc][ip]);
	if (optprop->mu[iv][lc]!=NULL)
	  if (optprop->mu[iv][lc][ip]!=NULL)
	    free(optprop->mu[iv][lc][ip]);
	if (optprop->phase[iv][lc]!=NULL)
	  if (optprop->phase[iv][lc][ip]!=NULL)
	    free(optprop->phase[iv][lc][ip]);
      }
      if (optprop->moment[iv][lc]!=NULL)
	free (optprop->moment[iv][lc]);
      if (optprop->ntheta[iv][lc]!=NULL)
	free (optprop->ntheta[iv][lc]);
      if (optprop->theta[iv][lc]!=NULL)
	free (optprop->theta[iv][lc]);
      if (optprop->mu[iv][lc]!=NULL)
	free (optprop->mu[iv][lc]);
      if (optprop->phase[iv][lc]!=NULL)
	free (optprop->phase[iv][lc]);
    }

    free (optprop->moment[iv]);
    free (optprop->ntheta[iv]);
    free (optprop->theta[iv]);
    free (optprop->mu[iv]);
    free (optprop->phase[iv]);
  }
  
  free (optprop->g1); 
  free (optprop->g2);
  free (optprop->ff);
  free (optprop->f);
  free (optprop->dscale);
  free (optprop->ssa);
  free (optprop->dtau);
  free (optprop->nmom);
  free (optprop->moment);
  free (optprop->ntheta);
  free (optprop->theta);
  free (optprop->mu);
  free (optprop->phase);

  return 0;
}


int calloc_ssprop_struct (ssprop_struct *out, int nlambda) 
{
  out->nlambda = nlambda;

  out->extinc = calloc (nlambda, sizeof(double *));
  out->albedo = calloc (nlambda, sizeof(double *));
  out->f   = calloc (nlambda, sizeof(double *));
  
  out->r0     = calloc (nlambda, sizeof(double));
  out->dr     = calloc (nlambda, sizeof(double));

  out->nreff  = calloc (nlambda, sizeof(size_t));
  out->nphamat = 1;

  out->reff   = calloc (nlambda, sizeof(double *)); 

  out->nleg   = calloc (nlambda, sizeof(int *));
  out->legen  = calloc (nlambda, sizeof(float ***));

  out->ntheta = calloc (nlambda, sizeof(int **));
  out->theta  = calloc (nlambda, sizeof(float ***));
  out->mu     = calloc (nlambda, sizeof(double ***));
  out->phase  = calloc (nlambda, sizeof(float ***));

  out->alloc          = 1;  /* basic memory (extinc[iv], ...) allocated */
  out->alloc_moments  = 0;  /* phase function moments not yet allocated */
  out->alloc_explicit = 0;  /* phase function values not yet allocated  */

  return 0;
}

int free_ssprop_struct (ssprop_struct *out) 
{
  int iv=0, ir=0, ip=0;
  /* ??? should we also free the fields if out->alloc_xxx is set ??? */

  if (out->alloc) {
    for (iv=0; iv<out->nlambda; iv++) {
      free (out->extinc[iv]);
      free (out->albedo[iv]);
      free (out->f[iv]);
      free (out->reff[iv]);

      if (out->alloc_moments) {
	for (ir=0; ir<out->nreff[iv]; ir++) {
	  if (out->legen[iv][ir]!=NULL) {
	    for (ip=0; ip<out->nphamat; ip++) {
	      if (out->legen[iv][ir][ip]!=NULL)
		free (out->legen[iv][ir][ip]);
	    }
	    free (out->legen[iv][ir]);
	  }
	}

	free (out->nleg[iv]);
	free (out->legen[iv]);
      }
      if (out->alloc_explicit) {
	if (out->ntheta[iv]!=NULL) {
	  for (ir=0; ir<out->nreff[iv]; ir++) {
	    if (out->ntheta[iv][ir]!=NULL) {
	      for (ip=0; ip<out->nphamat; ip++) {
		free (out->theta[iv][ir][ip]);
		free (out->mu[iv][ir][ip]);
		free (out->phase[iv][ir][ip]);
	      }
	      free (out->ntheta[iv][ir]);
	      free (out->theta[iv][ir]);
	      free (out->mu[iv][ir]);
	      free (out->phase[iv][ir]);
	    }
	  }

	  free (out->ntheta[iv]);
	  free (out->theta[iv]);
	  free (out->mu[iv]);
	  free (out->phase[iv]);
	}
      }
    }

    free (out->extinc);
    free (out->albedo);
    free (out->f);
    free (out->reff);

    free (out->r0);    
    free (out->dr);    

    free (out->nreff); 

    free (out->nleg);  
    free (out->legen); 

    free (out->ntheta);
    free (out->theta);
    free (out->mu);
    free (out->phase);
  }

  out->alloc = 0; 
  out->alloc_moments=0;
  out->alloc_explicit=0;

  return 0;
}



/******************************************/
/* allocate memory for a caoth_out_struct */
/******************************************/

int calloc_caoth_out ( caoth_out_struct *out,
		       char             *name,
		       char             *fullname,
		       int               nlambda,
		       int               nlev,
		       int               nphamat,
		       int               calloc_ssprop )
{
  int status=0;

  out->name = (char *) calloc (strlen(name)+1, sizeof (char));
  strcpy (out->name, name);
  out->fullname = (char *) calloc (strlen(fullname)+1, sizeof (char));
  strcpy (out->fullname, fullname);
  out->zd     = calloc (nlev, sizeof(float));
  out->nlev   = nlev;
  out->newsiz = 1;

  status = calloc_caoth_microphys_struct (&(out->microphys), out->nlev);
  if (status!=0) {
    fprintf (stderr, "Error %d allocating memory for cloud microphysics\n", status);
    return status;
  }

  status = calloc_optprop (&(out->optprop), nlambda, out->nlev, nphamat);
  if (status!=0) {
    fprintf (stderr, "Error %d allocating memory for cloud optical properties\n", status);
    return status;
  }

  if (calloc_ssprop) {
    status = calloc_ssprop_struct (&(out->ssprop), nlambda);
    if (status!=0) {
      fprintf (stderr, "Error %d allocating memory for cloud single scattering properties\n",
	       status);
      return status;
    }
  }

  return 0;
}


/**************************************/
/* free memory for a caoth_out_struct */
/**************************************/

int free_caoth_out (caoth_out_struct *out, int free_ssprop)
{

  free (out->zd);

  free_caoth_microphys_struct (&(out->microphys));
  
  //aky free_optprop (&(out->optprop));

  if (free_ssprop)
    free_ssprop_struct (&(out->ssprop));

  return 0;
}


int cp_optprop (optprop_struct *target,
		optprop_struct *source,
		int             nlev,
		int             alloc_moment,
		int             quiet)
{
  int iv=0, lev=0, k=0, ip=0;
  
  if (target->alloc && source->alloc) {
    
    /* check if number of wavelength bands equal */
    if (target->nlam != source->nlam) {
      fprintf (stderr, "Error: cp_optprop(), nlam differs between source and target\n");
      fprintf (stderr, "       source: %d,  target: %d\n", 
	       source->nlam, target->nlam);
      return -1;
    }

    for (lev=0; lev<nlev-1; lev++) {
      for (iv=0; iv<source->nlam; iv++) {
	target->dtau  [iv][lev] = source->dtau  [iv][lev];
	target->ssa   [iv][lev] = source->ssa   [iv][lev];
	target->f     [iv][lev] = source->f     [iv][lev];
        target->dscale[iv][lev] = source->dscale[iv][lev];

	target->g1    [iv][lev] = source->g1    [iv][lev];
	target->g2    [iv][lev] = source->g2    [iv][lev];
	target->ff    [iv][lev] = source->ff    [iv][lev];
	
	if (source->nmom[iv][lev] > 0) {
	  
	  /* allocate memory (a) if the user wants it; or    */
	  /* (b) if the number of moments in the target is 0 */
	  if (alloc_moment || (target->nmom[iv][lev]==0 && source->nmom[iv][lev]>0)) {
	    target->nmom[iv][lev] = source->nmom[iv][lev];
            /* Dimension of phase matrix elements */
            target->moment[iv][lev] = calloc (source->nphamat, sizeof(float *));
            for (ip=0; ip<source->nphamat; ip++)
              target->moment[iv][lev][ip] = calloc (source->nmom[iv][lev], sizeof(float));
          }

	  /* check if source and target compatible */
	  if (target->nmom[iv][lev] != source->nmom[iv][lev]) {
	    fprintf (stderr, "Error, number of moments in source (%d) and target (%d) different\n",
		     source->nmom[iv][lev], target->nmom[iv][lev]);
	    return -1;
	  }

          for (ip=0; ip<source->nphamat; ip++)
            for (k=0; k<source->nmom[iv][lev]; k++)
              target->moment[iv][lev][ip][k] = source->moment[iv][lev][ip][k];
        }

	/* XXX might be necessary! */
	/*
	   if (alloc_moment || (target->nphamat=0 && source->nphamat>0)) {
	   //(ntheta[iv][lev]==NULL && source->ntheta[iv][lev]!=NULL)) {
	   target->ntheta[iv][lev] = calloc (source->nphamat, sizeof(int));
	   target->theta [iv][lev] = calloc (source->nphamat, sizeof(float *));
	   target->mu    [iv][lev] = calloc (source->nphamat, sizeof(double *));
	   target->phase [iv][lev] = calloc (source->nphamat, sizeof(float *));
	   }
	*/

        if(source->ntheta[iv][lev]!=NULL){
          for (ip=0; ip<source->nphamat; ip++) {
            if (source->ntheta[iv][lev][ip] > 0) {
              /* allocate memory (a) if the user wants it; or    */
              /* (b) if the number of phases in the target is 0 */
              if (alloc_moment || (target->ntheta[iv][lev][ip]==0 && source->ntheta[iv][lev][ip]>0)) {
                /* Dimension of phase matrix elements */
                target->ntheta[iv][lev][ip] = source->ntheta[iv][lev][ip];
                target->theta [iv][lev][ip] = calloc (source->ntheta[iv][lev][ip], sizeof(float));
                target->mu    [iv][lev][ip] = calloc (source->ntheta[iv][lev][ip], sizeof(double));
                target->phase [iv][lev][ip] = calloc (source->ntheta[iv][lev][ip], sizeof(float));
                
              }
              
              /* check if source and target compatible */
              if (target->ntheta[iv][lev][ip] != source->ntheta[iv][lev][ip]) {
                fprintf (stderr, "Error, number of moments in source (%d) and target (%d) different\n",
                         source->ntheta[iv][lev][ip], target->ntheta[iv][lev][ip]);
                return -1;
              }
              
              for (k=0;k<source->ntheta[iv][lev][ip];k++) {
                target->theta [iv][lev][ip][k] = source->theta [iv][lev][ip][k];
                target->mu    [iv][lev][ip][k] = source->mu    [iv][lev][ip][k];
                target->phase [iv][lev][ip][k] = source->phase [iv][lev][ip][k];
              }
            }
	  }
	}
      } /* iv loop */
    } /* lev loop */
    target->nphamat=source->nphamat; 
  }
  else {
    /*
    if (!quiet) {

      if (!source->alloc) {
	fprintf (stderr, " *** Warning, not copying optical properties in cp_optprop() because\n");
	fprintf (stderr, " *** memory for source not allocated\n");
      }

      if (!target->alloc) {
	fprintf (stderr, " *** Warning, not copying optical properties in cp_optprop() because\n");
	fprintf (stderr, " *** memory for target not allocated\n");
      }
    }
    */
  }

  return 0;
}
    

/*****************************************************/
/* copy contents of a caoth_out_struct to another;   */
/* NO memory allocation is done here!!!              */
/*****************************************************/

int cp_caoth_out ( caoth_out_struct *target,
		   caoth_out_struct  source, 
		   int               copy_ssprop,
		   int               alloc_moment,
		   int               quiet )
{
  int iv=0, ir=0, itheta=0, k=0, ip=0, lev=0, status=0;

  if (target->nlev < source.nlev) {
    fprintf (stderr, "Error: nlev differs between source and target ");
    fprintf (stderr, "(line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    fprintf (stderr, "       source: %d,  target: %d\n", source.nlev, target->nlev);
    return -1;
  }

  target->newsiz = source.newsiz;

  /* copy microphysical properties */
  for (lev=0; lev<source.nlev; lev++) {
    target->zd[lev] = source.zd[lev];
    if (target->microphys.alloc) {
      target->microphys.lwc  [lev] = source.microphys.lwc  [lev];
      target->microphys.effr [lev] = source.microphys.effr [lev];

      target->microphys.lwc_layer [lev] = source.microphys.lwc_layer  [lev];
      target->microphys.effr_layer[lev] = source.microphys.effr_layer [lev];
    }
  }

  /* copy single scattering properties */
  if (copy_ssprop)  {

    if (!source.ssprop.alloc) {
      fprintf (stderr, "Error, no memory for single scattering properties allocated\n");
      fprintf (stderr, "       (line %d, function %s in %s), source structure!\n",
	       __LINE__, __func__, __FILE__);
      return -1;
    }

    if (!target->ssprop.alloc) {
      fprintf (stderr, "Error, no memory for single scattering properties allocated\n");
      fprintf (stderr, "       (line %d, function %s in %s), target structure!\n",
	       __LINE__, __func__, __FILE__);
      return -1;
    }

    target->ssprop.type = source.ssprop.type;

    if (source.ssprop.alloc_moments || source.ssprop.alloc_explicit) {
      target->ssprop.nphamat = source.ssprop.nphamat;
      
      for (iv=0; iv<source.optprop.nlam; iv++) {
	target->ssprop.nreff[iv] = source.ssprop.nreff[iv];
	target->ssprop.r0[iv]    = source.ssprop.r0[iv];
	target->ssprop.dr[iv]    = source.ssprop.dr[iv];

	/* allocate memory for moments if required */
	if (source.ssprop.alloc_moments && !target->ssprop.alloc_moments) {
	  target->ssprop.extinc[iv] = calloc (target->ssprop.nreff[iv], sizeof(double));
	  target->ssprop.albedo[iv] = calloc (target->ssprop.nreff[iv], sizeof(double));
	  target->ssprop.f[iv] = calloc (target->ssprop.nreff[iv], sizeof(double));
          
	  target->ssprop.nleg[iv]   = calloc (target->ssprop.nreff[iv], sizeof(int));
	  target->ssprop.legen[iv]  = calloc (target->ssprop.nreff[iv], sizeof(float *));

          for (ir=0; ir<target->ssprop.nreff[iv]; ir++) {
            target->ssprop.legen[iv][ir] = calloc (target->ssprop.nphamat, sizeof(float *));

            for (ip=0; ip<target->ssprop.nphamat; ip++)
              target->ssprop.legen[iv][ir][ip] = calloc (source.ssprop.nleg[iv][ir],
							 sizeof(float));
          }
        }

	/* allocate memory for phase function if required */

	/* the following "if" statements are a quick and dirty fix. A clean fix would */
	/* be to only copy the ssprop for certain wavelengths, i.e. where the interpolation */
	/* in lambda was actually done initially. To this end this fct needs nlambda_lower */
	/* and nlambda_upper in the output structure. iv should only run between these ... */
	if(source.ssprop.theta [iv]!=NULL) {
	  if(source.ssprop.theta [iv][0]!=NULL) {
	    if (source.ssprop.alloc_explicit && !target->ssprop.alloc_explicit) {
	      target->ssprop.extinc[iv] = calloc (target->ssprop.nreff[iv], sizeof(double));
	      target->ssprop.albedo[iv] = calloc (target->ssprop.nreff[iv], sizeof(double));
	      target->ssprop.f[iv] = calloc (target->ssprop.nreff[iv], sizeof(double));
          
	      target->ssprop.ntheta[iv]  = calloc (target->ssprop.nreff[iv], sizeof(int *));
	      target->ssprop.theta [iv]  = calloc (target->ssprop.nreff[iv], sizeof(float **));
	      target->ssprop.mu    [iv]  = calloc (target->ssprop.nreff[iv], sizeof(double **));
	      target->ssprop.phase [iv]  = calloc (target->ssprop.nreff[iv], sizeof(float **));

	      for (ir=0; ir<target->ssprop.nreff[iv]; ir++) {
		target->ssprop.ntheta[iv][ir] = calloc (target->ssprop.nphamat, sizeof(int));
		target->ssprop.theta [iv][ir] = calloc (target->ssprop.nphamat, sizeof(float *));
		target->ssprop.mu    [iv][ir] = calloc (target->ssprop.nphamat, sizeof(double*));
		target->ssprop.phase [iv][ir] = calloc (target->ssprop.nphamat, sizeof(float *));
	      }

	      for (ir=0; ir<target->ssprop.nreff[iv]; ir++)
		for (ip=0; ip<target->ssprop.nphamat; ip++) {
		  target->ssprop.theta [iv][ir][ip] = calloc (source.ssprop.ntheta[iv][ir][ip],
							      sizeof (float));
		  target->ssprop.mu    [iv][ir][ip] = calloc (source.ssprop.ntheta[iv][ir][ip],
							      sizeof (double));
		  target->ssprop.phase [iv][ir][ip] = calloc (source.ssprop.ntheta[iv][ir][ip],
							      sizeof(float));
		}
	    }
	  }
	}
      }

      target->ssprop.alloc_moments  = source.ssprop.alloc_moments; 
      target->ssprop.alloc_explicit = source.ssprop.alloc_explicit;

      /* finally copy source to target */
      for (iv=0; iv<source.optprop.nlam; iv++) {
	if (target->ssprop.alloc_moments) {
	  for (ir=0; ir<target->ssprop.nreff[iv]; ir++) {
	    target->ssprop.extinc [iv][ir] = source.ssprop.extinc [iv][ir];
	    target->ssprop.albedo [iv][ir] = source.ssprop.albedo [iv][ir];
            target->ssprop.f      [iv][ir] = source.ssprop.f      [iv][ir];
	    target->ssprop.nleg   [iv][ir] = source.ssprop.nleg   [iv][ir];
	    
	    for (k=0; k<target->ssprop.nleg[iv][ir]; k++)
	      for (ip=0; ip<target->ssprop.nphamat; ip++)
		target->ssprop.legen [iv][ir][ip][k] = source.ssprop.legen [iv][ir][ip][k];
	  }
	}


	/* the following "if" statements are a quick and dirty fix. A clean fix would */
	/* be to only copy the ssprop for certain wavelengths, i.e. where the interpolation */
	/* in lambda was actually done initially. To this end this fct needs nlambda_lower */
	/* and nlambda_upper in the output structure. iv should only run between these ... */
	if(source.ssprop.theta [iv]!=NULL) {
	  if(source.ssprop.theta [iv][0]!=NULL) {
	    if (target->ssprop.alloc_explicit) {
	      for (ir=0; ir<target->ssprop.nreff[iv]; ir++) {
		target->ssprop.extinc [iv][ir] = source.ssprop.extinc [iv][ir];
		target->ssprop.albedo [iv][ir] = source.ssprop.albedo [iv][ir];
		target->ssprop.f      [iv][ir] = source.ssprop.f      [iv][ir];

		for (ip=0; ip<target->ssprop.nphamat; ip++) {
		  target->ssprop.ntheta [iv][ir][ip] = source.ssprop.ntheta [iv][ir][ip];
	    
		  for (itheta=0; itheta<target->ssprop.ntheta[iv][ir][ip]; itheta++) {
		    target->ssprop.mu    [iv][ir][ip][itheta]
		      = source.ssprop.mu    [iv][ir][ip][itheta];
		    target->ssprop.theta [iv][ir][ip][itheta]
		      = source.ssprop.theta [iv][ir][ip][itheta];
		    target->ssprop.phase [iv][ir][ip][itheta]
		      = source.ssprop.phase [iv][ir][ip][itheta];
		  } /* end for itheta */
		} /* end for ip */
	      } /* end for ir */
	    } /* end if explicit */
	  }
	}
      } /* end for iv */
    } /* end if moments || explicit */
  } /* end if copy ssprop */
  

  /* copy optical properties (in this file cloud.c) */
  status = cp_optprop (&(target->optprop), &(source.optprop), source.nlev, alloc_moment, quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d copying cloud optical properties (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  
  return 0;
}


/*****************************************************/
/* copy contents of a cloud_out_struct to another;   */
/* NO memory allocation is done here!!!              */
/*****************************************************/

int cp_caoth3d_out ( caoth_out_struct   *target,
		     caoth3d_out_struct  source,
		     int                 copy_ssprop,
		     int                 alloc_moment,
		     int                 quiet,
		     int                 iipa,
		     int                 jipa )
{
  int lev=0;

  /* ulrike: replaced source.nlev by source.nlyr  */
  if (target->nlev < source.nlyr) {
    fprintf (stderr, "Error: nlev differs between source and target (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    fprintf (stderr, "       source: %d,  target: %d\n", source.nlyr, target->nlev);
    return -1;
  }
  
  /*   target->newsiz = source.newsiz; // ulrike: auskommentiert*/

  /* copy microphysical properties */
  for (lev=0; lev<source.nlyr; lev++) {
    target->zd[lev] = source.zd[source.nlyr-lev]; /* ulrike */
    
    if (target->microphys.alloc) {
      if (source.threed[source.nlyr-lev-1]) {

	/*ulrike - used for hu?! evtl -2 wenn segmentation fault*/
	target->microphys.lwc  [lev+1] = source.lwc  [source.nlyr-lev-1][iipa][jipa];
	target->microphys.effr [lev+1] = source.reff [source.nlyr-lev-1][iipa][jipa];

	/* used for mie */
	target->microphys.lwc_layer [lev] = source.lwc  [source.nlyr-lev-1][iipa][jipa];
	target->microphys.effr_layer[lev] = source.reff [source.nlyr-lev-1][iipa][jipa];
      }
    }
  }

  /* copy single scattering properties  ulrike: not necessary (?) */

  return 0;
}


/***********************************************************************************/
/* Function: interpolate_ssprop_in_reff                                   @62_30i@ */
/* Description:                                                                    */
/* Interpolate the Legendre coefficients which are given on a grid of              */
/* effective radii r0, r0 + dr, ..., r0 + (nreff-1) * dr to a given                */
/* radius reff; memory for the array of Legendre coefficients coeffc[]             */
/* is allocated automatically.                                                     */
/*                                                                                 */
/* CE: 24.06.2008: Modified interpolation method. Interpolation no longer linear   */
/* for all optical properties. Single scattering albedo is weighted with extinction*/
/* and Legendre coefficients are weighted with extinction * albedo.                */
/* Also included interpolation of delta scaling factor                             */
/*                                                                                 */
/* RPB: 07.05.2009: Completely rewritten for speed up and clearness                */
/*                  added interpolation of phase function                          */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: last mayor change by Robert Buras                                       */
/*                                                                        @i62_30@ */
/***********************************************************************************/

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
				int **ntheta_new )
{
  int k=0, ip=0, status=0;
  int lower=0, upper=0;
  int min=0, max=0;
  
  double reffl=0.0;
  double interp_optprop_weight_0=0.0;
  double interp_optprop_weight[2]={0.,0.};
  double norm=0.0;

  int **tmp_ntheta=NULL;
  float ***tmp_theta=NULL, ***tmp_phase=NULL;
  double ***tmp_mu=NULL;

  if (reff<r0 || reff> r0 + (double) (nreff-1) * dr) {
    fprintf (stderr, "Error, reff = %f out of range (%f - %f)\n",
	     reff, r0, r0 + (double) (nreff-1) * dr);
    return -1;
  }

  /* Indices of effective radius */
  lower = (int) ((reff - r0) / dr);
  upper = lower+1;
  
  reffl = r0 + (double) lower * dr;
  
  /* upper boundary value, trivial case */
  if (lower==nreff-1) {
    /* dummy interpolation */
    upper = lower;
    interp_optprop_weight_0 = 0.0;
   
    /*
     *ext = extinc[lower];
     *alb = albedo[lower];

     if (interpolate_moments) {
     *ncoeffc = nleg[nreff-1];
     *coeffc = calloc (nphamat, sizeof(float *));
     for (ip=0; ip<nphamat; ip++){
     (*coeffc)[ip]  = calloc (nleg[nreff-1], sizeof(float*));
     for (k=0; k<*ncoeffc; k++)
     (*coeffc)[ip][k] = legen[nreff-1][ip][k];
     }
     }    

     if (interpolate_explicit) {
     *ntheta_new = calloc (nphamat, sizeof(int));
     *theta_new  = calloc (nphamat, sizeof(float *));
     *mu_new     = calloc (nphamat, sizeof(double *));
     *phase_new  = calloc (nphamat, sizeof(float *));


     }

     *f   = fdelta[lower];
    
     return 0;
    */
  }
  else {
    /* else not upper boundary value, interpolation required */
  
    /* prepare interpolation factor interp_optprop_weight */
    interp_optprop_weight_0 = ( reff - reffl ) / dr;
  }


  /*Interpolate extinction linearly*/
  interp_optprop_weight[0] = 1. - interp_optprop_weight_0;
  interp_optprop_weight[1] =      interp_optprop_weight_0;
  *ext = interp_optprop_weight[0] * extinc[lower]
       + interp_optprop_weight[1] * extinc[upper];
  
  /* vanishing extinction makes interpolation of other variables obsolete */
  if ( *ext == 0. )
    return 0;

  /* Interpolate single scattering albedo */
  /*Single scattering albedo should be weighted with extinction, this way the 
    interpolation means that the optical properties are averaged */
  interp_optprop_weight[0] = ( 1. - interp_optprop_weight_0 ) * extinc[lower] / *ext;
  interp_optprop_weight[1] = (      interp_optprop_weight_0 ) * extinc[upper] / *ext;
  *alb = interp_optprop_weight[0] * albedo[lower]
       + interp_optprop_weight[1] * albedo[upper];
  
  /* vanishing albedo makes interpolation of other variables obsolete */
  if ( *alb == 0. )
    return 0;

  /* The moments and phases need to be weighted with alb*ext */ 

  if (interpolate_moments || interpolate_explicit) {
    interp_optprop_weight[0] = ( 1. - interp_optprop_weight_0 ) *
      albedo[lower] * extinc[lower] / ( *alb * *ext );
    interp_optprop_weight[1] = (      interp_optprop_weight_0 ) *
      albedo[upper] * extinc[upper] / ( *alb * *ext );
  }

  /* Interpolate moments */
  if (interpolate_moments) {
    min = (nleg[lower] < nleg[upper] ? nleg[lower] : nleg[upper]);
    max = (nleg[lower] > nleg[upper] ? nleg[lower] : nleg[upper]);

    *ncoeffc = max;

    *coeffc = calloc (nphamat, sizeof(float *));

    for (ip=0; ip<nphamat; ip++){
      (*coeffc)[ip]  = calloc (*ncoeffc, sizeof(float));

      for (k=0; k<min; k++)      
	(*coeffc)[ip][k] = interp_optprop_weight[0] * legen[lower][ip][k]
	                 + interp_optprop_weight[1] * legen[upper][ip][k];

      if (nleg[lower] < nleg[upper])
	for (k=min; k<*ncoeffc; k++)
	  (*coeffc)[ip][k] = interp_optprop_weight[1] * legen[upper][ip][k];
      else
	for (k=min; k<*ncoeffc; k++)
	  (*coeffc)[ip][k] = interp_optprop_weight[0] * legen[lower][ip][k];

    }
  }
  else
    *ncoeffc = 0;

  /* interpolate phase function */
  if (interpolate_explicit) {
    /* allocate new theta */
    *ntheta_new = calloc (nphamat, sizeof(int));
    *theta_new  = calloc (nphamat, sizeof(float *));
    *mu_new     = calloc (nphamat, sizeof(double *));
    *phase_new  = calloc (nphamat, sizeof(float *));

    /* allocate temporary theta */
    tmp_ntheta = calloc(nphamat, sizeof(int *));
    tmp_theta  = calloc(nphamat, sizeof(float **));
    tmp_mu     = calloc(nphamat, sizeof(double **));
    tmp_phase  = calloc(nphamat, sizeof(float **));

    for (ip=0; ip<nphamat; ip++) {
      tmp_ntheta[ip] = calloc(2, sizeof(int));
      tmp_theta [ip] = calloc(2, sizeof(float *));
      tmp_mu    [ip] = calloc(2, sizeof(double *));
      tmp_phase [ip] = calloc(2, sizeof(float *));

      tmp_ntheta[ip][0] = ntheta[lower][ip];
      tmp_ntheta[ip][1] = ntheta[upper][ip];

      tmp_theta [ip][0] = theta [lower][ip];
      tmp_theta [ip][1] = theta [upper][ip];

      tmp_mu    [ip][0] = mu    [lower][ip];
      tmp_mu    [ip][1] = mu    [upper][ip];

      tmp_phase [ip][0] = phase [lower][ip];
      tmp_phase [ip][1] = phase [upper][ip];
    }

    /* perform interpolation */
    status = sort_and_add_weighted_phase (2, interp_optprop_weight,
					  tmp_ntheta, tmp_theta, tmp_mu, tmp_phase,
					  &(*ntheta_new),
					  &(*theta_new ),
					  &(*mu_new    ),
					  &(*phase_new ),
					  nphamat, 0, 1 );

    if (status) {
      fprintf(stderr,"Error! something (%d) went wrong in sort_and_add_weighted_phase \n",status);
      return status;
    }

    for (ip=0; ip<nphamat; ip++) {
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
  else {
    *ntheta_new = calloc (nphamat, sizeof(int));
    for (ip=0; ip<nphamat; ip++)
      (*ntheta_new)[ip]=0;
  }

  /* Interpolate delta scaling factor */
  /* The delta scaling factor needs to be weighted with the UNSCALED alb*ext */
  interp_optprop_weight[0] = ( 1. - interp_optprop_weight_0 ) * extinc[lower] * albedo[lower] / ( 1. - fdelta[lower]);
  interp_optprop_weight[1] =        interp_optprop_weight_0   * extinc[upper] * albedo[upper] / ( 1. - fdelta[upper]);
  norm = interp_optprop_weight[1] + interp_optprop_weight[0];

  /* prevent division by zero */
  if ( norm == 0. )
    return 0;

  interp_optprop_weight[0] /= norm;
  interp_optprop_weight[1] /= norm;

  *f = interp_optprop_weight[0] * fdelta[lower]
     + interp_optprop_weight[1] * fdelta[upper];
  /* To avoid rounding errors for f=0.0 */
  if (*f<1e-8)
    *f=0.0;
  
  /*  debugging output 
      fprintf(stderr, "reff %g reffl %g dr %g \n", reff, reffl, dr);
      fprintf(stderr, "f_low %g, f_up %g \n", fdelta[lower], fdelta[upper]);
      fprintf(stderr, "ext_low %g, ext_up %g, alb_low %g, alb_up %g \n",
      extinc[lower], extinc[upper], albedo[lower], albedo[upper]);
      fprintf(stderr, "ext %g, alb %g, fdel %g \n", *ext, *alb, *f);
      fprintf(stderr, "unscaled ext %g, alb %g \n", ext_us, alb_us);
  */

  return 0;
}


/***********************************************************************************/
/* Function: read_caoth_file                                              @62_30i@ */
/* Description:                                                                    */
/* Read caoth data from file; microphysical properties are expected                */
/* in three columns:                                                               */
/*   column CLD_ZD:    z [km]                                                      */
/*   column CLD_LWC:   lwc [g / m3]                                                */
/*   column CLD_EFFR:  effective radius [micron]                                   */
/* Memory for microphysical and optical properties is allocated                    */
/* automatically. This function can be used for caoth                              */
/* because all output is written to caoth_out_struct *caoth                        */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

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
		      cf_out_struct    *cf )
{
  int lev=0, rows=0, min_columns=0, max_columns=0, status=0, lc=0;
  float **caoth_data=NULL;
    
  switch (source) {
  /* read caoth from one file */
  case CAOTH_FROM_1D:
  case CAOTH_FROM_IPA_FILES:

    status = ASCII_file2float ( filename,
				&rows, &max_columns, &min_columns, 
				&caoth_data );
    if (status!=0) {
      fprintf (stderr, "Error %d opening file %s\n", status, filename);
      return status;
    }
  
    if (max_columns!=min_columns) {
      fprintf (stderr, " !! ATTENTION !! Inconsistent number of columns\n");
      fprintf (stderr, "     min = %d, max =%d\n", min_columns, max_columns);
    }

    if (min_columns<3) {
      fprintf (stderr, "Error, found less than 3 columns in %s\n",
	       filename);
      return -1;
    }

    break;

  case CAOTH_FROM_ECMWF:
    /* read cloud from ECMWF file */
    status = read_ECMWF_clouds ( filename,
				 filename_ic_reff,
				 latitude,
				 longitude, 
				 UTC,
				 time_interpolate, 
				 caoth_name, 
				 reff_prop,
				 properties,
				 reff_fixed,
				 altitude,
				 press_atm,
				 cloud_overlap, 
				 verbose,
				 quiet,
				 &rows,
				 &max_columns,
				 &min_columns, 
				 &caoth_data,
				 cf );
    if (status)
      return fct_err_out ( status, "read_ECMWF_clouds", ERROR_POSITION );
    break;

  case CAOTH_FROM_ECHAM:
    status = read_ECHAM_clouds ( filename,
				 latitude,
				 longitude,
				 UTC,
				 time_interpolate, 
				 caoth_name,
				 z_model_layer, 
				 verbose,
				 quiet, 
				 &rows,
				 &max_columns,
				 &min_columns,
				 &caoth_data,
				 cf );
    if (status)
      return fct_err_out ( status, "read_ECHAM_clouds", ERROR_POSITION );
    break;

  default:  
    fprintf (stderr, "Error, unknown '%s source = %d' (line %d, function %s in %s)\n",
	     caoth_fullname, source, __LINE__, __func__, __FILE__);
    fprintf (stderr, "  This is a program bug, please report it to us. \n");
    return -1;
  }

  /* allocate memory for 1D caoth strucure */

  /* CE: set nphamat here to 6 to avoid valgrind errors with IPA files */
  /* nphamat from optprop files is not known at this place */
  if (optical) {
    status = calloc_caoth_out ( caoth,
				caoth_name,
				caoth_fullname,
				nlambda,
				rows,
				6,
				1 );
    if (status)
      return mem_err_out ("output_caoth struct ", ERROR_POSITION );
  }
  else {
    caoth->nlev   = rows;
    caoth->zd     = calloc (caoth->nlev, sizeof(float));
    caoth->newsiz = 1;

    status = calloc_caoth_microphys_struct ( &(caoth->microphys), caoth->nlev );
    if (status)
      return mem_err_out ( "caoth_microphys struct", ERROR_POSITION );
  }

  if (verbose) {
    fprintf (stderr, "  lc     zd      cwc          eff      \n");
    fprintf (stderr, "        [km]    [g/m3]    [micro meter]\n");
  }
  for (lev=0; lev<caoth->nlev; lev++) {
    caoth->zd            [lev] = caoth_data[lev][CLD_ZD];
    caoth->microphys.lwc [lev] = caoth_data[lev][CLD_LWC];
    caoth->microphys.effr[lev] = caoth_data[lev][CLD_EFFR];
    if (verbose)
      fprintf (stderr, " %3d %7.3f %e %e\n", lev, caoth->zd[lev],
	       caoth->microphys.lwc[lev], caoth->microphys.effr[lev]);
  }

  /* free temporary memory */
  ASCII_free_float(caoth_data, rows);
  
  /* newsiz is a flag, telling the Fortran function wcloud() if a */
  /* new liquid water content and/or effective radius are given   */
  /* (newsiz=1) or if the values haven't changed since the last   */
  /* call (newsiz=0).                                             */ 
  caoth->newsiz = 1;

  /* test if levels are sorted in descending order */
  for (lev=0; lev<caoth->nlev-1; lev++)
    if (caoth->zd[lev] <= caoth->zd[lev+1]) {
      fprintf (stderr, "Error, input grid in %s not sorted in descending order!\n", filename);
      fprintf (stderr, "z[%d] = %f, z[%d] = %f\n",
	       lev, caoth->zd[lev], lev+1, caoth->zd[lev+1]);
      return -1;
    }
  
  /* calculate layer microphysical properties. */
  /* needs to be done in exactly the same      */
  /* way as in wcloud.f                        */
  if (caoth->microphys.alloc) {
    if (layer) {
      for (lc=0; lc<caoth->nlev-1; lc++) {
	caoth->microphys.lwc_layer [lc] = caoth->microphys.lwc [lc+1];
	caoth->microphys.effr_layer[lc] = caoth->microphys.effr[lc+1];
      }
    }
    else {
      for (lc=0; lc<caoth->nlev-1; lc++) {
	if (caoth->microphys.lwc[lc+1]>0 && caoth->microphys.lwc[lc]>0) {
	  caoth->microphys.lwc_layer [lc] = 
	    0.5 * (caoth->microphys.lwc [lc+1] + caoth->microphys.lwc [lc]);
	  caoth->microphys.effr_layer[lc] = 
	    0.5 * (caoth->microphys.effr[lc+1] + caoth->microphys.effr[lc]);
	}
	else {
	  caoth->microphys.lwc_layer [lc] = 0.0;
	  caoth->microphys.effr_layer[lc] = 0.0;
	}
      }
    }
  }

  /* determine minimum and maximum effective radii */
  /* of cloudy layers; we assume that the layer    */
  /* properties are within this range, as they are */
  /* linearely interpolated from the level prop's  */

  caoth->microphys.effrmin=+FLT_MAX;
  caoth->microphys.effrmax=-FLT_MAX;
  for (lev=0; lev<caoth->nlev; lev++) {
    if (caoth->microphys.lwc [lev] > 0) {
      if (caoth->microphys.effr[lev] > caoth->microphys.effrmax)
	caoth->microphys.effrmax = caoth->microphys.effr[lev]; 
      if (caoth->microphys.effr[lev] < caoth->microphys.effrmin)
	caoth->microphys.effrmin = caoth->microphys.effr[lev]; 
    }
  }

  return 0;
}





/* Read caoth data from file using read_caoth_file() and convert */
/* microphysical properties to optical properties                */

int read_and_convert_caoth_file ( input_struct      input,
				  char             *filename,
				  wl_out_struct     wl_out,
				  float             altitude,
				  atm_out_struct    atm_out,
				  caoth_inp_struct  input_caoth,
				  /* output */
				  caoth_out_struct *output_caoth,
				  cf_out_struct    *cf_out,
				  int              *nmom )
{
  int iv=0, lev=0, nlyr=0, newkey=0, status=0;
  int nstring = (int) strlen(input.filename[FN_PATH]);

  void  F77_FUNC (wcloud, WCLOUD) (float *lambda_r, int *newsiz, int *nlyr,
				   char *filepath, int *nstring, float *lwc, float *wceffr, 
				   float *tmp_wc_dtau, float *tmp_wc_gg, float *tmp_wc_ssa,
				   float *zd, int *wclyr);

  /* read caoth file */
  status = read_caoth_file ( input_caoth.source,
			     input_caoth.name,
			     input_caoth.fullname,
			     filename,
			     input.filename[FN_IC_REFF],
			     input.latitude,
			     input.longitude, 
			     input.UTC,
			     input.atm.time_interpolate, 
			     input_caoth.reff_prop,
			     input_caoth.properties,
			     input_caoth.reff,
			     altitude,
			     atm_out.microphys.press[0][0], 
			     input.cloud_overlap, 
			     atm_out.z_model_layer,
			     1,
			     wl_out.nlambda_r,
			     input_caoth.layer,
			     input.verbose,
			     input.quiet,
			     output_caoth,
			     cf_out );
  if (status)
    return fct_err_out ( status, "read_caoth_file", ERROR_POSITION );

  /* default: assume HG caoth - flag needed by MYSTIC */
  output_caoth->cldproperties=CLD_OPTPROP;
  output_caoth->nonHG=0;

  switch (input_caoth.properties) {
  case PROP_HU:

    switch (input.ck_scheme) {
    case CK_FU:

      /* call Fu and Liou code by Fred Rose */
      
#if HAVE_FULIOU
      status = ckdfucld ( output_caoth->zd, 
			  output_caoth->microphys.effr_layer, 
			  output_caoth->microphys.lwc_layer, 
			  output_caoth->nlev,
			  output_caoth->optprop.dtau, 
			  output_caoth->optprop.g1, 
			  output_caoth->optprop.ssa);
      if (status)
	return fct_err_out ( status, "ckdfucld", ERROR_POSITION );
#else
      fprintf (stderr, "Error, Fu and Liou not supported!\n");
      return -1;
#endif
      
      break;
      
    default:   /* old uvspec wcloud.f for Hu and Stamnes */
      
      nlyr = output_caoth->nlev-1;

      /* wcloud calculates lwc_layer and reff_layer internally */
      for (iv=0; iv<wl_out.nlambda_r; iv++)
	F77_FUNC (wcloud, WCLOUD) ( &wl_out.lambda_r[iv],
				    &(output_caoth->newsiz),
				    &nlyr,
				    input.filename[FN_PATH],
				    &nstring,
				    output_caoth->microphys.lwc,
				    output_caoth->microphys.effr,
				    output_caoth->optprop.dtau[iv],
				    output_caoth->optprop.g1[iv],
				    output_caoth->optprop.ssa[iv],
				    output_caoth->zd,
				    &input_caoth.layer );

      break;
    }

    break;

  case PROP_FU:

    if (!input.quiet)
      fprintf (stderr, " ... using Fu et al. [1996/98] ice cloud properties\n");

    if (input_caoth.fu2yang) {

      /* convert from Fu (1996/98) to Key et al. (2002) effective radius */

      for (lev=0; lev<output_caoth->nlev; lev++) {
	output_caoth->microphys.effr[lev] /= (3.0*sqrt(3.0)/4.0);

	/* now used */
        if ( output_caoth->microphys.lwc_layer[lev]>0 && !input.quiet )
	  fprintf (stderr, "     yang2fu: scaling effective radius from %f to",
		   output_caoth->microphys.effr_layer[lev]);

	output_caoth->microphys.effr_layer[lev] /= (3.0*sqrt(3.0)/4.0);

        if ( output_caoth->microphys.lwc_layer[lev]>0 && !input.quiet )
	  fprintf (stderr, " %f\n", output_caoth->microphys.effr_layer[lev]);
      }

    }

    /* calculate caoth optical properties */
    switch (input.ck_scheme) {
    case CK_FU:
      
#if HAVE_FULIOU
      status = ckdfuice ( output_caoth->zd, 
			  output_caoth->microphys.effr_layer, 
			  output_caoth->microphys.lwc_layer, 
			  output_caoth->nlev,
			  output_caoth->optprop.dtau, 
			  output_caoth->optprop.g1, 
			  output_caoth->optprop.ssa, 
			  input_caoth.unscaled );
      if (status)
	return fct_err_out ( status, "ckdfuice", ERROR_POSITION );
#else
      fprintf (stderr, "Error, Fu and Liou not supported!\n");
      return -1;
#endif
      break;
      
    default:
      
      for (iv=wl_out.nlambda_rte_lower; iv<=wl_out.nlambda_rte_upper; iv++) {
	
	if (wl_out.lambda_r[iv]<4000.0) {
	  status = ic_fu96 ( wl_out.lambda_r[iv],
			     output_caoth->nlev-1,
			     input.filename[FN_PATH],
			     output_caoth->microphys.lwc_layer,
			     output_caoth->microphys.effr_layer, 
			     output_caoth->optprop.dtau[iv],
			     output_caoth->optprop.g1[iv],
			     output_caoth->optprop.ssa[iv],
			     output_caoth->optprop.f[iv],
			     output_caoth->zd,
			     input_caoth.layer,
			     input_caoth.unscaled);
	  if (status)
	    return fct_err_out ( status, "ic_fu96", ERROR_POSITION );
	}
	else {
	  status = ic_fu98 ( wl_out.lambda_r[iv],
			     output_caoth->nlev-1,
			     input.filename[FN_PATH],
			     output_caoth->microphys.lwc_layer,
			     output_caoth->microphys.effr_layer, 
			     output_caoth->optprop.dtau[iv],
			     output_caoth->optprop.g1[iv],
			     output_caoth->optprop.ssa[iv], 
			     output_caoth->zd,
			     input_caoth.layer);
	  if (status)
	    return fct_err_out ( status, "ic_fu98", ERROR_POSITION );
	}
      }
    }

    break;

  case PROP_KEY:
  case PROP_YANG:

    if (input_caoth.properties==PROP_KEY) {
      if (!input.quiet)
	fprintf (stderr, " ... using Key et al. [2002] ice cloud properties\n");

      newkey=0;
    }

    if (input_caoth.properties==PROP_YANG) {
      if (!input.quiet) {
	fprintf (stderr, " ... using Key et al. [2002] ice cloud properties below 3.4micron;\n");
	fprintf (stderr, " ... new data by Yang/Mayer above 3.4 micron.\n");
      }

      newkey=1;
    }

    for (iv=wl_out.nlambda_rte_lower; iv<=wl_out.nlambda_rte_upper; iv++) {

      status = ic_yang ( wl_out.lambda_r[iv],
			 output_caoth->nlev-1,
			 input_caoth.habit,
			 input.filename[FN_PATH],
			 output_caoth->microphys.lwc_layer,
			 output_caoth->microphys.effr_layer, 
			 output_caoth->optprop.dtau[iv], 
			 output_caoth->optprop.ff[iv],
			 output_caoth->optprop.g1[iv],
			 output_caoth->optprop.g2[iv], 
			 output_caoth->optprop.ssa[iv], 
			 output_caoth->zd,
			 input_caoth.layer,
			 newkey );
      if (status)
	return fct_err_out ( status, "ic_yang", ERROR_POSITION );
    }
    break;

  case PROP_ECHAM4:
    if (strcasecmp(output_caoth->name,"wc")==0) {
      for (iv=wl_out.nlambda_rte_lower; iv<=wl_out.nlambda_rte_upper; iv++) {
        status = wc_echam4 ( wl_out.lambda_r[iv],
			     output_caoth->nlev-1,
			     output_caoth->microphys.lwc,
			     output_caoth->microphys.effr,
			     output_caoth->optprop.dtau[iv],
			     output_caoth->optprop.g1[iv],
			     output_caoth->optprop.ssa[iv],
			     output_caoth->zd,
			     input_caoth.layer );
	if (status)
	  return fct_err_out ( status, "wc_echam4", ERROR_POSITION );
      }
    }
    else if (strcasecmp(output_caoth->name,"ic")==0) {
      for (iv=wl_out.nlambda_rte_lower; iv<=wl_out.nlambda_rte_upper; iv++) {
	status = ic_echam4 ( wl_out.lambda_r[iv],
			     output_caoth->nlev-1,
			     output_caoth->microphys.lwc,
			     output_caoth->microphys.effr,
			     output_caoth->optprop.dtau[iv],
			     output_caoth->optprop.g1[iv],
			     output_caoth->optprop.ssa[iv],
			     output_caoth->zd,
			     input_caoth.layer );
	if (status)
	  return fct_err_out ( status, "ic_echam4", ERROR_POSITION );
      }
    }
    else {
      fprintf(stderr,"Error! You can not use echam4 for profile %s\n",output_caoth->name);
      return -1;
    }
    break;  

  case PROP_MIE:
  case PROP_FILE:
  case PROP_IC_MIE:
  case PROP_BAUM:
  case PROP_BAUM_V36:
  case PROP_HEY: 
  case PROP_YANG2013:

    /* Read optical properties file */
    if (!(output_caoth->ssprop.alloc && output_caoth->ssprop.alloc_moments)) {

      /* we don't interpolate level properties to layer properties here */
      /* (could be done but requires some work)                         */
      if (!input_caoth.layer) {
	fprintf (stderr,
		 "Error, %s_properties%s requires optical properties defined per layer.\n",
		 input_caoth.name1, input_caoth.name2);
	fprintf (stderr, "Please specify %s_layer%s and retry if that is what you mean!\n",
		 input_caoth.name1, input_caoth.name2);
	return -1;
      }

      status = read_caoth_prop ( input_caoth.properties_filename, 
				 wl_out.lambda_r,
				 wl_out.nlambda_r, 
				 wl_out.nlambda_rte_lower,
				 wl_out.nlambda_rte_upper, 
				 input_caoth.interpolate,
				 0,
				 output_caoth,
				 input_caoth.properties,
				 NULL,
				 NULL,
				 0,
				 nmom,
				 input.rte.polradtran[POLRADTRAN_NSTOKES],
				 input.rte.nstr,
				 input.rte.solver,
				 input.rte.disort_icm,
				 input.verbose,
				 input.quiet );
      if (status!=0) {
	fprintf (stderr, "Error %d reading property file %s for %s\n", 
		 status, input_caoth.properties_filename, input_caoth.fullname);
	return status;
      }

    }

    /* non Henyey-Greenstein - flag needed by MYSTIC */
    output_caoth->nonHG=1;
    output_caoth->cldproperties=CLD_LWCREFF;

    break;

  case PROP_RAYTRACING:
    if (strncasecmp(output_caoth->name,"ic", 2)==0) {
      status = read_raytracing_file (input.filename[FN_RAYTRACING], &(output_caoth->raytracing_prop), &(output_caoth->n_raytracing_prop), input.quiet);

      if (status!=0) {
	fprintf (stderr, "Error %d reading %s\n", status, input.filename[FN_RAYTRACING]);
	return status;
      }
    }

    for (iv=wl_out.nlambda_rte_lower; iv<=wl_out.nlambda_rte_upper; iv++) {

      /* simple geometrical optics calculation of the extinction cross section */
      /* tau = 3 IWC / (2 rho_ice reff); geometrical optics cross section      */
      /* sigma_geo = 2 pi r**2 is twice the geometrical cross section          */

      status = ic_geom ( wl_out.lambda_r[iv],
			 output_caoth->nlev-1,
			 output_caoth->microphys.lwc_layer,
			 output_caoth->microphys.effr_layer, 
			 output_caoth->zd,
			 output_caoth->optprop.dtau[iv], 
			 output_caoth->optprop.ssa[iv]);
      
      if (status)
	return fct_err_out ( status, "ic_geom", ERROR_POSITION );
    }
    
    break;
    
    
  default:
    fprintf (stderr, "Error, unknown property %d for %s. This is not\n",
	     input_caoth.properties, input_caoth.fullname);
    fprintf (stderr, "supposed to happen and indicates a coding error.\n");

    return -1;
  }

  return 0;
}


/* convert microphysical and single scattering properties */
/* to optical properties at model layers                  */

/* ????? the following function used lwc[lc+1] and reff[lc+1] ????? */
/* ????? which was obviously wrong; replaced by lwc_layer[lc] ????? */
/* ????? and effr_layer[lc]                                   ????? */

int ssprop2optprop (caoth_out_struct *cld, int iv, int quiet)
{

  int lc=0, status=0;

  for (lc=0; lc<cld->nlev-1; lc++) {  /* count layers, 0 ... nlev-2 */
    
    /* first free memory, if already allocated */
    if (cld->optprop.moment[iv][lc] != NULL) {
      free (cld->optprop.moment[iv][lc]);
      cld->optprop.moment[iv][lc] = NULL;
    }

    if (cld->microphys.lwc_layer[lc] > 0) {
      /*CE: Interpolation must only be performed if nreff in optical properties
        file is at least 2 !!!*/
      if(cld->ssprop.nreff[iv] > 1){
        /* interpolate Legendre coefficients to specified effective radius; */
        /* memory for moments is allocated here in this function.           */
	status = interpolate_ssprop_in_reff (cld->microphys.effr_layer[lc],
					     cld->ssprop.nreff[iv],
					     cld->ssprop.r0[iv], 
					     cld->ssprop.dr[iv],
					     cld->ssprop.extinc[iv], 
					     cld->ssprop.albedo[iv],
					     cld->ssprop.f[iv],
					     cld->ssprop.alloc_moments,
					     cld->ssprop.legen[iv],  
					     cld->ssprop.nleg[iv],
					     cld->ssprop.alloc_explicit,
					     cld->ssprop.ntheta[iv],
					     cld->ssprop.theta[iv],
					     cld->ssprop.mu[iv],
					     cld->ssprop.phase[iv],
					     cld->ssprop.nphamat,
					     &(cld->optprop.dtau[iv][lc]), 
					     &(cld->optprop.ssa[iv][lc]),
					     &(cld->optprop.f[iv][lc]),
					     &(cld->optprop.moment[iv][lc]), 
					     &(cld->optprop.nmom[iv][lc]),
					     &(cld->optprop.phase[iv][lc]),
					     &(cld->optprop.theta[iv][lc]),
					     &(cld->optprop.mu[iv][lc]),
					     &(cld->optprop.ntheta[iv][lc]));
          
	if (status!=0) {
	  fprintf (stderr, "Error %d interpolating Legendre coefficients to reff = %f in layer %d.\n", 
		   status, cld->microphys.effr_layer[lc], lc);
	  return status;
	}
      }
      else{ /*no interpolation required if Nreff =1*/
        /* CE: This should probably be an error ???? */
	if(cld->microphys.effr_layer[lc] != cld->ssprop.r0[iv] && quiet!=1 ){
          fprintf(stderr, " *** Warning: The effective radius in the optical \n");
          fprintf(stderr, " ***          properties file %f does not correspond to the \n",
		  cld->ssprop.r0[iv]);
          fprintf(stderr, " ***          effective radius of the cloud %f in layer %d. \n",
		  cld->microphys.effr_layer[lc], lc );
        }
        /* ???CE Allocation not necessary here, all quantities are already allocated */
        /*
          cld->optprop.moment[iv][lc]  = calloc(cld->ssprop.nphamat, sizeof(float *)); 
          for (ip=0; ip<cld->ssprop.nphamat; ip++)
          cld->optprop.moment[iv][lc][ip]  = calloc(cld->ssprop.nleg[iv][0], sizeof(float));
          cld->optprop.ntheta[iv][lc] = calloc (cld->ssprop.nphamat, sizeof(int));
          cld->optprop.theta[iv][lc]  = calloc (cld->ssprop.nphamat, sizeof(float *));
          cld->optprop.mu[iv][lc]     = calloc (cld->ssprop.nphamat, sizeof(double *));
          cld->optprop.phase[iv][lc]  = calloc (cld->ssprop.nphamat, sizeof(float *));
        */
        
        cld->optprop.moment[iv][lc]=cld->ssprop.legen[iv][0];
        cld->optprop.nmom[iv][lc]=cld->ssprop.nleg[iv][0];
        cld->optprop.ssa[iv][lc]=cld->ssprop.albedo[iv][0];
        cld->optprop.f[iv][lc]=cld->ssprop.f[iv][0];
        cld->optprop.dtau[iv][lc]=cld->ssprop.extinc[iv][0];
	if(cld->ssprop.alloc_explicit) { /* bug fix by RPB on 14.10.2010: when using simple pmom netcdf file, this is not adequate */
	  cld->optprop.phase[iv][lc]=cld->ssprop.phase[iv][0];
	  cld->optprop.theta[iv][lc]=cld->ssprop.theta[iv][0];
	  cld->optprop.mu[iv][lc]=cld->ssprop.mu[iv][0];
	  cld->optprop.ntheta[iv][lc]=cld->ssprop.ntheta[iv][0];
	}
      }

      /* Number of phase matrix elements constant */
      cld->optprop.nphamat = cld->ssprop.nphamat;

      /* convert from extinction / lwc to optical thickness dtau */
      cld->optprop.dtau[iv][lc] *= cld->microphys.lwc_layer[lc] * (cld->zd[lc] - cld->zd[lc+1]);
      
      if (cld->optprop.nmom[iv][lc]>1)
        cld->optprop.g1[iv][lc] = cld->optprop.moment[iv][lc][0][1];
    }
    else {
      cld->optprop.dtau[iv][lc] = 0;
      cld->optprop.ssa [iv][lc] = 0;
      cld->optprop.f[iv][lc] = 0; 
      cld->optprop.g1  [iv][lc] = 0;
      
      /* reset Legendre coefficients */
      cld->optprop.nmom[iv][lc] = 0;

      /* reset phase function ??? */
      /* this caused problems when the array was not allocated! */
      /* if it should be necessary, find out in advance whether ntheta is fully allocated! */
      /*
	for (ip=0; ip < cld->ssprop.nphamat; ip++)
      	cld->optprop.ntheta[iv][lc][ip] = 0;
      */
    }
  }
  

  return 0;
}


/***********************************************************************************/
/* Function: aer_ssprop2optprop                                                    */
/* Description:                                                                    */
/* Convert microphysical and single scattering properties of aerosol particles     */
/* to optical properties at model layers                                           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: aerosol.c                                                                */
/* Known bugs:                                                                     */
/* Author: Claudia Emde , 22.02.2008                                               */    
/*                                                                                 */
/***********************************************************************************/

int aer_ssprop2optprop (aer_out_struct *aer, int iv, int i_aer, float *rh,
			       int nlambda, int nlambda_lower, int nphamat, int verbose)
{
  
  int lc=0, status=0;
  int i_hum=0; 
  float diff1=0.0, diff2=0.0;
  int i_hum_nearest=0;
  
  if (iv==nlambda_lower) {
    status = calloc_optprop(&aer->optprop_opac[i_aer], nlambda, aer->nlev, nphamat);
    if (status!=0) {
      fprintf (stderr, "Error %d allocating memory for aer struct\n", status);
      return status;
    }
  }

  if(verbose) {
    if (aer->ssprop.nreff[iv] > 1) {
      fprintf(stderr, " ... search for nearest relative humidity point in aerosol data, iv = %d \n", iv);
      fprintf(stderr, "     level    rh(atmos)   rh(aerofile) \n");
      fprintf(stderr, "     ----------------------------------\n");
    }
    else {
      fprintf(stderr, " ... aerosol does not swell, no relative humidity adjustment needed\n");
    }
  }

  for (lc=0; lc<aer->nlev-1; lc++) {  /* count layers, 0 ... nlev-2 */

    aer->optprop_opac[i_aer].moment[iv][lc] = NULL;

    aer->optprop_opac[i_aer].ntheta[iv][lc] = NULL;
    aer->optprop_opac[i_aer].theta [iv][lc] = NULL;
    aer->optprop_opac[i_aer].mu    [iv][lc] = NULL;
    aer->optprop_opac[i_aer].phase [iv][lc] = NULL;
    
    if (aer->massdens[i_aer][lc] > 0) {
      /* Search fo closest relative humidity point in database. */
      /* ???? It still needs to be checked, whether interpolation is the better option. ????*/
      /* *interpolate_ssprop_in_reff* only works for equidistant reff/hum grid, so it needs to */
      /* be modified if we want to do interpolation. */
      if(aer->ssprop.nreff[iv] > 1){ /* else hum 0 is the only point. In this case aerosols do not swell. */
        diff1 = fabs( 0.5*(rh[lc]+rh[lc+1]) - aer->ssprop.reff[iv][0] );
        for (i_hum=1; i_hum<aer->ssprop.nreff[iv]; i_hum++){
          diff2=fabs( 0.5*(rh[lc]+rh[lc+1]) - aer->ssprop.reff[iv][i_hum]);
          if( diff2 < diff1 ) {
            i_hum_nearest = i_hum;
            diff1=diff2;
          }
        }
      }
      
      if(verbose && aer->ssprop.nreff[iv] > 1 )
        fprintf(stderr, "  %8d     %9.5f    %9.5f \n", 
               lc, 0.5*(rh[lc]+rh[lc+1]), aer->ssprop.reff[iv][i_hum_nearest ] );
      
      if (aer->ssprop.alloc_moments) {
	aer->optprop_opac[i_aer].nmom[iv][lc]=aer->ssprop.nleg[iv][i_hum_nearest];
	aer->optprop_opac[i_aer].moment[iv][lc]=aer->ssprop.legen[iv][i_hum_nearest];
      }

      if (aer->ssprop.alloc_explicit) {
	aer->optprop_opac[i_aer].ntheta[iv][lc]=aer->ssprop.ntheta[iv][i_hum_nearest];
	aer->optprop_opac[i_aer].theta [iv][lc]=aer->ssprop.theta [iv][i_hum_nearest];
	aer->optprop_opac[i_aer].mu    [iv][lc]=aer->ssprop.mu    [iv][i_hum_nearest];
	aer->optprop_opac[i_aer].phase [iv][lc]=aer->ssprop.phase [iv][i_hum_nearest];
      }
      
      aer->optprop_opac[i_aer].ssa[iv][lc]=aer->ssprop.albedo[iv][i_hum_nearest];
      aer->optprop_opac[i_aer].dtau[iv][lc]=aer->ssprop.extinc[iv][i_hum_nearest];
      
      /* Number of phase matrix elements constant */
      aer->optprop_opac[i_aer].nphamat = aer->ssprop.nphamat;
      
      /* convert from extinction / mass density to optical thickness dtau */
      /* The profile is given in g/m^3, we have exactly the same as for clouds, the optical properties */
      /* are normalized to 1g/m^3 hence there is no conversion factior. */
      aer->optprop_opac[i_aer].dtau[iv][lc] *= aer->massdens[i_aer][lc] *
        (aer->zd[lc] - aer->zd[lc+1]);
      
      if (aer->optprop_opac[i_aer].nmom[iv][lc]>1)
        aer->optprop_opac[i_aer].g1[iv][lc] = aer->optprop_opac[i_aer].moment[iv][lc][0][1];
      
    }
    
    else {
      aer->optprop_opac[i_aer].dtau[iv][lc] = 0;
      aer->optprop_opac[i_aer].ssa [iv][lc] = 0;
      aer->optprop_opac[i_aer].g1  [iv][lc] = 0;
      
    /* reset Legendre coefficients */
      aer->optprop_opac[i_aer].nmom[iv][lc] = 0;
    }
  }

  return 0;
}


/*****************************************************/
/* read caoth properties from optical property files */
/*****************************************************/

int read_caoth_prop ( char             *filename,
		      float            *lambda_r,
		      int               nlambda_r, 
		      int               nlambda_lower,
		      int               nlambda_upper,
		      int               interpolate,
		      int               aerosol,
		      caoth_out_struct *caoth,
		      int               caoth_properties,
		      aer_out_struct   *aer,
		      float            *rh,
		      int               i_aer,
		      int              *nmom,
		      int               nstokes,
		      int               nstrmax,
		      int               solver,
		      int               disort_icm,
		      int               verbose,
		      int               quiet )
{
  int status=0, lc=0, iw=0, ifi=0, i=0, iv=0, iws=0, niw=0;
  int lambda_optprop_lower=0, lambda_optprop_upper=0;
  int type=0;

  ssprop_struct tmp_ssprop, *ssprop;
  float *tmp_lambda=NULL;
  int tmp_nlambda=0, nfiles=0;

  char momfilename[FILENAME_MAX]="";
  char type_string[8]="";
  char tempstr[20]="";

  glob_t momfiles;

  size_t *nlam=NULL;
  double **wavelen=NULL, *wave_arr=NULL, nre=0, nim=0;
  double wvln=0.0;
  int read_scattering_phase_function=0, *iwstart=NULL, *ifile=NULL, nlambda_read=0;

  /* define whether we rather want phases or moments or both BCA replace numbers by defines */
  switch(solver) {
  case SOLVER_FDISORT2:
  case SOLVER_DISORT:
    switch(disort_icm) {
    case DISORT_ICM_OFF:
      read_scattering_phase_function=READ_SPF_MOMENTS;
      break;
    case DISORT_ICM_PHASE:
      read_scattering_phase_function=READ_SPF_BOTH;
      break;
    case DISORT_ICM_MOMENTS:
      read_scattering_phase_function=READ_SPF_ALL_MOMENTS;
      break;
    default:
      fprintf (stderr, "Error, unknown disort_intcor %d in read_caoth_prop\n",
	       disort_icm);
      return -1;
    }
    break;
  case SOLVER_SSLIDAR:
    read_scattering_phase_function=READ_SPF_PHASES;
    break;
  case SOLVER_MONTECARLO:
    read_scattering_phase_function=READ_SPF_PHASES;
    break;
  case SOLVER_FDISORT1:
  case SOLVER_SDISORT:
  case SOLVER_FTWOSTR:
  case SOLVER_SOS:
  case SOLVER_POLRADTRAN:
  case SOLVER_SPSDISORT:
  case SOLVER_TZS:
  case SOLVER_SSS:
  case SOLVER_SSSI:
  case SOLVER_NULL:
  case SOLVER_RODENTS:
  case SOLVER_TWOSTREBE:
  case SOLVER_TWOMAXRND:
  case SOLVER_TWOSTR:
    read_scattering_phase_function=READ_SPF_MOMENTS;
    break;
  default:
     fprintf (stderr, "Error, unknown rte_solver %d in read_caoth_prop\n", solver);
     return -1;
  }

  if (aerosol) {
    sprintf (type_string,"aerosol");
    ssprop = &(aer->ssprop);
  }
  else {
    sprintf (type_string,"cloud");
    ssprop = &(caoth->ssprop);
  }

  if (verbose)
    fprintf (stderr, " ... reading %s optical properties from %s\n", type_string, filename);
  
  if (interpolate) {
    /* read optical properties from file */
    /* try netcdf first                  */

    if (aerosol)
      /* either *.mie.cdf (for spherical aerosols), or *.tmatrix.cdf for aspherical */
      sprintf (momfilename, "%s.*.cdf", filename);
    else {
      if (caoth_properties==PROP_IC_MIE || caoth_properties==PROP_MIE)
	sprintf (momfilename, "%s*.mie.cdf", filename);
      else if (caoth_properties==PROP_BAUM || caoth_properties==PROP_BAUM_V36)
	sprintf (momfilename, "%s*.baum.cdf", filename);
      else if (caoth_properties==PROP_HEY)
        sprintf (momfilename, "%s.sol.cdf", filename);
      else if (caoth_properties==PROP_YANG2013)
	sprintf (momfilename, "%s.*.cdf", filename);
      /* User defined */
      else{
	if (!quiet)
	  fprintf(stderr, " ... user defined optical properties file . \n");
	sprintf (momfilename, "%s*", filename);
      }
    }
    
    /* get list of filenames that match the pattern momfilename???.???.cdf */
    status = glob (momfilename, 0, NULL, &momfiles);
    nfiles = momfiles.gl_pathc; /* number of files found */
   
    if(nfiles==0) { /* if no file was found try it without searching for .*.cdf */
      sprintf (momfilename, "%s", filename);
      status = glob (momfilename, 0, NULL, &momfiles);
      nfiles = momfiles.gl_pathc; /* number of files found */
    }

    if (!quiet)
      fprintf (stderr, " ... reading %s optical properties from %s\n", type_string, momfilename);
    
    if (!aerosol) {
      /* if no netcdf files found, try ASCII */
      if (nfiles==0) {
	if (caoth_properties==PROP_IC_MIE)
	  sprintf (momfilename, "%s*.mie", filename);
	else if (caoth_properties==PROP_BAUM)
	  sprintf (momfilename, "%s*.baum", filename);
      
	if (!quiet)
	  fprintf (stderr, " ... no netcdf cloud optical property files found, trying ASCII instead!\n");
      
	status = glob (momfilename, 0, NULL, &momfiles);
	nfiles = momfiles.gl_pathc; /* number of files found */
      }
    }

    if (nfiles==0) {
      fprintf (stderr, "Error, found neither netcdf nor ASCII optical property files.\n");
      fprintf (stderr, "Check your %s properties and retry!\n", type_string);
      return -1;
    }

    
    /* remove extension .cdf if necessary   */
    /* because read_mie_table automatically */
    /* appends .cdf to the filename         */
    for (iw=0; iw<nfiles; iw++)
      if (strlen(momfiles.gl_pathv[iw])>=4)
	if (!strcmp(momfiles.gl_pathv[iw]+strlen(momfiles.gl_pathv[iw])-4, ".cdf"))
	  momfiles.gl_pathv[iw][strlen(momfiles.gl_pathv[iw])-4]=0;

    wavelen = calloc(nfiles, sizeof(double *));
    nlam    = calloc(nfiles, sizeof(size_t));
    iwstart = calloc(nfiles+1, sizeof(int));

    /* first determine required wavelength range */
    tmp_nlambda=0;
    for (ifi=0; ifi<nfiles; ifi++) {
      status = read_mie_table_lambda (momfiles.gl_pathv[ifi], &(nlam[ifi]),
				      &(wavelen[ifi]), quiet);
      
      if (status!=0) {
	fprintf (stderr, "Error %d reading wavelength from %s\n", 
		 status, momfiles.gl_pathv[ifi]);
	return status;
      }

      iwstart[ifi]=tmp_nlambda;
      tmp_nlambda+=nlam[ifi];
    }
    iwstart[ifi]=tmp_nlambda;

    /* copy wavelengths to one-dimensional array */
    tmp_lambda = calloc (tmp_nlambda, sizeof (float));
    ifile = calloc(tmp_nlambda, sizeof (int));

    i=0;
    for (ifi=0; ifi<nfiles; ifi++)
      for (iw=0; iw<nlam[ifi]; iw++) {
	tmp_lambda[i] = (float) wavelen[ifi][iw];
	ifile[i] = ifi;
	i++;
      }

    for (ifi=0; ifi<nfiles; ifi++)
      free(wavelen[ifi]);
    free(wavelen);

    /* check if wavelengths are sorted in ascending order;  */
    /* we could in principle sort them, but if they are not */
    /* sorted we assume that the user did something he did  */
    /* not really want.                                     */
    
    for (iw=0; iw<tmp_nlambda-1; iw++)
      if (tmp_lambda[iw+1] <= tmp_lambda[iw]) {
	fprintf (stderr, "Error, %s optical properties files not sorted by ascending wavelength\n", type_string);
	fprintf (stderr, "       %f nm >= %f nm\n", tmp_lambda[iw], tmp_lambda[iw+1]);
	return -1;
      }
    
    if (!quiet)
      fprintf (stderr, " ... found %d %s optical properties files covering the wavelength range %f - %f nm with a total of %d discrete wavelengths\n",
	       nfiles, type_string, tmp_lambda[0], tmp_lambda[tmp_nlambda-1],tmp_nlambda);
    
    /* check if wavelength range is sufficient */
    if (tmp_lambda[0]>lambda_r[nlambda_lower] || 
	tmp_lambda[tmp_nlambda-1]<lambda_r[nlambda_upper]) {
      fprintf (stderr, "Error, wavelength range covered by %s optical properties files (%f - %f nm)\n", 
	       type_string, tmp_lambda[0], tmp_lambda[tmp_nlambda-1]);
      fprintf (stderr, "does not cover the required range (%f - %f nm)\n",
	       lambda_r[nlambda_lower], lambda_r[nlambda_upper]);
      return -1;
    }

    /* determine required wavelength range */
    for (iw=0; iw<tmp_nlambda; iw++)
      if (tmp_lambda[iw] > lambda_r[nlambda_lower])
	break;
    
    lambda_optprop_lower=iw-1;
    
    for (iw=0; iw<tmp_nlambda; iw++)
      if (tmp_lambda[iw] >= lambda_r[nlambda_upper])
	break;
    
    lambda_optprop_upper=iw;

    /* CE: */
    /* There are numerical problems in interpolation of single scattering properties if*/
    /* the wavelength to be computed is exactly on the boarder of the wavelength range */
    /* of the optical properties files. Therefore wavelength interval is extended. */
    if (lambda_optprop_lower > 0)
      lambda_optprop_lower-=1;
    if (lambda_optprop_upper < tmp_nlambda-1)
      lambda_optprop_upper+=1;
    
    if (!quiet)
      fprintf (stderr, " ... selected %s files %d - %d, wavelengths %f - %f nm\n",
	       type_string, ifile[lambda_optprop_lower], ifile[lambda_optprop_upper],
	       tmp_lambda[lambda_optprop_lower], tmp_lambda[lambda_optprop_upper]);

    /* define number of wavelengths to be read */
    nlambda_read = lambda_optprop_upper - lambda_optprop_lower + 1;

    /* allocate memory; need only single scattering properties, */
    /* hence nlev is set to 0                                   */
    status = calloc_ssprop_struct (&tmp_ssprop, nlambda_read);
    if (status!=0) {
      fprintf (stderr, "Error %d allocating memory for %s ssprop struct\n",
	       status, type_string);
      return status;
    }

    wave_arr = calloc(nlambda_read, sizeof(double));

    /* read mie tables */
    /* from each file, "niw" wavelengths are read starting with wavelength "iws" */
    /* the array index of the first wavelength to be read is "iw"                */ 
    iw=0;
    for (ifi=ifile[lambda_optprop_lower]; ifi<=ifile[lambda_optprop_upper]; ifi++) {
      iws = lambda_optprop_lower - iwstart[ifi];
      if (iws < 0) iws=0;
      niw = ( iwstart[ifi+1]-1 < lambda_optprop_upper ? iwstart[ifi+1]-1 : lambda_optprop_upper )
	- iws - iwstart[ifi] + 1;

      status = read_mie_table (momfiles.gl_pathv[ifi],
			       (tmp_ssprop.r0+iw),
			       (tmp_ssprop.dr+iw),
			       (tmp_ssprop.nreff+iw),
                               (tmp_ssprop.reff+iw),
			       wave_arr+iw, &nre, &nim,
			       (tmp_ssprop.extinc+iw),
			       (tmp_ssprop.albedo+iw),
                               (tmp_ssprop.f+iw),
			       (tmp_ssprop.nleg+iw),
			       (tmp_ssprop.legen+iw),
                               &(tmp_ssprop.nphamat),
			       (tmp_ssprop.ntheta+iw),
			       (tmp_ssprop.theta+iw),
			       (tmp_ssprop.mu+iw),
			       (tmp_ssprop.phase+iw),
			       &(tmp_ssprop.alloc_moments),
			       &(tmp_ssprop.alloc_explicit),
			       aerosol,
                               nstokes,
			       nstrmax,
			       read_scattering_phase_function,
			       iws, niw,
			       quiet);
      
      if (status!=0) {
	fprintf (stderr, "Error %d reading %s\n", status, momfilename);
	return status;
      }

      if (verbose)
	fprintf (stderr, " ... read %3zd reff's from file %3d, %s, %11.4f - %11.4f nm\n", 
		 tmp_ssprop.nreff[iw], ifi, momfiles.gl_pathv[ifi], wave_arr[iw], wave_arr[iw+niw-1]);

      iw += niw;

    }
    
    if (!quiet)
      fprintf (stderr, " ... interpolating %s single scattering properties to internal wavelength grid\n", type_string);

    status = interpolate_ssprop_in_lambda (&tmp_ssprop, ssprop,
					   wave_arr, nlambda_read,
					   lambda_r, nlambda_r,
					   nlambda_lower, nlambda_upper);

    if (status!=0) {
      fprintf (stderr, "Error %d returned by interpolate_ssprop_in_lambda()\n", status);
      return status;
    }

    free(tmp_lambda);
    free(wave_arr);
  }
  else {
    
    if(caoth_properties==PROP_BAUM){
      fprintf (stderr, "Error: Optical properties for Baum parameterization must be \n");
      fprintf (stderr, "       in netcdf format (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    if (!quiet) {
      fprintf (stderr, "WARNING: You have not specified to interpolate the optical properties!\n");
      fprintf (stderr, "         Although this might work, you should check thoroughly whether\n");
      fprintf (stderr, "         the optical properties are correct. Use the option 'verbose'.\n");
      fprintf (stderr, "         In particular, you are now using the first wavelength from the\n");
      fprintf (stderr, "         optprop file for the first wavelength you specified in your\n");
      fprintf (stderr, "         input file, and the second for the second etc. Is this what\n");
      fprintf (stderr, "         you want?\n");
    }

    for (iv=0; iv<nlambda_r; iv++) {
      
      /* read optical properties from file */
      strcpy (momfilename, filename);
      
      if (iv+1<100)
	strcat (momfilename, "0");
      
      if (iv+1<10) 
	strcat (momfilename, "0");
      
      sprintf (tempstr, "%d", iv+1);
      
      strcat (momfilename, tempstr);
      strcat (momfilename, ".mie");
      
      status = read_mie_table (momfilename,
			       &(caoth->ssprop.r0[iv]),
			       &(caoth->ssprop.dr[iv]),
			       &(caoth->ssprop.nreff[iv]),
                               &(caoth->ssprop.reff[iv]),
			       &wvln, &nre, &nim,
			       &(caoth->ssprop.extinc[iv]),
			       &(caoth->ssprop.albedo[iv]),
                               &(caoth->ssprop.f[iv]),
			       &(caoth->ssprop.nleg[iv]),
			       &(caoth->ssprop.legen[iv]),
                               &(caoth->ssprop.nphamat),
			       &(caoth->ssprop.ntheta[iv]),
			       &(caoth->ssprop.theta[iv]),
			       &(caoth->ssprop.mu[iv]),
			       &(caoth->ssprop.phase[iv]),
			       &(caoth->ssprop.alloc_moments),
			       &(caoth->ssprop.alloc_explicit),
			       0,
                               nstokes,
			       nstrmax,
			       read_scattering_phase_function,
			       0, 1,
			       quiet);
      
      if (status!=0) {
	fprintf (stderr, "Error %d reading %s\n", status, momfilename);
	return status;
      }
    }
  }
 
  /* convert microphysical and single scattering properties */
  /* to optical properties  */ 
  /* CE: Changed this loop to go only from nlambda_lower to lambda_upper, the wavelengths to be calculated. */
  /* this doesn't break anything */
  /*for (iv=0; iv<nlambda_r; iv++) {*/
  for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {

    if (aerosol)
      status = aer_ssprop2optprop (aer, iv, i_aer, rh, nlambda_r, nlambda_lower,  tmp_ssprop.nphamat, verbose);
    else
      status = ssprop2optprop (caoth, iv, quiet);  /* in cloud.c, this file */
    
    if (status!=0) {
      fprintf (stderr, "Error %d calculating %s optical properties\n", status, type_string);
      return status;
    }
  }

  /* evtl put this into ssprop2optprop ? NO!!!, there only define optprops! */
  if (ssprop->alloc_moments) {
    if (ssprop->alloc_explicit)
      type = PHASE_HYBRID; /* read Legendre moments and phase function */
    else
      type = PHASE_MOMENTS; /* read Legendre moments */
  }
  else {
    if (ssprop->alloc_explicit)
      type = PHASE_EXPLICIT; /* read Legendre moments and phase function */
    else
      type = PHASE_NONE;
  }
  if (aerosol)
    aer->ssprop.type = type;
  else
    caoth->ssprop.type = type;


  /* increase number of moments if more found */
  if (aerosol) {
    for (iv=nlambda_lower; iv<nlambda_upper; iv++)
     
      for (lc=0; lc<aer->nlev-1; lc++){
       	if (aer->optprop_opac[i_aer].nmom[iv][lc]-1 > *nmom)
	  *nmom = aer->optprop_opac[i_aer].nmom[iv][lc]-1;
      }
  }
  else {
    for (iv=0; iv<nlambda_r; iv++)
      for (lc=0; lc<caoth->nlev-1; lc++)
	if (caoth->optprop.nmom[iv][lc]-1 > *nmom)
	  *nmom = caoth->optprop.nmom[iv][lc]-1;
  }
  
  /* free tmp_ssprop (ulrike/rpb 11.05.2010) */

  if (interpolate)
    free_ssprop_struct(&tmp_ssprop);
  
  return 0;
}


/* overwrite water cloud properties with user-defined properties */

/* CE: 24.06.2008: Modified this function so that all user defined 
   properties are applied to un-deltascaled optical properties. 
   This assures that for instances *ic_properties baum* and *ic_properties baum_detailed*
   give the same result when ic_set_tau is used.
   If delta scaling is not applied, the delta scaling factor f equals 0 and 
   all options give the same result as before the modificaltion of this function.
*/

int apply_user_defined_properties_to_caoth ( caoth_inp_struct  cldin, 
					     int               nlambda_r,
					     float            *lambda_r,
					     float             altitude,
					     /* output */
					     caoth_out_struct *cldout )
{
  int iv=0, lc=0, status=0;
  float scale_factor=0.0, sum=0.0;
  double *tautot = NULL, *lambda = NULL;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double tau550=0.0;
  double unscaled_gg=0.0, unscaled_ssa=0.0;


  for (iv=0; iv<nlambda_r; iv++) {
    
    /* Set the cloud asymmetry factor to a user-defined value */
    if (cldin.modify[MODIFY_VAR_GG][MODIFY_TYPE_SET] >= -1.0) {
      for (lc=0; lc<cldout->nlev-1; lc++) {
        
        /* If delta scaling is applied, for example in Fu parameterization, we need the 
           delta scaled asymmetry parameter in optprop_struct*/
	cldout->optprop.g1[iv][lc] = (cldin.modify[MODIFY_VAR_GG][MODIFY_TYPE_SET] - cldout->optprop.f[iv][lc])/
          (1.0 - cldout->optprop.f[iv][lc]) ;

        /* To assure that the asymmetry-parameter has the correct value, g2 of the Double-
           HG-function is set to 0. */
        cldout->optprop.g2[iv][lc] = 0.0;
        cldout->optprop.ff[iv][lc] = 1.0;

	/* We don't want a user-defined asymmetry parameter if an explicit  */
	/* phase function was defined (e.g. ic_properties mie, baum);       */
	/* it is not a good idea to overwrite the first moment because then */
	/* the phase function might look quite weird: The first moment      */
	/* defines the linear term in the Legendre expansion; changing that */
	/* results basically in a change of the slope of the phase function */
	/* which easily results in negative values. An alternative would be */
	/* to replace the explicit phase function by a simple HG, but is    */
	/* that what the user really wants?                                 */

        if (cldout->optprop.nmom[iv][lc] > 1) {
	  fprintf (stderr, "Error, it does not make sense to use an explicit phase function\n");
	  fprintf (stderr, "and then to manually over-write the asymmetry parameter.\n"); 
	  return -1;
	}
      }
    }
    
    /* Scale the cloud asymmetry factor with a user-defined value */
    if (cldin.modify[MODIFY_VAR_GG][MODIFY_TYPE_SCALE] >= 0.0) {
      for (lc=0; lc<cldout->nlev-1; lc++) {

        /* The un-deltascaled asymmetry parameter should be scaled. So first we 
           unscale the asymmetry parameter */
        unscaled_gg=cldout->optprop.g1[iv][lc]*(1.0 - cldout->optprop.f[iv][lc]) +
          cldout->optprop.f[iv][lc];
        
        /* Now we apply the scaling factor to the un-deltascaled assymetry parameter */
        unscaled_gg *= cldin.modify[MODIFY_VAR_GG][MODIFY_TYPE_SCALE];

        /* And deltascale the result */
	cldout->optprop.g1[iv][lc] = (unscaled_gg -  cldout->optprop.f[iv][lc])/
          (1.0 - cldout->optprop.f[iv][lc]);
	

        /* g2 is nonzero in case of Key/Yang  (not delta-scaled)*/
        cldout->optprop.g2[iv][lc] *=  cldin.modify[MODIFY_VAR_GG][MODIFY_TYPE_SCALE];
        
	/* check if we are still within the limits */
	if (cldout->optprop.g1[iv][lc]>1.0)
	  cldout->optprop.g1[iv][lc] = 1.0;
	
	if (cldout->optprop.g1[iv][lc]<-1.0)
	  cldout->optprop.g1[iv][lc] = -1.0;
	
        
	/* We don't want a user-defined asymmetry parameter if an explicit  */
	/* phase function was defined (e.g. ic_properties mie, baum);       */
	/* see comment above.                                               */
	if (cldout->optprop.nmom[iv][lc] > 1) {
	  fprintf (stderr, "Error, it does not make sense to use an explicit phase function\n");
	  fprintf (stderr, "and then to manually over-write the asymmetry parameter.\n"); 
	  return -1;
	}
      }
    }

    /* Set the cloud single scattering albedo to a user-defined value */
    if (cldin.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SET] >= 0.0)
      for (lc=0; lc<cldout->nlev-1; lc++)
        /* We assume that the user means the un-deltascaled ssa, so ssa is here
           delta-scaled */
        cldout->optprop.ssa[iv][lc] = (1.0 - cldout->optprop.f[iv][lc]) * cldin.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SET] /
          (1.0 - cldin.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SET] * cldout->optprop.f[iv][lc]) ;
    
    /* Scale the cloud single scattering albedo with a user-defined value */
    if (cldin.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SCALE] >= 0.0)
      for (lc=0; lc<cldout->nlev-1; lc++) {
        
        /* Here also we assume that the user wants to scale the un-deltascaled ssa. 
           So we unscale the ssa first */
        unscaled_ssa = cldout->optprop.ssa[iv][lc]/
          (1.0 + cldout->optprop.f[iv][lc] * ( cldout->optprop.ssa[iv][lc] - 1.0));
                                 
        /* Apply scaling factor */
        unscaled_ssa *= cldin.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SCALE];

        /* Delta scale the result */
	cldout->optprop.ssa[iv][lc] = (1.0 - cldout->optprop.f[iv][lc])*unscaled_ssa/
          (1.0 - unscaled_ssa* cldout->optprop.f[iv][lc]);
        
        if (cldout->optprop.ssa[iv][lc]>1.0)
	  cldout->optprop.ssa[iv][lc] = 1.0;
      }
  
  
    /* Set the cloud optical thickness to a constant user-defined value */
    /* (CE) Important: This needs to be done *after* ssa is modified !!!*/ 
    if (cldin.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SET] >= 0.0) {
      
      sum = 0.0;
      
      for (lc=0; lc<cldout->nlev-1; lc++) {
        
        unscaled_ssa = cldout->optprop.ssa[iv][lc]/
          (1.0 + cldout->optprop.f[iv][lc] * ( cldout->optprop.ssa[iv][lc] - 1.0));
        
        if (cldout->zd[lc] > altitude) {    /* layers with upper boundary above the surface */
	  if (cldout->zd[lc+1] > altitude){  /* layers completely above the surface          */
            
            /* Un-deltascale optical thickness. If no delta scaled properties are used 
               f is 0 and the optical thickness is not changed here */
            sum += cldout->optprop.dtau[iv][lc]/
              (1.0 - unscaled_ssa*cldout->optprop.f[iv][lc]) ;
          }
	  else                              /* upper boundary above, lower boundary below surface */
	    sum += cldout->optprop.dtau[iv][lc]/  
              (1.0 - unscaled_ssa*cldout->optprop.f[iv][lc]) * 
              (cldout->zd[lc]-altitude) / (cldout->zd[lc]-cldout->zd[lc+1]);
	}
      }
      
      scale_factor = 0.0;
      if (sum>0)
	scale_factor = cldin.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SET]/sum;
      
      for (lc=0; lc<cldout->nlev-1; lc++)
	cldout->optprop.dtau[iv][lc] *= scale_factor;
      
    } 

    /* Scale the cloud optical thickness (RPB)*/ 
    if (cldin.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE] >= 0.0)
      for (lc=0; lc<cldout->nlev-1; lc++)
	cldout->optprop.dtau[iv][lc] *= cldin.modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE];

  }


  /* Set the cloud optical thickness to a user-defined value at 550nm */
  if (cldin.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET] >= 0.0) {
    lambda = calloc (nlambda_r, sizeof(double));
    tautot = calloc (nlambda_r, sizeof(double));
    
    for (iv=0; iv<nlambda_r; iv++) {
      lambda[iv] = (double) lambda_r[iv];

      for (lc=0; lc<cldout->nlev-1; lc++) {

	unscaled_ssa = cldout->optprop.ssa[iv][lc]/
          (1.0 + cldout->optprop.f[iv][lc] * ( cldout->optprop.ssa[iv][lc] - 1.0));
        
        if (cldout->zd[lc] > altitude) {    /* layers with upper boundary above the surface */
	  if (cldout->zd[lc+1] > altitude)  /* layers completely above the surface          */
            
            /* un-deltascale optical thickness */
	    tautot[iv] += cldout->optprop.dtau[iv][lc] /
              (1.0 - unscaled_ssa *
               cldout->optprop.f[iv][lc]) ;
          else        /* upper boundary above, lower boundary below surface */
            tautot[iv] += cldout->optprop.dtau[iv][lc]/
              (1.0 - unscaled_ssa * cldout->optprop.f[iv][lc]) * 
	      (cldout->zd[lc]-altitude) / (cldout->zd[lc]-cldout->zd[lc+1]);
	}
      }
    }
    
    status = linear_coeffc (lambda, tautot, nlambda_r, &a0, &a1, &a2, &a3);
    if (status!=0) {
      fprintf (stderr, "Error %d calculating interpolation coefficients for cloud optical thickness.\n", 
	       status);
      return status;
    }

    status = calc_splined_value (550.0, &tau550, lambda, nlambda_r, a0, a1, a2, a3);
    if (status!=0) {
      fprintf (stderr, "Error %d interpolating cloud optical thickness to 550nm.\n", 
	       status);
      fprintf (stderr, "In the current implementation of 'wc_modify tau550' and \n");
      fprintf (stderr, "'ic_modify tau550', it is required\n");
      fprintf (stderr, "that the internal wavelength range includes 550nm.\n");
      fprintf (stderr, "The range of the internal wavelength grid is from %f to %f nm. \n", lambda[0], lambda[nlambda_r-1]); 
      fprintf (stderr, "Hope to change that in future!\n");
      return status;
    }
    
    if (tau550>0.0) {
      scale_factor = cldin.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET] / tau550;
      
      for (iv=0; iv<nlambda_r; iv++)
	for (lc=0; lc<cldout->nlev-1; lc++)
	  cldout->optprop.dtau[iv][lc] *= scale_factor;
    }


    free(a0); free(a1); free(a2); free(a3);
    free(lambda); free(tautot);
  }

  return 0;
}


/***********************************************************************************/
/* Function: ic_geom                                                      @62_30i@ */
/* Description: simple geometrical optics calculation of the extinction cross      */
/*              section; tau = 3 IWC / (2 rho_ice reff); geometrical optics cross  */
/*              section sigma_geo = 2 pi r**2 is 2x the geometrical cross section  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Bernhard Mayer                                                          */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int ic_geom (float wavelength, int nlyr,
	     float *iwc, float *reff, float *zd,
	     float *dtau, float *ssa)
{

  int lc=0, status=0;
  float rho_ice = 0.917;

  fprintf (stderr, "ICGEOM\n");
  
  for (lc=0; lc<nlyr; lc++)  {
    dtau[lc] = 3.0 * iwc[lc] / 2.0 / rho_ice / reff[lc] * (zd[lc]-zd[lc+1]) * 1000.0;
    ssa[lc]  = 1.0;

    if (dtau[lc]>0)
      fprintf (stderr, "DTAU=%f, IWC=%f, REFF=%f\n", dtau[lc], iwc[lc], reff[lc]);
  }

  return status;   /* if ok */
}


/*agf*/
int ic_yang (float wavelength, int nlyr, int habit,
	     char *path,
	     float *iwc, float *reff, 
	     float *dtau, 
	     float *f, float *g1, float *g2, float *ssa, 
	     float *zd, int layer, int newkey)
{
  float ext=0;
  int lc=0, status=0;
  char mipath[256]="";
  
  sprintf(mipath,"%s/ic/yang56/",path);

  if (layer==0) {
    fprintf (stderr, "Error, Yang et al. [2000] can only be used with cloud properties\n");
    fprintf (stderr, "defined per layer; please use ic_layer!\n");
    return -1;
  }

  for (lc=0; lc<nlyr; lc++) {
  
    /* calculate dtau[lc], g1[lc], g2[lc], f[lc], ssa[lc] from */
    /* iwc[lc], reff[lc], wavelength                           */

    if (iwc[lc]>0.0) {

      if (newkey) {
#if HAVE_YANG

	status = yang56 (wavelength/1000.0, reff[lc], habit, mipath, &ext, ssa+lc, g1+lc, g2+lc, f+lc, newkey);
	
	if (status!=0) {
	  fprintf (stderr, "Error %d returned by yang56()\n", status);
	  return status;
	}
#else
	fprintf (stderr, " ***********************************************************************\n");
	fprintf (stderr, " * You have built uvspec without Yang/Key/Mayer et al support and      *\n");
	fprintf (stderr, " * hence cannot use 'ic_properties yang'.                              *\n");
	fprintf (stderr, " ***********************************************************************\n");
	return -1;
#endif
      }
      else {

	status = yang56 (wavelength/1000.0, reff[lc], habit, mipath, &ext, ssa+lc, g1+lc, g2+lc, f+lc, newkey);
	
	if (status!=0) {
	  fprintf (stderr, "Error %d returned by yang56()\n", status);
	  return status;
	}
      }

      dtau[lc] = ext*iwc[lc]*(zd[lc]-zd[lc+1]);
    }
  }
  
  return 0;   /* if ok */
}




/*******************************************************************/
/* Ice cloud properties according to                               */
/* Q. Fu, An Accurate Parameterization of the Solar Radiative      */
/* Properties of Cirrus Clouds for Climate Models,                 */
/* Journal of Climate 9, 2058-2082, 1996.                          */
/*******************************************************************/

int ic_fu96 (float wavelength, int nlyr,
	     char *path,
	     float *iwc, float *reff, 
	     float *dtau, float *g, float *ssa, float *f,
	     float *zd, int layer, int unscaled)
{
  int i=0;
  int lc=0, status=0;
  char filename [FILENAME_MAX]="";

  static double *al1=NULL, *al2=NULL, *a0=NULL, *a1=NULL;
  static double *bl1=NULL, *bl2=NULL, *b0=NULL, *b1=NULL, *b2=NULL, *b3=NULL;
  static double *cl1=NULL, *cl2=NULL, *c0=NULL, *c1=NULL, *c2=NULL, *c3=NULL;
  static double *dl1=NULL, *dl2=NULL, *d0=NULL, *d1=NULL, *d2=NULL, *d3=NULL;

  static int na=0, nb=0, nc=0, nd=0;

  static int first=1;

  double D=0;
  double asy=0, ext=0, coalb=0, del=0;
  double fw=0;

  float *deff = calloc (nlyr+1, sizeof(float));

  if (layer==0) {
    fprintf (stderr, "Error, Fu [1996] can only be used with cloud properties\n");
    fprintf (stderr, "defined per layer. Please use ic_layer!\n");
    return -1;
  }

  if (first) {  /* read only once */
    first=0;
    
    /* read extinction data */
    sprintf (filename,"%s/ic/fu96/fu96.ext",path);
    status = read_4c_file (filename, &al1, &al2, &a0, &a1, &na);
    if (status!=0) {
      fprintf (stderr, "Error %d reading Fu [1996] data from %s\n", 
	       status, filename);
      return status;
    }
    
    /* read absorption data */
    sprintf (filename,"%s/ic/fu96/fu96.ssa",path);
    status = read_6c_file (filename, &bl1, &bl2, &b0, &b1, &b2, &b3, &nb);
    if (status!=0) {
      fprintf (stderr, "Error %d reading Fu [1996] data from %s\n", 
	       status, filename);
      return status;
    }
    
    /* read asymmetry parameter data */
    sprintf (filename,"%s/ic/fu96/fu96.asy",path);
    status = read_6c_file (filename, &cl1, &cl2, &c0, &c1, &c2, &c3, &nc);
    if (status!=0) {
      fprintf (stderr, "Error %d reading Fu [1996] data from %s\n", 
	       status, filename);
      return status;
    }

    /* read delta-fraction parameter data */
    sprintf (filename,"%s/ic/fu96/fu96.del",path);
    status = read_6c_file (filename, &dl1, &dl2, &d0, &d1, &d2, &d3, &nd);
    if (status!=0) {
      fprintf (stderr, "Error %d reading Fu [1996] data from %s\n", 
	       status, filename);
      return status;
    }

    /* check if wavelength grids are identical */
    if (na!=nb || na!=nc || na!=nd) {
      fprintf (stderr, "Error, something weird happened to the Fu [1996] data files.\n");
      fprintf (stderr, "Found different wavelength grids.\n");
      return -1;
    }

    for (i=0; i<na; i++)
      if (al1[i] != bl1[i] || al1[i] != cl1[i] || al1[i] != dl1[i]) {
	fprintf (stderr, "Error, something weird happened to the Fu [1996] data files.\n");
	fprintf (stderr, "Found different wavelength grids.\n");
	return -1;
      }
    
    /* convert from micron to nm */
    for (i=0; i<na; i++) {
      al1[i] *= 1000.0;
      al2[i] *= 1000.0;
    }
  }
  
  /* determine the tabulated wavelength interval containing the actual wavelength */
  if (wavelength < al1[0] || wavelength > al2[na-1]) {
    fprintf (stderr, "Error, wavelength %f nm not covered by Fu [1996]\n", wavelength);
    fprintf (stderr, "Allowed range is %.1f - %.1f nm\n", al1[0], al2[na-1]);
    return -1;
  }


  /* Fu et al. requires effective diameter */
  for (lc=0; lc<=nlyr; lc++)
    deff[lc] = 2.0 * reff[lc];

  for (i=0; i<na-1; i++)
    if (wavelength < al2[i])
      break;

  for (lc=0; lc<nlyr; lc++) {
    if (iwc[lc]>0) {
      /* calculate dtau[lc], g[lc], and ssa[lc] from */
      /* iwc[lc], deff[lc], wavelength               */
      
      D = deff[lc];

      if (D < MIN_DEFF_FU96 || D > MAX_DEFF_FU96) {
	fprintf (stderr, "Error, effective radius %f um not covered by Fu [1996]\n",D/2.);
	fprintf (stderr, "Allowed range is %7.3f - %7.3f um\n",
		 MIN_DEFF_FU96/2.0, MAX_DEFF_FU96/2.0);
	return -1;
      }
      
      /* extinction */
      ext = a0[i] + a1[i]/D;
      
      /* single scattering albedo */
      coalb = b0[i] + b1[i] * D + b2[i] * D * D + b3[i] * D * D * D;
      
      /* asymmetry parameter */
      asy = c0[i] + c1[i] * D + c2[i] * D * D + c3[i] * D * D * D;

      /* delta-fraction */
      del = d0[i] + d1[i] * D + d2[i] * D * D + d3[i] * D * D * D;
      
      dtau[lc] = ext*iwc[lc]*(zd[lc]-zd[lc+1])*1000.0;
      ssa [lc] = 1.0-coalb;
      g   [lc] = asy;
     

      /* Fu [1996], Eq. 3.8 */
      f[lc]  = 0.5 / ssa[lc] + del;
      fw     = f[lc] * ssa[lc];

      /* Fu [1996], pg. 2067, left, first paragraph */
      if (f[lc]>g[lc])
	f[lc] = g[lc];
      
      /* delta-scaling, Fu [1996], eqs. A.2a - A.2c */
      if (!unscaled) {
	dtau[lc] = (1.0 - fw) * dtau[lc];
	ssa [lc] = (1.0 - f[lc]) * ssa[lc] / (1.0 - fw);
	g   [lc] = (g[lc] - f[lc]) / (1.0 - f[lc]);
      }
      else 
        /* in this case we need to set the delta-scaling factor to 0 */
	/* because it otherwise might be applied to the data later   */
	/* which is not desired                                      */
	f[lc] = 0.0;

    }
  }
  
  free(deff);

  return 0;   /* if ok */
}

/*******************************************************************/
/* Ice cloud properties according to                               */
/* Q. Fu, P. Yang, and W.B. Sun, An Accurate Parameterization of   */
/* the Infrared Radiative Properties of Cirrus Clouds for Climate  */
/* Models, Journal of Climate 11, 2223-2237, 1998.                 */
/*******************************************************************/

int ic_fu98 (float wavelength, int nlyr,
	     char *path,
	     float *iwc, float *reff, 
	     float *dtau, float *g, float *ssa, 
	     float *zd, int layer)
{
  int i=0, lower=0, upper=0;
  int lc=0, status=0;
  char filename [FILENAME_MAX]="";

  static double *al=NULL, *a0=NULL, *a1=NULL, *a2=NULL;
  static double *bl=NULL, *b0=NULL, *b1=NULL, *b2=NULL, *b3=NULL;
  static double *cl=NULL, *c0=NULL, *c1=NULL, *c2=NULL, *c3=NULL;

  static int na=0, nb=0, nc=0;

  static int first=1;

  double D=0;
  double extlower=0, extupper=0, ext=0;
  double abslower=0, absupper=0, abs=0;
  double asylower=0, asyupper=0, asy=0;

  float *deff = calloc (nlyr+1, sizeof(float));

  if (layer==0) {
    fprintf (stderr, "Error, Fu et al. [1998] can only be used with cloud properties\n");
    fprintf (stderr, "defined per layer. Please use ic_layer!\n");
    return -1;
  }

  if (first) {  /* read only once */
    first=0;
    
    /* read extinction data */
    sprintf (filename,"%s/ic/fu98/fu98.ext",path);
    status = read_4c_file (filename, &al, &a0, &a1, &a2, &na);
    if (status!=0) {
      fprintf (stderr, "Error %d reading Fu et al. [1998] data from %s\n", 
	       status, filename);
      return status;
    }
    
    /* read absorption data */
    sprintf (filename,"%s/ic/fu98/fu98.abs",path);
    status = read_5c_file (filename, &bl, &b0, &b1, &b2, &b3, &nb);
    if (status!=0) {
      fprintf (stderr, "Error %d reading Fu et al. [1998] data from %s\n", 
	       status, filename);
      return status;
    }
    
    /* read asymmetry parameter data */
    sprintf (filename,"%s/ic/fu98/fu98.asy",path);
    status = read_5c_file (filename, &cl, &c0, &c1, &c2, &c3, &nc);
    if (status!=0) {
      fprintf (stderr, "Error %d reading Fu et al. [1998] data from %s\n", 
	       status, filename);
      return status;
    }

    /* check if wavelength grids are identical */
    if (na!=nb || na!=nc) {
      fprintf (stderr, "Error, something weird happened to the Fu et al. [1998] data files.\n");
      fprintf (stderr, "Found different wavelength grids.\n");
      return -1;
    }

    for (i=0; i<na; i++)
      if (al[i] != bl[i] || al[i] != cl[i]) {
	fprintf (stderr, "Error, something weird happened to the Fu et al. [1998] data files.\n");
	fprintf (stderr, "Found different wavelength grids.\n");
	return -1;
      }
    
    /* convert from micron to nm */
    for (i=0; i<na; i++)
      al[i] *= 1000.0;
    
  }

  
  /* determine tabulated wavelengths above and below the actual wavelength */
  if (wavelength < al[0] || wavelength > al[na-1]) {
    fprintf (stderr, "Error, wavelength %f nm not covered by Fu et al. [1998]\n", wavelength);
    fprintf (stderr, "Allowed range is %.1f - %.1f nm\n", al[0], al[na-1]);
    return -1;
  }


  /* Fu et al. requires effective diameter */
  for (lc=0; lc<=nlyr; lc++)
    deff[lc] = 2.0 * reff[lc];


  for (i=0; i<na-1; i++)
    if (wavelength < al[i])
      break;

  lower=i-1;
  upper=i;
    
  for (lc=0; lc<nlyr; lc++) {
  

    if (iwc[lc]>0) {
      /* calculate dtau[lc], g[lc], and ssa[lc] from */
      /* iwc[lc], deff[lc], wavelength               */
      
      D = deff[lc];
      
      if (D < MIN_DEFF_FU98 || D > MAX_DEFF_FU98) {
	fprintf (stderr, "Error, effective radius %f um not covered by Fu et al. [1998]\n", D/2.0);
	fprintf (stderr, "Allowed range is %6.2f - %6.2f um\n",
		 MIN_DEFF_FU98/2.0, MAX_DEFF_FU98/2.0);
	return -1;
      }

      /* extinction */
      extlower = a0[lower] + a1[lower]/D + a2[lower] / (D * D);
      extupper = a0[upper] + a1[upper]/D + a2[upper] / (D * D);
      ext = extlower + (wavelength - al[lower]) / (al[upper]-al[lower]) * (extupper-extlower);
      
      /* absorption */
      abslower = b0[lower] + b1[lower] * D + b2[lower] * D * D + b3[lower] * D * D * D;
      absupper = b0[upper] + b1[upper] * D + b2[upper] * D * D + b3[upper] * D * D * D;
      abs = (abslower + (wavelength - al[lower]) / (al[upper]-al[lower]) * (absupper-abslower)) / D;
      
      /* asymmetry parameter */
      asylower = c0[lower] + c1[lower] * D + c2[lower] * D * D + c3[lower] * D * D * D;
      asyupper = c0[upper] + c1[upper] * D + c2[upper] * D * D + c3[upper] * D * D * D;
      asy = asylower + (wavelength - al[lower]) / (al[upper]-al[lower]) * (asyupper-asylower);
      
      dtau[lc] = ext*iwc[lc]*(zd[lc]-zd[lc+1])*1000.0;
      ssa [lc] = (ext-abs)/ext;
      g   [lc] = asy;
    }
  }

  free(deff);
  
  return 0;   /* if ok */
}


/*******************************************************************/
/* ECHAM4 water cloud optical properties, see                      */
/* Roeckner et al., MPI Report 218, Table 2, page 26.              */
/*******************************************************************/

int wc_echam4 (float wavelength, int nlyr,
	       float *iwc, float *reff, 
	       float *dtau, float *g, float *ssa,
	       float *zd, int layer)
{
  int lc=0;
  double logreff=0, ext=0;

  /* check allowed wavelength range */
  if (wavelength < 200.0 || wavelength > 4900.0) {
    fprintf (stderr, "Error, wavelength %f out if range; allowed is 200 - 4000nm.\n",
	     wavelength);
    return -1;
  }

  if (layer==0) {
    fprintf (stderr, "Error, ECHAM4 can only be used with cloud properties\n");
    fprintf (stderr, "defined per layer; please use ic_layer!\n");
    return -1;
  }

  for (lc=0; lc<nlyr; lc++) {
    if (iwc[lc+1]>0) {
      /* calculate dtau[lc], g[lc], and ssa[lc] from */
      /* iwc[lc], deff[lc], wavelength               */
      
      logreff = log10(reff[lc+1]);

      /* optical thickness */
      if (wavelength<=680) 
	ext = 1.8706058417*pow(reff[lc+1],-1.0756364457);
      else
	ext = 1.9655460426*pow(reff[lc+1],-1.0778999732);

      dtau[lc] = ext*iwc[lc+1]*(zd[lc]-zd[lc+1])*1000.0;
      
      /* single scattering albedo */
      if (wavelength<=680) 
        ssa [lc] = 1.0;
      else
        ssa [lc] = 0.9854369057 
	    + 0.013584242533*logreff 
	    - 0.024856960461*logreff*logreff
	    + 0.0055918314369*logreff*logreff*logreff;

      /* asymmetry parameter */
      if (wavelength<=680) 
	g [lc] = 0.78756640717
	  + 0.10660598895*logreff
	  - 0.031012468401*logreff*logreff;
      else
	g [lc] = 0.79208639802 
	  -0.044930076174*logreff
	  + 0.18980672305*logreff*logreff
	  - 0.082590933352*logreff*logreff*logreff;
    }
  }

  return 0;
}


/*******************************************************************/
/* ECHAM4 ice cloud optical properties, see                        */
/* Roeckner et al., MPI Report 218, Table 2, page 26.              */
/*******************************************************************/

int ic_echam4 (float wavelength, int nlyr,
	       float *iwc, float *reff, 
	       float *dtau, float *g, float *ssa,
	       float *zd, int layer)
{
  int lc=0;
  double logreff=0, ext=0;
  static int first=1;

  /* check allowed wavelength range */
  if (wavelength < 200.0 || wavelength > 4900.0) {
    fprintf (stderr, "Error, wavelength %f out if range; allowed is 200 - 4000nm.\n",
	     wavelength);
    return -1;
  }

  if (layer==0) {
    fprintf (stderr, "Error, ECHAM4 can only be used with cloud properties\n");
    fprintf (stderr, "defined per layer; please use ic_layer!\n");
    return -1;
  }

  if (first) {
    fprintf (stderr, "*** WARNING: Check if the additional factor 0.91 to g should really be applied!\n");
    fprintf (stderr, "*** See Roeckner et al. 1996, MPI-Report 218, pg 26.\n");
    first=0;
  }

  for (lc=0; lc<nlyr; lc++) {
    if (iwc[lc+1]>0) {
      /* calculate dtau[lc], g[lc], and ssa[lc] from */
      /* iwc[lc], deff[lc], wavelength               */
      
      logreff = log10(reff[lc+1]);

      /* optical thickness */
      if (wavelength<=680) 
	ext = 1.9056067426*pow(reff[lc+1],-1.0318784654);
      else
	ext = 2.1666771102*pow(reff[lc+1],-1.0634702711);

      dtau[lc] = ext*iwc[lc+1]*(zd[lc]-zd[lc+1])*1000.0;
      
      /* single scattering albedo */
      if (wavelength<=680) 
        ssa [lc] = 1.0;
      else
        ssa [lc] = 0.98475089485
	    + 0.0053152066002*logreff 
	    - 0.0061150583857*logreff*logreff
	    - 0.0032775655896*logreff*logreff*logreff;

      /* asymmetry parameter */
      if (wavelength<=680) 
	g [lc] = 0.7700034985
	  + 0.19598466851*logreff
	  - 0.11836420885*logreff*logreff
	  + 0.025209205131*logreff*logreff*logreff;
      else
	g [lc] = 0.83631171237
	  - 0.19965998649*logreff
	  + 0.46130320487*logreff*logreff
	  - 0.29719270332*logreff*logreff*logreff
	  + 0.062554483594*logreff*logreff*logreff*logreff;

      /* multiply with 0.91 (see Table 2 of Roeckner et al. 1996, MPI report 218) */
      g[lc] *= 0.91;
    }
  }

  return 0;
}



int interpolate_ssprop_in_lambda (ssprop_struct *src, ssprop_struct *tgt,
				  double *src_lambda, int src_nlambda,
				  float *tgt_lambda, int tgt_nlambda,
				  int nlambda_lower, int nlambda_upper)
{
  int ivn=0, i=0, kmax=0, ip1=0;
  double interp_optprop_weight_0=0.0, interp_optprop_weight[2]={0.,0.}, value[2]={0.,0.}, norm=0.0;
  int ir=0, iv=0, k=0, ip=0, status=0;
  int maxleg=0;

  int **tmp_ntheta=NULL;
  float ***tmp_theta=NULL, ***tmp_phase=NULL;
  double ***tmp_mu=NULL;

  /* if extinction_weight is set the interpolation is weighted with extinction */
  int extinction_weight = 1;

  if (!(src->alloc_moments || src->alloc_explicit)) {
    fprintf (stderr, "Error, not possible to interpolate ssprop because no memory allocated\n");
    return -1;
  }

  /* check source optical properties,                        */
  /* if effective radius grid identical for all wavelengths; */
  /* otherwise cannot do interpolation                       */
  for (iv=0; iv<src_nlambda-1; iv++) {
    
    if (src->nreff[iv] != src->nreff[iv+1]) {
      fprintf (stderr, "Error, effective radius grids in optical property files not identical:\n"); 
      fprintf (stderr, "nreff[%d] = %zd, nreff[%d] = %zd\n", 
	       iv, src->nreff[iv], iv+1, src->nreff[iv+1]);
      return -1;
    }
    
    for (ir=0; ir<src->nreff[iv]; ir++){
      if( src->reff[iv][ir] != src->reff[iv+1][ir]) {
        fprintf (stderr, "Error, effective radius grids in optical property files not identical:\n"); 
        return -1;
      }
    }
  }		 
  
  status = calloc_ssprop_struct (tgt, tgt_nlambda);
  if (status!=0) {
    fprintf (stderr, "Error %d allocating memory for target cloud structure\n", status);
    return status;
  }

  tgt->nphamat     = src->nphamat;
  tgt->type        = src->type;

  /* copy effective radii */
  for (iv=0; iv<tgt_nlambda; iv++) {
    tgt->r0[iv]    = src->r0[0];
    tgt->dr[iv]    = src->dr[0];
    tgt->nreff[iv] = src->nreff[0];
  }

  for (iv=0; iv<tgt_nlambda; iv++) {
    tgt->extinc[iv] = calloc (tgt->nreff[iv], sizeof(double));
    tgt->albedo[iv] = calloc (tgt->nreff[iv], sizeof(double));
    tgt->f[iv]      = calloc (tgt->nreff[iv], sizeof(double));
    tgt->reff[iv]   = calloc (tgt->nreff[iv], sizeof(double));
    
    tgt->legen [iv] = calloc (tgt->nreff[iv], sizeof(float **));
    tgt->nleg  [iv] = calloc (tgt->nreff[iv], sizeof(int));

    tgt->ntheta[iv] = calloc (tgt->nreff[iv], sizeof(int **));
    tgt->theta [iv] = calloc (tgt->nreff[iv], sizeof(float **));
    tgt->mu    [iv] = calloc (tgt->nreff[iv], sizeof(double **));
    tgt->phase [iv] = calloc (tgt->nreff[iv], sizeof(float **));
    
    /* copy effective radius grid */
    /* Since the radius grid should be the same for all wavelength,
       it does not need to be interpolated. */
    for (ir=0; ir<tgt->nreff[iv]; ir++)
      tgt->reff[iv][ir] = src->reff[0][ir];
  }

  /* this part is completely new. It no longer uses arp_wvn2  */
  /* since that was extremely inefficient                     */
  /* the outermost loop is now the new wavelength index j     */
  /* where we initially locate the wavelength index i and     */
  /* calculate weight = ( xn[j] - x[i] ) / ( x[i+1] - x[i] )  */
  /* then all variables are interpolated trivially by using   */
  /* yn[j] = weight * y[i+1] + ( 1. - weight ) * y[i]         */
  /* extinction weighting with weight_beta                    */
  /* for phase function, first sort mus, then interpolate if  */
  /* necessary */

  /* test that initial wavelength grid is ascending */
  for (iv=0; iv<src_nlambda-1; iv++)
    if (src_lambda[iv]>=src_lambda[iv+1]) {
      fprintf (stderr, "wavelength not ascending %d %f %f\n",iv, src_lambda[iv], src_lambda[iv+1]);
      return -1;
    }

  for (ivn=nlambda_lower; ivn<=nlambda_upper; ivn++) {
    /* Chris: Hotfix here because of numerical error in the range 470-479nm */
    if ( tgt_lambda[ivn] < src_lambda[0]-1.E-10 || tgt_lambda[ivn] > src_lambda[src_nlambda-1]+1.E-10 ) {
      fprintf (stderr, "internal wavelength not within external wavelength grid %d %f [%f - %f]\n",
               ivn, tgt_lambda[ivn], src_lambda[0], src_lambda[src_nlambda-1]);
      return -1;
    }

    if (src_nlambda>1) {

      /* locate position within external wavelength grid */
      i = locate ( src_lambda, src_nlambda, (double) tgt_lambda[ivn] );
      if ( i == src_nlambda )
	i = src_nlambda-1;
      ip1 = i+1;

      /* calculate weight */
      interp_optprop_weight_0 = ( tgt_lambda[ivn] - src_lambda[i] ) / ( src_lambda[ip1] - src_lambda[i] );
    }
    else {
      /* special case src_nlambda==0 */
      /* this is a simple trick to handle this special case, it causes unneccessarry overhead, but simplifies code reading */
      i = 0;
      ip1 = i;
      interp_optprop_weight_0 = 0.;
    }

    /* interpolate extinction */
    for (ir=0; ir<src->nreff[0]; ir++) {
      interp_optprop_weight[0] = 1. - interp_optprop_weight_0;
      interp_optprop_weight[1] =      interp_optprop_weight_0;
      tgt->extinc[ivn][ir] = interp_optprop_weight[0] * src->extinc[i][ir] + interp_optprop_weight[1] * src->extinc[ip1][ir];

      /* vanishing extinction makes interpolation of other variables obsolete */
      if ( tgt->extinc[ivn][ir] == 0. )
	continue;

      /* weights weighted by extinction */
      if (extinction_weight) {
	interp_optprop_weight[0] *= src->extinc[i  ][ir] / tgt->extinc[ivn][ir];
	interp_optprop_weight[1] *= src->extinc[ip1][ir] / tgt->extinc[ivn][ir];
      }

      /* interpolate albedo */
      tgt->albedo[ivn][ir] = interp_optprop_weight[0] * src->albedo[i][ir] + interp_optprop_weight[1] * src->albedo[ip1][ir];

      /* vanishing albedo makes interpolation of other variables obsolete */
      if ( tgt->albedo[ivn][ir] == 0. )
	continue;

      /* weights weighted by albedo */
      if (extinction_weight) {
	interp_optprop_weight[0] *= src->albedo[i  ][ir] / tgt->albedo[ivn][ir];
	interp_optprop_weight[1] *= src->albedo[ip1][ir] / tgt->albedo[ivn][ir];
      }

      if (src->alloc_moments) {

	/* determine maximum number of Legendre coefficients */
	/* for this effective radius                         */
	maxleg=0;
	for (iv=0; iv<src_nlambda; iv++)
	  if (src->nleg[iv][ir]>maxleg)
	    maxleg=src->nleg[iv][ir];

	/* get index of maximum nonzero moment of phase function */
	/* NOTE: we could save memory if we set this for each phase matrix element independently XXX */
	/* NOTE2: we might speed this up using src->nleg XXX */
	kmax = 0;
	for (iv=i; iv <= ip1; iv++)
	  for (ip=0; ip < src->nphamat; ip++){
	    for (k=maxleg-1; k>=0; k--)  
	      if (src->legen[iv][ir][ip][k]!=0)
		break;
	    if (k > kmax)
	      kmax = k;
	  }
	
	tgt->nleg[ivn][ir] = kmax+1;
	
	/* allocate memory for moments */
	tgt->legen[ivn][ir] = calloc (src->nphamat, sizeof(float *));
	for (ip=0; ip<src->nphamat; ip++)
	  tgt->legen[ivn][ir][ip]  = calloc (tgt->nleg[ivn][ir], sizeof(float));

	/* interpolate moments */
	for (ip=0; ip<src->nphamat; ip++)
	  for (k=0; k<tgt->nleg[ivn][ir]; k++) {
	    if (src->nleg[i][ir] < k)
	      value[0] = 0.;
	    else
	      value[0] = src->legen[i][ir][ip][k];
	    if (src->nleg[ip1][ir] < k)
	      value[1] = 0.;
	    else
	      value[1] = src->legen[ip1][ir][ip][k];

	    tgt->legen[ivn][ir][ip][k] =
	      interp_optprop_weight[0] * value[0] + interp_optprop_weight[1] * value[1];
	  }

	tgt->alloc_moments = 1; /* memory has been allocated */
      }

      if (src->alloc_explicit) {

	/* allocate phase dimension */
	tgt->ntheta [ivn][ir] = calloc (tgt->nphamat, sizeof(int));
	tgt->theta  [ivn][ir] = calloc (tgt->nphamat, sizeof(float *));
	tgt->mu     [ivn][ir] = calloc (tgt->nphamat, sizeof(double *));
	tgt->phase  [ivn][ir] = calloc (tgt->nphamat, sizeof(float *));

	/* allocate temporary theta as large as can possibly get */
	tmp_ntheta = calloc(src->nphamat, sizeof(int *));
	tmp_theta = calloc(src->nphamat, sizeof(float **));
	tmp_mu = calloc(src->nphamat, sizeof(double **));
	tmp_phase = calloc(src->nphamat, sizeof(float **));

	for (ip=0; ip<src->nphamat; ip++) {
	  tmp_ntheta[ip] = calloc(2, sizeof(int));
	  tmp_theta [ip] = calloc(2, sizeof(float *));
	  tmp_mu    [ip] = calloc(2, sizeof(double *));
	  tmp_phase [ip] = calloc(2, sizeof(float *));

	  tmp_ntheta[ip][0] = src->ntheta[i  ][ir][ip];
	  tmp_ntheta[ip][1] = src->ntheta[ip1][ir][ip];

	  tmp_theta [ip][0] = src->theta [i  ][ir][ip];
	  tmp_theta [ip][1] = src->theta [ip1][ir][ip];

	  tmp_mu    [ip][0] = src->mu    [i  ][ir][ip];
	  tmp_mu    [ip][1] = src->mu    [ip1][ir][ip];

	  tmp_phase [ip][0] = src->phase [i  ][ir][ip];
	  tmp_phase [ip][1] = src->phase [ip1][ir][ip];
	}

	/* sort theta values and interpolate */
	status = sort_and_add_weighted_phase (2, interp_optprop_weight,
					      tmp_ntheta, tmp_theta, tmp_mu, tmp_phase,
					      &(tgt->ntheta [ivn][ir]),
					      &(tgt->theta[ivn][ir]),
					      &(tgt->mu[ivn][ir]),
					      &(tgt->phase[ivn][ir]),
					      tgt->nphamat, 0, 1 );

	if (status) {
	  fprintf(stderr,"Error! something (%d) went wrong in sort_and_add_weighted_phase \n",status);
	  return status;
	}

	for (ip=0; ip<src->nphamat; ip++) {
	  free(tmp_ntheta[ip]);
	  free(tmp_theta[ip]);
	  free(tmp_mu[ip]);
	  free(tmp_phase[ip]);
	}
	free(tmp_ntheta);
	free(tmp_theta);
	free(tmp_mu);
	free(tmp_phase);

	tgt->alloc_explicit = 1; /* memory has been allocated */

      }

      /* interpolate scaling factor */
      if (extinction_weight) {
	interp_optprop_weight[0] = ( 1. - interp_optprop_weight_0 ) * src->extinc[i  ][ir] * src->albedo[i  ][ir] / ( 1. - src->f[i  ][ir]);
	interp_optprop_weight[1] =        interp_optprop_weight_0   * src->extinc[ip1][ir] * src->albedo[ip1][ir] / ( 1. - src->f[ip1][ir]);
	norm = interp_optprop_weight[1] + interp_optprop_weight[0];

	/* prevent division by zero */
	if ( norm == 0. )
	  continue;

	interp_optprop_weight[0] /= norm;
	interp_optprop_weight[1] /= norm;
      }

      tgt->f[ivn][ir] = interp_optprop_weight[0] * src->f[i][ir] + interp_optprop_weight[1] * src->f[ip1][ir];

      /* To avoid rounding errors for f=0.0 */
      if (tgt->f[ivn][ir]<1e-8)
	tgt->f[ivn][ir]=0.0;
 
    } /* end loop ir */
  } /* end loop ivn */

  return 0;
}


int read_isccp_reflectivity (int type, float sza, float tau, char *path, int quiet,
			     float *ref)
{
  int i=0;
  int u=0, l=0;
  int utau=0, ltau=0, usza=0, lsza=0;
  double ref1=0, ref2=0;
  static int read[2]={0,0};
  static size_t ntau=0, nsza=0;
  
  static double *szas=NULL, *taus=NULL;
  static double **refs=NULL;

  char filename[FILENAME_MAX]="";

#if HAVE_LIBNETCDF
  /* try to open the CDF file */

  int ncid=0, status=0;

  int idd_nsza=0, idd_ntau=0;
  int id_sza=0, id_tau=0, id_ref=0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  size_t refstart[2] = {0,0};
  size_t refcount[2] = {0,0};
#endif
 

  if (!read[type]) {
    
    strcpy (filename, path);
    
    switch (type) {
    case ISCCP_WATER:
      strcat(filename, "wc/isccp/wc_reflectivity.cdf");
      break;
    case ISCCP_ICE:
      strcat(filename, "ic/isccp/ic_reflectivity.cdf");
      break;
    default:
      fprintf (stderr, "Error, unknown ISCCP cloud type\n");
      return -1;
    }

    /* read reflectivity file */

#if HAVE_LIBNETCDF

    /* open netcdf file */
    status = nc_open (filename, NC_NOWRITE, &ncid);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return -1;
    }

    /* read file */
    
    /* get dimension id for "sza" */
    status = nc_inq_dimid (ncid, "sza", &idd_nsza);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    /* get dimension length for "sza" */
    status = nc_inq_dimlen (ncid, idd_nsza, &nsza);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    /* get dimension id for "tau" */
    status = nc_inq_dimid (ncid, "tau", &idd_ntau);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    /* get dimension length for "tau" */
    status = nc_inq_dimlen (ncid, idd_ntau, &ntau);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }

    /* allocate memory for fields */
    szas    = calloc (nsza, sizeof(double));
    taus    = calloc (ntau, sizeof(double));
    
    refs    = calloc (nsza, sizeof(double *));
    for (i=0; i<nsza; i++)
      refs[i] = calloc (ntau, sizeof(double));
    
    /* get variable id for "sza" */
    status = nc_inq_varid (ncid, "sza", &id_sza);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
      }
    
    /* get variable id for "tau" */
    status = nc_inq_varid (ncid, "tau", &id_tau);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    /* get variable id for "ref" */
    status = nc_inq_varid (ncid, "ref", &id_ref);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    /* read "sza" */
    count[0] = nsza;
    status = nc_get_vara_double (ncid, id_sza, start, count, szas);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    /* read "tau" */
    count[0] = ntau;
    status = nc_get_vara_double (ncid, id_tau, start, count, taus);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    /* read reflectivities */
    refstart[1] = 0;  /* start with first reflectivity    */
    refcount[0] = 1;  /* read one reflectivity at a time  */
      
    for (i=0; i<nsza; i++) {
      
      refstart[0] = i;
      refcount[1] = ntau;
      
      status = nc_get_vara_double (ncid, id_ref, refstart, refcount, refs[i]);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading %s\n", status, filename);
	return status;
      }
    }
    
    /* close netcdf file */
    nc_close(ncid);

    if (!quiet)
      fprintf (stderr, " ... read ISCCP reflectivities from %s\n", filename);

#else
    fprintf (stderr, "Error, need netcdf for the sssi solver\n");
    fprintf (stderr, "Please install netcdf and re-compile!\n");
    return -1;
#endif

    /************************************************************/
    /* now we have the cloud reflectivity data stored in arrays */
    /* szas[isza], taus[itau], and refs[isza][itau];            */
    /* next step is the interpolation of the reflectivity to    */
    /* the user-defined (sza, tau)                              */
    /************************************************************/

    read[type] = 1;
  }  

  /* check if values out of range */
  if (tau<taus[0])  {
    fprintf (stderr, "Error, tau smaller than minimum in reflectivity table, %f\n", taus[0]);
    return -1;
  }

  if (tau>taus[ntau-1])  {
    fprintf (stderr, "Error, tau larger than maximum in reflectivity table, %f\n", taus[ntau-1]);
    return -1;
  }

  if (sza<szas[0])  {
    fprintf (stderr, "Error, sza smaller than minimum in reflectivity table, %f\n", szas[0]);
    return -1;
  }

  if (sza>szas[nsza-1])  {
    fprintf (stderr, "Error, sza larger than maximum in reflectivity table, %f\n", szas[nsza-1]);
    return -1;
  }

  /* now search for sza in sza array */
  l = 0;
  u = nsza-1;
  while (1==1) {
    if (u <= l+1)
      break;

    /* midpoint */
    i=(l+u)/2;

    if (sza==szas[i])
      u=l=i;
    else {
      if (sza<szas[i])
	u = i;
      else
	l = i;
    }
  }

  lsza=l;
  usza=u;

  /* now search for tau in tau array */
  l = 0;
  u = ntau-1;
  while (1==1) {
    if (u <= l+1)
      break;

    /* midpoint */
    i=(l+u)/2;

    if (tau==taus[i])
      u=l=i;
    else {
      if (tau<taus[i])
	u = i;
      else
	l = i;
    }
  }

  ltau=l;
  utau=u;

  /* l and u are the indices (upper and lower) of the */
  /* neighbouring data points; now do simply a linear */
  /* interpolation ...                                */

  /* first, interpolation in sza */
  if (lsza==usza) {
    ref1 =  refs[lsza][ltau];
    ref2 =  refs[lsza][utau];
  }
  else {
    ref1 =  refs[lsza][ltau] + (refs[usza][ltau]-refs[lsza][ltau]) / (szas[usza] - szas[lsza]) * (sza - szas[lsza]);
    ref2 =  refs[lsza][utau] + (refs[usza][utau]-refs[lsza][utau]) / (szas[usza] - szas[lsza]) * (sza - szas[lsza]);
  }

  /* second, interpolation in tau */
  if (ltau==utau)
    *ref = ref1;
  else
    *ref = ref1 + (ref2-ref1) / (taus[utau] - taus[ltau]) * (tau - taus[ltau]);
  
  return 0;
}


/***********************************************************************************/
/* Function: read_ECMWF_clouds                                                     */
/*                                                                                 */
/* Description:                                                                    */
/*  * read cloud data cloud water content (cwc), level heights (zd)                */ 
/*    and cloud fraction (cf) from ECMWF file                                      */
/*  * convert cwc (kg/kg) to cloud water density (cwd) (g/m3)                      */
/*  * parametrisation of the effective radius                                      */
/*  * check for inconsistent data (cf == 0) and cwc != 0                           */
/*    set cwc in this case to 0 (also for cloudoverlapp off) for sensitivity       */
/*    studies (this should not mask the effect of overlap change)                  */
/*                                                                                 */
/*                                                                                 */
/* Parameters (input):                                                             */
/*     filename          ECMWF data file                                           */
/*     file_ic_reff      data file with ice cloud effective radius data            */
/*     lat               latitude in deg N                                         */
/*     lon               longitude in deg E                                        */
/*     UTC               time in UTC                                               */
/*     time_interpolate  switch for time interpolation                             */
/*     cloud             either 'wc' or 'ic' indicate type of cloud                */
/*     reff_prop         id_number for effective radius parametrisation            */
/*     ic_prop           id_number for optical parametrisation                     */
/*     reff_fixed        fixed effective radius (if needed)                        */
/*     altitude          surface height                                            */
/*     press_atm         pressure profile of the atmosphere file                   */
/*     cloud_overlap     id_number for overlap schema                              */
/*     verbose           switch for feedback to user                               */
/*     quiet             switch for feedback to user                               */
/*                                                                                 */
/* Parameters (output):                                                            */
/*     row               number of rows of cloud data                              */
/*     max_columns       number of columns of cloud data                           */
/*     min_columns       number of columns of cloud data                           */
/*     cld_data          cld_data[*][0] == z     [km]                              */
/*                       cld_data[*][1] == cwd   [kg/kg]                           */
/*                       cld_data[*][0] == reff  [micro meter]                     */
/*     cf                cloud fraction structure                                  */
/*                                                                                 */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    cloud.c                                                               */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Jan 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int read_ECMWF_clouds ( char            *filename,
			char            *file_ic_reff,
			float            lat,
			float            lon,
			struct tm        UTC,
			int              time_interpolate, 
			char            *cloud, 
			int              reff_prop,
			int              ic_prop,
			float            reff_fixed,
			float            altitude,
			float           *press_atm,
			int              cloud_overlap, 
			int              verbose,
			int              quiet,
			int             *rows,
			int             *max_columns,
			int             *min_columns, 
			float         ***cld_data,
			cf_out_struct   *cf )
{
  int status=0;

#if HAVE_MYSTIC

#if HAVE_LIBNETCDF

  int ncid    = NOT_DEFINED_INTEGER;
  int ncid2   = NOT_DEFINED_INTEGER;
  int id_data = NOT_DEFINED_INTEGER;

  char T_name    [2]="";
  char CLWC_name [5]="";
  char CIWC_name [5]="";
  char CC_name   [4]="";
  char phase     [8]="";

  int lc = NOT_DEFINED_INTEGER; /* loop index */
  int t  = NOT_DEFINED_INTEGER; /* loop index */
  int nt = NOT_DEFINED_INTEGER;
 
  size_t nlat  = NOT_DEFINED_INTEGER;
  size_t nlon  = NOT_DEFINED_INTEGER;
  size_t nlev  = NOT_DEFINED_INTEGER;
  size_t nlay  = NOT_DEFINED_INTEGER;  

  size_t ilat  = NOT_DEFINED_INTEGER;   /* index for lat, lon, time in netCDF file */
  size_t ilon  = NOT_DEFINED_INTEGER;   /* index for lat, lon, time in netCDF file */
  size_t itime = NOT_DEFINED_INTEGER;   /* index for lat, lon, time in netCDF file */
  long ilat2, ilon2;

  int itime1 = NOT_DEFINED_INTEGER, itime2 = NOT_DEFINED_INTEGER;

  double *ECMWF_lat  = NULL;
  double *ECMWF_lon  = NULL;

  float dt = NOT_DEFINED_FLOAT;
  float dx = NOT_DEFINED_FLOAT;

  float **p_level = NULL;
  float **p_layer = NULL;
  float **T       = NULL;
  float  *T_level = NULL;
  float **CWC     = NULL; /* cloud water content, liquid or ice in kg/kg */

  float **CC      = NULL; /* cloud cover, water and ice together */

  float *p_int    = NULL; /* for integration of height */
  float *T_int    = NULL; /* for integration of height */
  float *zd_layer = NULL; /* layer heights (important for reff interpolation) */

  float *z        = NULL;
  float *reff     = NULL;
  float Deff = NOT_DEFINED_FLOAT;
  float D = NOT_DEFINED_FLOAT, L = NOT_DEFINED_FLOAT, ou2fu = NOT_DEFINED_FLOAT;
  float *z_reff   = NULL;
  int   nreff = NOT_DEFINED_INTEGER;
  float *dens_air = NULL; /* number density of air moleculs in moleculs per cm^3 */
  float factor    = NOT_DEFINED_FLOAT; /* unit conversion factor */

  float c0=326.3,c1=12.42,c2=0.197,c3=0.0012;

  /* float min_cf = 0.0; */

  float T_in_C=-1;

  //20120816ak first_layers is done is not used, commented
  //  int first_layer_done=FALSE;

  float fu2yang = 3.0*sqrt(3.0)/4.0;

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) { 
    fprintf (stderr, "Error '%s' opening netCDF file '%s' (line %d, function %s in %s)\n", 
                      nc_strerror(status), filename, __LINE__, __func__, __FILE__);
    return -abs(status);
  }

  /* check format */
  if      ( (status = nc_inq_varid (ncid, "t", &id_data)) == NC_NOERR ) {
    strcpy (T_name,    "t");
    strcpy (CLWC_name, "clwc");
    strcpy (CIWC_name, "ciwc");
    strcpy (CC_name,   "cc");
  }
  else if ( (status = nc_inq_varid (ncid, "T", &id_data)) == NC_NOERR ) {
    strcpy (T_name,    "T");
    strcpy (CLWC_name, "CLWC");
    strcpy (CIWC_name, "CIWC");
    strcpy (CC_name,   "CC");
  }
  else {
    fprintf (stderr, "Error '%s' while getting id for temperature from '%s' \n", nc_strerror(status), filename);
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -1;
  }

  /* read latitude array */
  alloc_and_read_netCDF_1D_double(ncid, "lat", &nlat, "lat", &(ECMWF_lat));

  /* search correct latitude index */
  ilat = get_grid_index( lat, ECMWF_lat, nlat, FALSE);
  if (ilat < 0) {
    fprintf (stderr, "Error -1 finding index for lat=%5.2f in %s (line %d, function %s in %s)\n", lat, filename, __LINE__, __func__, __FILE__);
    return -1;
  }

  /* read longitude */
  alloc_and_read_netCDF_1D_double(ncid,"lon", &nlon, "lon", &(ECMWF_lon));

  /* search correct latitude index */
  ilon = get_grid_index( lon, ECMWF_lon, nlon, TRUE);
  if (ilon < 0) {
    fprintf (stderr, "Error -2 finding index for lon=%5.2f in %s (line %d, function %s in %s)\n", lon, filename, __LINE__, __func__, __FILE__);
    return -2;
  }

  if (verbose)
    fprintf (stderr, " *** ECMWF %s cloud data at lat = %5.2f (%5.2f), lon = %5.2f (%5.2f) \n", 
                     cloud, ECMWF_lat[ilat], lat, ECMWF_lon[ilon], lon);

  /* get time index */
  status = get_time_index (ncid, UTC, time_interpolate, 
                           &(nt), &(itime1), &(itime2), &(dt),
                           verbose, quiet);
  if (status != 0) {
    fprintf (stderr, "Error %d, during get_time_index (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }
  /* alloc nt timesteps for pressure */
  if ((p_level = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for pressure (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -10;
  }
  /* alloc nt timesteps for pressure */
  if ((p_layer = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for pressure (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -11;
  }
  /* alloc nt timesteps for temperature */
  if ((T = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for temperature (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -12;
  }
  /* alloc nt timesteps for temperature */
  if ((CWC = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for Cloud Water Content (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -13;
  }
  /* alloc nt timesteps for temperature */
  if ((CC = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for Cloud Fraction (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -14;
  }

  /* read nt (= 1 or 2) time steps */
  for (t=0;t<=nt-1;t++) {

    if (t == 0)
      itime = itime1;
    if (t == 1)
      itime = itime2;

    /* alloc and read pressure */
    /* requires that SP, hyai, hybi, hyam, and hybm are in the ECMWF netCDF file */
    alloc_and_read_ECMWF_netCDF_pressure (ncid, &(p_level[t]), &(nlev), &(p_layer[t]), &(nlay), itime, ilat, ilon, verbose);

    /* nlev == nlay, as we throw the uppermost level (p=0.0) away, in order to combine the ECMWF data with background data */
    /* !!! not any more !!! */

    /* allocate and read temperature */
    status = alloc_and_read_netCDF_column_float(ncid, T_name, &(T[t]), nlay, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading temperature '%s' from netCDF file %s (line %d, function %s in %s)\n", 
                        status, T_name, filename, __LINE__, __func__, __FILE__);
     return status;
    }

    /* allocate and read cloud water */
    if      (strncasecmp(cloud,"wc",2)==0) {
      strcpy (phase,  "liquid"); 
      /* allocate and read CLWC */
      status = alloc_and_read_netCDF_column_float(ncid, CLWC_name, &CWC[t], nlay, itime, ilat, ilon, verbose);
      if (status != 0) {
        fprintf (stderr, "Error %d reading Cloud Liquid Water Content '%s' from netCDF file %s (line %d, function %s in %s)\n", 
                         status, CLWC_name, filename, __LINE__, __func__, __FILE__);
        return status;
      }
    }
    else if (strncasecmp(cloud,"ic",2)==0) {
      strcpy (phase,  "ice"); 
      /* allocate and read CIWC */
      status = alloc_and_read_netCDF_column_float(ncid, CIWC_name, &CWC[t], nlay, itime, ilat, ilon, verbose);
      if (status != 0) {
        fprintf (stderr, "Error %d reading Cloud Ice Water Content '%s' from netCDF file %s (line %d, function %s in %s)\n", 
                         status, CIWC_name, filename, __LINE__, __func__, __FILE__);
        return status;
      }
    }
    else {
      fprintf(stderr,"Unknown cloud type '%s' (line %d, function %s in %s)\n", cloud, __LINE__, __func__, __FILE__ );
      return -1; 
    }

    /* allocate and read CloudCover == CloudFraction */
    /* if ( cloud_overlap != CLOUD_OVERLAP_OFF ) { */
      status = alloc_and_read_netCDF_column_float(ncid, CC_name, &CC[t], nlay, itime, ilat, ilon, verbose);
      if (status != 0) {
        fprintf (stderr, "Error %d reading Cloud Cover '%s' from netCDF file %s (line %d, function %s in %s)\n", 
                         status, CC_name, filename, __LINE__, __func__, __FILE__);
        return status;
      }
    /*}*/
  }
  nc_close (ncid);

  if (nt == 1) {
    /* no time interpolation needed, do nothing */
  }
  else {
    /* write time interpolated data into the zero'th entry */
    for (lc=0; lc<nlev; lc++)
      p_level[0][lc] = (1.0-dt)* p_level[0][lc]  +  dt* p_level[1][lc];
    for (lc=0; lc<nlay; lc++) {
      p_layer[0][lc] = (1.0-dt)* p_layer[0][lc]  +  dt* p_layer[1][lc];
      T      [0][lc] = (1.0-dt)* T      [0][lc]  +  dt* T      [1][lc];
      CWC    [0][lc] = (1.0-dt)* CWC    [0][lc]  +  dt* CWC    [1][lc];
      /* if ( cloud_overlap != CLOUD_OVERLAP_OFF ) */
      CC     [0][lc] = (1.0-dt)* CC     [0][lc]  +  dt* CC     [1][lc];
    }
  }

/*   /\* additional verbose output data on layers *\/ */
/*   fprintf (stderr, " *** ADDITIONAL VERBOSE OUTPUT: ECMWF clouds on layers \n"); */
/*   fprintf (stderr,"  lc   p_level    p_layer     p_level      T         CWC    \n"); */
/*   for (lc=0; lc<nlay; lc++) { */
/*     fprintf (stderr," %3d %10.4f %10.4f %10.4f %10.5f %12.6e", */
/*       lc+1, p_level[0][lc+1], p_layer[0][lc], p_level[0][lc], T[0][lc], CWC[0][lc]); */
/*     if ( cloud_overlap != CLOUD_OVERLAP_OFF )     fprintf (stderr,", CC=%12.6e", CC[0][lc]); */
/*     fprintf (stderr,"\n"); */
/*   } */

  /* allocate temperature data on levels */
  if ((T_level = calloc (nlev, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for T_layer (line %d, function %s in %s)\n", 
                    __LINE__, __func__, __FILE__);
    return -30;
  }

  /* inter- and extra-polate from LAYER to LEVEL values for of T on the log(p/p0) grid */
  /* interpolation */
  for (lc=1; lc<nlev-1; lc++) {
    dx = ( log ( p_level[0][lc] / p_level[0][nlev-1]) - log ( p_layer[0][lc-1] / p_level[0][nlev-1])) /
         ( log ( p_layer[0][lc] / p_level[0][nlev-1]) - log ( p_layer[0][lc-1] / p_level[0][nlev-1]) );
    T_level[lc] = (1.0-dx) *  T[0][lc-1] + dx *  T[0][lc];
  }
  /* extrapolate last layer midpoint to surface */
  dx = ( log ( p_level[0][nlev-1] / p_level[0][nlev-1]) - log ( p_layer[0][nlay-1] / p_level[0][nlev-1]) ) / 
       ( log ( p_layer[0][nlay-1] / p_level[0][nlev-1]) - log ( p_layer[0][nlay-2] / p_level[0][nlev-1]) );
  T_level[nlev-1] = (1.0+dx) *  T[0][nlay-1] - dx *  T[0][nlay-2];
  /* assume constant temperature between last layer midpoint and TOA */
  T_level[  0   ] = T [0][0];

  /* if uppermost pressure level is 0, than replace it by uppermost background atmosphere pressure */
  if ( p_level[0][0] == 0.0 )
    p_level[0][0] = press_atm[0];

  /* calculate z from T and p using hydrostatic equation (in ancillary.c) */
  status = calculate_z_from_p_and_T (p_level[0], T_level, &(z), nlev, altitude, verbose);
  if (status!=0) {
    fprintf (stderr, "Error %d calculating z-grid of ECMWF cloud data (line %d, function %s in %s) \n", 
                     status, __LINE__, __func__, __FILE__);
    return status;
  }

/*   for (lc=0; lc<nlev; lc++) */
/*     fprintf (stderr," ### %3d cloud p=%e T=%e, z=%e\n", lc, p_level[0][lc], T_level[lc], z[lc] ); */

  if ( z[0] > 120.0 ) /* not nice, assumes standard atm. */
    z[0] = 120.0;

  /* /\* additional verbose output data on levels *\/ */
  /* fprintf (stderr, " ### ECMWF p and T on levels for integration of z ### \n"); */
  /* for (lc=0; lc<nlev; lc++) { */
  /*   fprintf (stderr,"lc=%3d, p=%10.4f, T=%10.5f, z=%10.5f \n", lc, p_level[0][lc], T_level[lc], z[lc]); */
  /* } */

  /* calculate zd_layer for interpolation of effective radius */
  /* alloc nlay+1 timesteps for pressure */
  if ((p_int = calloc (nlay+1, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for p_int (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -10;
  }
  /* alloc nlay+1 timesteps for pressure */
  if ((T_int = calloc (nlay+1, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for T_int (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -10;
  }

  /* copy data */
  for (lc=0; lc<nlay; lc++) {
    p_int[lc] = p_layer[0][lc];
    T_int[lc] = T[0][lc];
  }
  p_int[nlay] = p_level[0][nlev-1];
  T_int[nlay] = T_level[nlev-1];

  /* calculate z from T and p using hydrostatic equation (in ancillary.c) */
  status = calculate_z_from_p_and_T (p_int, T_int, &(zd_layer), nlay+1, altitude, verbose);
  if (status!=0) {
    fprintf (stderr, "Error %d calculating zd_layer-grid of ECMWF cloud data (line %d, function %s in %s) \n", 
                     status, __LINE__, __func__, __FILE__);
    return status;
  }

  free(T_int);
  free(p_int);

  /* allocate and read cloud water */
  if      (strncasecmp(cloud,"wc",2)==0) {

    nreff=4;

    /* alloc 2 levels for height of cloud droplet effective radius */
    if ((z_reff = calloc (nreff, sizeof (float *))) == NULL) {
      fprintf (stderr,"Error, allocating memory for 'z_reff' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -17;
    }
    /* alloc 2 levels for cloud droplet effective radius */
    if ((reff = calloc (nreff, sizeof (float *))) == NULL) {
      fprintf (stderr,"Error, allocating memory for 'reff' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -16;
    }

    reff[0]=45.0; /* 45 micro meter at top of atmosphere */
    reff[1]=45.0; /* 45 micro meter at uppermost layer midpoint */
    reff[2]=10.0; /* 10 micro meter at surface + 10m */
    reff[3]=10.0; /* 10 micro meter at surface */

    if (z[0] > zd_layer[0]) {     /* if top of atmosphere is above 72.0km (top of atmosphere of the ECMWF), than ... */
      z_reff[0] = z[0];        /* reff = 45 from uppermost level */
      z_reff[1] = zd_layer[0]; /*           to top of atm. of the ECMWF */
    }
    else {
      z_reff[0] = zd_layer[0]+0.01; /* dummy value */
      z_reff[1] = zd_layer[0];      /* reff = 45 in 72 km height */
    }
    /* z_reff[1] = zd_layer[0]; */
    z_reff[2] = z[nlev-1] + 0.010; /* lowermost layer (10m) reff = 10 micro m */
    z_reff[3] = z[nlev-1];

    /* linear interpolation of the effective radius to the z-scale of the cloud data */
    status = interpolate_profile (z_reff, &(reff), nreff, zd_layer, nlay, INTERP_METHOD_LINEAR, quiet);  /* !! ATTENTION POINTER ARITHMETIC !! */
    if (status!=0) {                                                                                          /* last entry is lowermost level */
      fprintf (stderr, "Error %d during interpolate_profile of r_eff (line %d, function %s in %s) \n", 
               status, __LINE__, __func__, __FILE__);
      return status;
    }

    /* for (lc=0; lc<nlay; lc++) */
    /*   fprintf (stderr, " ### z = %8.4f ... %8.4f, reff = %8.5f  \n", z[lc], z[lc+1], reff[lc] ); */

    free(z_reff);
  }
  else if (strncasecmp(cloud,"ic",2)==0) {
    switch (reff_prop) {
    case REFF_FIXED:

      if (verbose) fprintf(stderr," ... fixed r_eff_ice = %f\n", reff_fixed);

      /* alloc 2 levels for height of cloud droplet effective radius */
      if ((z_reff = calloc (2, sizeof (float *))) == NULL) {
        fprintf (stderr,"Error, allocating memory for 'z_reff' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -17;
      }
      /* alloc 2 levels for cloud droplet effective radius */
      if ((reff = calloc (2, sizeof (float *))) == NULL) {
        fprintf (stderr,"Error, allocating memory for 'reff' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -16;
      }

      reff[0] = reff_fixed;
      reff[1] = reff_fixed;

      /* height scale for effective radius range from surface to top of atmosphere */
      z_reff[0]=z[0];
      z_reff[1]=z[nlev-1];

      /* linear interpolation of the effective radius to the z-scale of the cloud data */
      status = interpolate_profile (z_reff, &(reff), 2, zd_layer, nlay, INTERP_METHOD_LINEAR, quiet);
      if (status!=0) {
        fprintf (stderr, "Error %d during interpolate_profile of r_eff (line %d, function %s in %s) \n", 
                 status, __LINE__, __func__, __FILE__);
        return status;
      }

      free(z_reff);

      break;
    case REFF_OU:
      if (verbose) fprintf(stderr," ... parametrisation for r_eff_ice of Ou and Liou, Atm. Res. 1995, 'Ice microphysics ...'\n");

      /* alloc nlev levels for cloud droplet effective radius */
      if ((reff = calloc (nlay, sizeof (float *))) == NULL) {
        fprintf (stderr,"Error, allocating memory for 'reff' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -16;
      }

      /* if (verbose) { */
      /*   fprintf(stderr,"-------------------------------------------------------------------------------------------\n"); */
      /*   fprintf(stderr,"  lc    z[km]    z[km]     T[C]   CWC[kg/kg]    r_eff(Ou)  r_eff(Fu)   ou2fu       CC      \n"); */
      /*   fprintf(stderr,"-------------------------------------------------------------------------------------------\n"); */
      /* } */

      for (lc=0; lc<nlay; lc++) {
        T_in_C = T[0][lc] - 273.16;
        /* parametrisation valid between -60 to -20 degree Celsius */
        if ( -60.0 <= T_in_C && T_in_C <= -20.0 )
          Deff = c0 + c1*T_in_C + c2*T_in_C*T_in_C + c3*T_in_C*T_in_C*T_in_C;
        else if ( T_in_C < -60.0 )
          Deff = c0 + c1*(-60.0) + c2*(-60.0)*(-60.0) + c3*(-60.0)*(-60.0)*(-60.0);    /* below -60 C -> deff = const */
        else if ( T_in_C > -20.0 )
          Deff = c0 + c1*(-20.0) + c2*(-20.0)*(-20.0) + c3*(-20.0)*(-20.0)*(-20.0);    /* above -20 C -> deff = const */
        
        /* relationship between Deff and reff from Wyser, 1998, JC, The effective Radius in Ice Clouds, equ. (34) */
        /* reff[lc] = -2.2054 + 0.56383*Deff + 5.6416E-3*Deff*Deff - 3.0954E-5*Deff*Deff*Deff + 1.2601E-7*Deff*Deff*Deff*Deff; */
        /* fprintf(stderr," ... effective size d_eff[%3d] = %e, r_eff =%e, T=%7.2f C (%7.2f K) \n", lc, Deff, reff[lc], T_in_C, T[0][lc]); */
        reff[lc] = 0.5 * Deff;

        /* convert to YANG definition of the effective radius */
        D = 2*reff[lc]/1000.0;    /* width in micro meter */
        L = pow(D/0.185,1/0.53);  /* length in micro meter, see Heymmsfield, 1972, JAS 29, 1358 - 1366 */
                                  /* hmmmm only valid for L>0.3 mm == D>0.98 mue and bullet crystals */
        ou2fu = D*L / ( D*L + sqrt(3.)/4.*D*D ); /* !!! THIS IS AN APPROXIMATION !!! NEEDS TO BE IMPROVED !!! */

        reff[lc] *= ou2fu*fu2yang;

        /*  if (verbose) */
        /*    fprintf(stderr," %3d  %8.4f %8.4f  %7.3f  %e  %9.4f  %9.4f  %8.6f  %e \n", */
        /*                     lc, z[lc+1], z[lc], T_in_C, CWC[0][lc], reff[lc]/ou2fu, reff[lc], ou2fu, CC[0][lc]); */
      }
      break;
    case REFF_FILE:
      if (verbose) fprintf (stderr, " ... read r_eff_ice from file %s\n", file_ic_reff );

      if (nt == 1) {

        status = get_all_netCDF_indices (file_ic_reff, lat, lon,
                                         &(ncid), &(ilat2), &(nlat), &(ilon2), &(nlon), 
                                         &(ECMWF_lat), &(ECMWF_lon),
                                         UTC, time_interpolate,
                                         &(nt), &(itime1), &(itime2), &(dt),
                                         verbose, quiet);

        /* open netcdf file */
        status = nc_open (file_ic_reff, NC_NOWRITE, &ncid2);
        if (status!=NC_NOERR) { 
          fprintf (stderr, "Error '%s' opening netCDF file '%s' (line %d, function %s in %s)\n", 
                   nc_strerror(status), file_ic_reff, __LINE__, __func__, __FILE__);
          return -abs(status);
        }

        /* allocate and read effective radius for ice cloud cristals */
        status = alloc_and_read_netCDF_column_float(ncid2, "reff_ice", &reff, nlay, itime, ilat2, ilon2, verbose);
        if (status != 0) {
          fprintf (stderr, "Error %d reading Ice Cloud effective radius '%s' from netCDF file\n   %s (line %d, function %s in %s)\n", 
                   status, "reff_ice", file_ic_reff, __LINE__, __func__, __FILE__);
          return status;
        }

        /* close netCDF file */
        nc_close (ncid2);
      }
      else {
        fprintf (stderr, "Error, ECMWF_ic_reff file not implemented for 2 time steps \n");
        return -1;
      }
      break;
    default:
      fprintf (stderr, "Error, unknown reff_property %d (line %d, function %s in %s)\n", reff_prop, __LINE__, __func__, __FILE__);
      return -1;
    }

    free(zd_layer);


    /* check reff limits of the parametrisations */ 
    switch (ic_prop) {
    case PROP_FU:
      if (verbose) fprintf (stderr, " ... check, if effective radius is inside limits of the Fu parameterisation [%f,%f] \n", 9.316, 64.799*fu2yang );
      for (lc=0; lc<nlay; lc++) {
        /* if (verbose) fprintf (stderr, "     check: ic_reff[%3d] = %f \n", lc, reff[lc] ); */
        /* reff_max for Fu */
        if (reff[lc] > (64.799*fu2yang)) {
          if (!quiet && CWC[0][lc] != 0.0) {
            fprintf (stderr, " *** Warning, reff_ice[%3d] = %f outside Fu parameterisation [%f,%f] \n", lc, reff[lc], 9.316, 64.799*fu2yang );
            fprintf (stderr, "    reff_ice[%3d] = %f, IWC = %e ->", lc, reff[lc], CWC[0][lc] );
          }
          /* extinction coefficient ~ IWC / reff = const -> IWC' = IWC reff / reff' */
          CWC[0][lc] = CWC[0][lc] * (64.799*fu2yang) / reff[lc];        
          reff[lc]   = 64.799*fu2yang;
          if (!quiet && CWC[0][lc] != 0.0) fprintf (stderr, " changed to reff_ice = %f IWC = %e \n", reff[lc], CWC[0][lc] );
        }
        if (reff[lc] < (9.316*fu2yang)) {
          if (!quiet && CWC[0][lc] != 0.0) {
            fprintf (stderr, " *** Warning, reff_ice[%3d] = %f outside Fu parameterisation [%f,%f] \n", lc, reff[lc], 9.316, 64.799*fu2yang);
            fprintf (stderr, "    reff_ice[%3d] = %f, IWC = %e \n", lc, reff[lc], CWC[0][lc] );
          }
          CWC[0][lc] = CWC[0][lc] * (9.316*fu2yang) / reff[lc];       
          reff[lc]   = 9.316*fu2yang;
          if (!quiet && CWC[0][lc] != 0.0) fprintf (stderr, " changed to reff_ice = %f, IWC = %e \n", reff[lc], CWC[0][lc] );
        }
      }
      break;
    case PROP_FILE:
    case PROP_KEY:
    case PROP_YANG:
    case PROP_ECHAM4:
    case PROP_IC_MIE:
      /* no check here */
      break;
    case PROP_BAUM:
    case PROP_HEY:
    case PROP_YANG2013:
      if (verbose) fprintf (stderr, " ... check, if effective radius is inside limits of the Baum parametrisation [%f,%f] \n", 5.0, 90.0 );
      for (lc=0; lc<nlay; lc++) {
        if (reff[lc] > 90.0) {
          if (!quiet && CWC[0][lc] != 0.0) {
            fprintf (stderr, " *** Warning, reff_ice[%3d] = %f outside Baum parametrisation [%f,%f] \n", lc, reff[lc], 5.0, 90.0 );
            fprintf (stderr, "    reff_ice[%3d] = %f, IWC = %e \n", lc, reff[lc], CWC[0][lc] );
          }
          CWC[0][lc] = CWC[0][lc] * 90.0 / reff[lc];     
          reff[lc] = 90.0 ;
          if (!quiet && CWC[0][lc] != 0.0) fprintf (stderr, " changed to reff_ice = %f, IWC = %e \n", reff[lc], CWC[0][lc] );
        }
        if (reff[lc] < 5.0 ) {
          if (!quiet && CWC[0][lc] != 0.0) {
            fprintf (stderr, " *** Warning, reff_ice[%3d] = %f outside Baum parametrisation [%f,%f] \n", lc, reff[lc], 5.0, 90.0 );
            fprintf (stderr, "    reff_ice[%3d] = %f, IWC = %e \n", lc, reff[lc], CWC[0][lc] );
          }
          CWC[0][lc] = CWC[0][lc] * 5.0 / reff[lc];     
          reff[lc]   = 5.0;
          if (!quiet && CWC[0][lc] != 0.0) fprintf (stderr, " changed to reff_ice = %f, IWC = %e \n", reff[lc], CWC[0][lc] );
        }
      }
      break;
    case PROP_BAUM_V36:
      if (verbose) fprintf (stderr, " ... check, if effective radius is inside limits of the Baum (v3.6) parametrisation [%f,%f] \n", 5.0, 60.0 );
      for (lc=0; lc<nlay; lc++) {
        if (reff[lc] > 60.0) {
          if (!quiet && CWC[0][lc] != 0.0) {
            fprintf (stderr, " *** Warning, reff_ice[%3d] = %f outside Baum (v 3.6) parametrisation [%f,%f] \n", lc, reff[lc], 5.0, 60.0 );
            fprintf (stderr, "    reff_ice[%3d] = %f, IWC = %e \n", lc, reff[lc], CWC[0][lc] );
          }
          CWC[0][lc] = CWC[0][lc] * 60.0 / reff[lc];     
          reff[lc] = 60.0 ;
          if (!quiet && CWC[0][lc] != 0.0) fprintf (stderr, " changed to reff_ice = %f, IWC = %e \n", reff[lc], CWC[0][lc] );
        }
        if (reff[lc] < 5.0 ) {
          if (!quiet && CWC[0][lc] != 0.0) {
            fprintf (stderr, " *** Warning, reff_ice[%3d] = %f outside Baum parametrisation [%f,%f] \n", lc, reff[lc], 5.0, 60.0 );
            fprintf (stderr, "    reff_ice[%3d] = %f, IWC = %e \n", lc, reff[lc], CWC[0][lc] );
          }
          CWC[0][lc] = CWC[0][lc] * 5.0 / reff[lc];     
          reff[lc]   = 5.0;
          if (!quiet && CWC[0][lc] != 0.0) fprintf (stderr, " changed to reff_ice = %f, IWC = %e \n", reff[lc], CWC[0][lc] );
        }
      }
    default:
      fprintf (stderr, "Error, unknown ice cloud property %d. This is not\n", ic_prop);
      fprintf (stderr, "       supposed to happen and indicates a coding error.\n");
      return -1;
    }
  }
  else {
    fprintf(stderr,"Unknown cloud type '%s' (line %d, function %s in %s)\n", cloud, __LINE__, __func__, __FILE__ );
    return -1; 
  }


  /* write numbers to the output structures */
  (*rows)        = nlev;
  (*max_columns) = 3;
  (*min_columns) = 3;
  /* if ( cloud_overlap != CLOUD_OVERLAP_OFF ) */
    cf->nlev     = nlev;
  /* else */ 
  /*  cf->nlev     = 0; */

  /* allocate array, which will be returned */
  if ((status = ASCII_calloc_float (cld_data,(*rows),3)) != 0) {
    fprintf (stderr," Error allocating cld_data (line %d, function %s in %s)\n",__LINE__, __func__, __FILE__);   
    return status;
  }

  factor = 1000.0 * 1.e+6 * 1.e-3 * MOL_MASS_AIR / AVOGADRO;
     /* (kg water)/(kg air) (ECMWF) -> (g water) / (m3 air) (libRadtran) */
     /* 1000. == convert cloud water content  from (kg/kg -> g/kg)       */
     /* 1.e+6 == convert air number density   from (1/cm3 -> 1/m3)       */
     /* MOL_MASS_AIR == weight in g / mol particles                      */
     /* AVOGADRO == Avogadro constant == 1 mol Teilchen                  */
     /* 1.e-3 == convert air mass   density   from (   g  ->  kg )       */

  /* write uppermost level as a boundary for a possible cloud in the first layer */
  (*cld_data)[0][0] = z[0];
  (*cld_data)[0][1] = 0.0;
  (*cld_data)[0][2] = 0.0;

  /* if (verbose) { */
  /*   fprintf (stderr,"#-----------------------------------------------------------------------------\n"); */
  /*   fprintf (stderr,"# lc    z[km]    z[km]    CWC[kg/kg]  rho_air[g/m3]   CWD[g/m3]     reff[mue] \n"); */
  /*   fprintf (stderr,"#-----------------------------------------------------------------------------\n"); */
  /* } */

  /* alloc nlev timesteps for dens_air */
  if ((dens_air = calloc (nlev, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'dens_air' (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -10;
  }

  for (lc=1; lc<nlev; lc++) {
    (*cld_data)[lc][0] = z[lc];                                           /* z(level) */
    dens_air[lc-1] = p_layer[0][lc-1] * 1E-04 / (BOLTZMANN * T[0][lc-1]); /* INDEX SHIFT FROM LAYER TO LEVEL PROP */ 
                                                                          /* dens_air is number concentration of air moleculs per cm^3 */
    (*cld_data)[lc][1] = CWC[0][lc-1] * factor * dens_air[lc-1];          /* INDEX SHIFT FROM LAYER TO LEVEL PROP */  
                                                                          /* cld_data[g/m3] = CWC[kg/kg](layer) * factor * dens_air(layer) */
    (*cld_data)[lc][2] = reff[lc-1];                                

    /* if (verbose) { */
    /*   fprintf (stderr," %3d  %8.4f  %8.4f  %e  %11.5e  %e  %e\n",  */
    /*                     lc, z[lc], z[lc-1], CWC[0][lc-1], dens_air[lc-1]*factor, (*cld_data)[lc][1], (*cld_data)[lc][2]); */
    /* } */
  }

  free(reff);

  /*if ( cloud_overlap != CLOUD_OVERLAP_OFF ) {*/

    if ( ( cf->cf = calloc(cf->nlev, sizeof(double)) ) == NULL) {
      fprintf (stderr,"Error, allocating memory for cloud fraction (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -10;
    }
    if ( ( cf->zd = calloc(cf->nlev, sizeof(double)) ) == NULL) {
      fprintf (stderr,"Error, allocating memory for cloud fraction (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -10;
    }

    /* uppermost layer boundary */
    cf->zd[0] = z[0];
    cf->cf[0] = 0.0;

    for (lc=1; lc<nlev; lc++) {
      cf->zd[lc] = z[lc];
      cf->cf[lc] = CC[0][lc-1]; /* shifted by one, as we write a layer property to a z-level */
      /* fprintf(stderr," %3d  zd = %8.5f, cf->cf=%e, cf->cf=%e  \n", lc, z[lc], CC[0][lc-1], cf->cf[lc]); */
    }

    /* ATTENTION NO FINAL DECISION HERE JET WHAT TO DO, IF CWC != 0 AND CF == 0 */
    /* THROW AWAY CWC OR ASSUME A RESONABE VALUE FOR CF ??? */

/*     /\* check consistence of input data: cloud fraction == 0  <=> cloud water content == 0 *\/ */
/*     for (lc=0; lc<nlev; lc++) */
/*       if ( cf->cf[lc] == 0.0 && (*cld_data)[lc][1] != 0.0 ) { */
/*         min_cf = 10.0; */
/*         if ( lc !=      0 ) if (cf->cf[lc-1] != 0.0)                           min_cf = cf->cf[lc-1]; */
/*         if ( lc != nlev-1 ) if (cf->cf[lc+1] != 0.0 && cf->cf[lc+1] <  min_cf) min_cf = cf->cf[lc+1]; */
/*         if ( min_cf == 10.0) min_cf = 0.001; */
/*         fprintf (stderr, " *** Warning, cloud fraction[%3d]=%e, but %s cloud water density = %e g/m3 \n", */
/*                            lc, cf->cf[lc], phase, (*cld_data)[lc][1] ); */
/*         fprintf (stderr, "     Set cloud fraction to %f \n", min_cf ); */
/*         cf->cf[lc] = min_cf ; */
/*       } */

    for (lc=0; lc<nlev; lc++) {
      if ( cf->cf[lc] != 0.0 ) {
        /* fprintf (stderr, " read_ECMWF_clouds: cf[%3d] = %f,  no scaling any more here \n", lc, cf->cf[lc]); */
        if ( cloud_overlap == CLOUD_OVERLAP_OFF ) cf->cf[lc] = 1.0; 
      }
      else if ( (*cld_data)[lc][1] != 0.0 ) { 
        if (!quiet) {
          fprintf (stderr, " *** Warning, cloud fraction[%3d]=%e, but %s cloud water density = %e g/m3 \n", 
                           lc, cf->cf[lc], phase, (*cld_data)[lc][1] );
          fprintf (stderr, "     Set cloud %s water content to 0.0 g/m3 \n", phase);
        }
        (*cld_data)[lc][1] = 0.0;
      }
    }

    /*}*/


  free(z);

  /* verbose output */ 
  if (verbose) {
    //20120816ak first_layers is done is not used, commented
    //    first_layer_done=FALSE;

    if ( cloud_overlap != CLOUD_OVERLAP_OFF ) {
      fprintf (stderr, "#-------------------------------------------------------------------------------------------------------\n");
      fprintf (stderr, "# lc |    z     |    z    |  CWD_cloud |  r_eff   |   rho_air  |  CWC_cloud   |     CC    |  CWC_mean  |\n");
      fprintf (stderr, "#    |   [km]   |   [km]  |   [g/m3]   |  [mue]   |   [g/m3]   |   [kg/kg]    |           |   [kg/kg]  |\n");
      fprintf (stderr, "#-------------------------------------------------------------------------------------------------------\n");
    }
    else {
      fprintf (stderr, "#------------------------------------------------------------------------------\n");
      fprintf (stderr, "# lc |    z     |    z    |  CWD_cloud |  r_eff   |   rho_air  |  CWC_cloud   |\n");
      fprintf (stderr, "#    |   [km]   |   [km]  |   [g/m3]   |  [mue]   |   [g/m3]   |   [kg/kg]    |\n");
      fprintf (stderr, "#------------------------------------------------------------------------------\n");
    }

    for (lc=1; lc<nlev; lc++) {
      if ( lc < nlay-1 ) { /* check also CWC one level below (lc+1) */
        /*if ( (*cld_data)[lc][1] != 0.0 || (*cld_data)[lc+1][1] != 0.0 ) { */ /* || first_layer_done == TRUE */
           fprintf (stderr, "%5d %9.5f %9.5f  %e %8.4f  %e  %e ",
                    lc, (*cld_data)[lc][0], (*cld_data)[lc-1][0], (*cld_data)[lc][1], (*cld_data)[lc][2], dens_air[lc-1]*factor,
                    (*cld_data)[lc][1]/(dens_air[lc-1]*factor));

           if ( cloud_overlap != CLOUD_OVERLAP_OFF )
             fprintf (stderr, " %e %e ", cf->cf[lc], (*cld_data)[lc][1] / (dens_air[lc-1]*factor) * cf->cf[lc] );

           fprintf (stderr, "\n");
           /* fprintf (stderr,"mass dens air [kg/m3] = %e, water content [g/kg] = %e \n",  */
           /*          dens_air[lc] * 1.e+6 * 1.e-3 * air_mol_mass / AVOGADRO, CWC[lc]*1000.0); */
	   //20120816ak first_layers is done is not used, commented
	   //           first_layer_done=TRUE;
        /*}*/
      }
      else { /* lowermost layer */
        /*/if ( (*cld_data)[lc][1] != 0.0 ) { */  /* || first_layer_done == TRUE */
           fprintf (stderr, "%5d %9.5f %9.5f  %e %8.4f  %e  %e ",
                    lc, (*cld_data)[lc][0], (*cld_data)[lc-1][0], (*cld_data)[lc][1], (*cld_data)[lc][2], dens_air[lc-1]*factor,
                    (*cld_data)[lc][1]/(dens_air[lc-1]*factor));
           if ( cloud_overlap != CLOUD_OVERLAP_OFF )
             fprintf (stderr, " %e %e ", cf->cf[lc], (*cld_data)[lc][1] / (dens_air[lc-1]*factor) * cf->cf[lc] );
           fprintf (stderr, "\n");
	   //20120816ak first_layers is done is not used, commented
	   //           first_layer_done=TRUE;
        /*}*/
      }
    }
  }

  /* for (lc=0; lc<nlev; lc++) { */
  /*   fprintf (stderr, "read_ECMWF_clouds: %5d %9.5f %10.4f   %e %8.4f ", */
  /*            lc, (*cld_data)[lc][0], p_layer[0][lc-1], (*cld_data)[lc][1], (*cld_data)[lc][2] ); */
  /*   if ( cloud_overlap != CLOUD_OVERLAP_OFF ) */
  /*     fprintf (stderr, " %11.7f ", cf->cf[lc] ); */
  /*   fprintf (stderr, "\n"); */
  /* } */

  for (t=0;t<=nt-1;t++) {
    free(p_level[t]);
    free(p_layer[t]);
    free(T[t]);
    free(CWC[t]);
    free(CC[t]);
  }
  free(T_level);
  free(p_level);
  free(p_layer);
  free(dens_air);
  free(T);
  free(CWC);
  free(CC);

  /* need to shift by one layer as cloud fractions are automatically considered layer properties */
  for (lc=0;lc<cf->nlev-1;lc++) {
    cf->cf[lc] = cf->cf[lc+1];
  }

#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use the ECMWF data file option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ******************************************************************\n");
  return -1;
#endif

#else
  fprintf (stderr, " Error, ECMWF cloud input is not available in this verion \n");
#endif

  return status;

}

/********************************************************************************************************/

int read_ECHAM_clouds ( char            *filename,
			float            lat,
			float            lon,
			struct tm        UTC,
			int              time_interpolate, 
			char            *cloud,
			float           *zd,
			int              verbose,
			int              quiet,
			int             *rows,
			int             *max_columns,
			int             *min_columns,
			float         ***cld_data,
			cf_out_struct   *cf )
{
  int status=0;

#if HAVE_LIBNETCDF

  int ncid=NOT_DEFINED_INTEGER;

  int lc = NOT_DEFINED_INTEGER; /* loop index */
  int t  = NOT_DEFINED_INTEGER; /* loop index */
  int nt = NOT_DEFINED_INTEGER;
 
  size_t nlat  = NOT_DEFINED_INTEGER;
  size_t nlon  = NOT_DEFINED_INTEGER;
  size_t nlev  = 39; /* constant for ECHAM data files; must be changed for 
                        ECHAM5 */
  
  size_t ilat  = NOT_DEFINED_INTEGER;  
  size_t ilon  = NOT_DEFINED_INTEGER;  
  size_t itime = NOT_DEFINED_INTEGER;  

  int itime1 = NOT_DEFINED_INTEGER, itime2 = NOT_DEFINED_INTEGER;

  double *ECHAM_lat  = NULL;
  double *ECHAM_lon  = NULL;
  
  float dt = NOT_DEFINED_FLOAT;
  
  float **CWP     = NULL; /* cloud water content, liquid or ice */
  float *CWC      = NULL;
  float **reff    = NULL;
  float **splitfactor  = NULL;
  float **totcc3D = NULL;
  
  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, filename);
    return status;
  }

  /* read latitude array */
  alloc_and_read_netCDF_1D_double(ncid,"lat", &nlat, "lat", &(ECHAM_lat));

  /* search correct latitude index */
  ilat = get_grid_index( lat, ECHAM_lat, nlat, FALSE);
  if (ilat < 0) {
    fprintf (stderr, "Error -1 finding index for lat=%5.2f in %s\n", lat, filename);
    return -1;
  }
  
  /* read longitude */
  alloc_and_read_netCDF_1D_double(ncid,"lon", &nlon, "lon", &(ECHAM_lon));
  
  /* search correct latitude index */
  ilon = get_grid_index( lon, ECHAM_lon, nlon, TRUE);
  if (ilon < 0) {
    fprintf (stderr, "Error -2 finding index for lon=%5.2f in %s\n", lon, filename);
    return -2;
  }

  /* get time index */
  status = get_time_index (ncid, UTC, time_interpolate,
                           &(nt), &(itime1), &(itime2), &(dt),
                           verbose, quiet);
  if (status != 0) {
    fprintf (stderr, "Error %d, during get_time_index in get_float_from_netCDF_map (ancillary.c)\n", status);
    return status;
  }

  /* alloc nt timesteps for cloud cover profile */
  if ((totcc3D = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for pressure in read_ECHAM_clouds() (in cloud.c)\n");
    return -10;
  }
  
  /* alloc nt timesteps for cloud water content */
  if ((CWP = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for cloud water path in read_ECHAM_clouds() (in cloud.c)\n");
    return -13;
  }
  
  if ((splitfactor = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for splitfactor in read_ECHAM_clouds() (in cloud.c)\n");
    return -11;
  }
  
  /* alloc nt timesteps for temperature */
  if ((reff = calloc (nt, sizeof (float *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for effective radius in read_ECHAM_clouds() (in cloud.c)\n");
    return -13;
  }
  
  /* read nt (= 1 or 2) time steps */
  for (t=0;t<=nt-1;t++) {
    
    if (t == 0)
      itime = itime1;
    if (t == 1)
      itime = itime2;
    
    status = alloc_and_read_netCDF_column_float(ncid, "totcc3D", &totcc3D[t], nlev, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading cloud cover from netCDF file %s\n", status, filename);
      return status;
    }
    
  status = alloc_and_read_netCDF_column_float(ncid, "iwp", &CWP[t], nlev, itime, ilat, ilon, verbose);
  if (status != 0) {
      fprintf (stderr, "Error %d reading CWC from netCDF file %s\n", status, filename);
      return status;
    }
    
    status = alloc_and_read_netCDF_column_float(ncid, "splitfactor", &splitfactor[t], nlev, itime, ilat, ilon, verbose);
    if (status != 0) {
      fprintf (stderr, "Error %d reading splitfactor from netCDF file %s\n", status, filename);
      return status;
    }

    if (strncasecmp(cloud,"wc",2)==0) {
      status = alloc_and_read_netCDF_column_float(ncid, "reffwat", &(reff[t]), nlev, itime, ilat, ilon, verbose);
      if (status != 0) {
        fprintf (stderr, "Error %d reading effective radius from netCDF file %s\n", status, filename);
        return status;
      } 
    }
    else if (strncasecmp(cloud,"ic",2)==0) {
      status = alloc_and_read_netCDF_column_float(ncid, "reffice", &(reff[t]), nlev, itime, ilat, ilon, verbose);
      if (status != 0) {
        fprintf (stderr, "Error %d reading effective radius from netCDF file %s\n", status, filename);
        return status;
      }
    }
    else {
      fprintf(stderr,"Unknown 'cloud option %s' in read_ECHAM_clouds \n", cloud);
      return -1; 
    }
  }
  
  /* time interpolation only if nt==2, for nt==1 do nothing */
  if (nt == 2) {
    /* write time interpolated data into the zero'th entry */
    for (lc=0; lc<nlev; lc++) {
      totcc3D[0][lc] = (1.0-dt)*totcc3D[0][lc] + dt*totcc3D[1][lc];
      CWP[0][lc] = (1.0-dt)*CWP[0][lc] + dt* CWP[1][lc];
      splitfactor[0][lc] = (1.0-dt) * splitfactor[0][lc] + 
        dt * splitfactor[1][lc];
      reff[0][lc] = (1.0-dt)*reff[0][lc] + dt*reff[1][lc];
    }
  }
 
  (*rows)        = nlev;
  (*max_columns) = 3;
  (*min_columns) = 3;
  
  if ((status = ASCII_calloc_float (cld_data,(*rows),3)) != 0) {
    fprintf (stderr," Error allocating cld_data in read_ECHAM_cloud (cloud.c) \n");   
    return status;
  }

  cf->cf = calloc(nlev, sizeof(double));
  cf->zd = calloc(nlev, sizeof(double));
  cf->zd[0]= zd[0];
  (*cld_data)[0][0] = zd[0];

  for (lc=1; lc<nlev; lc++){
    /* Cloud fraction, here already converted to layer property */
    cf->cf[lc] = totcc3D[0][lc+1];
   
    cf->zd[lc] = zd[lc];
    /* Altitude levels */
    (*cld_data)[lc][0] = zd[lc];
    /* Effective radius*/ 
    (*cld_data)[lc][2] = reff[0][lc];
    if (strncasecmp(cloud,"wc",2)==0) /*LWC*/ 
      (*cld_data)[lc][1] = (CWP[0][lc] * (1.0-splitfactor[0][lc]))/
        (zd[lc-1] - zd[lc]) / 1000.0;
    else  /* IWC*/
        (*cld_data)[lc][1] = (CWP[0][lc]*splitfactor[0][lc])/
          (zd[lc-1] - zd[lc])/1000.0;
    
  }
  cf->nlev=nlev;
 

  /* verbose output */ 
  if (verbose) {

    fprintf (stderr, " *** ECHAM %s cloud data at lat = %5.2f, lon = %5.2f \n", cloud, lat, lon);
    fprintf (stderr, "-----------------------------------------------\n");
    fprintf (stderr, "  lc |   z   | cloud water |   r_eff   |   cf  \n");
    fprintf (stderr, "     |  [km] |    [g/m3]   | [micro m] |       \n");
    fprintf (stderr, "-----------------------------------------------\n");

    for (lc=0; lc<nlev; lc++) {
      fprintf (stderr, "%5d %7.3f  %e  %7.2f  %7.2f \n",
               lc, (*cld_data)[lc][0], (*cld_data)[lc][1],
               (*cld_data)[lc][2], cf->cf[lc-1]);
    }
  }
  
  for (t=0;t<=nt-1;t++) {
    free(totcc3D[t]);
    free(reff[t]);
    free(splitfactor[t]);
    free(CWP[t]);
  }
  
  free(CWC);

#else
  fprintf (stderr, " ******************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot       *\n");
  fprintf (stderr, " * use the ECHAM data file option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ******************************************************************\n");
  return -1;
#endif

  return status;

}


/* copy cloud fraction structure */
int copy_cloud_fraction ( cf_out_struct *target, cf_out_struct source, int alloc )
{

  int status=0;
  int lc=0;

  if (alloc) {
    status = alloc_cloud_fraction( target, source.nlev );
    if (status!=0) {
      fprintf (stderr, "Error: Allocation of cf->cf (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__);
      return -1;
    }
  }

  for (lc=0;lc<target->nlev;lc++) {
    target->cf[lc] = source.cf[lc];
    target->zd[lc] = source.zd[lc];
  }

  /* we do not copy total cloud fraction here */
  /* as this function is used to copy ipa columns with cloud cover 1 inside the column */
  /* to the general cp structure with cloud cover somewhere between 0 and 1, */
  /* which we like to have for the output user */

  return status;

}

int alloc_cloud_fraction (cf_out_struct *cf, int nlev)
{
  int status=0;

  cf->nlev=nlev;

  if ((cf->cf = calloc(nlev, sizeof(float)))==NULL) {
    fprintf (stderr, "Error: Allocation of cf->cf (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__);
    return -1;
  }
  if ((cf->zd = calloc(nlev, sizeof(float)))==NULL) {
    fprintf (stderr, "Error: Allocation of cf->zd (line %d, function %s in %s) \n", __LINE__, __func__, __FILE__);
    return -1;
  }

  return status;

}

/***********************************************************************************/
/* Function: caoth_prop_switch                                            @62_30i@ */
/* Description: this subroutine is called once for wc (standard) and once          */
/*              for wctipa(-structure) to calculate tau for tipa dir; the          */
/*              subroutine contains all switches for caoth properties              */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Ulrike Wissmeier                                                        */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int caoth_prop_switch ( input_struct      input,
			caoth_inp_struct  input_caoth,
			wl_out_struct     wl_out,
			int               iv,
			/* Output */
			caoth_out_struct *output_caoth )
{
  int    status=0, nlyr=0, lev=0, newkey=0;
  int    nstring = (int) strlen(input.filename[FN_PATH]);
  
  /*ulrike: wcloud-function was copied from cloud.c; declaration of a subroutine */
  void  F77_FUNC (wcloud, WCLOUD) (float *lambda_r, int *newsiz, int *nlyr,
				   char *filepath, int *nstring, float *lwc, float *wceffr, 
				   float *tmp_wc_dtau, float *tmp_wc_gg, float *tmp_wc_ssa,
				   float *zd, int *wclyr);

  switch (input_caoth.properties) {
  case PROP_HU:
    if (input.ipa3d) { /* We are not sure whether this if statement is correct; ignoring for now */
      switch (input.ck_scheme) {
      case CK_FU:
	/* call Fu and Liou code by Fred Rose, ulrike: copied from cloud.c */
#if HAVE_FULIOU
	status = ckdfucld ( output_caoth->zd, 
			    output_caoth->microphys.effr_layer, 
			    output_caoth->microphys.lwc_layer, 
			    output_caoth->nlev,
			    output_caoth->optprop.dtau, 
			    output_caoth->optprop.g1, 
			    output_caoth->optprop.ssa );
	if (status)
	  return fct_err_out ( status, "ckdfucld", ERROR_POSITION );
#else
	fprintf (stderr, "Error, Fu and Liou not supported!\n");
	return -1;
#endif
      break;
      
    default:   /* old uvspec wcloud.f for Hu and Stamnes */
              
	nlyr = output_caoth->nlev-1;

	/* ulrike: set optprop to zero before calling wcloud! */
	for (lev=0; lev<nlyr; lev++) {
	  output_caoth->optprop.dtau[iv][lev] = 0.0;
	  output_caoth->optprop.g1[iv][lev]   = 0.0;
	  output_caoth->optprop.ssa[iv][lev]  = 0.0;
	}

	/* wcloud calculates lwc_layer and effr_layer internally */
	F77_FUNC (wcloud, WCLOUD) ( &(wl_out.lambda_r[iv]),
				    &(output_caoth->newsiz),
				    &nlyr,
				    input.filename[FN_PATH],
				    &nstring, 
				    output_caoth->microphys.lwc,
				    output_caoth->microphys.effr, 
				    output_caoth->optprop.dtau[iv], 
				    output_caoth->optprop.g1[iv],
				    output_caoth->optprop.ssa[iv], 
				    output_caoth->zd,
				    &input_caoth.layer );
	break;
      }
    } /* end if(input.ipa3d) */

    break;

  case PROP_ECHAM4:
  case PROP_EXPLICIT:
    if (input.ipa3d) {
      fprintf (stderr, "Error, wc_properties echam4 (and explicit) do not work with ipa_3d\n");
      fprintf (stderr, "       -  use mie or hu!\n");
      
      return -1;
      break;
    }
    break;

  case PROP_MIE:
  case PROP_FILE:
  case PROP_IC_MIE:
  case PROP_BAUM:
  case PROP_BAUM_V36:  
  case PROP_HEY:
  case PROP_YANG2013:

    /* recalculate the optical properties for the */
    /* new microphysics profile                   */
    status = ssprop2optprop (output_caoth, iv, input.quiet);
    if (status)
      return fct_err_out ( status, "ssprop2optprop", ERROR_POSITION );
    break;

  case PROP_KEY:
  case PROP_YANG:
    if (input.ipa3d) {
      if (input_caoth.properties==PROP_KEY) {
	if (!input.quiet)
	  fprintf (stderr, " ... using Key et al. [2002] ice cloud properties\n");

	newkey=0;
      }

      if (input_caoth.properties==PROP_YANG) {
	if (!input.quiet) {
	  fprintf (stderr, " ... using Key et al. [2002] ice cloud properties below 3.4micron;\n");
	  fprintf (stderr, " ... new data by Yang/Mayer above 3.4 micron.\n");
	}

	newkey=1;
      }

      for (iv=wl_out.nlambda_rte_lower; iv<=wl_out.nlambda_rte_upper; iv++) {
	status = ic_yang ( wl_out.lambda_r[iv],
			   output_caoth->nlev-1,
			   input_caoth.habit,
			   input.filename[FN_PATH],
			   output_caoth->microphys.lwc_layer,
			   output_caoth->microphys.effr_layer, 
			   output_caoth->optprop.dtau[iv], 
			   output_caoth->optprop.ff[iv],
			   output_caoth->optprop.g1[iv],
			   output_caoth->optprop.g2[iv], 
			   output_caoth->optprop.ssa[iv], 
			   output_caoth->zd,
			   input_caoth.layer,
			   newkey );
	if (status)
	  return fct_err_out ( status, "ic_yang", ERROR_POSITION );
      }
    } /*end if (input.ipa3d) */
    
    break;

  case PROP_FU:
    if (input.ipa3d) {
      nlyr = output_caoth->nlev-1;
	/* ulr: scaling of effr was done already in cloud3d.c */
	/* if (input.ic.fu2yang) // convert from Fu (1996/98) to Key et al. (2002) effective radius */
	/* ulr: set output for all height levels to zero */
	for (lev=0; lev<nlyr; lev++) {
	  output_caoth->optprop.dtau[iv][lev] = 0.0;
	  output_caoth->optprop.g1[iv][lev]   = 0.0;
	  output_caoth->optprop.ssa[iv][lev]  = 0.0;
	  output_caoth->optprop.f[iv][lev]    = 0.0;
	  if(output_caoth->microphys.lwc_layer[lev]==0)
	    output_caoth->optprop.ff[iv][lev]  = 0.0;
	}
      
      /* calculate caoth optical properties */
      switch (input.ck_scheme) {
      case CK_FU:
        
#if HAVE_FULIOU
        status = ckdfuice ( output_caoth->zd,
			    output_caoth->microphys.effr_layer, 
			    output_caoth->microphys.lwc_layer, 
			    output_caoth->nlev,
			    output_caoth->optprop.dtau, 
			    output_caoth->optprop.g1, 
			    output_caoth->optprop.ssa, 
			    input_caoth.unscaled );
	if (status)
	  return fct_err_out ( status, "ckdfuice", ERROR_POSITION );
#else
        fprintf (stderr, "Error, Fu and Liou not supported!\n");
        return -1;
#endif
	break;
	       
      default:
	/* ulrike: removed loop over wavelength since there is a lambda-loop already */

	if (wl_out.lambda_r[iv]<4000.0) {
	  status = ic_fu96 ( wl_out.lambda_r[iv],
			     output_caoth->nlev-1,
			     input.filename[FN_PATH],
			     output_caoth->microphys.lwc_layer,
			     output_caoth->microphys.effr_layer, 
			     output_caoth->optprop.dtau[iv],
			     output_caoth->optprop.g1[iv],
			     output_caoth->optprop.ssa[iv],
			     output_caoth->optprop.f[iv],
			     output_caoth->zd,
			     input_caoth.layer,
			     input_caoth.unscaled );
	  if (status)
	    return fct_err_out ( status, "ic_fu96", ERROR_POSITION );
	}
	else {
	  status = ic_fu98 ( wl_out.lambda_r[iv],
			     output_caoth->nlev-1,
			     input.filename[FN_PATH],
			     output_caoth->microphys.lwc_layer,
			     output_caoth->microphys.effr_layer, 
			     output_caoth->optprop.dtau[iv],
			     output_caoth->optprop.g1[iv],
			     output_caoth->optprop.ssa[iv], 
			     output_caoth->zd,
			     input_caoth.layer );
	  if (status)
	    return fct_err_out ( status, "ic_fu98", ERROR_POSITION );
	}

	for (lev=0; lev<nlyr; lev++)
	  if (output_caoth->optprop.dtau[iv][lev]>0.0)
	    output_caoth->optprop.ff[iv][lev]  = 1; /* ulrike: eigentlich sollte
							 ff==0 sein, aber im falle
							 von fu muss ff=1 sein*/
      } /*ulrike: switch (input.ck_scheme)*/
    } /* ulrike: end of if: input.ipa3d*/
    break; /* ulrike: end of case fu */
  default:
    fprintf (stderr, "Error, unknown property %d for %s. This is not\n",
	     input_caoth.properties, input_caoth.fullname);
    fprintf (stderr, "supposed to happen and indicates a coding error.\n");

    return -1;
    break;
  }
  
  return 0;
}


/***********************************************************************************/
/* Function: read_raytracing_file                                         @62_30i@ */
/* Description: read 5-column file with ice cloud description for rayrtracing;     */
/*              the five columns are                                               */
/*              1. name (ASCII label)                                              */
/*              2. fraction                                                        */
/*              3. oriented fraction                                               */
/*              4. angular distribution width                                      */
/*              5. degrees of freedom (n = n axes fixed)                           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Bernhard Mayer                                                          */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int read_raytracing_file (char *filename, crystal_prop_struct **raytracing_prop, int *n_raytracing_prop, int quiet)
{
  int i=0, status=0;
  int rows=0, min_columns=0, max_columns=0, max_length=0;
  char ***string=NULL, *dummy=NULL;
  
  status = ASCII_checkfile (filename, 
			    &rows,
			    &min_columns,
			    &max_columns,
			    &max_length);
    
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  if (max_columns<5) {
    fprintf (stderr, "Error, found less than four columns in %s\n", 
	     filename);
    return -1;
  }
  
    
  /* allocate memory */
  status = ASCII_calloc_string (&string,
				rows,
				max_columns,
				max_length);
    
  if (status!=0) {
    fprintf (stderr, "Error %d allocating memory\n", status);
    return status;
  }
    
  
  /* read file to string array */
  status = ASCII_readfile (filename, string);
  
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  
  /* allocate memory for Raytracing properties */
  *raytracing_prop = calloc(rows, sizeof(crystal_prop_struct));
  
  for (i=0; i<rows; i++) {
    (*raytracing_prop)[i].name = calloc (strlen(string[i][0])+1, sizeof(char));
    strcpy ((*raytracing_prop)[i].name, string[i][0]);

    (*raytracing_prop)[i].fraction = strtod (string[i][1], &dummy);
    (*raytracing_prop)[i].oriented_fraction = strtod (string[i][2], &dummy);
    (*raytracing_prop)[i].angdist_width = strtod (string[i][3], &dummy);
    (*raytracing_prop)[i].orientation_dof = strtol (string[i][4], NULL, 0);


    fprintf (stderr, "%s %f %f %f %d\n", 
	     (*raytracing_prop)[i].name,
	     (*raytracing_prop)[i].fraction,
	     (*raytracing_prop)[i].oriented_fraction,
	     (*raytracing_prop)[i].angdist_width,
	     (*raytracing_prop)[i].orientation_dof);
  }

  /* free memory */
  (void) ASCII_free_string (string, rows, max_columns);
  
  if (!quiet)
    fprintf (stderr, " ... read %d data points from %s\n", 
	     rows, filename);
  
  *n_raytracing_prop=rows;
  
  return 0;
}
