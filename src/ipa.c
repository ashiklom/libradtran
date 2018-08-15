/*--------------------------------------------------------------------
 * $Id: ipa.c 3276 2017-07-04 14:16:36Z Claudia.Emde $
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

#include "uvspec.h"
#include "ascii.h"
#include "errors.h"

#define EPSILON 1E-6

int setup_ipa ( input_struct   input,
		output_struct *output )
{
  int isp=0, isipa=0;

  /* Reading of ipa input files etc. is done in setup_wcloud */

  for (isp=0; isp<input.n_caoth; isp++)
    if (input.caoth[isp].ipa)
      isipa=1;

  if ( isipa == 0 ) {
    /* no independent pixel */
    output->nipa = 1;
    output->ipaweight = calloc (1, sizeof(double));
    output->ipaweight[0] = 1;
  }

  output->niipa=1; /* ulrike */
  output->njipa=1; /* ulrike */

  return 0;   /* if o.k. */
}


/*********************************/
/* compare two altitude grids    */
/*********************************/

int compare_levels_ipa_df ( float *z1,
			    int    n1,
			    float *z2,
			    int    n2 )
{
  int kc=0;

  if (n1!=n2) {
    fprintf (stderr, "Error, number of layers differing (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    fprintf (stderr, "%d in 3D, %d in 1D\n", n1, n2);
    return -1;
  }

  for (kc=0; kc<=n1; kc++)
    if (fabs (z1[kc]/1000.0 - z2[n1-kc]) > EPSILON * z1[kc]) {
      fprintf (stderr, "Error, difference at level %d, altitude %g vs. %g (line %d, function %s in %s)\n",	
	       kc, z1[kc], z2[n1-kc], __LINE__, __func__, __FILE__);
      return -1;
    }

  return 0;  /* if o.k. */
}


/***********************************************************************************/
/* Function: read_caoth_ipa                                               @62_30i@ */
/* Description:                                                                    */
/*  read ipa file for caoth                                                        */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras (RPB)                                                      */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int read_caoth_ipa ( input_struct       input,
		     caoth_inp_struct   input_caoth,
		     float              altitude,
		     wl_out_struct      wl_out,
		     atm_out_struct     atm_out,
		     /* output */
		     int               *nipa,
		     float            **ipaweight,
		     caoth_out_struct  *output_caoth,
		     caoth_out_struct **output_caoth_ipa,
		     cf_out_struct     *cf_out )
{
  int ipa=0, min_columns=0, max_columns=0, max_length=0, status=0, maxlev=0, maxic=0;
  char ***string=NULL;
  char *dummy=NULL;

  float tmpweight=0;
  int tmpnipa=0;
  int first=0;

  if (!input.quiet)
    fprintf (stderr, " ... reading ipa data from %s for %s\n",
	     input_caoth.filename, input_caoth.fullname);
      
  status = ASCII_checkfile ( input_caoth.filename, &tmpnipa,
			     &min_columns, &max_columns, &max_length );
  if (status)
    return fct_err_out ( status, "ASCII_checkfile", ERROR_POSITION );

  if (min_columns<2) {
    fprintf (stderr, "Error, need at least 2 columns in independent pixel file %s\n",
	     input_caoth.filename);
    return -1;
  }

  if ( *nipa == 0 ) {
    /* first call to read_caoth_ipa, initialize */
    *nipa = tmpnipa;
    *ipaweight = calloc (*nipa, sizeof(float));
    if (*ipaweight==NULL)
      return mem_err_out ( "ipaweight", ERROR_POSITION );
    first=1;
  }
  else {
    /* check if caoth have same number of ipa columns */
    if ( tmpnipa != *nipa ) {
      fprintf (stderr, "Error. %s has different number of ipa columns\n",
	       input_caoth.fullname);
      fprintf (stderr, "than former, %d and %d.\n", tmpnipa, *nipa);
      return -1;
    }
  }

  /* read ASCII file */
  status = ASCII_calloc_string ( &string, *nipa, max_columns, FILENAME_MAX+1 );
  if (status)
    return mem_err_out ( "string", ERROR_POSITION );
  
  status = ASCII_readfile ( input_caoth.filename, string );
  if (status)
    return fct_err_out ( status, "ASCII_readfile", ERROR_POSITION );

  /* allocate memory for IPA arrays */
  *output_caoth_ipa = calloc (*nipa, sizeof(caoth_out_struct));
  if (*output_caoth_ipa==NULL)
    return mem_err_out ( "output_caoth_ipa", ERROR_POSITION );

  for (ipa=0; ipa<*nipa; ipa++) {
    tmpweight = (float) strtod (string[ipa][1], &dummy);
    if (first==1)
      (*ipaweight)[ipa] = tmpweight;
    else {
      if ((*ipaweight)[ipa] != tmpweight) {
	fprintf (stderr, "Error, %s has not same weight as former:\n",
		 input_caoth.fullname);
	fprintf (stderr, "       %f and %f in row %d\n",
		 (*ipaweight)[ipa], tmpweight, ipa);
	return -1;
      }
    }

    if (!input.quiet)
      fprintf (stderr, " ... reading %s\n", string[ipa][0]);

    if (strlen(input_caoth.properties_filename) > 0) {
      status = read_caoth_file ( input_caoth.source,
				 input_caoth.name,
				 input_caoth.fullname,
				 string[ipa][0],
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
				 0,
				 wl_out.nlambda_r,
				 input_caoth.layer,
				 input.verbose,
				 input.quiet,
				 &((*output_caoth_ipa)[ipa]),
				 cf_out );
      if (status)
	return fct_err_out ( status, "read_caoth_file", ERROR_POSITION );
    }
    else {
      status = read_and_convert_caoth_file (input,
					    string[ipa][0],
					    wl_out,
					    altitude,
					    atm_out,
					    input_caoth,
					    &((*output_caoth_ipa)[ipa]),
					    cf_out,
					    &(atm_out.nmom) );
      if (status)
	return fct_err_out ( status, "read_and_convert_caoth_file", ERROR_POSITION );
    }
  }

  /* free string */
  (void) ASCII_free_string (string, *nipa, max_columns);
  
  /* allocate memory for general wc structure where the single scattering properties */
  /* are stored, and where the independent columns are copied later                  */
  maxlev=0;
  maxic=0;
  for (ipa=0; ipa<*nipa; ipa++) {
    if (maxlev < (*output_caoth_ipa)[ipa].nlev) {
      maxlev = (*output_caoth_ipa)[ipa].nlev;
      maxic = ipa;
    }
  }
  
  /* ????? */
  /* CE dirty fix here, for some reason nphamat is not set here */
  if( output_caoth->optprop.nphamat==0)
    output_caoth->optprop.nphamat=6;
  /* ????? */
  
  /* allocate memory */
  status = calloc_caoth_out ( output_caoth,
			      input_caoth.name,
			      input_caoth.fullname,
			      wl_out.nlambda_r,
			      maxlev, 
			      output_caoth->optprop.nphamat,
			      1);
  if (status)
    return mem_err_out ( "output_caoth struct ", ERROR_POSITION );

  /* Read optical properties file if available */
  if (strlen(input_caoth.properties_filename) > 0) {
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
			       &(atm_out.nmom),
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

    if (input.verbose)
      fprintf (stderr, " ... read optical properties from %s for %s\n",
	       input_caoth.properties_filename, input_caoth.fullname);
  }

  /* copy independent pixel microphysical properties of the pixel with
     the most levels to output_caoth */
  status = cp_caoth_out (output_caoth, (*output_caoth_ipa)[maxic], 0, 0, input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d copying ipa %s to %s ",
	     status, input_caoth.fullname, input_caoth.fullname );
    fprintf (stderr, "(line %d, function %s in %s)\n", 
	     __LINE__, __func__, __FILE__ );
    return status;
  }

  if (input.verbose)
    fprintf (stderr, " ... read %d independent columns from %s for %s\n",
	     *nipa, input_caoth.filename, input_caoth.fullname);

  return 0;
}
