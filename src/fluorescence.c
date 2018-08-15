/*--------------------------------------------------------------------
 * $Id: fluorescence.c 2986 2013-10-03 13:40:35Z svn-kylling $
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

/************************************/
/* prototypes of internal functions */
/************************************/

static int read_fluorescence (char *filename, flu_out_struct *flu_out); 

/**************************************************************/
/* Setup surface fluorescence.                                      */
/**************************************************************/

int setup_fluorescence (input_struct input, output_struct *output)
{
  int status=0, iv=0, linear=0, i=0;
  double fluorescence = NOT_DEFINED_DOUBLE;
  double *tmp=NULL;

  char *filename=NULL; 

  char function_name[]="setup_fluorescence";
  char file_name[]="fluorescence.c";
  
  /* copy fluorescence input to output structure*/
  output->flu.source        = input.flu.source;
  output->flu.fluorescence  = input.flu.fluorescence;
  
  /* allocate requiered arrays */
  output->flu.fluorescence_r = (double *) calloc (output->wl.nlambda_r, sizeof(double));
  
  switch(input.flu.source) {
  case FLUORESCENCE_CONSTANT:
    
    fluorescence = input.flu.fluorescence;
    
    /* security check of fluorescence range */
    if ( fluorescence < 0.0  ) {
      fprintf (stderr, "\nError, fluorescence value = %f is out of range [ 0, 1 ] !!! \n", fluorescence );
      fprintf (stderr,   "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    }

    if (input.verbose)
      fprintf (stderr, "     wavelength indenpendent fluorescence of %f \n", fluorescence);
    
    for (iv=0;iv<output->wl.nlambda_r;iv++)
      output->flu.fluorescence_r[iv] = fluorescence;
    break;

  case FLUORESCENCE_FROM_FLUORESCENCE_FILE:

    filename = calloc (FILENAME_MAX+1, sizeof(char));
    if (input.flu.source == FLUORESCENCE_FROM_FLUORESCENCE_FILE)
      /* copy fluorescence file name */
      strcpy (filename, input.filename[FN_FLUORESCENCE]);
    else {
      fprintf (stderr, "Error, unknown fluorescence source %d in %s (%s) \n", input.flu.source, function_name, file_name);
      return -1;
    }

    status = read_fluorescence (filename, &output->flu);
    if (status != 0) {
      fprintf (stderr, "Error reading fluorescence file '%s' \n", input.filename[FN_FLUORESCENCE]);
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return status;
    }
    
    /**** Interpolate fluorescence to internal wavelength grid ****/
    
    linear = 1;  /**** Use linear interpolation for the fluorescence, otherwise
		       negative fluorescences may occurr AKy, ****/
    
    tmp = (double *) calloc (output->wl.nlambda_r, sizeof(double));
    for (i=0;i<output->wl.nlambda_r;i++) 
      tmp[i] = (double) output->wl.lambda_r[i];

    status = arb_wvn2_double (output->flu.file_nlambda_fluorescence, output->flu.file_lambda_fluorescence,
		       output->flu.file_fluorescence,
		       output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
		       tmp, output->flu.fluorescence_r,
		       linear);

    free(tmp);
    if (status != 0) {
      fprintf (stderr, "Error %d returned by arb_wvn2()\n", status);
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return status;
    }
    free(filename);
    break;
    
  default:
    return 0;
  }

   /**** Check interpolated fluorescence values ****/

   for (iv=0; iv<output->wl.nlambda_r; iv++) {
     /*      if (input.verbose) { */
     /*        fprintf (stderr, "lambda   fluorescence"); */
     /*        fprintf (stderr, "%10.3f %9.5f\n",output->wl.lambda_r[iv], output->alb.fluorescence_r[iv]); */
     /*      } */
     if ( output->flu.fluorescence_r[iv] < 0.0)  {
       fprintf (stderr, "Fluorescence interpolated from %s out of range:\n", input.filename[FN_FLUORESCENCE]);
       fprintf (stderr, "wavelength: %f, fluorescence %f\n",
		output->wl.lambda_r[iv], output->flu.fluorescence_r[iv]);
       fprintf (stderr, "Modify fluorescence file %s\n", input.filename[FN_FLUORESCENCE]);
       return -1;
     }
   }
   return 0;
}





/**************************************************************/
/* Read 2-column fluorescence file.                                 */
/**************************************************************/

static int read_fluorescence (char *filename, flu_out_struct *flu_out)
{
  int status=0;
  
  status = read_2c_file (filename,
			 &(flu_out->file_lambda_fluorescence),
			 &(flu_out->file_fluorescence),
			 &(flu_out->file_nlambda_fluorescence));

  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  return 0;
}
