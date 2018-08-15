/*--------------------------------------------------------------------
 * $Id: vac2air.c 2623 2011-12-23 10:52:38Z robert.buras $
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
#ifndef _AIX
#include <getopt.h>
#else
#include <unistd.h>
#endif
#include <math.h>
#include "ascii.h"
#include "numeric.h"
#include "wavelength.h"

#define PROGRAM "VAC2AIR"
#define VERSION "0.8b"
#define DATE    "May 28, 1998"



/*************************************************************/
/* print usage information                                   */
/*************************************************************/

static void print_usage (char *filename)
{
  fprintf (stderr, "%s %s - Vaccuum to air and vice versa conversion \n\n", PROGRAM, VERSION);
  fprintf (stderr, "written by  Bernhard Mayer,\n");
  fprintf (stderr, "            National Center for Atmospheric Research (NCAR)\n");
  fprintf (stderr, "            P.O. Box 3000, Boulder, CO 80307, USA\n");
  fprintf (stderr, "Version %s finished %s\n\n", VERSION, DATE);
  fprintf (stderr, "Be aware that this program is a beta version! Please report any\n");
  fprintf (stderr, "kind of error to the author, including the error message and the\n");
  fprintf (stderr, "database being processed when the error occured. Thanks!\n\n");

  fprintf (stderr, "USAGE: %s [options] <input file>\n\n",
	   filename);
  fprintf (stderr, "%s will shift wavelength dependent quantities from vacuum to air\n", PROGRAM);
  fprintf (stderr, "or vice versa. The first column of the input file must be the wavelength\n");
  fprintf (stderr, "in nm. The data will be shifted and then interpolated back to the\n");
  fprintf (stderr, "original wavelength scale. Note that during that process, at least one\n");
  fprintf (stderr, "data point will be lost. The output will be written to stdout.\n");
  fprintf (stderr, "A formula from the 1997/98 'CRC handbook of Chemistry and Physics'\n");
  fprintf (stderr, "is used, according to which the wavelength in air corresponds to\n");
  fprintf (stderr, "a temperature of +15 deg C and a pressure of 1013.25 mbar.\n\n");
  fprintf (stderr, " command line options:\n");
  fprintf (stderr, "  -h  show this page\n");
  fprintf (stderr, "  -r  convert from air to vacuum (default: vacuum to air)\n");
  fprintf (stderr, "  -l  use linear interpolation for the data (default: spline)\n");
}
  


/*************************************************************/
/* parse command line options                                */
/*************************************************************/

static int get_options (int argc, char **argv,
			char *programname,
			char *infilename,
			char *reverse,
			char *linear)
{
  int c=0;

  /* default settings */
  *reverse = 0;
  *linear  = 0;


  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  /* get command line arguments */
  /* The following did not work on my linux box, GETOPT_ANY must
     be supplied from somewhere, but where, configure?. 
     A. Kylling 13.06.98 */
  /*  #ifndef _AIX
  optmode = GETOPT_ANY;
  #endif */

  while ((c=getopt (argc, argv, "hrl")) != EOF)  {
    switch(c)  {
    case 'r':
      *reverse = 1;
      break;
    case 'l':
      *linear = 1;
      break;
    case 'h':  /* help */
      print_usage (programname);
      return (-1);
      break;
    case '?':
      print_usage (programname);
      return (-1);
      break;
    default:
      print_usage (programname);
      return (-1);
    }
  }
      

  /* check number of remaining command line arguments */
  if (argc - optind != 1)  {
    print_usage (programname);
    return -1;
  }
  

  /* save command line arguments */
  strncpy (infilename,    argv[optind+0], FILENAME_MAX);

  return 0;  /* if o.k. */
}





int main(int argc, char **argv) 
{
  int i=0, j=0, status=0;

  char programname   [FILENAME_MAX]="";
  char infilename    [FILENAME_MAX]="";

  double *lambda=NULL, *irradiance=NULL, *irradiance_shifted=NULL;

  int rows=0, min_columns=0, max_columns=0;

  double **value=NULL, **shifted_value=NULL;

  char reverse=0, linear=0;



  /* get command line options */
  status=get_options (argc, argv, programname, infilename, &reverse, &linear);

  if (status!=0)
    return status;
  

  /* read input file */
  status = ASCII_file2double (infilename,   
			      &rows, &max_columns, &min_columns, 
			      &value);

  if (status!=0) {
    fprintf (stderr, "error %d reading input file %s\n", status, infilename);
    return status;
  }

  if (min_columns<2) {
    fprintf (stderr, "error, input file must have at least two columns\n");
    return -1;
  }


  if (min_columns != max_columns) {
    fprintf (stderr, " WARNING! %s is not a rectangular matrix; only the first\n", infilename);
    fprintf (stderr, " %d columns will be used for the analysis\n", min_columns);
  }

  /* extract wavelength column to array lambda[] */
  lambda = ASCII_column(value, rows, 0);


  /* allocate memory for the shifted array */
  status = ASCII_calloc_double(&shifted_value, rows, min_columns);
  
  if (status!=0) { 
    fprintf (stderr, " error allocating memory for shifted array\n");
    return status;
  }

  /* copy wavelength column to shifted array */
  for (i=0; i<rows; i++)
    shifted_value[i][0] = lambda[i];


  for (i=1; i<min_columns; i++) {

    /* copy i-th column to irradiance[] */
    irradiance = ASCII_column(value, rows, i);
    

    /* shift i-th column */
    status = vac2air (lambda, irradiance, rows,
		      reverse, linear,
		      &irradiance_shifted);
    
    if (status!=0) {
      fprintf(stderr, "error shifting %dth column\n", i+1);
      return status;
    }

				    
    /* copy shifted data to result array */
    for (j=0; j<rows; j++)
      shifted_value[j][i] = irradiance_shifted[j];


    /* free memory */
    free(irradiance);
    free(irradiance_shifted);
  }


  /* write data to stdout */
  for (i=0; i<rows; i++) {
    for (j=0; j<min_columns; j++)
      fprintf (stdout, "%g ", shifted_value[i][j]);
    
    fprintf (stdout, "\n");
  }

  /* free memory */
  ASCII_free_double (value, rows);
  ASCII_free_double (shifted_value, rows);
  free (lambda);

  return 0;   /* if o.k. */
}
