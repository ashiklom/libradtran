/*--------------------------------------------------------------------
 * $Id: conv.c 2623 2011-12-23 10:52:38Z robert.buras $
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

/*  @15c@                

@code{conv} convolutes a spectrum with a given filter function.

The different options to @code{conv} are displayed when executing:

@example
conv -h
@end example

    @c15@ */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ascii.h"
#include "numeric.h"


#define PROGRAM "CONV"
#define VERSION "1.99b"


/*************************************************************/
/* print usage information                                   */
/*************************************************************/

static void print_usage (char *filename)
{
  fprintf (stderr, "%s %s - convolute ASCII data\n\n", PROGRAM, VERSION);
  fprintf (stderr, "written by  Bernhard Mayer,\n");
  fprintf (stderr, "            DLR, eMail bernhard.mayer@dlr.de\n");
  fprintf (stderr, "Version %s finished March 18, 1997\n\n", VERSION);
  fprintf (stderr, "Be aware that this program is an beta version! Please report any\n");
  fprintf (stderr, "kind of error to the author, including the error message and the\n");
  fprintf (stderr, "database being processed when the error occured. Thanks!\n\n");
  fprintf (stderr, "USAGE: %s <input file> <convolution file>\n\n", filename);
  fprintf (stderr, "%s will convolute ASCII files with the given convolution\n", PROGRAM);
  fprintf (stderr, "function. For this purpose, the ASCII file is interpolated LINEARELY\n");
  fprintf (stderr, "to the resolution of the convolution file, which MUST be given\n");
  fprintf (stderr, "in equidistant wavelength steps. Output is written to stdout.\n");

}

int conv(char *specfilename, char *convfilename, int *rows, int* columns, 
	 double **x, double ***ys) {
  int max_columns=0, min_columns=0;
  int status=0;
  int k=0, n_conv=0;
  double **data=NULL;
  double *y_spec=NULL, *x_conv=NULL, *y_conv=NULL;

  /*  fprintf (stderr, " ... reading file %s ...\n", specfilename);*/ 

  /* read file */ 
  status = ASCII_file2double (specfilename,  
			      rows, &max_columns, &min_columns,  
			      &data); 
  if (status!=0) { 
    fprintf (stderr, "Error %d opening file %s\n", status, specfilename); 
    return status; 
  } 
  
  /*  fprintf (stderr, " ... counted %d data points\n", rows);*/
  
  if (max_columns!=min_columns) { 
    fprintf (stderr, " !! ATTENTION !! Inconsistent number of columns\n"); 
    fprintf (stderr, "     min = %d, max =%d\n", min_columns, max_columns); 
  } 
  
  if (min_columns<2) { 
    fprintf (stderr, "Error, too few columns in %s\n", specfilename); 
    return 0; 
  } 
  
  /*  fprintf (stderr, " ... reading file %s ...\n", convfilename);*/ 
  status = read_2c_file (convfilename, &x_conv, &y_conv, &n_conv); 
  
  if (status!=0)  { 
    fprintf (stderr, "Error %d reading file %s\n", status, convfilename); 
    return status; 
  } 
  
 /*  fprintf (stderr, " ... convoluting ...\n");*/ 

  /* allocate memory for result array */ 
   *ys = calloc (min_columns-1, sizeof(double *)); 

   /* do convolution */ 
  *x = ASCII_column (data, *rows, 0); 
  for (k=1; k<min_columns; k++) { 
    y_spec = ASCII_column (data, *rows, k); 
    
    if ((status = int_convolute (*x, y_spec, *rows, 
 				 x_conv, y_conv, n_conv, 
 				 &((*ys)[k-1]), 0)) < 0)  { 
    
      switch (status)  { 
      case SPEC_NOT_EQUIDISTANT: 
 	fprintf (stderr, "Error, wavelengths in %s not equidistant!\n", specfilename); 
 	break; 
      case CONV_NOT_EQUIDISTANT: 
 	fprintf (stderr, "Error, wavelengths in %s not equidistant!\n", convfilename); 
 	break; 
      case CONV_NOT_CENTERED: 
 	fprintf (stderr, "Error, no center defined in %s\n", convfilename); 
 	break; 
      case SPEC_CONV_DIFFERENT: 
 	fprintf (stderr, "Error, different wavelength steps in %s and %s!\n",  
 		 specfilename, convfilename); 
 	break; 
      default: 
 	fprintf (stderr, "Error, status %d returned by conv()\n", status); 
 	break; 
      } 
    
      return status; 
    } 

  } 

  free(y_spec); 
  free (x_conv);
  free (y_conv);
    
  *columns = min_columns;
  return 0;
}

int main(int argc, char **argv)
{  
  char programname  [FILENAME_MAX+200] = "";
  char specfilename [FILENAME_MAX+200] = "";
  char convfilename [FILENAME_MAX+200] = "";

  double *x=NULL;
  double **ys=NULL;

  int i=0, k=0, status=0;
  int rows=0, columns=0;

  strcpy (programname, argv[0]);

  /* check for command line arguments */
  if (argc != 3)  {
    print_usage (argv[0]);
    return (-1);
  }

  /* get filenames from command line */
  strcpy (specfilename, argv[1]);
  strcpy (convfilename, argv[2]);

  status = conv(specfilename, convfilename, &rows, &columns, &x, &ys);

  if (status!=0) { 
    fprintf (stderr, "Error, status %d returned by conv\n", status); 
    return status; 
  } 
    
  for (i=0; i<rows; i++)  {
    fprintf (stdout, "%14.8e  ", x[i]);
    for (k=0; k<columns-1; k++) 
      fprintf (stdout, "%14.8e  ", ys[k][i]);
    fprintf (stdout, "\n");
  }

  /* free memory  */
  free(x); 
  for (k=0; k<columns-1; k++)
    free (ys[k]);
  
  free (ys);

  return 0;
}

