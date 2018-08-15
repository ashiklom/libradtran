/*--------------------------------------------------------------------
 * $Id: spline.c 2623 2011-12-23 10:52:38Z robert.buras $
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

/*  @14c@                

@code{spline} interpolates discrete data points using natural cubic
splines or linear interpolation. The x-values in the first column must 
be in ascending order. 

The different options to @code{spline} are displayed when executing:

@example
spline -h
@end example

    @c14@ */


#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h> 
#include <math.h> 
#include <getopt.h>
#include "ascii.h"
#include "numeric.h"

#define PROGRAM "SPLINE"
#define VERSION "3.99d"



static void print_usage (char *name)
{
  fprintf (stderr, "Usage: %s [OPTION] FILE\n", name);
  fprintf (stderr, "%s interpolates between given daxta points using\n", PROGRAM);
  fprintf (stderr,"natural cubic splines (default) or linear interpolation (option -l).\n\n");
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "  -b FLOAT      first value to be interpolated\n");
  fprintf (stderr, "  -s FLOAT      step width for equidistant interpolation\n");
  fprintf (stderr, "  -w FLOAT      weighting factor for approximating spline\n");
  fprintf (stderr, "  -x FILE       file containing x-values to be interpolated\n");
  fprintf (stderr, "  -l            linear interpolation\n");
  fprintf (stderr, "  -q            be quiet!\n\n");
  fprintf (stderr, "Input-file FILE: ASCII-file, blank delimited\n");
  fprintf (stderr, "  x1  y1  z1  ...\n");
  fprintf (stderr, "  x2  y2  z2  ...\n\n");
  fprintf (stderr, "Example: %s test.prn -b 2.5 -s 0.5  calculates:\n", name);
  fprintf (stderr, "  2.5  y(x=2.5)\n");
  fprintf (stderr, "  3.0  y(x=3.0) ...\n");
}




int main (int argc, char **argv)
{
  char filename[FILENAME_MAX+200] = "", xfilename[FILENAME_MAX+200]="";
  int rows=0;
  int newnumber=0;
  int status=0;
  int linear=0;
  int i=0, k=0, m=0;
  double *new_x=NULL, *new_y=NULL;
  double weighting=-1;
  double start=0, step=0;
  double *x=NULL, *y=NULL, *w=NULL, *x_user=NULL;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  
  int max_columns=0, min_columns=0;
  double **data=NULL;
  double **result=NULL;

  double ynew=0;
  int c=0;
  char *dummy=NULL;
  int user_start=0, user_step=0, user_x=0;

  int quiet=0;

  char allocated=0;

/*  fprintf (stderr, "SPLINE %s, 1995 by Bernhard Mayer\n", VERSION);
  fprintf (stderr, "! beta ! beta ! beta ! beta ! beta ! beta !\n\n"); */
  
  while ((c=getopt (argc, argv, "b:s:w:x:lcq")) != EOF)  {
    switch(c)  {
    case 'w':   
      weighting = strtod (optarg, &dummy);
      break;
    case 'b': 
      user_start = 1;
      start = strtod (optarg, &dummy);
      break;
    case 's':
      user_step = 1;
      step = strtod (optarg, &dummy);
      break;
    case 'x':
      user_x = 1;
      strcpy (xfilename, optarg);
      break;
    case 'l':
      linear = 1;
      break;
    case 'q':
      quiet = 1;
      break;
    case '?':
      break;
      return (-1);
    default:
      print_usage (argv[0]);
      return (-1);
    }
  }
      

  /* check number of remaining command line arguments */
  if (argc - optind != 1)  {
    print_usage (argv[0]);
    return -1;
  }
  

  strcpy (filename, argv[optind]);
	  

  /* read input file */
  if (!quiet)
    fprintf (stderr, " ... reading data from %s\n", filename);

  /* read file */
  status = ASCII_file2double (filename, &rows, &max_columns, &min_columns, &data);
  if (status!=0) {
    fprintf (stderr, "Error %d opening file %s\n", status, filename);
    return status;
  }

  if (rows<=0) {
    fprintf (stderr, "Error, empty file. Exiting\n");
    return -1;
  }

  if (!quiet)
    fprintf (stderr, " ... counted %d data points\n", rows);

  if (max_columns!=min_columns) {
    fprintf (stderr, " !! ATTENTION !! Inconsistent number of columns\n");
    fprintf (stderr, "     min = %d, max =%d\n", min_columns, max_columns);
  }
  
  if (min_columns<2) {
    fprintf (stderr, "Error, too few columns in %s\n", filename);
    return 0;
  }

  /* create array of weighting factors */
  w = (double *) calloc (rows, sizeof(double));
  
  for (i=0; i<rows; i++)
    w[i] = weighting;
    

  if (user_x)  {
    if (!quiet)
      fprintf (stderr, " ... reading interpolation grid from %s\n", xfilename);
    
    /* read file with user x values */
    status = read_1c_file (xfilename, &x_user, &newnumber);
    if (status!=0)  {
      fprintf (stderr, "Error reading file %s\n", xfilename);
      return (-1);
    }

    /* allocate memory for result array */
    result = calloc (newnumber, sizeof(double));
    for (m=0; m<newnumber; m++)
      result[m] = calloc (min_columns, sizeof(double));
  }

  for (k=1; k<min_columns; k++) {
    x = ASCII_column (data, rows, 0);
    y = ASCII_column (data, rows, k);
 
  
    /* do interpolation */
    if (user_x)  {
      if (!linear)  {   /* spline interpolation */
	if (weighting<0)  {    /* interpolating spline */
	  
	  /* calculate interpolating spline coefficients */
	  if (!quiet)
	    fprintf (stderr, " ... spline interpolation\n");
	  
	  status = spline_coeffc (x, y, rows, &a0, &a1, &a2, &a3);
	  if (status!=0)  {
	    fprintf (stderr, "sorry cannot do spline interpolation\n");
	    fprintf (stderr, "spline_coeffc() returned status %d\n", status);
	  }
	}
	else  {                /* approximating spline */
	  
	  if (!quiet)
	    fprintf (stderr, " ... spline approximation\n");
	  
	  /* calculate approximating spline coefficients */
	  status = appspl_coeffc (x, y, w, rows, &a0, &a1, &a2, &a3);
	  if (status!=0)  {
	    fprintf (stderr, "sorry cannot do spline approximation\n");
	    fprintf (stderr, "appspl_coeffc() returned status %d\n", status);
	  }
	}	
      }
      else  {   /* linear interpolation */
	
	if (!quiet)
	  fprintf (stderr, " ... linear interpolation\n");
	
	status = linear_coeffc (x, y, rows, &a0, &a1, &a2, &a3);
	if (status!=0)  {
	  fprintf (stderr, "sorry cannot do linear interpolation\n");
	  fprintf (stderr, "linear_coeffc() returned status %d\n", status);
	}
      }
      
      for (i=0; i<newnumber; i++)  {
	status = calc_splined_value (x_user[i], &ynew, x, rows, a0, a1, a2, a3);
	
	/* copy data to result array */
	if (status==0) {
	  result[i][0] = x_user[i];
	  result[i][k] = ynew;
	}
	else {
	  result[i][0] = x_user[i];
	  result[i][k] = NAN;
	}
      } 
    }
    else  {   /* equidistant steps */
      
      /* set first value to be interpolated to default */
      if (!user_start)
	start=ceil(x[0]);
      
      /* set step width to default 1.0 */
      if (!user_step)
	step=1.0;
      
      /* do spline interpolation */
      if (!linear)  {   /* spline interpolation */
	if (weighting<0)  {     /* interpolating spline */
	  if (!quiet)
	    fprintf (stderr, " ... spline interpolation\n");
	  if ((status = spline (x, y, rows, start, step, 
				&newnumber, &new_x, &new_y)) != 0)  {
	    fprintf (stderr, "sorry, cannot do interpolation!\n");
	    fprintf (stderr, "spline() returned status %d\n", status);
	    return (status);
	  }
	}
	else  {                 /* approximating spline */
	  if (!quiet)
	    fprintf (stderr, " ... spline approximation\n");

	  if ((status = appspl(x, y, w, rows, start, step, 
			       &newnumber, &new_x, &new_y)) != 0)  {
	    fprintf (stderr, "sorry, cannot do interpolation!\n");
	    fprintf (stderr, "appspl() returned status %d\n", status);
	    return (status);
	  }
	}
      }
      else  {   /* linear interpolation */
	if (!quiet)
	  fprintf (stderr, " ... linear interpolation\n");

	if ((status = linear_eqd (x, y, rows, start, step, 
				  &newnumber, &new_x, &new_y)) != 0)  {
	  fprintf (stderr, "sorry, cannot do interpolation!\n");
	  fprintf (stderr, "linear_eqd() returned status %d\n", status);
	  return (status);
	}
      }
      
      
      /* allocate memory for result array */
      if (!allocated) {
	allocated=1;

	/* allocate memory for result array */
	result = calloc (newnumber, sizeof(double));

	for (m=0; m<newnumber; m++)
	  result[m] = calloc (min_columns, sizeof(double));
      }


      /* copying data to result array */
      for (i=0; i<newnumber; i++) {
	result[i][0] = new_x[i];
	result[i][k] = new_y[i];
      }
      
      /* free memory of interpolated arrays */
      free(new_x);
      free(new_y);
    }

    free(x);
    free(y);
  }

  /* write new values to stdout */
  if (!quiet)
    fprintf (stderr, " ... writing % d data pairs to stdout\n", newnumber);
  
  for (i=0; i<newnumber; i++) {
    fprintf (stdout, "%20.11g ", result[i][0]);
    for (m=1; m<min_columns; m++)
      fprintf (stdout, " %20.11g", result[i][m]);
    fprintf (stdout, "\n");
  }

  ASCII_free_double (data, rows);

  for (m=0; m<newnumber; m++)
    free(result[m]);
  free(result);
  
  
  /* free memory */
  free(w);

  if (user_x)
    free (x_user);
  
  return 0;
}
