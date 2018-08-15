/*--------------------------------------------------------------------
 * $Id: integrate.c 2623 2011-12-23 10:52:38Z robert.buras $
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
  
/*  @13c@                

@code{integrate} calculates the integral between limits x_min and x_max
by interpolating the data points (x[i], y[i]) with natural cubic
splines or linear interpolation. x_min and x_max are the minimum
and maximum values of the first column in the input_file. The 
x-values in the first column must be in ascending order.

The different options to @code{integrate} are displayed when executing:

@example
integrate -h
@end example

    @c13@ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ascii.h"
#include "numeric.h"
#include "getopt.h"

#define PROGRAM "integrate"
#define VERSION "1.0"



static void print_usage (char *name)
{
  fprintf (stderr, "\nUsage:   %s input_file\n\n", name);
  fprintf (stderr, "%s calculates the integral between limits x_min and x_max\n", PROGRAM);
  fprintf (stderr, "by interpolating the data points (x[i], y[i]) with natural cubic\n");
  fprintf (stderr, "splines (default) or linear interpolation. x_min and x_max are the minimum\n");
  fprintf (stderr, "and maximum values of the first column in the input_file. The \n");
  fprintf (stderr, "x-values in the first column must be in ascending order.\n\n");
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "  -l            linear interpolation between data points\n");
  fprintf (stderr, "  -p            turn off printing of messages\n");
  fprintf (stderr, "  -a <a>        left integration limit (optional)\n");
  fprintf (stderr, "  -b <b>        right integration limit (optional)\n");
  fprintf (stderr, "  -q            be quiet!\n\n");
  fprintf (stderr, "Input-file: ASCII-file, blank delimited\n");
  fprintf (stderr, "              x1  y1  z1  ...\n");
  fprintf (stderr, "              x2  y2  z2  ...\n\n");
  fprintf (stderr, "For example  %s -p ../test/integrate_test.in  will create:\n", PROGRAM);
  fprintf (stderr, "               5.000000e-01    1.000000e+00    5.000000e-01\n");
}




int main (int argc, char **argv)
{
  char filename[FILENAME_MAX+200] = "";
  int print=1, adef=0, bdef=0;

  char *dummy=NULL;
  int k=0, rows=0, linear=0, status=0;
  double *x=NULL, *y=NULL;
  
  int max_columns=0, min_columns=0;
  double **data=NULL;

  double a=0, b=0, ans=0;
  int c=0;
  
  while ((c=getopt (argc, argv, "plc?ha:b:q")) != EOF)  {
    switch(c)  {
    case 'l':
      linear = 1;
      break;
    case 'h':
      print_usage (argv[0]);
      return (-1);
      break;
    case 'a':
      a = strtod (optarg, &dummy);
      adef=1;
      break;
    case 'b':
      b = strtod (optarg, &dummy);
      bdef=1;
      break;
    case 'q':  /* don't print messages */  
    case 'p':  /* obsolete Syntax, undocumented but probably already in use by somebody */
      print = 0;
      break;
    case '?':
      print_usage (argv[0]);
      return (-1);
      break;
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
  if (print) {
    if (linear)
      fprintf (stderr, "Using linear interpolation\n");
    else
      fprintf (stderr, "Using cubic spline interpolation\n");
    
    fprintf (stderr, "Results reported for all columns\n");
    fprintf (stderr, " ... reading data from file \"%s\" ...\n", filename);
  }

  /* read file */
  status = ASCII_file2double (filename, &rows, &max_columns, &min_columns, &data);
  if (status!=0) {
    fprintf (stderr, "Error %d opening file %s\n", status, filename);
    return status;
  }

  if (print)
    fprintf (stderr, " ... counted %d data points\n", rows);
  
  if (rows == 0) { 
    fprintf (stderr, "%d data points, aborting.\n", rows);
    status = 99;
    ans = 0.0;
    fprintf (stdout, "%13.6e\n", ans);
    return status;
  }

  if (max_columns!=min_columns) {
    fprintf (stderr, " !! ATTENTION !! Inconsistent number of columns\n");
    fprintf (stderr, "     min = %d, max =%d\n", min_columns, max_columns);
  }
  
  if (min_columns<2) {
    fprintf (stderr, "Error, too few columns in %s\n", filename);
    return 0;
  }

  for (k=1; k<min_columns; k++) {

    x = ASCII_column (data, rows, 0);
    y = ASCII_column (data, rows, k);

    if (!adef)
      a = x[0];
    if (!bdef)
      b = x[rows-1];

    if (a<x[0]||a>x[rows-1]) {
      fprintf (stderr, "Error, a = %g out of range [%g, %g]!\n", a, x[0], x[rows-1]);
      return -1;
    }

    if (b<x[0]||b>x[rows-1]) {
      fprintf (stderr, "Error, b = %g out of range [%g, %g]!\n", b, x[0], x[rows-1]);
      return -1;
    }
    
    if (a>b) {
      fprintf (stderr, "Error, left boundary larger than right boundary!\n");
      return -1;
    }

    
    if (linear)
      status = integrate_linear (x, y, rows, a, b, &ans);
    else
      status = integrate_spline (x, y, rows, a, b, &ans);

    if (status!=0) {
      fprintf (stderr, "Error %d integrating\n", status);
      return status;
    }
    
    if (print)
      if (k==1)
	fprintf (stderr, "Results:\n");
    
    if (k > 1) 
      fprintf(stdout, "   ");
    
    fprintf(stdout, "%13.6e", ans);
  }
  
  fprintf(stdout, "\n");
  
  ASCII_free_double (data, rows);
  free(x);
  free(y);
  
  return 0;
}
