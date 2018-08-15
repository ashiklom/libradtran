/*--------------------------------------------------------------------
 * $Id: plkavg.c 2623 2011-12-23 10:52:38Z robert.buras $
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

/*  @19c@                

@code{plkavg} calculates the integrated Planck radiance [W / (m2 sterad)] 
for a given temperature and wavenumber interval.

The different options to @code{plkavg} are displayed when executing:

@example
plkavg -h
@end example

    @c19@ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "ascii.h"
#include "f77-uscore.h"
#include "solver.h"


#define PROGRAM "PLKAVG"
#define VERSION "0.80a"
#define DATE    "March 1, 2002"


static void print_usage (char *name)
{
  fprintf (stderr, "usage:   %s  <lower wvn> <upper wvn> <T>\n", name);
  fprintf (stderr, "%s calculates the integrated Planck radiance [W/(m2 sterad)]\n", PROGRAM);
  fprintf (stderr, "for a given temperature [K] and wavenumber interval [cm-1].\n\n");
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "  -h             Display this message.\n");
  fprintf (stderr, "  -f <filename>  Read input from three-column file.\n");
  fprintf (stderr, "  -g <filename>  Read wavenumbers from two-column file;\n");
  fprintf (stderr, "                 need to define temperature as well.\n");
  fprintf (stderr, "  -q             Be quiet.\n\n");
} 


/*************************************************************/
/* parse command line options                                */
/*************************************************************/

static int get_options (int argc, char **argv,
			char *programname,
			float *wvnlo,
			float *wvnhi,
			float *T,
			char *threefilename,
			char *twofilename,
			int *quiet)
{
  int c=0;
  char *dummy=NULL;

  /* reset parameters */
  strcpy (programname, "");
  strcpy (threefilename, "");
  strcpy (twofilename, "");
  *wvnlo = -1;
  *wvnhi = -1;
  *T     = -1;
  *quiet = 0;

  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  
  while ((c=getopt (argc, argv, "hf:g:q")) != EOF)  {
    switch(c)  {
    case 'f':
      strcpy (threefilename, optarg);
      break;
    case 'g':
      strcpy (twofilename, optarg);
      break;
    case 'q':
      *quiet=1;
      break;
    case 'h':
    case '?':
    default:
      print_usage (argv[0]);
      return (-1);
    }
  }
      

  /* check number of remaining command line arguments */

  if (strlen(twofilename)>0 && strlen(threefilename)>0) {
    print_usage (argv[0]);
    return -1;
  }

  if (strlen(threefilename)==0 && strlen(twofilename)==0)  {
    if (argc - optind != 3)  {
      print_usage (argv[0]);
      return -1;
    }

    *wvnlo = strtod (argv[optind+0], &dummy);
    *wvnhi = strtod (argv[optind+1], &dummy);
    *T     = strtod (argv[optind+2], &dummy);
  }

  if (strlen(threefilename)>0) {
    if (argc - optind != 0) {
      print_usage (argv[0]);
      return -1;
    }
  }

  if (strlen(twofilename)>0) {
    if (argc - optind != 1) {
      print_usage (argv[0]);
      return -1;
    }

    *T = strtod (argv[optind], &dummy);
  }



  return 0;  /* if o.k. */
}





int main (int argc, char **argv)
{
  int i=0, n=0, status=0, quiet=0;
  float r=0;

  float wvnlo=0, wvnhi=0, T=0;
  float *awvnlo=NULL, *awvnhi=NULL, *aT=NULL;

  char programname[FILENAME_MAX]="", threefilename[FILENAME_MAX]="",
    twofilename[FILENAME_MAX]="";


  /* get command line options */
  status = get_options (argc, argv,
			programname, &wvnlo, &wvnhi, &T, 
			threefilename, twofilename, &quiet);
  
  if (status!=0)
    return status;
  
  
  if (strlen(threefilename)>0) {
    if (!quiet)
      fprintf (stderr, " ... reading wavenumbers and temperature from %s\n", threefilename);

    status = read_3c_file_float (threefilename, 
				 &awvnlo, &awvnhi, &aT, &n);
    if (status!=0) {
      fprintf (stderr, "Error %d reading three column from %s\n", status, threefilename);
      return status;
    }
    
    for (i=0; i<n; i++) {
      F77_FUNC (cplkavg, CPLKAVG) (&(awvnlo[i]), &(awvnhi[i]), &(aT[i]), &r);
      fprintf (stdout, "%g %g %g %g\n", awvnlo[i], awvnhi[i], aT[i], r); 
    }
  }
  else {
    if (strlen(twofilename)>0) {
      if (!quiet)
	fprintf (stderr, " ... reading wavenumbers from %s\n", twofilename);

      status = read_2c_file_float (twofilename, 
				   &awvnlo, &awvnhi, &n);
      if (status!=0) {
	fprintf (stderr, "Error %d reading two column from %s\n", status, twofilename);
	return status;
      }
    
      for (i=0; i<n; i++) {
	F77_FUNC (cplkavg, CPLKAVG) (&(awvnlo[i]), &(awvnhi[i]), &T, &r);
	fprintf (stdout, "%g %g %g %g\n", awvnlo[i], awvnhi[i], T, r);
      }
    }
    else {
      if (!quiet)
	fprintf (stderr, " ... integrated Planck radiance, %g - %g cm-1, %g K:\n",
		 wvnlo, wvnhi, T);
      
      F77_FUNC (cplkavg, CPLKAVG) (&wvnlo, &wvnhi, &T, &r);
      fprintf (stdout, "%g\n", r);
    }
  }

  free(awvnlo); free(awvnhi); free(aT);

  return 0;
}
