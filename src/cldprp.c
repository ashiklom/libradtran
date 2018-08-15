/*--------------------------------------------------------------------
 * $Id: cldprp.c 2859 2013-02-25 17:04:24Z robert.buras $
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

/*  @18c@                

@code{cldprp} calculates wavelength-dependent cloud properties 
using one of several parameterizations.

The different options to @code{cldprp} are displayed when executing:

@example
cldprp -h
@end example

    @c18@ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "ascii.h"
#include "numeric.h"
#include "ckdfu.h"
#include "cloud.h"
#include "f77-uscore.h"

#if HAVE_KEY56
#include "yang56.h"
#endif



#define PROGRAM "CLDPRP"
#define DATE    "November 23, 2001"


#define HU   1
#define KEY  2
#define YANG 3
#define FU96 4
#define FU98 5


static void print_usage (char *name)
{
  fprintf (stderr, "usage:   %s  input_file \n", name);
  fprintf (stderr, "%s calculates optical properties of water clouds\n", PROGRAM);
  fprintf (stderr, "according to Hu and Stamnes [1993], and ice clouds according to\n");
  fprintf (stderr, "Key et al. [2002] / Yang et al. [2002], Fu [1996] or Fu et al. [1998].\n\n");
  fprintf (stderr, "The wavelengths covered are:\n");
  fprintf (stderr, "Hu and Stamnes [1993]: 0.290-150 um, water droplets.\n");
  fprintf (stderr, "Fu [1996]            : 0.250-4.99 um, hexagonal ice crystals.\n");
  fprintf (stderr, "Fu et al. [1998]     : 3.969-100 um, hexagonal ice crystals.\n");
  fprintf (stderr, "Key et al. [2002]    : 0.2-5.0 um, various ice crystals.\n");
  fprintf (stderr, "Yang et al. [2000]   : 0.2-100.0 um, various ice crystals.\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Required arguments:\n");
  fprintf (stderr, "  The input file is a three-column file with wavelength [nm],\n");
  fprintf (stderr, "  liquid (ice) water content [g/m3], and effective droplet radius [um];\n");
  fprintf (stderr, "  for the Yang et al parameterization, a 4th column may be specified\n");
  fprintf (stderr, "  with the particle habit. The following particle habits are recognised\n");
  fprintf (stderr, "  (enter the single digit number representing the wanted particle habit\n");
  fprintf (stderr, "   in the 4th column if using Yang et al.)\n");
  fprintf (stderr, "         0 : Solid column\n");
  fprintf (stderr, "         1 : Hollow column\n");
  fprintf (stderr, "         2 : Rough aggregate\n");
  fprintf (stderr, "         3 : Rosette 4\n");
  fprintf (stderr, "         4 : Rosette 6\n");
  fprintf (stderr, "         5 : Plate\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "(the default parameterization is Hu and Stamnes[1993])\n\n");
  fprintf (stderr, "  -h            Show this help page\n");
  fprintf (stderr, "  -f            Use Fu et al. [1998] ice cloud parameterizations\n");
  fprintf (stderr, "  -g            Use Fu  [1996] ice cloud parameterizations\n");
  fprintf (stderr, "  -k            Use Key et al. [2002] ice cloud parameterization\n");
  fprintf (stderr, "  -y            Use Yang et al. [2000] ice cloud parameterization\n");
  fprintf (stderr, "  -u            Provide unscaled properties rather than scaled\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "The output is six-columns: wavelength [nm],\n");
  fprintf (stderr, "liquid (ice) water content [g/m3], effective droplet radius [um],\n");
  fprintf (stderr, "the extinction coefficient [1/km], the asymmetry parameter,\n");
  fprintf (stderr, "and the single scattering albedo. For the Key et al. / Yang et al. parameterizations,\n");
  fprintf (stderr, "g1, g2, and f (double Henyey-Greenstein, Key et al. [2002], equation 12)\n");
  fprintf (stderr, "are printed as well. For the Fu [1996] parameterization the\n"); 
  fprintf (stderr, "forward deltafraction is output as the seventh column.\n");
  
}




/*************************************************************/
/* parse command line options                                */
/*************************************************************/

static int get_options (int argc, char **argv,
			char *programname,
			char *infilename,
			int *properties,
			int *unscaled)
{
  int c=0;

  /* reset parameters */
  strcpy (infilename, "");
  *properties=HU;
  *unscaled=0;

  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  

  while ((c=getopt (argc, argv, "fghkyu")) != EOF)  {
    switch(c)  {
    case 'f':  /* use Fu et al [1998] */
      *properties = FU98;
      break;
    case 'g':  /* use Fu [1996] */
      *properties = FU96;
      break;
    case 'y':  /* use Yang et al [2000] */
      *properties = YANG;
      break;
    case 'k':  /* use Key et al [2002]  */
      *properties = KEY;
      break;
    case 'u':  /* unscaled properties   */
      *unscaled = 1;
      break;
    case 'h':  /* help */
    case '?':
    default:
      print_usage (programname);
      return (-1);
    }
  }
      
  /* check number of remaining command line arguments */
  if (argc - optind != 1)  {
    print_usage (argv[0]);
    return -1;
  }

  strcpy (infilename, argv[optind]);

  return 0;  /* if o.k. */
}





int main (int argc, char **argv)
{
  int i=0, newkey=0, status=0;
  char path[FILENAME_MAX]="";
  char yangpath[FILENAME_MAX]="";
  char fupath[FILENAME_MAX]="";
  int npath;

  float *wvl=NULL, *lwc=NULL, *reff=NULL, *habit=NULL;
  int hab=0;

  int n=0;

  int properties=0, unscaled=0;

  void  F77_FUNC (wcloud, WCLOUD) (float *lambda_r, int *newsiz, int *nlyr,
			  char *filepath, int *nstring, float *wccon, float *wceffr, 
			  float *dtau, float *gg, float *ssa, float *zd,
			  int *wclyr);


  char programname[FILENAME_MAX], infilename[FILENAME_MAX];

  int newsiz=1, nlyr=1, iclayer=1, wclayer=1, nohabit=0;
  float tau=0, gg=0, ssa=0, g1=0, g2=0, ff=0, f=0;

  float zd   [2] = {1.0, 0.0};
  float effr [2] = {0.0, 0.0};
  float wccon[2] = {0.0, 0.0};
  float iwc  [2] = {0.0, 0.0};

  /* get command line options */
  status = get_options (argc, argv,
			programname, infilename, &properties, &unscaled);
  
  if (status!=0)
    return status;
  

  fprintf (stderr, " ... reading data from %s\n", infilename);


  /* read wavelength file */
  status = read_4c_file_float (infilename, &wvl, &lwc, &reff, &habit, &n);

  if (status!=0) {
    nohabit=1;

    if (properties==YANG || properties==KEY)
      fprintf (stderr, " ... warning, no habit defined, assuming 0\n");

    /* try to read 3 columns */
    status = read_3c_file_float (infilename, &wvl, &lwc, &reff, &n);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, infilename);
      return status;
    }
  }


  strcpy (path, "../data/");
  strcpy (yangpath, "../data/ic/yang56/");
  strcpy (fupath, "../data");


  npath = strlen(path);

  for (i=0; i<n; i++) {
    effr [0] = reff[i]; effr [1] = reff[i];
    wccon[0] = lwc [i]; wccon[1] = lwc [i];
    iwc  [0] = lwc [i]; iwc  [1] = lwc [i];

    if (!nohabit)
      hab = (int) (habit[i]+0.5);
    else 
      hab = 0;

    switch(properties) {
      
    case FU96:
      status = ic_fu96 (wvl[i], nlyr, fupath, iwc, effr, 
			&tau, &gg, &ssa, &f, zd, iclayer, unscaled);
      
      if (status!=0) {
	fprintf (stderr, "Error %d returned by ic_fu96()\n", status);
	return status;
      }
      
      fprintf (stdout, "%g %g %g %g %g %g %g\n", 
	       wvl[i], lwc[i], reff[i], tau, gg, ssa, f);
      
      break;
      
    case FU98:
      status = ic_fu98 (wvl[i], nlyr, fupath, iwc, effr, 
			&tau, &gg, &ssa, zd, iclayer);
      
      if (status!=0) {
	fprintf (stderr, "Error %d returned by ic_fu98()\n", status);
	return status;
      }
      
      fprintf (stdout, "%g %g %g %g %g %g\n", 
	       wvl[i], lwc[i], reff[i], tau, gg, ssa);
      
      break;
      
    case HU:
      F77_FUNC (wcloud, WCLOUD) (&(wvl[i]), &newsiz, &nlyr,
			path, &npath, wccon, effr, 
			&tau, &gg, &ssa, 
			zd, &wclayer);
      

      
      fprintf (stdout, "%g %g %g %g %g %g\n", wvl[i], lwc[i], reff[i], tau, gg, ssa);
      
      break;
      
    case KEY:
    case YANG:

      if (properties==KEY)
	newkey=0;

      if (properties==YANG)
	newkey=1;

      if (newkey) {
#if HAVE_YANG
	status = yang56 (wvl[i]/1000.0, reff[i], hab, yangpath, 
			 &tau, &ssa, &g1, &g2, &ff, newkey);
	
	if (status!=0) {
	  fprintf (stderr, "Error %d returned by yang56()\n", status);
	  return status;
	}
#else
	fprintf (stderr, " ***********************************************************************\n");
	fprintf (stderr, " * You have built uvspec without Yang/Key/Mayer support and hence      *\n");
	fprintf (stderr, " * cannot use 'ic_properties yang'.                                    *\n");
	fprintf (stderr, " ***********************************************************************\n");
	return -1;
#endif
      }
      else {

	status = yang56 (wvl[i]/1000.0, reff[i], hab, yangpath, 
			 &tau, &ssa, &g1, &g2, &ff, newkey);
	
	if (status!=0) {
	  fprintf (stderr, "Error %d returned by yang56()\n", status);
	  return status;
	}
      }

      fprintf (stdout, "%g %g %g %g %g %g %g %g %g\n", wvl[i], lwc[i], reff[i], 
	       tau*lwc[i], g1*ff+(1.0-ff)*g2, ssa,
	       g1, g2, ff);

      break;

    default:
      fprintf (stderr, "Error, unknown cloud properties %d", properties);
      return status;
    }
  }


  return 0;
}
