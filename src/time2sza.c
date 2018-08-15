/*--------------------------------------------------------------------
 * $Id: time2sza.c 2623 2011-12-23 10:52:38Z robert.buras $
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
#include <math.h>

#include "getopt.h"
#include "f77-uscore.h"

/*************************************************************/
/* print usage information                                   */
/*************************************************************/

static void usage()
{
  fprintf (stderr, "\nConverts start and stop times for a scan into\n");
  fprintf (stderr, "solar zenith angles as a function of wavelength.\n");
  fprintf (stderr, "All options must be set.\n\n");
  fprintf (stderr, "No error checking yet, so be careful.\n\n");

  fprintf (stderr, "Usage: time2sza [-hdlbstuvw]\n");
  fprintf (stderr, " -d day of year (Jan 1 is day 1)\n");
  fprintf (stderr, " -b latitude (degrees positive north from Equator)\n");
  fprintf (stderr, " -l longitude (degrees positive west from Greenwich)\n");
  fprintf (stderr, " -s start_time (decimal hours of Greenwich time)\n");
  fprintf (stderr, " -t end_time (decimal hours of Greenwich time)\n");
  fprintf (stderr, " -u start_wavelength (nanometers)\n");
  fprintf (stderr, " -v end_wavelength (nanometers)\n");
  fprintf (stderr, " -w step_wavelength (nanometers)\n");
  fprintf (stderr, " -h    Print this message.\n");  
}

int main(int argc, char **argv) 
{

  /*  extern doublereal F77_FUNC (cozena, COZENA) ();
  doublereal sza=0.0;*/
  double F77_FUNC (cozena, COZENA) (float *day, float *hour, float *dlat, float *dlon);
  double sza=0.0;

  int c=0;
  float day=0.0, hour=0.0, latitude=0.0, longitude=0.0, start_time=0.0,
        end_time=0.0, start_wvn=0.0, end_wvn=0.0, step_wvn=0.0;

  float lambda, time_step;
  float pi = 3.1415926;
  int i, n_lambda=0;

  while ((c=getopt (argc, argv, "hd:l:b:s:t:u:v:w:")) != EOF)  {
    switch(c)  {
    case 'b': 
      latitude = atof(optarg);
      break;
    case 'd': 
      day = atof(optarg);
      break;
    case 'l': 
      longitude = atof(optarg);
      break;
    case 's': 
      start_time = atof(optarg);
      break;
    case 't': 
      end_time = atof(optarg);
      break;
    case 'u': 
      start_wvn = atof(optarg);
      break;
    case 'v': 
      end_wvn = atof(optarg);
      break;
    case 'w': 
      step_wvn = atof(optarg);
      break;
    case 'h': 
      usage();
      return (-1);
      break;
    default:
      usage();
      return (-1);
    }
  }
  n_lambda = (int) ((end_wvn-start_wvn)/step_wvn) + 1.000001;
  time_step = (end_time-start_time) / (n_lambda-1);
  for (i=0;i<n_lambda;i++) {
    hour = start_time + i*time_step;
    lambda = start_wvn + i*step_wvn;
    sza = acos(F77_FUNC (cozena, COZENA) (&day, &hour, &latitude, &longitude))*180./pi;
    printf("%8.3f %7.4f\n",lambda, sza);
  }
  return 0;
}


