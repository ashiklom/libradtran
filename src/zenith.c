/*--------------------------------------------------------------------
 * $Id: zenith.c 2623 2011-12-23 10:52:38Z robert.buras $
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
#include <getopt.h>
#include <string.h> 
#include <math.h> 
#include "sun.h"
#include "sunpos.h"


#define PROGRAM "ZENITH"
#define VERSION "1.99g"
#define DATE    "July 6, 2004"


/*************************************************************/
/* print usage information                                   */
/*************************************************************/

/*
  Documentation in latex for the User's Guide. Extracted by sdoc2.awk 
  and input in the doc/tools.tex.
*/
/*
  <lpdoc>
  \subsection{Solar zenith and azimuth angle - \codeidx{zenith}\index{solar zenith angle}}

  The \code{zenith} tool calculates the solar zenith and azimuth angle for
  a given time and location. Output is to stdout and is self-explanatory (unless 
  the \code{-q} option is used).

  The solar zenith and azimuth angles are calculated using the algorithm of 
  \citet{BlancoMuriel2001}. If the \code{-S} option is invoked the 
  \citet{Spencer1971} algorithm is used.

  The \code{zenith} tool is invoked by
  \begin{Verbatim}[fontsize=\footnotesize]
    zenith [options] <day> <month> <hour> <min> [sec]
  \end{Verbatim}
   where the various options are
   \begin{description}
      \item[-a] $<$latitude$>$  Latitude  (North positive)
      \item[-o] $<$longitude$>$ Longitude (West positive)
      \item[-s] $<$std. long$>$ Standard Longitude (West positive)
                 this is the longitude to which the time zone refers
                 (-15 deg for central Europe, corresponds to UTC+1).
      \item[-l] $<$location$>$  Instead of \code{-a}, \code{-o} and \code{-s} define a location.
                 possible locations are ifu, dlrop.
      \item[-y] $<$yyyy$>$      year; not used if -S specified, default: 2003.
      \item[-S] Use the Spencer algorithm.
      \item[-e] Calculate eccentricity.
      \item[-t] $<$UTC + x$>$   Time zone; e.g. -t2 means UTC + 2.
      \item[-q] Be quiet.
      \item[-h] Print help message.
   \end{description}

   The options below apply if the solar zenith angle is wanted as a
   function of wavelength. This is useful for simulation of scanning 
   spectroradiometer measurements. Output is two columns with wavelength 
   and solar zenith angle. All options must be specfied. However $<$hour$>$
    and $<$min$>$ should not be specified. To avoid too much output use 
   the \code{-q} option.

   \begin{description}
      \item[-B] start\_time (decimal hours of Greenwich time)
      \item[-E] end\_time (decimal hours of Greenwich time)
      \item[-u] start\_wavelength (nanometers)
      \item[-v] end\_wavelength (nanometers)
      \item[-w] step\_wavelength (nanometers)
   \end{description}


The following invocation of \code{zenith} calculates the solar
zenith and azimuth angles at the time and location of the writing
of this text
\begin{Verbatim}[fontsize=\footnotesize]
     zenith -a 62.462052 -o -6.303358 -s -15 -y 2010 4 3 9 35
\end{Verbatim}

  </lpdoc>

*/

static void print_usage (char *filename)
{
  fprintf (stderr, "%s %s - calculate solar zenith and azimuth\n\n", PROGRAM, VERSION);
  fprintf (stderr, "written by  Bernhard Mayer,\n");
  fprintf (stderr, "            DLR, eMail bernhard.mayer@dlr.de\n");
  fprintf (stderr, "modified by  Arve Kylling to include Blanco-Muriel algorithm\n");
  fprintf (stderr, "Version %s finished %s\n\n", VERSION, DATE);
  fprintf (stderr, "Be aware that this program is a beta version! Please report any\n");
  fprintf (stderr, "kind of error to the author, including the error message and the\n");
  fprintf (stderr, "database being processed when the error occured. Thanks!\n\n");
  fprintf (stderr, "USAGE: %s [options] <day> <month> <hour> <min> [sec]\n\n", filename);
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "  -h             display help message\n");
  fprintf (stderr, "  -a <latitude>  Latitude  (North positive)\n");
  fprintf (stderr, "  -o <longitude> Longitude (West positive !!!)\n");
  fprintf (stderr, "  -s <std. long> Standard Longitude (West positive !!!);\n");
  fprintf (stderr, "                 this is the longitude to which the time zone refers\n");
  fprintf (stderr, "                 (-15 deg for central Europe, corresponds to UTC+1)\n");
  fprintf (stderr, "  -l <location>  Instead of \"-a -o -s\" define a location;\n");
  fprintf (stderr, "                 possible locations are ifu, dlrop.\n");
  fprintf (stderr, "  -y <yyyy>      year; not used if -S specified, default: 2003\n");
  fprintf (stderr, "  -e             Calculate eccentricity\n");
  fprintf (stderr, "  -q             Be quiet\n");
  fprintf (stderr, "  -t <UTC + x>   Time zone; e.g. -t2 means UTC + 2\n");
  fprintf (stderr, "  -S             Use Spencer [1971] algorithm instead of default\n\n");
  fprintf (stderr, " \n");
  fprintf (stderr, "The following options apply if the solar zenith angle is wanted as a \n");
  fprintf (stderr, "function of wavelength. Useful for simulation of scanning spectroradiometer \n");
  fprintf (stderr, "measurements. Output is two columns with wavelength and solar zenith angle.\n");
  fprintf (stderr, "All options must be specfied. However <hour> and <min> should not be \n");
  fprintf (stderr, "specified. To avoid too much output use -q option.\n");
  fprintf (stderr, " \n");
  fprintf (stderr, "  -B start_time (decimal hours of Greenwich time)\n");
  fprintf (stderr, "  -E end_time (decimal hours of Greenwich time)\n");
  fprintf (stderr, "  -u start_wavelength (nanometers)\n");
  fprintf (stderr, "  -v end_wavelength (nanometers)\n");
  fprintf (stderr, "  -w step_wavelength (nanometers)\n");
  fprintf (stderr, " \n");
  fprintf (stderr, "The default algorithm is taken from Blanco-Muriel et al., Solar \n");
  fprintf (stderr, "Energy, 70, 431-441, 2001.\n");
}



/*************************************************************/
/* parse command line options                                */
/*************************************************************/

static int get_options (int argc, char **argv,
			char *programname,
			int *doy,
			int *print,
			int *time_std,
			int *month, 
			int *day,
			int *hour, 
			int *min, 
			int *sec,
			int *year,
			double *latitude,
			double *longitude,
			double *long_std,
			int *stdlong,
			int *ecc,
			int *Spencer,
			char *locstr,
			int *wvlsza,
			double *start_time,
			double *end_time,
			double *start_wvl,
			double *end_wvl,
			double *step_wvl)
{
  int c=0;
  char *dummy=NULL;
  int stdtime=0;

  char lat=0, lon=0;

  *day   = 0;
  *month = 0;
  *hour  = 0;
  *min   = 0;
  *sec   = 0;


  /* reset parameters to default values */
  *latitude  = -999;
  *longitude = -999;
  *long_std  = 0;  /* Greenwich, UTC */
  strcpy (locstr, "");
  *ecc       = 0;

  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  while ((c=getopt (argc, argv, "ha:o:pqs:l:y:et:SB:E:u:v:w:")) != EOF)  {
    switch(c)  {
    case 'h':  /* help */
      print_usage (programname);
      return (-1);
      break;
    case 'a':  /* latitude  */  
      *latitude  = (double) strtod (optarg, &dummy);
      lat=1;
      break;   
    case 'o':  /* longitude */  
      *longitude = (double) strtod (optarg, &dummy);
      lon=1;
      break;   
    case 'p':  /* no print messages */  
    case 'q':  /* no print messages */  
      *print = 0;
      break;   
    case 's':  /* standard longitude */  
      *long_std = (double) strtod (optarg, &dummy);
      *stdlong=1;
      break;   
    case 'y':  /* year */  
      *year = (int) strtod (optarg, &dummy);
      break;   
    case 'l':
      strcpy (locstr, optarg);
      break;
    case 'e':  /* eccentricity */  
      *ecc = 1;
      break;   
    case 't':  /* time zone */  
      *long_std = (double) strtod (optarg, &dummy) * -15.0;
      stdtime=1;
      break;   
    case 'S':  /* Use Spencer sza algorithm */  
      *Spencer = 1;
      break;   
    case 'B':  /* Start time */  
      *start_time = (double) strtod (optarg, &dummy);
      (*wvlsza)++;
      break;   
    case 'E':  /* End time */  
      *end_time = (double) strtod (optarg, &dummy);
      (*wvlsza)++;
      break;   
    case 'u':  /* Start wavelength */  
      *start_wvl = (double) strtod (optarg, &dummy);
      (*wvlsza)++;
      break;   
    case 'v':  /* End wavelength */  
      *end_wvl = (double) strtod (optarg, &dummy);
      (*wvlsza)++;
      break;   
    case 'w':  /* Wavelength step */  
      *step_wvl = (double) strtod (optarg, &dummy);
      (*wvlsza)++;
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
  if (argc - optind != 2 && *wvlsza > 0 && *wvlsza!=5) {
    print_usage (programname);
    return -1;
  }
  if (argc - optind != 4 && argc - optind != 5 && *wvlsza==0) {
    print_usage (programname);
    return -1;
  }
       
  /* check arguments */
  if ((!lat&!lon) && strlen(locstr)==0) {
    print_usage (programname);
    fprintf (stderr, "\n** You need to specify a location, either by latitude and longitude");
    fprintf (stderr, "\n** or by location name!!\n");
    
    return -1;
  }

  if ((lat||lon) && strlen(locstr)>0) {
    print_usage (programname);
    fprintf (stderr, "\n** It does not make sense to define both location and\n");
    fprintf (stderr,   "** explicit latitude/longitude!!\n");

    return -1;
  }

  if (strlen(locstr)==0 && ((lat&&!lon) || (!lat&&lon))) {
    print_usage (programname);
    fprintf (stderr, "\n** You need to define both latitude and longitude!!\n");
    return -1;
  }

  /* check coordinates range if latitude and longitude have been */
  /* defined explicitely                                         */
  if (strlen(locstr)==0 && (*latitude<-90 || *latitude>90))  {
    fprintf (stderr, "\n** Latitude out of range **\n");
    return -1;
  }

  if (strlen(locstr)==0 && (*longitude<-180 || *longitude>180))  {
    fprintf (stderr, "\n** Longitude out of range **\n");
    return -1;
  }

  if (*long_std<-180 || *long_std>180)  {
    fprintf (stderr, "\n** Standard Longitude out of range **\n");
    return -1;
  }

  if (strlen(locstr)==0 && !*stdlong && !stdtime) {
    *long_std=0;  /* UTC, yeeeeees */
    fprintf (stderr, " ... standard longitude not specified - using default %3.1f (UTC)\n", *long_std);
  }

  if (*stdlong && stdtime) {
    print_usage (programname);
    fprintf (stderr, "\n** It does not make sense to define both standard longitude\n");
    fprintf (stderr,   "** and time zone!!\n");

    return -1;
  }

  if (*wvlsza != 5 && *wvlsza>0) {
    print_usage (programname);
    fprintf (stderr, "\n** All of the -B, -E, -u, -v and -w options must be specified.\n");

    return -1;
  }

  if (stdtime)    /* already converted time zone to standard longitude */
    *stdlong=1;   /* when reading time zone; the two variables are     */
                  /* only required for the above tests                 */
                 
  if (!*Spencer && *year <= 0) {
    *year = 2003; /* Default year. 2003 is chosen for no specific reason, */
                  /* although it was the year I onsigthed                 */
                  /* "Solo Con I Tuoi Pedali" - AK                        */
    if (*print) {
      fprintf (stderr, " ... year not specified - using default %04d\n", *year);
      fprintf (stderr, " ... (for more accurate results specify year by using -y yyyy)\n");
    }
  }

  /* get command line arguments */
  *day   = (int) strtod (argv[optind+0], &dummy);
  *month = (int) strtod (argv[optind+1], &dummy);
  if (*wvlsza == 0) {
    *hour  = (int) strtod (argv[optind+2], &dummy);
    *min   = (int) strtod (argv[optind+3], &dummy);
  }
  if (argc - optind == 5)
    *sec = strtod (argv[optind+4], &dummy);
  else
    *sec = 0; 


  /* calculate standard time [sec] and check range */
  if ((*time_std = *hour*3600 + *min*60 + *sec) > 86400)  {
    print_usage (programname);
    fprintf (stderr, "\n** Time out of range **\n");
    return -1;
  }

  /* calculate day of year */
  if ((*doy = day_of_year(*day, *month)) < 0)  {
    print_usage (programname);
    fprintf (stderr, "\n** Date out of range **\n");
    return -1;
  }

  return 0;  /* if o.k. */
}


static void dectime2hhmmss(double dectime, int *hour, int *min, int *sec)
{
  char s[10] = "";
  char *dummy=NULL;

  sprintf(s, "%2f", dectime);
  *hour = (int) strtod (s, &dummy);

  sprintf(s, "%2f", (dectime-*hour)*60);
  *min  = (int) strtod (s, &dummy);

  sprintf(s, "%2f", ((dectime-*hour)*60-*min)*60);
  *sec  = (int) strtod (s, &dummy);

}


int main(int argc, char ** argv)
{
  char programname[FILENAME_MAX]="";
  char locstr[255]="";

  int doy=0, std=0;

  double zenith    =  0.0;
  double azimuth   =  0.0;
  
  double latitude=0, longitude=0, long_std=0, tmp_long_std=0, utcplus=0;

  struct cTime udtTime;
  struct cLocation udtLocation;
  struct cSunCoordinates udtSunCoordinates;

  int print =1;
  int status=0;

  int time_std=0, time_lat=0;

  int hour_lat  = 0;
  int min_lat   = 0;
  int sec_lat   = 0;

  int year  = 0;
  int month = 0;
  int day   = 0;
  int hour  = 0;
  int min   = 0;
  int sec   = 0;

  int s_time = 0;
  int t_step = 0;

  char dummy[10] = "";

  int ecc=0;
  double eccent=0;

  int wvlsza=0;
  double start_time=0;
  double end_time=0;
  double start_wvl=0;
  double end_wvl=0;
  double step_wvl=0;
  double tmptime=0;

  int Spencer=0;

  float lambda;
  int i, n_lambda=0;


  /* get command line options */
  status = get_options (argc, argv,
			programname,
			&doy, &print, &time_std, &month, &day, &hour, &min, &sec, &year,
			&latitude, &longitude, &long_std, &std, &ecc, &Spencer, locstr,
			&wvlsza, &start_time, &end_time, &start_wvl, &end_wvl, &step_wvl);

  if (status!=0)
    return status;

  if (strlen(locstr)>0) {
    if (print)
      fprintf (stderr, " ... using location %s\n", locstr);

    status = location (locstr, &latitude, &longitude, &tmp_long_std);
    if (status!=0) {
      fprintf (stderr, "Error %d calculating coordinates of location %s\n", 
	       status, locstr);
      return status;
    }

    /* if no time zone specified */
    if (!std) 
      long_std = tmp_long_std;
  }

  if (wvlsza==0) {
    n_lambda = 1;
    s_time = time_std;
  }
  else {
    n_lambda = (int) ((end_wvl-start_wvl)/step_wvl) + 1.000001;
    t_step = (int) 3600*(end_time-start_time) / (n_lambda-1);
    s_time = (int) 3600*start_time;
  }

  for (i=0;i<n_lambda;i++) {
    time_std = s_time + i*t_step;
    lambda = start_wvl + i*step_wvl;

    /* calculate solar zenith angle */
    if (Spencer) {
      zenith  = solar_zenith  (time_std, doy, latitude, longitude, long_std);
      azimuth = solar_azimuth (time_std, doy, latitude, longitude, long_std);
      
    }
    else {
      
      if (wvlsza==5) {
	tmptime = (double) time_std/3600;
	dectime2hhmmss(tmptime, &hour, &min, &sec);
      }

      udtTime.iYear    = year;
      udtTime.iMonth   = month;
      udtTime.iDay     = day;
      udtTime.dHours   = hour;
      udtTime.dMinutes = min;
      udtTime.dSeconds = sec;
      
      udtLocation.dLatitude   = latitude;
      udtLocation.dLongitude  = -(longitude-long_std);
      
      sunpos (udtTime, udtLocation, &udtSunCoordinates);
      zenith  = udtSunCoordinates.dZenithAngle;
      azimuth = udtSunCoordinates.dAzimuth-180;
      
    }
    
    while (azimuth>=360.0)
      azimuth-=360.0;

    while (azimuth<0.0)
      azimuth+=360.0;


    /* calculate LOCAL APPARENT TIME */
    time_lat = LAT (time_std, doy, longitude, long_std);
    /* correct negative LATs and LATs larger than 24:00 */
    /* to get human-readable output                     */
    while (time_lat<0)
      time_lat += 86400;
    
    while (time_lat>86400)
      time_lat -= 86400;
    
    hour_lat = time_lat / 3600;
    min_lat  = (time_lat - hour_lat*3600) / 60;
    sec_lat  = (time_lat - hour_lat*3600 - min_lat*60);
    
    if (ecc)
      eccent = eccentricity (doy);
    
    if (print) {
      if (Spencer)
	fprintf (stderr, " ... using Spencer (1971) to calculate solar zenith and azimuth\n");
      else
	fprintf (stderr, " ... using Blanco-Muriel et al. (2001) to calculate solar zenith and azimuth\n");
      
      fprintf (stderr, "\n Latitude    %9.4f degree",   latitude);
      fprintf (stderr, "\n Longitude   %9.4f degree",   longitude);
      fprintf (stderr, "\n Std. Long.  %9.4f degree\n", long_std);
      
      utcplus = -long_std/15.0+1e-6;
      if (utcplus>=0)
	fprintf (stderr, " Time Zone:              UTC+%.1f hours\n", utcplus);
      else 
	fprintf (stderr, " Time Zone:              UTC-%.1f hours\n", -utcplus);
      
      fprintf (stderr, " Standard Time:          %s\n", time2str(dummy, hour,min,sec));
      fprintf (stderr, " Local Apparent Time:    %s\n", time2str(dummy, hour_lat,min_lat,sec_lat));
      fprintf (stderr, " Solar Zenith  Angle:  %8.3f degree\n", zenith);
      fprintf (stderr, " Solar Azimuth Angle:  %8.3f degree\n", azimuth);
      if (ecc)
	fprintf (stderr, " Eccentricity:           %6.4f\n", eccent);
      
      fprintf (stderr, "\n ... output to stdout: time, zenith, azimuth");
      if (ecc)
	fprintf (stderr, ", eccentricity\n");
      else 
	fprintf (stderr, "\n");
    }
    
    if (wvlsza>0) {
      fprintf (stdout, "%8.3f %7.4f",lambda, zenith);
    }
    else {
      fprintf (stdout, "%s %9.4f %9.4f", time2str(dummy, hour, min, sec), zenith, azimuth);
    }
    
    if (ecc)
      fprintf (stdout, " %9.6f\n", eccent);
    else 
      fprintf (stdout, "\n");
    
  }
  return 0;   /* if o.k. */
}
