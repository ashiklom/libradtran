/*--------------------------------------------------------------------
 * $Id: noon.c 2623 2011-12-23 10:52:38Z robert.buras $
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
#include <getopt.h>
#include "sun.h"
#include "sunpos.h"


#define IFU_LATITUDE   47.48
#define IFU_LONGITUDE -11.07
#define IFU_LONG_STD  -15.00


/*************************************************************/
/* print usage information                                   */
/*************************************************************/
/*
  Documentation in latex for the User's Guide. Extracted by sdoc2.awk 
  and input in the doc/tools.tex.
*/
/*
  <lpdoc>
  \subsection{Local noon time - \codeidx{noon}}

  The \code{noon} tool calculates the local noon time given a location 
  in terms of longitude and latitude or a location name using the \code{-l} option.
  Output is to stdout and is self-explanatory.

  The local noon time is calculated using the algorithm of \citet{BlancoMuriel2001}. 
  If the \code{-S} option is invoked the \citet{Spencer1971} algorithm is used.

  The \code{noon} tool is invoked by
  \begin{Verbatim}[fontsize=\footnotesize]
    noon [options] <day> <month>
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
   \item[-S]             Use the Spencer algorithm.
   \item[-h] Print help message.
\end{description}

The following invocation of \code{noon} calculates the noon time at the 
home location of one of the {\sl libRadtran} developers for his wedding date.
\begin{Verbatim}[fontsize=\footnotesize]
     noon -a 62.462052 -o -6.303358 -s -15 -y 1992 29 2
\end{Verbatim}


  </lpdoc>

*/

static void print_usage (char *filename)
{
  fprintf (stderr, "NOON 1.99b - calculate local noon time\n\n");
  fprintf (stderr, "written by  Bernhard Mayer,\n");
  fprintf (stderr, "            DLR, eMail bernhard.mayer@dlr.de\n");
  fprintf (stderr, "modified by  Arve Kylling to use Blanco-Muriel algorithm\n");
  fprintf (stderr, "Version 1.99b finished October 17, 2003\n\n");
  fprintf (stderr, "Be aware that this program is a beta version! Please report any\n");
  fprintf (stderr, "kind of error to the author, including the error message and the\n");
  fprintf (stderr, "database being processed when the error occured. Thanks!\n\n");
  fprintf (stderr, "USAGE: %s [options] <day> <month>\n\n", filename);
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
  fprintf (stderr, "  -S             Use the Spencer algorithm\n\n");
  fprintf (stderr, "If the coordinates are omitted, the location of Garmisch-Partenkirchen\n");
  fprintf (stderr, "is used for the calculation:\n");
  fprintf (stderr, "  Latitude 47.48, Longitude -11.07, Standard Longitude -15.00.\n");
}



/*************************************************************/
/* parse command line options                                */
/*************************************************************/

static int get_options (int argc, char **argv,
			char *programname,
			int *day,
			int *month,
			int *year,
			int *print,
			int *Spencer,
			char *locstr,
			double *latitude,
			double *longitude,
			double *long_std,
			int *stdlong)
{
  int c=0;
  char *dummy=NULL;
  int stdtime=0;

  char lat=0, lon=0;

  *day   = 0;
  *month = 0;


  /* reset parameters to default values */
  *latitude  = IFU_LATITUDE;
  *longitude = IFU_LONGITUDE;
  *long_std  = IFU_LONG_STD;
  
  strcpy (locstr, "");

  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  /* get command line arguments */
  /* The following line is not needed with the getopt.{c,h} I have */
  /* Arve Kylling, Aug 9, 1997 */
  /* optmode = GETOPT_ANY; */

  while ((c=getopt (argc, argv, "ha:o:ps:l:y:t:S")) != EOF)  {
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
    case 't':  /* time zone */  
      *long_std = (double) strtod (optarg, &dummy) * -15.0;
      stdtime=1;
      break;   
    case 'S':  /* Use Spencer sza algorithm */  
      *Spencer = 1;
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
  if (argc - optind != 2)  {
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
  *day    = strtod (argv[optind+0], &dummy);
  *month  = strtod (argv[optind+1], &dummy);

  return 0;  /* if o.k. */
}






int main(int argc, char ** argv)
{
  char programname[FILENAME_MAX]="";
  char locstr[255]="";

  double latitude=0, longitude=0, long_std=0, tmp_long_std=0, utcplus=0;
  double zenith=0;

  int doy=0, year=0, month=0, day=0;
  int hour=0, min=0, sec=0;
  int stdlong=0;

  int time=0;
  int shift=0;

  struct cTime udtTime;
  struct cLocation udtLocation;
  struct cSunCoordinates udtSunCoordinates;

  int print =1;
  int status=0;

  int Spencer=0;

  char dummy[10] = "";


  /* get command line options */
  status = get_options (argc, argv,
			programname, &day, &month, &year,
			&print, &Spencer, locstr, &latitude, &longitude, &long_std, &stdlong);

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
    if (!stdlong) 
      long_std = tmp_long_std;
  }

  /* calculate day of year */
  if ((doy = day_of_year(day, month)) < 0)  {
    print_usage (programname);
    fprintf (stderr, "\n** Date out of range **\n");
    return -1;
  }

  /* calculate LOCAL APPARENT TIME for midnight */
  /* time is calculated in UTC here!            */     
  /* standard longitude is corrected later.     */   
  shift = LAT (0, doy, longitude, 0);
  
  /* calculate local noon time */
  time = 43200 - shift;

  /* calculate solar zenith angle */
  if (Spencer)
    zenith  = solar_zenith  (time, doy, latitude, longitude, 0);
  else {
    
    hour = time / 3600;
    min  = (time - hour*3600) / 60;
    sec  = (time - hour*3600 - min*60);

    udtTime.iYear    = year;
    udtTime.iMonth   = month;
    udtTime.iDay     = day;
    udtTime.dHours   = hour;
    udtTime.dMinutes = min;
    udtTime.dSeconds = sec;
    
    udtLocation.dLatitude   = latitude;
    udtLocation.dLongitude  = -longitude;
    
    sunpos (udtTime, udtLocation, &udtSunCoordinates);
    zenith  = udtSunCoordinates.dZenithAngle;
  }
    

  /* correct time for standard longitude or time zone */
  time -= (int) ((long_std/15.0)*3600.0+0.5);
  if (time<0)
    time+=86400;
  while (time>86400) 
    time-=86400;
  hour = time / 3600;
  min  = (time - hour*3600) / 60;
  sec  = (time - hour*3600 - min*60);

  if (print) {
    fprintf (stderr, "Coordinates:\n");
    fprintf (stderr, " Latitude   %9.4f degree\n", latitude);
    fprintf (stderr, " Longitude  %9.4f degree\n", longitude);
    fprintf (stderr, " Std. Long. %9.4f degree\n", long_std);
    
    fprintf (stderr, " Local Noon Time:      %s ", 
	     time2str (dummy, hour,min,sec));

    utcplus = -long_std/15.0+1e-6;
    if (utcplus>=0)
      fprintf (stderr, "UTC+%.1f hours\n", utcplus);
    else 
      fprintf (stderr, "UTC-%.1f hours\n", -utcplus);
    
    fprintf (stderr, " Solar Zenith Angle: %9.4f degree\n\n", zenith);
    fprintf (stderr, " ... printing data to stdout\n");
    fprintf (stderr, " <local noon time>  <local noon time>  <solar zenith angle>\n");
    fprintf (stderr, "    [hh:mm:ss]          [hours]              [degrees]\n");
  }
  fprintf (stdout, "%s  %g  %g\n", time2str (dummy, hour,min,sec), time/3600.0, zenith);


  return 0;
}
