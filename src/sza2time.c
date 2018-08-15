/*--------------------------------------------------------------------
 * $Id: sza2time.c 2623 2011-12-23 10:52:38Z robert.buras $
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "sun.h"


#define PROGRAM "SZA2TIME"
#define VERSION "1.99d"
#define DATE    "July 6, 2004"


/*************************************************************/
/* print usage information                                   */
/*************************************************************/

static void print_usage (char *filename)
{
  fprintf (stderr, "%s %s - calculate standard time from solar zenith angle\n\n", PROGRAM, VERSION);
  fprintf (stderr, "written by  Bernhard Mayer,\n");
  fprintf (stderr, "            DLR, eMail bernhard.mayer@dlr.de\n");
  fprintf (stderr, "Version %s finished %s\n\n", VERSION, DATE);
  fprintf (stderr, "Be aware that this program is a beta version! Please report any\n");
  fprintf (stderr, "kind of error to the author, including the error message and the\n");
  fprintf (stderr, "database being processed when the error occured. Thanks!\n\n");
  fprintf (stderr, "USAGE: %s [options] <day> <month> <zenith>\n\n", filename);
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "  -h             display help message\n");
  fprintf (stderr, "  -a <latitude>  Latitude  (North positive)\n");
  fprintf (stderr, "  -o <longitude> Longitude (West positive !!!)\n");
  fprintf (stderr, "  -s <std. long> Standard Longitude (West positive !!!);\n");
  fprintf (stderr, "                 this is the longitude to which the time zone refers\n");
  fprintf (stderr, "                 (-15 deg for central Europe, corresponds to UTC+1)\n");
  fprintf (stderr, "  -l <location>  Instead of \"-a -o -s\" define a location;\n");
  fprintf (stderr, "                 possible locations are ifu, dlrop.\n");
  fprintf (stderr, "  -q             Be quiet\n");
  fprintf (stderr, "  -t <UTC + x>   Time zone; e.g. -t2 means UTC + 2\n");
  fprintf (stderr, "The algorithm is taken from Spencer (1971),\n");
  fprintf (stderr, "as described in Iqbal, An Introduction to Solar Radiation\n");
}



/*************************************************************/
/* parse command line options                                */
/*************************************************************/

static int get_options (int argc, char **argv,
			char *programname,
			int *doy,
			int *print,
			double *zenith,
			double *latitude,
			double *longitude,
			double *long_std,
			int *stdlong,
			char *locstr)
{
  int c=0;
  char *dummy=NULL;
  int stdtime=0;

  int day   = 0;
  int month = 0;

  char lat=0, lon=0;


  /* reset parameters to default values */
  *latitude  = -999;
  *longitude = -999;
  *long_std  = 0;  /* Greenwich, UTC */
  strcpy (locstr, "");


  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  /* get command line arguments */
  /* The following line is not needed with the getopt.{c,h} I have */
  /* Arve Kylling, Jul 25, 1997 */
  /* optmode = GETOPT_ANY; */

  while ((c=getopt (argc, argv, "ha:o:pqs:l:t:")) != EOF)  {
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
    case 'l':
      strcpy (locstr, optarg);
      break;
    case 's':  /* standard longitude */  
      *long_std = (double) strtod (optarg, &dummy);
      *stdlong=1;
      break;   
    case 't':  /* time zone */  
      *long_std = (double) strtod (optarg, &dummy) * -15.0;
      stdtime=1;
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
  if (argc - optind != 3)  {
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

  *stdlong=1;  /* already converted time zone to standard longitude */


  /* get command line arguments */
  day    = strtod (argv[optind+0], &dummy);
  month  = strtod (argv[optind+1], &dummy);
  *zenith = strtod (argv[optind+2], &dummy);

  if (*zenith<0 || *zenith>180)  {
    fprintf (stderr, "\n** Solar zenith angle out of range **\n");
    return -1;
  }

  /* calculate day of year */
  if ((*doy = day_of_year(day, month)) < 0)  {
    print_usage (programname);
    fprintf (stderr, "\n** Date out of range **\n");
    return -1;
  }


  return 0;  /* if o.k. */
}






int main(int argc, char ** argv)
{
  char programname[FILENAME_MAX]="";
  char locstr[255]="";

  double latitude  = 0;
  double longitude = 0;
  double long_std  = 0;
  double zenith    = 0;

  double hour1=0, hour2=0, tmp_long_std=0, utcplus=0;

  int doy=0, std=0;
  int hour=0, min=0, sec=0;

  int time1=0, time2=0;

  int status=0;
  int print =1;

  char dummy[9] = "";

  
  /* get command line options */
  status = get_options (argc, argv,
			programname,
			&doy, &print, &zenith,
			&latitude, &longitude, &long_std, &std, locstr);

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

  /* calculate time */
  if (zenith2time (doy, zenith, 
		   latitude, longitude, long_std, 
		   &time1, &time2) != 0)  {
    fprintf (stderr, "** Sorry, the specified zenith angle does not occur on this day! **\n");
    fprintf (stdout, "  %g  %g  %g\n", zenith, NAN, NAN);
    return (-1);
  }


  if (print) {
    fprintf (stderr, "\n  Latitude:             %7.2f degree",   latitude);
    fprintf (stderr, "\n  Longitude:            %7.2f degree",   longitude);
    fprintf (stderr, "\n  Standard Longitude:   %7.2f degree\n", long_std);

    utcplus = -long_std/15.0+1e-6;
    if (utcplus>=0)
      fprintf (stderr, "  Time Zone:            UTC+%.1f hours\n", utcplus);
    else 
      fprintf (stderr, "  Time Zone:            UTC-%.1f hours\n", -utcplus);

    fprintf (stderr, "  Solar Zenith  Angle:   %f deg\n\n", zenith);
  }

  hour1 = (double) time1 / 3600.0;
  hour2 = (double) time2 / 3600.0;

  hour   = time1 / 3600;
  time1 %= 3600;
  min    = time1 / 60;
  time1 %= 60;
  sec    = time1;

  if (print)
    fprintf (stderr, "  Standard time #1:      %s MEZ\n", time2str(dummy, hour, min, sec));

  hour   = time2 / 3600;
  time2 %= 3600;
  min    = time2 / 60;
  time2 %= 60;
  sec    = time2;

  if (print) {
    fprintf (stderr, "  Standard time #2:      %s MEZ\n", time2str(dummy, hour, min, sec));

    fprintf (stderr, "\n  Printing to stdout: zenith angle, morning time, afternoon time\n");
  }
  fprintf (stdout, "  %g  %g  %g\n", 
	   zenith, hour1, hour2); 

  return 0;
}
