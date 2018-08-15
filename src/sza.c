/*--------------------------------------------------------------------
 * $Id: sza.c 3231 2016-07-07 16:12:44Z bernhard.mayer $
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

#include "uvspec.h"
#include "ascii.h"
#include "f77-uscore.h"
#include <sun.h>
#include "solver.h"
#include "cdisort.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/************************************/
/* prototypes of internal functions */
/************************************/

static int read_sza (char* sza_filename, atm_inp_struct *atm_inp, atm_out_struct *atm_out);
static int get_sza_from_time_and_location (float lat, float lon, struct tm UTC, 
                                           float delta_time_start, float delta_time_end, float *sunshine_fraction,
                                           float *sza, float *phi0, int quiet, int verbose);

void F77_FUNC (qgaust, QGAUST) (int *NSTR, float *SZA, float *UMU0, int *RESULT);

/**************************************************************/
/* Setup solar zenith and azimuth.                            */
/**************************************************************/

int setup_sza (input_struct input, output_struct *output)
{
  int status=0, linear=0, iv=0;
  float umu0;
  float sza_tmp  = NOT_DEFINED_FLOAT;
  float phi0_tmp = NOT_DEFINED_FLOAT;

  char function_name[]="setup_sza";
  char file_name[]="sza.c";

  output->atm.sza_r  = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->atm.phi0_r = (float *) calloc (output->wl.nlambda_r, sizeof(float));
  output->sunshine_fraction = 1.0;

  if (input.source == SRC_SOLAR) {
    switch (input.atm.sza_source) {
    case SZA_BY_TIME_AND_LOCATION:

      /* if not all nessesary informations are given by the user   */
      /* this function will return the default values:             */
      /*            sza_tmp = 0, phi_tmp = 0                       */
      status = get_sza_from_time_and_location (input.latitude, input.longitude, input.UTC, 
                                               input.delta_time_start, input.delta_time_end, &(output->sunshine_fraction), 
                                               &(sza_tmp), &(phi0_tmp),
                                               input.quiet, input.verbose);

      if (status != 0) {
        fprintf (stderr, "Error %d returned by get_sza_from_time_and_location (line %d, function %s in %s) \n", 
                          status, __LINE__, __func__, __FILE__);
        return status;
      }

      for (iv=0;iv<output->wl.nlambda_r;iv++) {
        output->atm.sza_r [iv] = sza_tmp;    
        output->atm.phi0_r[iv] = phi0_tmp;    
      }

      if (input.verbose)
        fprintf (stderr, "     solar zenith angle from lat/lon/time, sza = %7.3f \n", output->atm.sza_r [0]); 

      break;

    case SZA_DIRECT_INPUT:

      for (iv=0;iv<output->wl.nlambda_r;iv++) {
        output->atm.sza_r [iv] = input.atm.sza;    
        output->atm.phi0_r[iv] = input.atm.phi0;    
      }

      if (input.verbose)
        fprintf (stderr, "     constant solar zenith angle, sza = %7.3f \n", output->atm.sza_r [0]); 

      break;
    case SZA_FROM_SZA_FILE:

      status = read_sza(input.filename[FN_SZA], &input.atm, &output->atm);

      if (status != 0) {
        fprintf (stderr, "Error %d reading sza file: %s\n", status, input.filename[FN_SZA]);
        return status;
      }

      /* Change last wavelength in sza file if noninteger wavelengths sza_file. */
      /* The user may ask for wl_end=420.25, which gives transmittance 421, and */
      /* which later give problems when interpolating sza to umu0, sighhh.      */ 

      if (output->atm.file_lambda_sza[output->atm.file_nlambda_sza-1] <
          output->wl.lambda_r[output->wl.nlambda_r-1])
        output->atm.file_lambda_sza[output->atm.file_nlambda_sza-1] =
	  output->wl.lambda_r[output->wl.nlambda_r-1];
    
      /**** Interpolate sza to internal wavelength grid ****/
    
      linear = 1;  /* changed to linear, BM March 29, 2001 */
      status = arb_wvn (output->atm.file_nlambda_sza, output->atm.file_lambda_sza,
		        output->atm.file_sza, 
		        output->wl.nlambda_r, output->wl.lambda_r, output->atm.sza_r, 
		        linear, 0);

      if (status!=0) {
        fprintf (stderr, "Error %d interpolating solar zenith angle\n", status);
        return status;
      }
    
      status = arb_wvn (output->atm.file_nlambda_sza, output->atm.file_lambda_sza,
		        output->atm.file_phi0, 
		        output->wl.nlambda_r, output->wl.lambda_r, output->atm.phi0_r, 
		        linear, 0);

      if (status!=0) {
        fprintf (stderr, "Error %d interpolating solar zenith angle\n", status);
        return status;
      }

      if (input.verbose)
        fprintf (stderr, "     wavelength dependent solar zenith angle, sza[iv=0] = %7.3f \n", output->atm.sza_r [0]); 

      break;
    
    case SZA_ECHAM:
      /* get the zenith angle from an ECHAM netcdf file */
      status = get_number_from_netCDF_map (input.latitude, input.longitude, input.UTC, input.atm.time_interpolate,
                                           input.filename[FN_SZA], &(sza_tmp), TYPE_FLOAT,
                                           "sza", input.verbose, input.quiet);
      if (status != 0) {
        fprintf (stderr, "Error %d returned by get_number_from_netCDF_map \n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }
    
      sza_tmp = acos(sza_tmp)*180.0/PI;
    
      /* Setup costant (not wavelength dependent) surface albedo */
      for (iv=0;iv< output->wl.nlambda_r;iv++){
        output->atm.sza_r[iv]=sza_tmp;
        /* If ECHAM za is given, phi is set to 0, assuming that you do only twostream 
           flux calculations where phi0 does not matter */
        output->atm.phi0_r[iv]=0.0;
      }
    
      if (input.verbose)
        fprintf (stderr, "ECHAM: solar zenith angle %g \n",
                 sza_tmp);

      break;
    
    default:
      fprintf (stderr, "\nError, determining solar zenith angle, wrong sza_source %d!\n",input.atm.sza_source);
      return -1;
    }
  

    /* test if one of the sza's coincides with a disort computational angle */
    switch (input.rte.solver) {
    case SOLVER_FDISORT1:
    case SOLVER_FDISORT2:
      status = F77_FUNC (dcheck, DCHECK)(&output->atm.nlyr, 
                                         &output->atm.nzout, &input.rte.nstr, 
                                         &input.rte.numu, &input.rte.nphi,
                                         &input.optimize_fortran, 
                                         &input.optimize_delta); 
      if (status!=0) {
        fprintf (stderr, "Error %d returned by dcheck in %s (%s)\n", 
                 status, function_name, file_name);
        return status;
      }	
      
      
      for (iv=0;iv<output->wl.nlambda_r;iv++) {
        umu0 = cos(output->atm.sza_r[iv]*PI/180.0);
        F77_FUNC  (qgaust, QGAUST) (&input.rte.nstr, &(output->atm.sza_r[iv]), &umu0, &status);
      }
      break;

    case SOLVER_DISORT:
      for (iv=0;iv<output->wl.nlambda_r;iv++) {
        umu0 = cos(output->atm.sza_r[iv]*PI/180.0);
	/*  c__gaussian_quadrature_test tests if umu0 is close to any cmu. If too close        */
        /*  it dithers umu0 and subsequently changes output->atm.sza_r.                        */
	status = c_gaussian_quadrature_test(input.rte.nstr, &output->atm.sza_r[iv],  umu0);
      }
      break;

    case SOLVER_SDISORT:
    case SOLVER_SPSDISORT:
    case SOLVER_MONTECARLO:
    case SOLVER_TWOSTR:
    case SOLVER_FTWOSTR:
    case SOLVER_RODENTS:
    case SOLVER_TWOSTREBE:
    case SOLVER_TWOMAXRND:
    case SOLVER_SOS:
    case SOLVER_POLRADTRAN:
    case SOLVER_TZS:
    case SOLVER_SSS:
    case SOLVER_SSSI:
    case SOLVER_SSLIDAR:
    case SOLVER_NULL:
      break;
    default:
      fprintf (stderr, "Error, unknown RTE solver %d\n", input.rte.solver);
      break;
    }
  }

  return 0;
}



/**************************************************************/
/* Read solar zenith angles from file.                        */
/**************************************************************/

static int read_sza (char* sza_filename, atm_inp_struct *atm_inp, atm_out_struct *atm_out) 
{
  int rows=0, min_columns=0, max_columns=0, status=0;
  float **sza_data=NULL;
  
  status = ASCII_file2float (sza_filename, 
			     &rows, &max_columns, &min_columns, 
			     &sza_data);
  if (status!=0) {
    fprintf (stderr, "Error %d opening file %s\n", status, sza_filename);
    return status;
  }

  if (max_columns!=min_columns) {
    fprintf (stderr, " !! ATTENTION !! Inconsistent number of columns\n");
    fprintf (stderr, "     min = %d, max =%d\n", min_columns, max_columns);
  }  

  if (min_columns<2) {
    fprintf (stderr, "Error, too few columns in %s\n", sza_filename);
    return -1;
  }
  
  atm_out->file_nlambda_sza = rows;
  
  atm_out->file_lambda_sza = ASCII_column_float(sza_data, atm_out->file_nlambda_sza, 0);
  atm_out->file_sza        = ASCII_column_float(sza_data, atm_out->file_nlambda_sza, 1);

  if (min_columns == 3)
    atm_out->file_phi0 = ASCII_column_float(sza_data, atm_out->file_nlambda_sza, 2);
  else		       
    atm_out->file_phi0 = calloc (atm_out->file_nlambda_sza, sizeof(float));
  
  ASCII_free_float (sza_data, rows);

  return 0;
}


static int get_sza_from_time_and_location (float lat, float lon, struct tm UTC, 
                                           float delta_time_start, float delta_time_end, float *sunshine_fraction,  
                                           float *sza, float *phi0, int quiet, int verbose)
{

  int status = 0;
  int Spencer=FALSE;  /* default is Blanco-Muriel et al. (2001) */
  int dt = NOT_DEFINED_INTEGER;
  int t = NOT_DEFINED_INTEGER;
  float cos_theta     = 0.0;
  float last_phi0     = 0.0;
  float sum_cos_theta = 0.0;
  float sum_phi0      = 0.0;
  int n_sun = 0;
  int n_time = 0;

  int time_std=0;

  struct cTime {
    int iYear;
    int iMonth;
    int iDay;
    double dHours;
    double dMinutes;
    double dSeconds;
  };

  struct cLocation {
    double dLongitude;
    double dLatitude;
  };

  struct cSunCoordinates {
    double dZenithAngle;
    double dAzimuth;
  };

  struct cTime udtTime;
  struct cLocation udtLocation;
  struct cSunCoordinates udtSunCoordinates;

  void sunpos(struct cTime udtTime, struct  cLocation udtLocation, struct cSunCoordinates *udtSunCoordinates);

  char *time_str;
  int hour_lat = 0, min_lat = 0, sec_lat = 0, time_lat = 0;
  double utcplus  = 0.0;
  double lon_std = 0.0;
  float lon_tmp = NOT_DEFINED_FLOAT;
  float EPSILON = 1.E-6;
  char dummy[10] = "";

  int seconds_per_day = 86400;

  /* calculate solar zenith angle */
  if (verbose)
    fprintf (stderr, " ... calculating sza with time, latitude and longitude \n");

  if ( fabs(lon - NOT_DEFINED_FLOAT) < EPSILON || fabs(lat - NOT_DEFINED_FLOAT) < EPSILON ) {
    if (verbose) {
      fprintf (stderr, "     latitude or longitude not defined, use default value for position of the sun \n");
      fprintf (stderr, "     sza = 0.0 , phi0 = 0.0 \n");
    }
    (*sza)  = 0.0;
    (*phi0) = 0.0;
    return 0;
  }
  else {

    if ( UTC.tm_hour == NOT_DEFINED_INTEGER ) {
      if (verbose) {
        fprintf (stderr, "     time not defined, use default value for position of the sun \n");
        fprintf (stderr, "     sza = 0.0 , phi0 = 0.0 \n");
      }
      (*sza)  = 0.0;
      (*phi0) = 0.0;
      return 0;
    }

    /* move longitude to range -180.0 ... 180.0 */
    if      ( lon < - 180.0 ) 
      lon_tmp = lon + 360.0;
    else if (lon >  + 180.0 )
      lon_tmp = lon - 360.0;
    else
      lon_tmp = lon;

    /* check range of lat and lon */
    if ( lon_tmp < -180.0 || lon_tmp > +180.0 || lat < -90 || lat > 90 ) {
      fprintf (stderr, "Error, latitude or longitude out of range for solar zenith angle calculations \n");
      fprintf (stderr, "       latitude = %8.3f degree, longitude = %8.3f degree \n", lat, lon_tmp);
      return -1;
    }

    /* calculate standard time [sec] and check range */
    if ((time_std = UTC.tm_hour*3600 + UTC.tm_min*60 + floor(UTC.tm_sec)) > seconds_per_day)  {
      fprintf (stderr, "Error, calculating time in s, (in sza.c)\n");
      fprintf (stderr, "    Time out of range **\n");
      return -2;
    }


    for (dt=(int)delta_time_start;dt<=(int)delta_time_end;dt=dt+60) {

      if ( Spencer ) {
        if (verbose && dt == -(int)delta_time_start )
          fprintf (stderr, "     using Spencer (1971) to calculate solar zenith and azimuth\n");
        /* attention, Spencer defines western longitudes as positive, therefor (-longitude) */
        t=time_std+dt;
        if (        t        < 0 ) t = t + seconds_per_day;
        if ( seconds_per_day < t ) t = t - seconds_per_day;
        (*sza)  = solar_zenith  (t, UTC.tm_yday, lat, -lon_tmp, lon_std);
        (*phi0) = solar_azimuth (t, UTC.tm_yday, lat, -lon_tmp, lon_std);  
      }
      else {
        if (verbose && dt == -(int)delta_time_start )
          fprintf (stderr, "     using Blanco-Muriel et al. (2001) to calculate solar zenith and azimuth\n");
        /*  if ( wvlsza==5) { */
        /*    tmptime = (double) time_std/3600.0; */
        /*    dectime2hhmmss(tmptime, &hour, &min, &sec); */
        /*  } */
        
        udtTime.iYear    = UTC.tm_year + 1900; /* C-standard is years after 1900*/
        udtTime.iMonth   = UTC.tm_mon  + 1;    /* C-standard is month after January [0 ... 11] */
        udtTime.iDay     = UTC.tm_mday;
        udtTime.dHours   = UTC.tm_hour;
        udtTime.dMinutes = UTC.tm_min;
        udtTime.dSeconds = UTC.tm_sec + dt;    /* all the time shift is placed in the second variable */
        
        udtLocation.dLatitude   = lat;
        udtLocation.dLongitude  = lon_tmp - lon_std;
        
        sunpos( udtTime, udtLocation, &udtSunCoordinates);
        (*sza)  = udtSunCoordinates.dZenithAngle;
        (*phi0) = udtSunCoordinates.dAzimuth-180.0;    
      }

      /* take care of periodicity of the solar azimuth */
      if (dt==(int)delta_time_start)
        last_phi0=(*phi0);
      else {
        if ( 350.0 < fabs((*phi0) - last_phi0) ) { /* jump over the periodic boundary */
          if ( (*phi0) - last_phi0 < 0.0 )
            (*phi0) = (*phi0) + 360.0;
          else 
            (*phi0) = (*phi0) - 360.0;
        }
        last_phi0 = (*phi0);
      }

      cos_theta = cos(PI/180.*(*sza));
      if ( cos_theta > 0.0 ) {
        sum_cos_theta += cos_theta;
        sum_phi0      += cos_theta * (*phi0);
        n_sun++;
      }
      n_time++;

    } /* for (dt=(int)delta_time_start;dt<=(int)delta_time_end;dt=dt+60) { */

    /* get averaged position of the sun */
    if ( n_sun != 0 ) {
      (*sza)  = 180.0 / PI * acos ( sum_cos_theta / n_sun );
      (*phi0) = sum_phi0 / sum_cos_theta;
    }
    else {
      /* all sun positions below horizon */
      (*sza)  = 180.0;
      (*phi0) =   0.0;
    }
    /* fraction of the total time, where the sun is shining */
    if (delta_time_start != 0.0 || delta_time_end != 0.0)          /* only for input of time intervals */
      (*sunshine_fraction) = ((float) n_sun) / ((float) n_time);
    else 
      (*sunshine_fraction) = 1.0;                                  /* otherwise original result        */

    if ((*phi0) <   0.) (*phi0) += 360.0;
    if ((*phi0) > 360.) (*phi0) -= 360.0;

    if (verbose) {
      fprintf (stderr,   "     Latitude:             %7.2f degree",   lat);
      fprintf (stderr, "\n     Longitude:            %7.2f degree",   lon_tmp);
      fprintf (stderr, "\n     Standard Longitude:   %7.2f degree\n", lon_std);
      
      utcplus = - lon_std / 15.0 + 1e-6;  /* 1h == 15 degrees longitude */
      if (utcplus>=0)
	fprintf (stderr, "     Time Zone:            UTC+%.1f hours\n", utcplus);
      else 
	fprintf (stderr, "     Time Zone:            UTC-%.1f hours\n", -utcplus);
      
      time_str = asctime( &(UTC) );
      fprintf (stderr, "     Standard Time:        %s", time_str);

      /* calculate LOCAL APPARENT TIME */
      time_lat = LAT (time_std, UTC.tm_yday, -lon_tmp, lon_std);
      /* correct negative LATs and LATs larger than 24:00 */
      /* to get human-readable output                     */
      while (time_lat<0)
        time_lat += seconds_per_day;    
      while (time_lat>seconds_per_day)
        time_lat -= seconds_per_day;
    
      hour_lat = time_lat / 3600;                           /* 3600 == s -> hours, attention INTEGER arithmetic */
      min_lat  = (time_lat - hour_lat*3600) / 60;           /*   60 == s -> min  */
      sec_lat  = (time_lat - hour_lat*3600 - min_lat*60);   

      fprintf (stderr, "     Local Apparent Time:             %s\n", time2str(dummy, hour_lat,min_lat,sec_lat) );
      if (delta_time_start > 0.0 || delta_time_end > 0.0 )
        fprintf (stderr, "     time_interval to average: [%6.1f min,%6.1f min]\n", delta_time_start/60.0, delta_time_end/60.0 );
      fprintf (stderr, "     Solar Zenith  Angle:  %8.3f degree\n", (*sza));
      fprintf (stderr, "     Solar Azimuth Angle:  %8.3f degree\n", (*phi0));
      if (delta_time_start > 0.0 || delta_time_end > 0.0 )
        fprintf (stderr, "     Sunshine fraction:    %8.3f per cent \n", (*sunshine_fraction) * 100.0);
    }
  }

  return status;

}



