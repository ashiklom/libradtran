/*--------------------------------------------------------------------
 * $Id: sslidar.c 2623 2011-12-23 10:52:38Z robert.buras $
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

#include "sslidar.h"
#include "locate.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

int ss_lidar (/* input: atmosphere */
	      int nlyr,           /* nlyr */
	      float *z,           /* = z-levels */
	      float altitude,     /* = altitude at which surface is */
	      float *dtau,        /* = dtauc in uvspec */
	      float *omega_0,     /* = ssalb */
	      double **phase_back,/* = phase function in backward direction (theta=180degrees) */
	      double albedo,      /* = surface albedo */
	      /* input: lidar */  
	      float lambda,       /* laser wavelength */
	      double E_0,         /* laser shot energy */
	      float position,     /* lidar position */
	      double umu,         /* cos(zenith angle) */
	      int nranges,        /* nranges */ 
	      double range_width, /* range bin width */
	      float *ranges,      /* particular ranges */
	      double efficiency,  /* efficiency */
	      double area,        /* detector area */
	      int polarisation,   /* whether or not polarisation */
	      /* OUTPUT / RESULT */
	      float *nphot,       /* number of photons detected per range bin */
	      float *nphot_q,     /* Q polarisation of nphot */
	      float *lidar_ratio  /* lidar ratio for each range bin */
	      )
{
  const double c_heisenberg = 6.62606896e-34; /* BCA unify all constants in all files */
  const double c_speed_of_light=2.99792458e8;

  double nphot_shot=0.0, performance=0.0,
    beta=0.0, alpha_dr=0.0, E_phot=0.0;
  float *z_loc=NULL;
  float zact=0.0, zact_last=0.0, range_border=0.0;
  double *alpha=NULL;
  double *lidar_ratio_atm=NULL, *lidar_ratio_q_atm=NULL;
  int ir=0, ia=0, upward=0, ia_last=0, iaa=0;
  int hit_surface=0;

  /* z [km], ranges [km], area [m^2], position [km], lambda [nm] */

  z_loc=calloc((size_t) nlyr+1, sizeof(float));
  for (ia=0;ia<nlyr;ia++)
    z_loc[ia] = z[ia]+altitude;
  z_loc[ia]=altitude;

  /* calculate alpha for all layers */
  alpha = calloc(nlyr, sizeof(double));
  for (ia=0; ia<nlyr; ia++)
    alpha[ia] = dtau[ia] / ( z_loc[ia]-z_loc[ia+1] );

  /* calculate photon energy */
  E_phot = c_heisenberg * c_speed_of_light / lambda * 1e9; /* conversion wavelength nm to m */

  /* number of photons per laser shot */
  nphot_shot = E_0 / E_phot;

  /* performance contains all constants of lidar equation */
  performance = nphot_shot * efficiency * area * range_width;

  lidar_ratio_atm   = calloc(nlyr, sizeof(double));
  lidar_ratio_q_atm = calloc(nlyr, sizeof(double));

  /* sum up phase_backs if polarisation is on; this assumed Q=1 for laser */
  if (polarisation)
    for (ia=0; ia<nlyr; ia++) {
      lidar_ratio_atm   [ia] = 4.*PI / ( phase_back [ia][0] + phase_back [ia][1] ) / omega_0[ia];
      lidar_ratio_q_atm [ia] = 4.*PI / ( phase_back [ia][4] - phase_back [ia][1] ) / omega_0[ia];
    }
  else
    for (ia=0; ia<nlyr; ia++)
      lidar_ratio_atm   [ia] = 4. * PI / phase_back [ia][0] / omega_0[ia];

  /* whether lidar is pointing upward or downward */
  upward = (umu<0.0);

  /* initialize position of last range bin evaluation */
  zact_last = position;
  ia_last   = flocate ( z_loc, nlyr, zact_last );
  alpha_dr = 0.0;

  for (ir=0; ir<nranges+1; ir++) { /* +1 because we want to go to end
				      of last range to check whether
				      surface is contained */

    if (ir<nranges)
      zact = position - umu * ranges[ir];
    else
      /* special case, calculate alpha_dr till end of last range bin */
      zact = position - umu * (ranges[ir-1] + 0.5 * range_width );
    ia = flocate ( z_loc, nlyr+1, zact );
    if (ia <= nlyr && ia >= 0) {

      if (ia == ia_last || umu==0.0) {
	/* center of range bin in same atmospheric layer as center of
	   last range bin */
	alpha_dr += ( zact_last - zact ) / umu * alpha[ia];
	if (ia==nlyr)
	  /* range is below ground, therefore hits surface; up to here
	     alpha_dr was integrated down to the surface */
	  hit_surface=1;

	if (ir==nranges) {
	  /* we check whether last range bin contains surface, then exit loop */
	  if (zact==altitude)
	    hit_surface=1;
	  break;
	}
      }
      else {
	/* sum up extinction in atmospheric layer where last range bin
	   was evaluated */
	alpha_dr += ( zact_last - z_loc[ia_last+1-upward] ) / umu * alpha[ia_last];
	/* sum up extinction in atmospheric layers in between */
	if (upward)
	  for (iaa=ia_last-1;iaa>ia;iaa--) {
	    alpha_dr += ( z_loc[iaa+1] - z_loc[iaa] ) / umu * alpha[iaa];
	  }
	else
	  for (iaa=ia_last+1;iaa<ia;iaa++) {
	    alpha_dr += - ( z_loc[iaa+1] - z_loc[iaa] ) / umu * alpha[iaa];
	  }
	if (ia==nlyr) {
	  /* range is below ground, therefore hits surface; up to here
	     alpha_dr was integrated down to the surface */
	  hit_surface=1;
	  /* ia==nlyr is below ground, take lowest atmospheric layer */
	  ia--;
	}
	if (ir==nranges) {
	  /* we check whether last range bin contains surface, then
	     exit loop */
	  if (zact==altitude)
	    hit_surface=1;
	  break;
	}

	if (!hit_surface)
	  /* sum up extinction in atmospheric layer where current
	     range bin will be evaluated */
	  alpha_dr += ( z_loc[ia+upward] - zact ) / umu * alpha[ia];
      } /* end else */

      /* evaluate signal */
      lidar_ratio [ir] = lidar_ratio_atm[ia];
      beta = alpha[ia] / lidar_ratio[ir];
      nphot[ir] = performance / ( ranges[ir] * ranges[ir] * 1e6 ) /* conversion km to m */
	* beta * exp ( - 2.0 * alpha_dr );

      /* in case we calculate polarisation, calculate Q part */
      if (polarisation) {
	beta = alpha[ia] / lidar_ratio_q_atm[ia];
	nphot_q[ir] = performance / ( ranges[ir] * ranges[ir] * 1e6 ) /* conversion km to m */
	  * beta * exp ( - 2.0 * alpha_dr );
      }

      if (hit_surface)
	/* range is below ground, stop loop; up to here alpha_dr was
	   integrated down to the surface */
	break;

      zact_last=zact;
      ia_last=ia;

      if (ir==nranges)
	break;
    }
    else
      break;
  } /* end for */

  /* surface */
  /* we assume equidistant range bins! */
  if (hit_surface) {
    zact = ( position - altitude ) / umu;
    ir = flocate ( ranges, nranges, zact );
    range_border = ranges[ir] + 0.5 * range_width;
    if ( zact > range_border )
      /* surface is in next range bin */
      ir++;
    else
      /* surface is in this range bin, reduce range_border */
      range_border -= range_width;
    /* scale atmospheric backscatter of current range bin, since only
       part of the range bin contains atmosphere */
    nphot[ir] *= ( zact - range_border ) / range_width;
    /* performance is different for surface reflection */
    performance = nphot_shot * efficiency * area;
    /* lambertian surface reflection */
    beta = albedo * 4. * umu / ( 4. * PI );
    nphot[ir] += performance / ( zact * zact * 1e6 ) /* conversion km to m */
	* beta * exp ( - 2.0 * alpha_dr );
  }
  
  /* outside of atmosphere */
  for (ir=ir+1; ir<nranges; ir++) {
    lidar_ratio[ir]=0.0;
    nphot[ir]=0.0;
  }

  free(alpha);
  free(z_loc);
  free(lidar_ratio_atm);
  free(lidar_ratio_q_atm);

  return 0;
}
