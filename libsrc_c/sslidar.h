/*--------------------------------------------------------------------
 * $Id: sslidar.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __sslidar_h
#define __sslidar_h

int ss_lidar (/* input: atmosphere */
	      int nlyr,           /* nlyr */
	      float *z,           /* = z-levels, including the zout-level */
	      float altitude,     /* = altitude at which surface is */
	      float *dtau,        /* = dtauc in uvspec */
	      float *omega_0,     /* = ssalb */
	      double **phase_back,/* = phase function in backward direction (theta=180degrees) */
	      double albedo,      /* = surface albedo */
	      /* input: lidar */
	      float lambda,       /* laser wavelength */
	      double E_0,         /* laser shot energy */
	      float zout,         /* lidar position */
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
	      );

#endif /* _SSLIDAR_H */

