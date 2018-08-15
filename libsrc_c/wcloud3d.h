/*--------------------------------------------------------------------
 * $Id: wcloud3d.h 3155 2015-08-10 10:48:25Z robert.buras $
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

#ifndef __wcloud3d_h
#define __wcloud3d_h

#if defined (__cplusplus)
extern "C" {
#endif
  
  
#include "mystic.h"

  /* prototypes */
int read_3D_caoth ( char      *filename, 
		    int       *Nx,
		    int       *Ny,
		    int       *Nz,
		    double    *delX,
		    double    *delY,
		    float    **z, 
		    float  ****lwc,
		    float  ****reff,
		    float  ****ext,
		    float  ****g1, 
		    float  ****g2, 
		    float  ****ff, 
		    float  ****f, 
		    float  ****ssa, 
		    float  ****dscale, 
		    float     *rmin,
		    float     *rmax,
		    int       *cldproperties,
		    int      **threed,
		    int        ixmin,
		    int        ixmax,
		    int        iymin,
		    int        iymax,
		    int        quiet );
  
int read_2D_albedo ( char       *filename, 
		     int       *Nx,
		     int       *Ny,
		     double    *delX,
		     double    *delY,
		     double  ***albedo,
		     int        ixmin,
		     int        ixmax,
		     int        iymin,
		     int        iymax,
		     int        quiet );

int read_2D_umu (char     *filename, 
		 int      *Nx,
		 int      *Ny,
		 double   *delX,
		 double   *delY,
		 double ***umu,
		 double ***phi,
		 int       ixmin,
		 int       ixmax,
		 int       iymin,
		 int       iymax,
		 int       quiet);

int read_2D_surface_labels (char *filename, 
			    char **label_lib, int nlabel_lib,
			    int *Nx, int *Ny,
			    double *delX, double *delY,
			    unsigned char ***label, int **nlabel,
			    int *therewereCaMs,
			    int quiet);

int read_2D_rossli (char *filename, 
		    int *Nx, int *Ny,
		    double *delX, double *delY,
		    rossli_brdf_spec ***rossli,
		    int *therewereCaMs,
		    int hotspot, int isAmbralsFile,
		    int ixmin, int ixmax, int iymin, int iymax,
		    int quiet);

int read_2D_elevation (char *filename, 
		       int *Nx, int *Ny,
		       double *delX, double *delY,
		       double ***elevation,
		       int quiet);


#if defined (__cplusplus)
}
#endif

#endif






