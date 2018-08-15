/*--------------------------------------------------------------------
 * $Id: sos.h 2626 2011-12-28 15:22:24Z svn-kylling $
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

#ifndef __sos_h
#define __sos_h

#if defined (__cplusplus)
extern "C" {
#endif

/* prototypes */
int sos (int nlyr, int newgeo,  int nstr, int nscat,
	 float albedo,
	 float radius, float *zd, float *ssalb, float **pmom,
	 float *dtauc, float zenang, int ntau, 
	 int numu, float *umu, float *utau,
	 float *rfldir, float *rfldn,  float *flup,
	 float *uavgso, float *uavgdn, float *uavgup,
	 float **u0u);

float chapmanfac(int lc, int *nfac, float **fac, float *dtauc );
void trans( int nlyr, float *chtau, float *trs);
void trans_double( int nlyr, double *chtau, double *trs);
int geocorfac(float **fac, int *nfac, float z_lay, int nlyr,
	      float *zd, float zenang, float r );


#if defined (__cplusplus)
}
#endif

#endif






